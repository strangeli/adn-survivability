#using Pkg
#Pkg.activate(".")

using Distributed
using ClusterManagers

SLURM = true

if length(ARGS) > 0
	N_workers = parse(Int, ARGS[1]) - 1
else
	N_workers = 1
end

if SLURM
    using ClusterManagers
	addprocs(SlurmManager(N_workers))
else
	addprocs(N_workers)
end
println(nprocs(), " process(es)")

@everywhere dir = @__DIR__
##

@everywhere begin
    using Pkg
    #Pkg.activate(".")

    using PowerDynamics
    using OrdinaryDiffEq
    using Distributions
	using Measurements

    using Random
    Random.seed!(42)
end

## for plotting/saving. don't need them everywhere
begin
    using Dates
    using CSV
    using DataFrames
    #using Query
    #using StatsPlots
end

## set per unit MV
@everywhere begin
    const base_power = 1E6 # 1MW
    const base_voltage = 20E3 # 20kV
    const base_current = base_power / base_voltage # 50A
    const base_admittance = base_power / base_voltage^2 # 0.0025Ω^-1
    const ω = 2 * π * 50.0 # 314.1593rad/s
    # per unit HV
    const base_voltage_HV = 110E3 # 110kV
    const base_admittance_HV = base_power / base_voltage_HV^2 # 8.264462809917356e-5
    const base_current_HV = base_power / base_voltage_HV
end

## some additional functions that might develop into new PD.jl features
@everywhere begin
    # ADN control layer
    include("$dir/control.jl")
    # custom types
    include("$dir/DGUnit.jl")
    include("$dir/OLTC.jl")
    # load actual data
    include("$dir/cigre_static.jl")
    # failure model
    include("$dir/short_circuit.jl")
    # static solution
    include("$dir/power_flow.jl")
end

##@everywhere Pkg.status()

@everywhere function FRT_curve(t; t_error=0.1, Usoll=1.)
    FRT_lower_limit = 0.9 *  Usoll
    FRT_upper_limit = 1.1 *  Usoll

    if t_error > t
        FRT_lower_limit = 0.9 *  Usoll
        FRT_upper_limit = 1.1 *  Usoll
    end

    #TODO: more elegant way for this? (smooth step not needed since this function runs on a ready solution object)
    # ]0, 0.1]
    if t > t_error && t <=  0.1 + t_error
        FRT_lower_limit = 0.15 *  Usoll
        FRT_upper_limit = 1.25 *  Usoll
    # ]0.1, 0.15]
    elseif t > 0.1 + t_error && t <= 0.15 + t_error
        FRT_lower_limit = 0.15 *  Usoll
        FRT_upper_limit = 1.2 * Usoll
    # ]0.15, 3]
    elseif t > 0.15 + t_error && t <= 3. + t_error
        FRT_lower_limit = ((0.85 - 0.15) * Usoll /(3. - 0.15)) * (t - t_error) + (0.15 * Usoll  - ((0.85 - 0.15) * Usoll /(3. - 0.15)) * 0.15)  # 43/380
        FRT_upper_limit = 1.2 * Usoll
    # ]3, 5]
    elseif t > 3. + t_error && t <= 5. + t_error
        FRT_lower_limit = 0.85 *  Usoll
        FRT_upper_limit = 1.2 * Usoll
    # ]5, 60]
    elseif t > 5. + t_error && t <= 60. + t_error
        FRT_lower_limit = 0.85 *  Usoll
        FRT_upper_limit = 1.15 *  Usoll
    # ]60, Inf]
    elseif t > 60. + t_error
        FRT_upper_limit = 1.1 *  Usoll
        FRT_lower_limit = 0.9 *  Usoll
    end

    (FRT_lower_limit, FRT_upper_limit)
end

@everywhere function CIGRE_prob(
	S_total, # S_total
    P_ref, #Pref
	Q_ref; # Qref
    nsc_node = 4,
    t_duration = 0.15,
    resistance = 10.0,
    quad = 0.,
    verbose = false,
    get_cpg = false,
	tspan = (0.0, 0.5), # total duration of a simulation
	initial_state = nothing,
)
    ## simulation parameter
    t_fault = 0.1 # onset of node short circuit
    Znsc = complex(resistance, 0.0) * base_admittance #  fixed by Emocosy  Ω
    @assert t_fault + t_duration < last(tspan)

    nsc = NodeShortCircuitAlt(;
        node_number = nsc_node,
        Y = (1 / Znsc),
        tspan_fault = (t_fault, t_fault + t_duration),
    )

    ## solve power flow
    busses_static, lines, T, elist, Zs, Yshs = CIGRE_static()
    pg_static = PowerGrid(busses_static, lines)
    power_flow = pf_sol(pg_static, initial_guess(pg_static), nothing)

    ## construct ADN
    busses = copy(busses_static)
    DG_locs = 2:12
    for i in DG_locs
        S_bkgrnd = zero(im)
        try
            S_bkgrnd = busses_static[i].S
        catch
            S_bkgrnd = complex(busses_static[i].P, busses_static[i].Q)
        end
        busses[i] = DGUnit(;
            K_pll = 1632.993, #Hz/pu
            T_int = 2.0,
            K_PT1 = 1.0,
            T_PT1 = 1e-8,
            K_FRT = 2.0,
            I_max = 1.0, # pu
            P_limit = 1.0, #pu
            Q_limit = 1.0, #pu
            Vref = 1.0, # pu
            Vdead = 0.1, # pu
            S_pq = V -> S_bkgrnd * (quad * V^2 + 1 - quad),
            Y_n = 0.0,
        )
    end
    pg = PowerGrid(busses, lines)
    cpg = ADN(pg, DGUnit, P_ref, Q_ref)

	if isnothing(initial_state)
	    ## find operation point
	    icguess_pf = initial_guess(cpg, power_flow[:, :u])
	    op = find_steady_state(cpg, icguess_pf)

	    verbose ? check_operationpoint(cpg, op) : nothing
	else
		op = State(cpg, initial_state)
	end


    ## define callbacks for nsc
    nsc_powergrid = nsc(cpg)
    nsc_ode = rhs(nsc_powergrid)
    ode = rhs(cpg)

    function errorState(integrator)
        sol1 = integrator.sol
        verbose ? println("searching fault state") : nothing
        #x2 = find_valid_initial_condition(nsc_powergrid, sol1[end]) # Jump the state to be valid for the new system.
        integrator.f = nsc_ode
        DiffEqBase.initialize_dae!(integrator, BrownFullBasicInit())
        verbose ?
        println("found fault state\n", State(cpg, integrator.u)[:, :v]) :
        nothing
        #integrator.u = x2
        #u_modified!(integrator,true)
    end

    function regularState(integrator)
        sol2 = integrator.sol
        verbose ? println("searching post-fault state") : nothing
        #x3 = find_valid_initial_condition(powergrid, sol2[end]) # Jump the state to be valid for the new system.
        integrator.f = ode
        DiffEqBase.initialize_dae!(integrator, BrownFullBasicInit())
        verbose ?
        println("found post-fault state\n", State(cpg, integrator.u)[:, :v]) :
        nothing
        #integrator.u = x3
        #u_modified!(integrator,true)
    end

    # discrete callback requires tstop
    cb1 = DiscreteCallback(
        ((u, t, integrator) -> t in nsc.tspan_fault[1]),
        errorState,
    )
    cb2 = DiscreteCallback(
        ((u, t, integrator) -> t in nsc.tspan_fault[2]),
        regularState,
    )

    # return new ODE problem
    if get_cpg
        return cpg,
        ODEProblem(
            ode,
            op.vec,
            tspan,
           callback = CallbackSet(cb1, cb2),
           tstops = [t_fault, t_fault + t_duration],
        )
    else
        return ODEProblem(
            ode,
            op.vec,
            tspan,
           callback = CallbackSet(cb1, cb2),
           tstops = [t_fault, t_fault + t_duration],
        )
    end
end


## check prob func
#
# @time cpg, prob = CIGRE_prob(;
#     P_ref = 24.0,
#     Q_ref = 5.0,
#     nsc_node = 4,
#     t_duration = 0.15,
#     resistance = 3.5,
#     quad = 0.,
#     get_cpg = true,
#     verbose = true,
# )
#
# @time _sol = solve(prob, Rodas5()) #ode23s = Rosenbrock23
#
# sol = PowerGridSolution(_sol, cpg)
#
# # plot(_sol, legend=false)
# # final_state = sol(prob.tspan[2])
# # check_operationpoint(cpg, final_state)
# # op = find_steady_state(cpg)
#
# using Plots
#
# T = range(prob.tspan...; step=0.00001)
#     plot(
#         T,
#         zeros(length(T)),
#         c = :gray,
#         legend = false,
#         ribbon = (-first.(FRT_curve.(T)),
#                     last.(FRT_curve.(T))
#                 ),
#     )
#     plot!(T, sol(T, 2:12, :v)', legend = false)
# #xlims!(0.09, 0.26)
# #plot!(sol.dqsol.t, first.(FRT_curve.(sol.dqsol.t)), c=:black)
# #plot!(sol.dqsol.t, last.(FRT_curve.(sol.dqsol.t)), c=:black)
#
# plot(sol, 2:12, :φ, c = :blue, legend = false, title="PLL")
#     plot!(sol, 2:12, :θ, c = :green, legend = false)
#
# plot(sol, 2:12, :P_g, legend = false)
#     plot!(sol, 2:12, :P_int, legend = false)
#     vline!(prob.kwargs[:tstops], c = :black, ls = :dash)
#
# plot(sol, 2:12, :Q_g, legend = false)
#     plot!(sol, 2:12, :Q_int, legend = false)
#     vline!(prob.kwargs[:tstops], c = :black, ls = :dash)
#
# # control error
# Δ = [
#     first(cpg.controller(u, nothing, t))[1:2]
#     for (t, u) in zip(sol.dqsol.t, sol.dqsol.u)
# ]
#     plot(
#     sol.dqsol.t,
#     first.(Δ),
#     label = "ΔP",
#     xlabel = "t [s]",
#     ylabel = "error [pu]",
#     )
#     plot!(sol.dqsol.t, last.(Δ), label = "ΔQ")
#
# # FRT currents
# DG = cpg.nodes[2] # take one, they are all identitical
#     plot(
#         sol.dqsol.t,
#         FRTcurrent.(sol(sol.dqsol.t, 2:12, :v), DG.K_FRT, DG.Vref, DG.Vdead)',
#         ylabel = "iqplus",
#         legend = false,
#     )
#
# # Dynamic Load
# DG = cpg.nodes[2] # take one, they are all identitical
#     plot(
#         sol.dqsol.t,
#         real.(DG.S_pq.(sol(sol.dqsol.t, 2:12, :v)))',
#         ylabel = "Load P",
#         legend = false,
#     )
#
# plot(
#     sol.dqsol.t,
#     imag.(DG.S_pq.(sol(sol.dqsol.t, 2:12, :v)))',
#     ylabel = "Load Q",
#     legend = false,
# )

## setup EnsembleProblem step
@everywhere S_total = 	complex(24.373141067954748, 6.115974820490637) * 1e6 / base_power
@everywhere P_ref(t) = t > 0.25 ? 1.1 * real(S_total) : real(S_total)
@everywhere  Q_ref(t) = t > 0.25 ? 1.1 * imag(S_total) : imag(S_total)


@everywhere cpg_step, prob_step = CIGRE_prob(
	S_total,
	t -> real(S_total),
	t -> imag(S_total);
    nsc_node = 4,
    t_duration = 0.,
    resistance = 3.5,
    quad = 0.,
    get_cpg = true,
    verbose = false,
	tspan = (0.,50.),
)

ode = remake(prob_step, callback = nothing)
dqsol = solve(ode, Rodas4())
sol = PowerGridSolution(dqsol, cpg_step)
Δ = [dqsol.prob.f.f.cpg.controller(u, nothing, t) |> first for (t, u) in zip(dqsol.t, dqsol.u)]
using Plots

plot(sol, 2:12, :P_g, legend=false)
plot(dqsol.t, first.(Δ))

# debug
# begin
# 	op = prob_step.u0
# 	pg = PowerGrid(cpg_step)
# 	step_rhs = ADN(pg, DGUnit, P_ref, Q_ref) |> rhs
# 	ode = ODEProblem(step_rhs, op, (0., 5.))
# 	sol = solve(ode, Rodas4())
# 	volt = abs.(complex.(sol[Int.(sol.prob.f.f.cpg.controller.args[6]),:] , sol[Int.(sol.prob.f.f.cpg.controller.args[7]),:]))'
# 	using Plots
# 	plot(sol.t, volt, label = false, xlabel="t [s]", ylabel="v [pu]")
# end
# rest
@everywhere begin
    # load flow sol
	# S_total = complex(24.373141067954748, 6.115974820490637) * 1e6 / base_power
	NT_step =  1_000_000
    Pdist = Uniform(-15., 15.)
	Qdist = Uniform(-10., 10.)
	#Wdist = Uniform(-π, π) #Uniform(angle(S_total)-π, angle(S_total)+π)
	P = real(S_total) .+ rand(Pdist, NT_step)
	Q = imag(S_total) .+ rand(Qdist, NT_step)
	# factor = rand(NT_step) # 2 factors for P and Q
	# Pfactor
	# Qfactor
end

@everywhere function prob_func_step(prob, i, repeat)
	Preff(t) = t > 0.25 ? P[i] : real(S_total)
	Qreff(t) = t > 0.25 ? Q[i] : imag(S_total)
	#op = prob_step.u0
	pg = PowerGrid(cpg_step)
	step_rhs = ADN(pg, DGUnit, Preff, Qreff) |> rhs
	return remake(prob_step, f = step_rhs, callback = nothing)
end


# TODO: check that :DTmin errors also don't survive
@everywhere begin
    r_idx = ones(Int, length(cpg_step.nodes))
    for n in 2:length(cpg_step.nodes)
        r_idx[n] = r_idx[n-1] + dimension(cpg_step.nodes[n-1])
    end

    i_idx = ones(Int, length(cpg_step.nodes)) .+ r_idx

    function observer_step(sol, i)
		state = sol[end]
        v = sqrt.(state[r_idx].^2 .+ state[i_idx].^2)
		#println(v)
		in_bounds_end = all(0.9 .< v .< 1.1) # change
		# this does not work unfortunately?
		Δ = sol.prob.f.f.cpg.controller(state, nothing, last(sol.t)) |> first
		ΔP = Δ[1]
		ΔQ = Δ[2]
        return ([P[i], Q[i], in_bounds_end, sol.retcode, minimum(v), maximum(v), ΔP, ΔQ, state], false)
    end
end

mcprob_step = EnsembleProblem(
    prob_step;
    prob_func = prob_func_step,
    output_func = observer_step,
    u_init = [],
)

sim_step = solve(mcprob_step, Rodas4(), EnsembleDistributed(), trajectories=NT_step)
@show sim_step.elapsedTime

result_step = DataFrame(
    P = [u[1] for u in sim_step.u],
    Q = [u[2] for u in sim_step.u],
	step_surv = [u[3] for u in sim_step.u],
	step_success = categorical([u[4] for u in sim_step.u]),
    step_min_v = [u[5] for u in sim_step.u],
	step_max_v = [u[6] for u in sim_step.u],
	ΔP = [u[7] for u in sim_step.u],
	ΔQ = [u[7] for u in sim_step.u],
)

	# check bounds detection
	#using StatsPlots
	#@df step_result scatter(:min_v, :max_v, group=:inbounds)
	#@df step_result scatter(:P, :Q, group=:inbounds)


_valid_idx_surv = findall(result_step.step_surv .& (result_step.step_success.== :Success))

error_P = abs.(result_step.ΔP ./ result_step.P) .<= 0.05
error_Q = abs.(result_step.ΔQ ./ result_step.Q) .<= 0.05

_valid_idx_error = findall(error_P .& error_Q)

valid_P_surv = [u[1] for u in sim_step.u][_valid_idx_surv]
valid_P_error = [u[1] for u in sim_step.u][_valid_idx_error]

valid_Q_surv = [u[2] for u in sim_step.u][_valid_idx_surv]
valid_Q_error = [u[2] for u in sim_step.u][_valid_idx_error]


_valid_idx_error_surv = findall(error_P .& error_Q .& result_step.step_surv .& (result_step.step_success.== :Success))
valid_states_after_step = [last(u) for u in sim_step.u][_valid_idx_error]

@everywhere valid_states_after_step = $valid_states_after_step
@everywhere valid_P = $valid_P_error
@everywhere valid_Q = $valid_Q_error


@everywhere cpg, prob = CIGRE_prob(
	S_total,
	t->real(S_total), # function
	t->imag(S_total); # function
    nsc_node = 4,
    t_duration = 0.15,
    resistance = 3.5,
    quad = 0.,
    get_cpg = true,
    verbose = true,
)



@everywhere begin
    Tdist = Normal(150e-3, 10e-3)
    Rdist = Uniform(3., 10.)

    NT = length(valid_states_after_step)#1000 #40_000

    T = rand(Tdist, NT)
    R = rand(Rdist, NT)
end

@everywhere function prob_func(prob, i, repeat)
    return CIGRE_prob(
		S_total,
        t ->  valid_P[i],
        t ->  valid_Q[i];
        nsc_node = 4,
        t_duration = T[i],
        resistance = R[i],
        quad = 0.,
        get_cpg = false,
        verbose = false,
		initial_state = valid_states_after_step[i],
    )
end

# TODO: check that :DTmin errors also don't survive
@everywhere begin
    r_idx = ones(Int, length(cpg.nodes))
    for n in 2:length(cpg.nodes)
        r_idx[n] = r_idx[n-1] + dimension(cpg.nodes[n-1])
    end
    i_idx = ones(Int, length(cpg.nodes)) .+ r_idx

    t_eval = range(prob.tspan...; step=0.001)
    lower_lim = first.(FRT_curve.(t_eval))
    upper_lim = last.(FRT_curve.(t_eval))

    function observer(sol, i)
        v = sqrt.(sol(t_eval, idxs=r_idx).^2 .+ sol(t_eval, idxs=i_idx).^2)
        in_bounds = [all(lower_lim[i] .< v[:, i] .< upper_lim[i]) for i in 1:length(t_eval)]
        surv_idx = findfirst(.!in_bounds)
        surv_time = isnothing(surv_idx) ? t_eval[end] : t_eval[surv_idx]
        survived = all(in_bounds) # for all time steps
        return ([survived, surv_time, sum(v[2:12,end])/11, sum(minimum(v[2:12,:], dims=2))/11, sol.retcode], false)
    end
end

mcprob = EnsembleProblem(prob; prob_func = prob_func, output_func = observer, u_init=[])

sim = solve(mcprob, Rodas4(), EnsembleDistributed(), trajectories=NT)
@show sim.elapsedTime

result = DataFrame(
	R = R,
	T = T,
	P = valid_P,
	Q = valid_Q,
	nsc_surv = first.(sim.u),
	nsc_surv_t = [u[2] for u in sim.u],
	nsc_fin_v = [u[3] for u in sim.u],
	nsc_min_v = [u[4] for u in sim.u],
	nsc_success = categorical(last.(sim.u)),
)

now = unix2datetime(time())
CSV.write("$dir/mc_res_$(Dates.format(now, "e_d_u_yy-H_M_S"))_step_nsc_node_4_quad0_error.csv", result)
