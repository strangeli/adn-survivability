#using Pkg
#Pkg.activate(".")

using Distributed
using ClusterManagers

slurmc = true

if length(ARGS) > 0
	N_workers = parse(Int, ARGS[1]) - 1
else
	N_workers = 1
end

if slurmc
    using ClusterManagers
	addprocs(SlurmManager(N_workers))
else
	addprocs(N_workers)
end
println(nprocs(), " process(es)")

@everywhere dir = @__DIR__

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

# for plotting/saving. don't need them everywhere
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

@everywhere Pkg.status()

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

@everywhere function CIGRE_prob(;
    P_ref = 0.0,
    Q_ref = 0.0,
    nsc_node = 9,
    t_duration = 0.15,
    resistance = 10.0,
    quad = 1.0,
    verbose = false,
    get_cpg = false,
)
    ## simulation parameter
    t_fault = 0.1 # onset of node short circuit
    Znsc = complex(resistance, 0.0) * base_admittance #  fixed by Emocosy  Ω
    tspan = (0.0, 0.5) # total duration of a simulation
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
    cpg = ADN(pg, DGUnit, t -> P_ref, t -> Q_ref)

    ## find operation point
    icguess_pf = initial_guess(cpg, power_flow[:, :u])
    op = find_steady_state(cpg, icguess_pf)

    verbose ? check_operationpoint(cpg, op) : nothing

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
#     nsc_node = 9,
#     t_duration = 0.15,
#     resistance = 3.5,
#     quad = 0,
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

## setup EnsembleProblem

@everywhere cpg, prob = CIGRE_prob(;
    P_ref = 24.0,
    Q_ref = 5.0,
    nsc_node = 9,
    t_duration = 0.15,
    resistance = 3.5,
    quad = 1,
    get_cpg = true,
    verbose = true,
)

@everywhere begin
    # load flow sol
    S_total = complex(24.373141067954748, 6.115974820490637) * 1e6 / base_power

    Sdist = Uniform(0, 2abs(S_total))
    Wdist = Uniform(0, 2π ) #Uniform(angle(S_total)-π, angle(S_total)+π)
    Tdist = Normal(150e-3, 10e-3)
    Rdist = Uniform(3., 10.)

    NT =  40_000 #1000

    S = rand(Sdist, NT)
    W = rand(Wdist, NT)
    T = rand(Tdist, NT)
    R = rand(Rdist, NT)
end

@everywhere function prob_func(prob, i, repeat)
    return CIGRE_prob(;
        P_ref = S[i]*cos(W[i]),
        Q_ref = S[i]*sin(W[i]),
        nsc_node = 9,
        t_duration = T[i],
        resistance = R[i],
        quad = 1,
        get_cpg = false,
        verbose = false,
    )
end

# naive binomial proportion
@everywhere function bp(sample)
    # add one failure + one success
    N = length(sample) + 2
    m = (1 + sum(sample)) / N
    m±sqrt(m * (1 - m) / N)
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
    S = S,
    W = W,
    T = T,
    R = R,
    surv = first.(sim.u),
    surv_t = [u[2] for u in sim.u],
    fin_v = [u[3] for u in sim.u],
    min_v = [u[4] for u in sim.u],
    success = categorical(last.(sim.u)),
)

result[!, :P] = result.S.*cos.(result.W)
result[!, :Q] = result.S.*sin.(result.W)

now = unix2datetime(time())
CSV.write("$dir/results/mc_res_$(Dates.format(now, "e_d_u_yy-H_M_S")).csv", result)
