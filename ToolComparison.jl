begin
    using Pkg
    Pkg.activate(".")

    using PowerDynamics
    using OrdinaryDiffEq
    using CSV
    using DataFrames
    using Plots
    gr()
end

##

dir = @__DIR__

## set per unit MV
begin
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
#include("PDpatches.jl")
begin
    # ADN control layer (needs ControlledPowerGrid from PDpatches)
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
    # load pre-defined test cases
    include("$dir/mwe.jl")
end


## setup system

# define events
t_step = 1.0 # step in the power flow reference
t_fault = 3.0 # onset of node short circuit
t_duration = 0.15 # duration of node short circuit

# # measured on HV side of trafo
# S_total = complex(24.40048903027117, 9.39761634666317) * 1e6 / base_power
# P_ref(t) = t > t_step ? 24.0257 : real(S_total)
# Q_ref(t) = t > t_step ? 8.0843 : imag(S_total)

# measured on LV side of trafo (reference value for controller taken to be positive)
S_total = complex(24.373141067954748, 6.115974820490637) * 1e6 / base_power
P_ref(t) = t > t_step ? 24. : real(S_total)
Q_ref(t) = t > t_step ? 5. : imag(S_total)

## load data from powerfactory for comparison

powerfactory_step = CSV.File("$dir/../results_power_factory/Step_change.csv", header = true, normalizenames = true) |> DataFrame
powerfactory10ohm = CSV.File("$dir/../results_power_factory/KS_10_Ohm.csv", header = true, normalizenames = true) |> DataFrame
#powerfactory2p5ohm = CSV.File("$dir/../results_power_factory/KS_2_5_Ohm.csv", header = true, normalizenames = true) |> DataFrame
powerfactory3p5ohm = CSV.File("$dir/../results_power_factory/KS_3_5_Ohm.csv", header = true, normalizenames = true) |> DataFrame

## access line power flow
using NetworkDynamics: GetGD

function get_trafo_flow(rhs, x)
    gd = rhs.f.open_loop_dyn(x.vec, nothing, 0., GetGD)
    #gs = ode.f.open_loop_dyn(GetGS)
    trafo_currents = first(gd.e)
    HV_to_LV = - complex(trafo_currents[1], trafo_currents[2]) * base_current_HV
    LV_to_HV = complex(trafo_currents[3], trafo_currents[4]) * base_current
    uHV = complex(gd.v[1][1], gd.v[1][2]) * base_voltage_HV
    uLV = complex(gd.v[2][1], gd.v[2][2]) * base_voltage
    power_to_HV = uHV * conj(HV_to_LV)
    power_from_HV = uLV * conj(LV_to_HV)
    absHV = abs(HV_to_LV)
    absLV = abs(LV_to_HV)
    println("$(abs(uHV)* 1e-3) kV \t\t\t\t\t | \t $(real(power_to_HV)* 1e-6) MW \t\t ()() \t $(real(power_from_HV)* 1e-6) MW \t | \t $(abs(uLV)* 1e-3) kV")
    println("$(abs(uHV)/base_voltage_HV) pu \t\t\t\t\t\t | \t $(imag(power_to_HV)* 1e-6) Mvar \t ()() \t $(imag(power_from_HV)* 1e-6) Mvar \t | \t $(abs(uLV)/base_voltage) pu")
    println("$(rad2deg(angle(uHV))) ° \t | \t $(absHV*1e-3/sqrt(3)) kA \t ()() \t $(absLV*1e-3/sqrt(3)) kA \t\t | \t $(rad2deg(angle(uLV))) °")
    return nothing #power_to_HV, power_from_HV, absHV, absLV
end


## check dynamics

cpg, power_flow = testbed_CIGRE(DGUnit; P_ref = t -> P_ref(0.0), Q_ref = t -> Q_ref(0.0), quad=1)
icguess_pf = initial_guess(cpg, power_flow.V .* exp.(power_flow.A .* 1im))
op = find_steady_state(cpg, icguess_pf)
ode = rhs(cpg)

cpg_step, _ = testbed_CIGRE(DGUnit; P_ref = t -> P_ref(t_step+1), Q_ref = t -> Q_ref(t_step+1), quad=1)
opstep = find_steady_state(cpg_step, icguess_pf)
ode_step = rhs(cpg_step)

# print power factory style
get_trafo_flow(ode, op)
get_trafo_flow(ode_step, opstep)

# map(x-> abs(complex(x[1], x[2]))*base_current*1e-3/sqrt(3), gd.e)

## simulate step

cpg, power_flow = testbed_CIGRE(DGUnit; P_ref = P_ref, Q_ref = Q_ref, quad=1)
ode = rhs(cpg)
perturb = Perturbation(9, :θ, Dec(0.))
u0 = perturb(op)

prob = ODEProblem(ode, u0.vec, (0.0, t_fault))

_sol = solve(prob, Rosenbrock23(), d_discontinuities = [t_step])
sol = PowerGridSolution(_sol, cpg)


plot(powerfactory_step.b_tnow_in_s, Array(powerfactory_step[:, 2:12]), label=false , ls=:dot, c=:black, lw=1.5)
    hline!(opstep[2:12, :v], c=:gray, label=false)
    plot!(sol, 2:12, :v, c=:red, label = false, xlabel="t [s]", ylabel="v [pu]")
    xlims!(prob.tspan...)
    ylims!(0.95, 1.05)
    #savefig("$dir/plots/step_change_voltage.pdf")


scatter(Array(powerfactory_step[1, 2:12]), op[2:12,:v], label="before", xlabel="v PF [pu]", ylabel="v PD [pu]", legend=:topleft)
    scatter!(Array(powerfactory_step[end, 2:12]), opstep[2:12,:v], label="after")
    plot!(0.95:0.001:1.03, 0.95:0.001:1.03, c=:gray, label=false)
    #savefig("$dir/plots/step_change_op.pdf")

Δ = [
    first(cpg.controller(u, nothing, t))[1:2]
    for (t, u) in zip(sol.dqsol.t, sol.dqsol.u)
]

plot(sol.dqsol.t, first.(Δ), label="ΔP", xlabel="t [s]", ylabel="error [pu]")
plot!(sol.dqsol.t, last.(Δ), label="ΔQ")

##  node short circuit

cpg, power_flow = testbed_CIGRE(DGUnit; P_ref = P_ref, Q_ref = Q_ref, quad=1)
icguess_pf = initial_guess(cpg, power_flow.V .* exp.(power_flow.A .* 1im))
op = find_steady_state(cpg, icguess_pf)

Ron = 10 # 2.5 Ω
Znsc = complex(Ron, 0.0) * base_admittance #  fixed by Emocosy
nsc_node = 9 #(8) in CIGRE - fixed by Emocosy

nsc = NodeShortCircuitAlt(;
    node_number = nsc_node,
    Y = (1 / Znsc),
    tspan_fault = (t_fault, t_fault + t_duration),
)

@time sol10 = simulate(nsc, cpg, op.vec, (0.0, 4.0); additional_tstops=t_step)

tarray = range(sol10.dqsol.t[1], sol10.dqsol.t[end], length=4001)
df10 = sol10(tarray, 2:12, :v) |> transpose |> DataFrame
names!(df10, [Symbol("MV_$a") for a in 1:11])
df10.t = tarray
CSV.write("$dir/KS_10ohm_julia.csv", df10)

Ron = 3.5 # 2.5 Ω
Znsc = complex(Ron, 0.0) * base_admittance #  fixed by Emocosy
nsc_node = 9 #(8) in CIGRE - fixed by Emocosy

nsc = NodeShortCircuitAlt(;
    node_number = nsc_node,
    Y = (1 / Znsc),
    tspan_fault = (t_fault, t_fault + t_duration),
)

@time sol3p5 = simulate(nsc, cpg, op.vec, (0.0, 4.0); additional_tstops=t_step)

tarray = range(sol3p5.dqsol.t[1], sol10.dqsol.t[end], length=4001)
df3p5 = sol3p5(tarray, 2:12, :v) |> transpose |> DataFrame
names!(df3p5, [Symbol("MV_$a") for a in 1:11])
df3p5.t = tarray
CSV.write("$dir/KS_3p5ohm_julia.csv", df3p5)


begin
    p1 = hline(opstep[2:12, :v], c=:gray)
    plot!(p1, sol10, 2:12, :v, c=:red, legend = false, xlabel="t [s]", ylabel="v [pu]")
    plot!(p1, powerfactory10ohm.b_tnow_in_s, Array(powerfactory10ohm[:, 2:12]), ls=:dot, c=:black)
    ylims!(0.94, 1.02)
    xlims!(0.75, 1.7)

    p2 = plot(sol10, 2:12, :v, c=:red, legend = false, xlabel="t [s]", ylabel="v [pu]")
    plot!(p2, powerfactory10ohm.b_tnow_in_s, Array(powerfactory10ohm[:, 2:12]), ls=:dot, c=:black)
    ylims!(0, 1.05)
    xlims!(2.95, 3.2)

    p3 = hline(opstep[2:12, :v], c=:gray)
    plot!(p3, sol3p5, 2:12, :v, c=:red, legend = false, xlabel="t [s]", ylabel="v [pu]")
    plot!(p3, powerfactory3p5ohm.b_tnow_in_s, Array(powerfactory3p5ohm[:, 2:12]), ls=:dot, c=:black)
    ylims!(0.94, 1.02)
    xlims!(0.75, 1.7)

    p4 = plot(sol3p5, 2:12, :v, c=:red, legend = false, xlabel="t [s]", ylabel="v [pu]")
    plot!(p4, powerfactory3p5ohm.b_tnow_in_s, Array(powerfactory3p5ohm[:, 2:12]), ls=:dot, c=:black)
    ylims!(0, 1.05)
    xlims!(2.95, 3.2)

    plot(p1, p2, p3, p4, layout = @layout [a b; c d])

    #hline!(op_fault[2:12, :v], c=:gray, label=false)
    #savefig("$dir/plots/KS.pdf")
end

Δ = [
    first(cpg.controller(u, nothing, t))[1:2]
    for (t, u) in zip(sol.dqsol.t, sol.dqsol.u)
]

plot(sol.dqsol.t, first.(Δ), label="ΔP", xlabel="t [s]", ylabel="error [pu]")
    plot!(sol.dqsol.t, last.(Δ), label="ΔQ")

plot(sol, 2:12, :θ, label = false, xlabel="t [s]", ylabel="PLL [rad]")

plot(sol, 2:12, :P_int, label = false, xlabel="t [s]", ylabel="Pint [pu]")
plot(sol, 2:12, :P_g, label = false, xlabel="t [s]", ylabel="Pg [pu]")

plot(sol, 2:12, :iqp, label = false, xlabel="t [s]", ylabel="Iq+ [pu]")
















## old experiments. here be dragons




##
fcpg = nsc(cpg)

fault_op = find_operationpoint(fcpg, icguess_pf)


# TODO idea: solve PF for voltages
# TODO better: solve constraints first and then other variables, compar DE.jl

nsc_op = State(fcpg, find_valid_initial_condition(fcpg, op.vec))

function stepwise_init(cpg, x, p = nothing, t0 = 0.0)
    rpg = rhs(cpg)
    M = rpg.mass_matrix
    algebraic_vars = [all(iszero, c) for c in eachcol(M)]
    diff_vars = .!algebraic_vars
    algebraic_eqs = [all(iszero, r) for r in eachrow(M)]
    diff_eqs = .!algebraic_eqs

    out = copy(x)

    dview = @view out[diff_vars]
    # fix differential vars y and solve g(x, y)=0
    function solve_with_fixed_diff(du, u)
        u[diff_vars] .= dview
        rpg(du, u, p, t0)
        du .= algebraic_eqs .* du
    end

    res = nlsolve(
        solve_with_fixed_diff,
        x,
        autodiff = :forward,
        show_trace = false,
        method = :trust_region,
    )
    out[algebraic_vars] .= res.zero[algebraic_vars]
    println(res.x_converged, res.f_converged, res.residual_norm)
    return out
end



rand = State(fcpg, 10 * randn(systemsize(fcpg)))
ic = State(fcpg, stepwise_init(fcpg, op.vec .+ rand.vec))
ic2 = State(fcpg, relax_alg(fcpg, op.vec))

@show op[2:12, :v]
@show fault_op[2:12, :v]
@show nsc_op[2:12, :v]

@show op[2:12, :P_int]
@show nsc_op[2:12, :P_int]

timespan = (0.0, 4.0)
problem = ODEProblem{true}(rhs(cpg), op_ss_pf.vec, timespan)
integrator = init(
    problem,
    Rodas5(autodiff = false),
    #force_dtmin = true,
    #save_everystep = false,
    tstops = [t_fault, t_fault + t_duration],
)
# step to failure
step!(integrator, t_fault, true)
check_error(integrator)
pre_failure = State(cpg, copy(integrator.u))
@show pre_failure[2:12, :v]
@show pre_failure[2:12, :iqp]
integrator.f = rhs(nsc(cpg))
DiffEqBase.initialize_dae!(integrator, BrownFullBasicInit())
valid_post_failure = State(cpg, copy(integrator.u))
@show valid_post_failure[2:12, :v]
@show valid_post_failure[2:12, :iqp]
small_step = 0.00001
step!(integrator, small_step, true)
check_error(integrator)
end_of_failure = State(cpg, copy(integrator.u))
@show end_of_failure[2:12, :v]
@show end_of_failure[2:12, :iqp]




timespan = (t_fault, t_fault+t_duration)
x2 = find_valid_initial_condition(nsc(cpg), op.vec)

dx = similar(x2)
rpg(dx, x2, nothing, t_fault)
@show dx
findall(abs.(dx) .> 0.1)
findall(isnan.(dx))

problem = ODEProblem{true}(rhs(nsc(cpg)), x2, timespan)
sol = solve(
    problem,
    ABDF2(),
    #Rodas4(autodiff = true),
    #force_dtmin = true,
)

using NLsolve
using LinearAlgebra
using LineSearches



## plots

plot(sol, 2:12, :v, label = false)
hline!(fault_op[2:12, :v], c = :gray, label = "op")
vline!(
    [t_step, t_fault, t_fault + t_duration],
    label = "fault",
    c = :black,
    ls = :dot,
)

plot(sol, 2:12, :φ, label = false)
plot!(sol, 2:12, :θ, ls = :dot, c = :black, label = "PLL")
hline!(op[:, :φ], c = :gray, label = "op")
#vline!([t_fault, t_fault + t_duration], label = "fault", c = :black)
#ylims!(-0.01, 0.01)

plot(sol, 2:12, :P_g)
plot!(sol, 2:12, :Q_g)
hline!(op[2:12, :P_g], label = "op")
hline!(op[2:12, :Q_g], label = "op")
vline!([t_fault, t_fault + t_duration], label = "fault", c = :black, ls = :dot)

plot(sol, 2:12, :P_int)
plot!(sol, 2:12, :Q_int)
hline!(op[2:12, :P_int], label = "op")
hline!(op[2:2, :Q_int], label = "op")
vline!([t_fault, t_fault + t_duration], label = "fault", c = :black, ls = :dot)

Δ = [
    cpg.controller(u, nothing, t) |> first
    for (t, u) in zip(sol.dqsol.t, sol.dqsol.u)
]

plot(sol.dqsol.t, first.(Δ))
plot!(sol.dqsol.t, last.(Δ))


plot(sol, 2:12, :iqp)


inc = (last(sol.dqsol.prob.tspan) - first(sol.dqsol.prob.tspan)) / 100.0
t_array = range(
    first(sol.dqsol.prob.tspan);
    stop = last(sol.dqsol.prob.tspan),
    step = inc,
)

# currents
node = nsc_node
# current set point
I_r = [
    complex(sol(t, node, :P_g), -sol(t, node, :Q_g)) / sol(t, node, :v)
    for t in t_array
]
plot(t_array, abs.(I_r), label = "I_r_$node")
