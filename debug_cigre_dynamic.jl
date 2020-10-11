## TODO Heavily work in progress

using Pkg
Pkg.activate(".")

using PowerDynamics
using OrdinaryDiffEq
using CSV
using Plots
plotly()

##

dir = @__DIR__

## set per unit MV
const base_power = 1E6 # 1MW
const base_voltage = 20E3 # 20kV
const base_current = base_power / base_voltage # 50A
const base_admittance = base_power / base_voltage^2 # 0.0025Ω^-1
const ω = 2 * π * 50.0 # 314.1593rad/s
# per unit HV
const base_voltage_HV = 110E3 # 110kV
const base_admittance_HV = base_power / base_voltage_HV^2 # 8.264462809917356e-5
const base_current_HV = base_power / base_voltage_HV

## some additional functions that might develop into new PD.jl features
#include("PDpatches.jl")
# ADN control layer (needs ControlledPowerGrid from PDpatches)
include("$dir/control.jl")
# custom types
include("$dir/DGUnit.jl")
include("$dir/OLTC.jl")
include("$dir/RLLine.jl")
# load actual data
include("$dir/cigre_static.jl")
# failure model
include("$dir/short_circuit.jl")
# static solution
include("$dir/power_flow.jl")
# load pre-defined test cases
include("$dir/mwe.jl")
#
# import PowerDynamics: systemsize
# @views systemsize(pg::PowerGrid) = sum(map(n -> dimension(n), pg.nodes)) + sum(map(n -> dimension(n), pg.lines))
# @views systemsize(pg::ControlledPowerGrid) = systemsize(PowerGrid(pg))



## setup system

# define events
t_step = 1.0 # step in the power flow reference
t_fault = 3.0 # onset of node short circuit
t_duration = 0.15 # duration of node short circuit

S_total = complex(24.40048903027117, 9.39761634666317) * 1e6 / base_power
P_ref(t) = t > t_step ? 0.9 * real(S_total) : real(S_total)
Q_ref(t) = t > t_step ? 0.9 * imag(S_total) : imag(S_total)
#P_ref(t) = - real(S_total)
#Q_ref(t) = - imag(S_total)

# edge data
# ode = rhs(pg)
# ode.f.graph_data.e

## find operating point
cpg_static, power_flow = testbed_CIGRE(DGUnit; P_ref = t -> P_ref(0.0), Q_ref = t -> Q_ref(0.0))
cpg = testbed_CIGRE_dynamic(DGUnit; P_ref = t -> P_ref(0.0), Q_ref = t -> Q_ref(0.0))

L = systemsize(cpg_static) + 4length(cpg.lines)

## init from load flow sol

icguess_pf = initial_guess(cpg_static, power_flow.V .* exp.(power_flow.A .* 1im))
op_static = find_steady_state(cpg_static, icguess_pf)
icguess_full = ones(L)
icguess_full[1:length(icguess_pf)] .= op_static.vec

for (l, line) in enumerate(cpg.lines)
    if line isa OLTC
        D=4
    else
        D = dimension(line)
    end
    offset = length(icguess_pf) + (l-1) * D
    Y = PiModel(line)
    volt = [
        complex(power_flow.real[line.from], power_flow.imag[line.from]),
        complex(power_flow.real[line.to], power_flow.imag[line.to]),
    ]
    current = Y * volt
    c_xy = [real(first(current)); imag(first(current)); -real(last(current)); -imag(last(current))]
    icguess_full[offset+1:offset+D] .= c_xy
end

ode = rhs(cpg)

dx = similar(icguess_full)
cpg.controller(icguess_full, nothing, 0.)
ode(dx, icguess_full, nothing, 0.)
dx

op_prob = ODEProblem(ode, icguess_full, (0.0, 1.0))
op = solve(
    SteadyStateProblem(op_prob),
    DynamicSS(Rodas4(); abstol = 1e-8, reltol = 1e-6, tspan = Inf),
    initialize_alg=BrownFullBasicInit(),
)

prob = ODEProblem(ode, op.vec, (0.0, 2.0))

_sol = solve(prob, Rodas4(autodiff=true), initialize_alg=BrownFullBasicInit(), d_discontinuities = [t_step])






## check dynamics
ode = rhs(cpg)
perturb = Perturbation(nsc_node, :θ, Dec(0.))
u0 = perturb(op)

prob = ODEProblem(ode, op.vec, (0.0, 2.0))

_sol = solve(prob, Rodas4(autodiff=true), initialize_alg=BrownFullBasicInit(), d_discontinuities = [t_step])
sol = PowerGridSolution(_sol, cpg)

plot(_sol)
plot(sol, 2:12, :v, label = false)
# hline(op_prestep[2:12, :v], c = :gray, label = "op")
# hline!(op_poststep[2:12, :v], c = :gray, label = "op")



## node short circuit

cpg, power_flow = testbed_CIGRE(DGUnit; P_ref = P_ref, Q_ref = Q_ref)

Ron = 4. # 2.5 Ω
Znsc = complex(Ron, 0.0) * base_admittance #  fixed by Emocosy
nsc_node = 9 #(8) in CIGRE - fixed by Emocosy

nsc = NodeShortCircuitAlt(;
    node_number = nsc_node,
    Y = (1 / Znsc),
    tspan_fault = (t_fault, t_fault + t_duration),
)

@time sol = simulate(nsc, cpg, op_ss_pf.vec, (0.0, 4.0))


plot(sol, 2:12, :v, label = false, xlabel="t [s]", ylabel="v [pu]")
plot(sol, 2:12, :θ, label = false, xlabel="t [s]", ylabel="PLL [rad]")
plot(sol, 2:12, :P_g, label = false, xlabel="t [s]", ylabel="Pg [pu]")
plot(sol, 2:12, :iqp, label = false, xlabel="t [s]", ylabel="Iq+ [pu]")

Δ = [
    cpg.controller(u, nothing, t) |> first
    for (t, u) in zip(sol.dqsol.t, sol.dqsol.u)
]

plot(sol.dqsol.t, first.(Δ), label="ΔP")
plot!(sol.dqsol.t, last.(Δ), label="ΔQ")


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
