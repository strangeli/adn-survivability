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

## import special functions

# ADN control layer
include("$dir/control.jl")
# custom types
include("$dir/DGUnit.jl")
include("$dir/OLTC.jl")
include("$dir/ExternalGrid.jl")
# load actual data
include("$dir/cigre_static.jl")
# load pre-defined test cases
include("$dir/mwe.jl")
# failure model
include("$dir/short_circuit.jl")
# power flow solution
include("$dir/power_flow.jl")


## setup system

# define events
t_step = 1.
t_fault = 3. # onset of node short circuit
t_duration = 0.15 # duration of node short circuit

Ron = 2.5 # Ω
Znsc = complex(Ron, 0.0) * base_admittance #  fixed by Emocosy
nsc_node = 2 # sc on low-voltage side

S_total = complex(-0.49998009701576795, -0.49758893704171214)
# for S_pq = - complex(0.5, 0.5):
#S_total = complex(0.5000200975012079,0.5024116164039469)
P_ref(t) = t > t_step ? 0.9 * real(S_total) : real(S_total)
Q_ref(t) = t > t_step ? 0.9 * imag(S_total) : imag(S_total)
#P_ref(t) = real(S_total)
#Q_ref(t) = imag(S_total)

## load data from powerfactory for comparison

powerfactory =
    CSV.File("$dir/../results_power_factory/FRT_Vac_id_iq.csv", header = true, normalizenames = true)

## setup system

begin
    quad = 0
    DG = DGUnit(;
        K_pll = 1632.993, #Hz/pu
        T_int = 2.,
        K_PT1 = 1.0,
        T_PT1 = 1e-8,
        K_FRT = 2.0,
        I_max = 1.0, # pu
        P_limit = 1.0, #pu
        Q_limit = 1.0, #pu
        Vref = 1.0, # pu
        Vdead = 0.1, #pu
        S_pq = V -> complex(0.5, 0.5)*(quad * V^2 + 1 - quad),
        Y_n = 0.0,
    )
end

# Minimalbeispiel -> S=complex(0.5, 0.5) positive !!
PQ = PQAlgebraic(
    P = -0.5,
    Q = -0.5,
    Y_n = 0.,
)


## find operating point
pg = testbed(PQ; P_ref = t->P_ref(0), Q_ref = t->Q_ref(0))
ic_guess = initial_guess(pg, complex(1., 0))
ic_pf = pf_sol(PowerGrid(pg), ic_guess, nothing)

ic_pf[1,:s] # -> S_total

cpg = testbed(DG; P_ref = t->P_ref(0), Q_ref = t->Q_ref(0))

# set all initial voltages to slack value
op = find_steady_state(cpg, initial_guess(cpg, ic_pf[:, :u]))
check_operationpoint(cpg, op)

# ode = rhs(cpg)
# du = similar(op.vec)
# p = cpg.controller(op.vec, nothing, 0.)
# ode(du, op.vec, nothing, 0.)
# du

@show FRTcurrent(op[2, :v], DG.K_FRT, DG.Vref, DG.Vdead)
@show VoltageDependence(DG.S_pq, op[2, :v], 1., 0., 0.)

@show op[2, :P_int]
@show op[2, :Q_int]

@show op[2, :P_g]
@show op[2, :Q_g]

@show op[2, :v] # 1.00241 pu
@show op[2, :φ] .|> rad2deg # 0.14 deg
@show op[1, :s] # -499.98 kW, -497.61 kvar
@show op[2, :s] # 500.0 kW, 500.0 kvar
# verketteter Strom / sqrt(3)
@show op[1, :iabs] * base_current_HV * 1e-3 / sqrt(3) # 0.004 kA
@show op[2, :iabs] * base_current * 1e-3 / sqrt(3) # 0.020 kA

## check dynamics

cpg = testbed(DG; P_ref = P_ref, Q_ref = Q_ref)

perturb = Perturbation(2, :θ, Dec(10.))
u0 = perturb(op)

prob = ODEProblem(rhs(cpg), u0.vec, (0.0, 2.))

_sol = solve(prob, Rodas4(), tstops=t_step)
sol = PowerGridSolution(_sol, cpg)

## node short circuit

cpg = testbed(DG; P_ref = P_ref, Q_ref = Q_ref)
nsc = NodeShortCircuitAlt(;
    node_number = nsc_node,
    Y = 1/Znsc,
    tspan_fault = (t_fault, t_fault + t_duration),
)
@time sol = simulate(nsc, cpg, op.vec, (0., 4.); additional_tstops=t_step)


## plots

# TODO: bugs in solution object, e.g. (sol, 2, :v), (sol, 2, :φ)
# (sol, 2, :iabs) triggers call to rhs and hangs
# --> compare with indexing sol

plot(sol, :, :v, label = false)#, yscale=:log10)
plot!(2.5 .+ powerfactory.b_tnow_in_s, powerfactory.s_Vac_ist, label = "PF", ms = 2, ls=:dash)
hline!(op[:, :v], c=:gray, label = "op")
#vline!([t_fault, t_fault + t_duration], label = "fault", c = :black, ls=:dot)

Δ = [
    cpg.controller(u, nothing, t) |> first
    for (t, u) in zip(sol.dqsol.t, sol.dqsol.u)
]

plot(sol.dqsol.t, first.(Δ))
plot!(sol.dqsol.t, last.(Δ))

plot(sol.dqsol.t, [FRTcurrent(sol(t, 2, :v), DG.K_FRT, DG.Vref, DG.Vdead) for t in sol.dqsol.t], ylabel="iqplus")

plot(sol, :, :φ, label = false)
plot!(sol, 2, :θ, c = :black, label = "PLL")
hline!(op[:, :φ], c = :gray, label = "op")
#vline!([t_fault, t_fault + t_duration], label = "fault", c = :black)
#ylims!(-0.01, 0.01)

plot(sol, 2, :P_g)
plot!(sol, 2, :Q_g)
hline!(op[2:2, :P_g], label = "op")
hline!(op[2:2, :Q_g], label = "op")
vline!([t_fault, t_fault + t_duration], label = "fault", c = :black, ls=:dot)

plot(sol, 2, :P_int)
plot!(sol, 2, :Q_int)
hline!(op[2:2, :P_int], label = "op")
hline!(op[2:2, :Q_int], label = "op")
vline!([t_fault, t_fault + t_duration], label = "fault", c = :black, ls=:dot)




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
