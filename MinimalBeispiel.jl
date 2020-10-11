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


## define OLTC

include("$dir/OLTC.jl")

## define DGUnit

import Base: @__doc__
using PowerDynamics: @DynamicNode
import PowerDynamics:
    construct_vertex, dimension, symbolsof, showdefinition, AbstractNode
using LinearAlgebra: Diagonal
using NetworkDynamics: ODEVertex

function Limit(x, l, u)
    if x < l
        return l
    elseif x > u
        return u
    else
        return x
    end
end



begin
    @__doc__ struct StaticGeneratorObs <: AbstractNode
        K_pll
        S_ref
        K_FRT
        I_max
        Vref
        Vdead
        Y_n
    end
    StaticGeneratorObs(; K_pll, S_ref, K_FRT, I_max, Vref, Vdead, Y_n) =
        StaticGeneratorObs(K_pll, S_ref, K_FRT, I_max, Vref, Vdead, Y_n)
    function construct_vertex(par::StaticGeneratorObs)
        K_pll = par.K_pll
        S_ref = par.S_ref
        K_FRT = par.K_FRT
        I_max = par.I_max
        Vref = par.Vref
        Vdead = par.Vdead
        Y_n = par.Y_n
        function rhs!(dx, x, e_s, e_d, p, t)
            u = complex(x[1], x[2])
            # current from incoming/outgoing edges, nodal shunt, PQ load background
            i = total_current(e_s, e_d) + Y_n * u
            Vamp = abs(u)

            # this is the only dynamic variable
            θ = x[3]

            # to keep track of the internal currents for debugging ...
            i_d_r = x[4] # behind FRT
            i_q_r = x[5] # behind FRT
            i_d_r_b = x[6] # before FRT
            i_q_r_b = x[7] # before FRT
            iqp = x[8]

            # PLL
            dθ = -K_pll * Vamp * sin(θ - angle(u))

            # FRT
            # additional reactive current during fault
            if Vamp > Vref + Vdead
                iqplus = K_FRT * (Vamp - (Vref + Vdead)) # pos for over_voltage
            elseif Vamp < Vref - Vdead
                iqplus = K_FRT * (Vamp - (Vref - Vdead)) # neg for under_voltage
            else
                iqplus = 0.0
            end

            # current set point
            I_r = conj(S_ref) / Vamp

            # save values for plotting
            di_d_r_b = i_d_r_b - real(I_r)
            di_q_r_b = i_q_r_b - imag(I_r)
            diqp = iqp - iqplus

            # FRT limiter
            if Vamp > Vref + Vdead || Vamp < Vref - Vdead
                iq_fault = Limit(imag(I_r) + iqplus, -I_max, I_max)
                id_lim = sqrt(I_max^2 - iq_fault^2)
                id_fault = Limit(real(I_r), 0.0, id_lim) #-id_lim, id_lim)
                I_r_limited = complex(id_fault, iq_fault)
            else
                id_normal = Limit(real(I_r), 0.0, I_max) #-I_max, I_max)
                iq_lim = sqrt(I_max^2 - id_normal^2)
                iq_normal = Limit(imag(I_r), -iq_lim, iq_lim)
                I_r_limited = complex(id_normal, iq_normal)
            end

            # save values for plotting
            di_d_r = i_d_r - real(I_r_limited)
            di_q_r = i_q_r - imag(I_r_limited)

            du = i - exp(θ * 1im) * I_r_limited # = 0 constraint

            try
                dx[1] = real(du)
                dx[2] = imag(du)
                dx[3] = dθ
                dx[4] = di_d_r
                dx[5] = di_q_r
                dx[6] = di_d_r_b
                dx[7] = di_q_r_b
                dx[8] = diqp
                return nothing
            catch e
                if typeof(e) === UndefVarError
                    throw(NodeDynamicsError("you need to provide $(e.var)"))
                else
                    throw(e)
                end
            end
        end
        ODEVertex(
            f! = rhs!,
            dim = 8,
            mass_matrix = Diagonal([0, 0, 1, 0, 0, 0, 0, 0]),
            sym = Symbol[
                :u_r,
                :u_i,
                :θ,
                :i_d_r,
                :i_q_r,
                :i_d_r_b,
                :i_q_r_b,
                :iqp,
            ],
        )
    end
    symbolsof(::StaticGeneratorObs) = begin
        [:u_r, :u_i, :θ, :i_d_r, :i_q_r, :i_d_r_b, :i_q_r_b, :iqp]
    end
    dimension(::StaticGeneratorObs) = begin
        8
    end
end


## setup system

# define events
t_fault = 0.5 # onset of node short circuit
t_duration = 0.15 # duration of node short circuit

Ron = 2.5 # Ω
Znsc = complex(Ron, 0.0) * base_admittance #  fixed by Emocosy
nsc_node = 2 # sc on low-voltage side

## setup system
SL = SlackAlgebraic(; U = complex(1.0, 0.0))

DG = StaticGeneratorObs(;
    K_pll = 1632.993, #Hz/pu
    S_ref = complex(0.5, 0.5), # pu
    K_FRT = 2.0,
    I_max = 1.0, # pu
    Vref = 1.0, # pu
    Vdead = 0.1, #pu
    Y_n = 0., #t ->(t_fault < t < t_fault + t_duration ? 1 / Znsc : 0.0),
)

T = OLTC(
    from = 1,
    to = 2,
    Uos = 110E3 / base_voltage_HV, # Bemessungsspannung Oberspannungsseite in kV
    Uus = 20E3 / base_voltage,# Bemessungsspannung Unterspannungsseite in kV
    k = 0, # Kennzahl der Schaltgruppe, Bsp: YD5-> 5
    ssp = 0.0, # Stufenschalterposition
    stufe_us = 0.625, # Stufung pro Stufenschalterstellung in %
    Sr = 25E6 / base_power, #Bemessungsscheinleistung in MVA
    uk = 12.0, # Kurzschlussspannung in %
    Pvk = 25E3 / base_power, # Kupferverluste in MW
    Pvl = 0.0, # Eisenverluste in kW
    iLeer = 0.0, # Leerlaufstrom in %
)

pg = PowerGrid([SL, DG], [T])


## load data from powerfactory for comparison

powerfactory =
    CSV.File("$dir/../results_power_factory/FRT_Vac_id_iq.csv", header = true, normalizenames = true)

## find operating point

# set all initial voltages to slack value
op = find_operationpoint(pg)

## compare power flow
@show op[2, :v] # 1.00241 pu
@show op[2, :φ] .|> rad2deg # 0.14 deg

@show op[1, :s] * 1e3 # -499.98 kW, -497.61 kvar
@show op[2, :s] * 1e3 # 500.0 kW, 500.0 kvar

# verketteter Strom / sqrt(3)
@show op[1, :iabs] * base_current_HV * 1e-3 / sqrt(3) # 0.004 kA
@show op[2, :iabs] * base_current * 1e-3 / sqrt(3) # 0.020 kA

@show op[2, :i_d_r] # 0.4988 p.u.
@show op[2, :i_q_r] # -0.4988 p.u.

# since FRT is not active, the currents should be identical
@assert op[2, :i_d_r_b] == op[2, :i_d_r]
@assert op[2, :i_q_r_b] == op[2, :i_q_r]

# iq_plus is zero
@assert op[2, :iqp] == 0

## simulate short circuit

# ode = rhs(pg)
# problem = ODEProblem(ode, op.vec, (0.0, 1.0))
# _sol = solve(
#     problem,
#     Rodas4(),
#     force_dtmin = true,
#     tstops = [t_fault, t_fault + t_duration],
# )
# sol = PowerGridSolution(_sol, pg)

# experimental version:
# include("$dir/short_circuit.jl")
#
nsc = NodeShortCircuit(;
    node_number = nsc_node,
    Y = 1/Znsc,
    tspan_fault = (t_fault, t_fault + t_duration),
)

sol = simulate(nsc, op, (0., 1.))



## plot result

scatter(powerfactory.b_tnow_in_s, powerfactory.s_Vac_ist, label = "PF", ms = 2)
plot!(sol, 2:2, :v, label = "PD")
xlims!(0.49, 0.66)
xlabel!("t [s]")
ylabel!("V_ac [pu]")

scatter(powerfactory.b_tnow_in_s, powerfactory.s_id_ref, label = "PF", ms = 2)
plot!(sol, 2, :i_d_r, label = "PD")
xlims!(0.49, 0.66)
xlabel!("t [s]")
ylabel!("i_d_ref [pu]")

scatter(powerfactory.b_tnow_in_s, powerfactory.s_iq_ref, label = "PF", ms = 2)
plot!(sol, 2, :i_q_r, label = "PD")
xlims!(0.49, 0.66)
xlabel!("t [s]")
ylabel!("i_q_ref [pu]")

# check FRT function
scatter(
    powerfactory.b_tnow_in_s,
    powerfactory.s_iq_ref,
    label = "i_q_ref PF",
    ms = 2,
)
plot!(sol, 2, :i_q_r, label = "FRT+limit")
plot!(sol, 2, :i_q_r_b, label = "before FRT")
plot!(sol, 2, :iqp, label = "iqplus")
xlims!(0.49, 0.66)
xlabel!("t [s]")
ylabel!("q currents [pu]")

scatter(powerfactory.b_tnow_in_s, powerfactory.s_id_ref, label = "PF", ms = 2)
plot!(sol, 2, :i_d_r, label = "FRT+limit")
plot!(sol, 2, :i_d_r_b, label = "before FRT")
xlims!(0.49, 0.66)
xlabel!("t [s]")
ylabel!("d currents [pu]")

## alternative model to compare short circuit power provided by slack bus

# PQ = PQAlgebraic(; S=complex(0., 0.), Y_n = 0.)
# pg = PowerGrid([SL, PQ], [T])

## to use the NodeShortCircuit, T must be converted to a PiModelLine

# # For ssp=0, replace it with PiModelLine
# # assume t_km = 1, t_mk real
# Y = PiModel(T) #./ admittance_m
# ue = sqrt(-Y[2, 2] / Y[1, 1])
# y = Y[1, 2] / ue
# ys = -Y[1, 1] - y
#
# L = PiModelLine(;
#     from = T.from,
#     to = T.to,
#     y = y,
#     y_shunt_km = ys,
#     y_shunt_mk = ys,
# )
# @assert iszero(T.ssp) && all(PiModel(L) .== PiModel(T))
