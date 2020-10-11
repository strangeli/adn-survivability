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
    @__doc__  struct StaticGeneratorObs <: AbstractNode
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


function VoltageDependence(P, V, V0, a, b)
    c = 1. - a - b
    return P * ( a * (V/V0)^2 + b * (V/V0) + c)
end

function FRTcurrent(V, K_FRT, Vref, Vdead)
    if V > Vref + Vdead
        iqplus = K_FRT * (V - (Vref + Vdead)) # pos for over_voltage
    elseif V < Vref - Vdead
        iqplus = - K_FRT * ((Vref - Vdead) - V) # neg for under_voltage
    else
        iqplus = 0.0
    end
    return iqplus
end


@doc """
```Julia
DGUnit(;I_r)
```
Pr_l_max = 2.0 # Holms Model: 10e12 / base_power
Rr_l_min = -2.0 # Holms Model: -10e12 / base_power
T = 10.0 # unit: s
K = 1.0 # unit: [y]/[u] = [P]/[P] = 1
Pmax = 2. # in this Limiter it is set to 200 MW in Holms Model!!
Pmin = -2. # in this limiter it is set to 0MW in Holms Model!
Qmax = 2. # in this Limiter it is set to 60 MW
Qmin = -2. # in this limiter it is set to -60 MW
Vac_ref = 1 * 20 * sqrt(2 / 3) # phase-ground = phase-phase * sqrt(2/3) in kV
V_dead = 0.1 * Vac_ref
k_FRT = 2.0 # p.u. or A/V or kA/kV stays the same as long as imax and Vac, Vac_ref and V_dead are all p.u or A,V or kA/kV
imax = 50 / sqrt(3) / base_current # assumption: 1 MVA nominal power of DGunit at 20kV and SB = 25 MVA
K_pll = 0.1 * base_voltage # theta in rad, Kpll before in rad/(Vs) now * 1000 V/kV --> rad/(kVs)
Y_shunt = 0.4 # shunt admittance for e.g. node short circuit
"""
begin
    @__doc__ struct DGUnit <: AbstractNode
        K_pll
        T_int
        K_PT1
        T_PT1
        K_FRT
        I_max
        P_limit
        Q_limit
        Vref
        Vdead
        S_pq
        Y_n
    end
    DGUnit(; K_pll, T_int, K_PT1, T_PT1, K_FRT, I_max, P_limit, Q_limit, Vref, Vdead, S_pq, Y_n) =
        DGUnit(K_pll, T_int, K_PT1, T_PT1, K_FRT, I_max, P_limit, Q_limit, Vref, Vdead, S_pq, Y_n)
    function construct_vertex(par::DGUnit)
        K_pll = par.K_pll
        T_int = par.T_int
        K_PT1 = par.K_PT1
        T_PT1 = par.T_PT1
        K_FRT = par.K_FRT
        I_max = par.I_max
        P_limit = par.P_limit
        Q_limit = par.Q_limit
        Vref = par.Vref
        Vdead = par.Vdead
        S_pq = par.S_pq
        Y_n = par.Y_n
        function rhs!(dx, x, e_s, e_d, p, t)
            u = complex(x[1], x[2])
            Vamp = abs(u)

            # this are the dynamic variables
            θ = x[3]
            P_int = x[4]
            Q_int = x[5]
            P_g = x[6]
            Q_g = x[7]

            # current from incoming/outgoing edges, nodal shunt, PQ load background
            i = total_current(e_s, e_d) + Y_n * u - conj(S_pq(Vamp)) / conj(u)

            # power flow tracking
            ΔP, ΔQ, global_over_voltage, global_under_voltage = p

            # detect error state
            over_voltage = Vamp > Vref + Vdead
            under_voltage = Vamp < Vref - Vdead
            # smooth_fault_state(x) = SmoothStep(x, Vref + Vdead, Vref - Vdead; order = 100)

            # don't integrate the error during faults
            if over_voltage | under_voltage | global_over_voltage | global_under_voltage
                dP_int = 0.
                dQ_int = 0.
                dP_g = 0.
                dQ_g = 0.
            else
                dP_int = -ΔP  #* smooth_fault_state(Vamp)
                dQ_int = -ΔQ  #* smooth_fault_state(Vamp)
                dP_g = (K_PT1 * P_int  - P_g)
                dQ_g = (K_PT1 * Q_int  - Q_g)
            end

            # dP_int = -ΔP
            # dQ_int = -ΔQ

            # dP_g = (K_PT1 * P_int / 2.0 - P_g) / T_PT1
            # dQ_g = (K_PT1 * Q_int / 2.0 - Q_g) / T_PT1

            # generator limits
            P_g = Limit(P_g, 0., P_limit) # 0.?
            Q_g = Limit(Q_g, -Q_limit, Q_limit)


            # PLL
            dθ = -K_pll * Vamp * sin(θ - angle(u))

            # FRT
            # additional reactive current during fault iq + iqplus
            iqplus = FRTcurrent(Vamp, K_FRT, Vref, Vdead)

            # current set point
            I_r = conj(complex(P_g, Q_g)) / Vamp

            # FRT limiter
            if over_voltage | under_voltage
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

            du = i - exp(θ * 1im) * I_r_limited  # = 0 constraint

            try
                dx[1] = real(du)
                dx[2] = imag(du)
                dx[3] = dθ
                dx[4] = dP_int
                dx[5] = dQ_int
                dx[6] = dP_g
                dx[7] = dQ_g
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
            dim = 7,
            mass_matrix = Diagonal([0, 0, 1, T_int, T_int, T_PT1, T_PT1]),
            sym = Symbol[
                :u_r,
                :u_i,
                :θ,
                :P_int,
                :Q_int,
                :P_g,
                :Q_g,
            ],
        )
    end
    symbolsof(::DGUnit) = begin
        [:u_r, :u_i, :θ, :P_int, :Q_int, :P_g, :Q_g]
    end
    dimension(::DGUnit) = begin
        7
    end
end
