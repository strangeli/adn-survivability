function testbed_chain(
    device,
    S,
    U = complex(1.0, 0.0),
    Y = complex(270.0, -380.0),
    P_ref = t -> -2,
    Q_ref = t -> 0.0,
)
    busses = [SlackAlgebraic(; U = U), device, PQAlgebraic(; S = S)]
    lines = [
        StaticLine(; from = 1, to = 2, Y = Y),
        StaticLine(; from = 2, to = 3, Y = Y),
    ]
    pg = PowerGrid(busses, lines)
    ADN(pg, typeof(device), P_ref, Q_ref)
end

function testbed(
    device,
    U = complex(1.0, 0.0);
    P_ref = t -> -2,
    Q_ref = t -> 0.0,
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
    busses = [SlackAlgebraic(; U = U), device]
    lines = [T,]
    pg = PowerGrid(busses, lines)
    ADN(pg, typeof(device), P_ref, Q_ref)
end

function testbed_externalgrid(
    device,
    U = complex(1.0, 0.0);
    P_ref = t -> -2,
    Q_ref = t -> 0.0,
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
    G = ExtDampedGrid(; H=1., D=1., Γ=1., Ω=ω, V=U)
    busses = [G, device]
    lines = [T,]
    pg = PowerGrid(busses, lines)
    ADN(pg, typeof(device), P_ref, Q_ref)
end

function testbed_dline(
    device,
    U = complex(1.0, 0.0);
    P_ref = t -> -2,
    Q_ref = t -> 0.0,
)
    Z = complex(1.41282, 2.0191199999999996)
    L = RLLine(;
        from = 1,
        to = 2,
        R = real(Z) * base_admittance,
        L = (imag(Z) / ω) * base_admittance,
        ω0 = ω,
    )
    busses = [SlackAlgebraic(; U = U), device]
    lines = [L]
    pg = PowerGrid(busses, lines)
    ADN(pg, typeof(device), P_ref, Q_ref)
end

function testbed_line(
    device,
    U = complex(1.0, 0.0),
    Y = complex(270.0, -380.0);
    P_ref = t -> -2,
    Q_ref = t -> 0.0,
)
    busses = [SlackAlgebraic(; U = U), device, device]
    lines = [
        StaticLine(; from = 1, to = 2, Y = Y),
        StaticLine(; from = 2, to = 3, Y = Y),
    ]
    pg = PowerGrid(busses, lines)
    ADN(pg, typeof(device), P_ref, Q_ref)
end

# # all values in p.u.
# device = DGUnitPLLPQTracking(;
#     K_pll = 0.1, #Hz/pu
#     K_PT1 = 1.0, # unit: [y]/[u] = [P]/[P] = 1
#     T_PT1 = 10.0, # unit: s
#     S_pq = complex(-1.0, -0.25), # pu
#     Y_n = 0.0,
# )

function testbed_CIGRE(
    DG_unit_type;
    idx = 1:12,
    P_ref = t -> -2,
    Q_ref = t -> 0.0,
    quad = 0,
)
    busses_static, lines, T, elist, Zs, Yshs = CIGRE_static()
    dir = @__DIR__
    power_flow = pf_sol(
        busses_static,
        lines;
        matfile = "$dir/LFCigre.mat",
        outfile = nothing,
    )
    busses = copy(busses_static)
    DG_locs = 2:12
    for i in DG_locs
        S_bkgrnd = zero(im)
        try
            S_bkgrnd = busses_static[i].S
        catch
            S_bkgrnd = complex(busses_static[i].P, busses_static[i].Q)
        end
        busses[i] = eval(Symbol(DG_unit_type))(;
            K_pll = 1632.993, #Hz/pu
            T_int = 2.,
            K_PT1 = 1.0,
            T_PT1 = 1e-8,
            K_FRT = 2.0,
            I_max = 1.0, # pu
            P_limit = 1.0, #pu
            Q_limit = 1.0, #pu
            Vref = 1.0, # pu
            Vdead = 0.1, # pu
            S_pq = V -> S_bkgrnd * (quad * (V/power_flow.V[i])^2 + 1 - quad),
            Y_n = 0.0,
        )
    end
    pg = PowerGrid(busses[idx], lines[1:length(idx)-1])
    ADN(pg, DG_unit_type, P_ref, Q_ref), power_flow
end





function testbed_CIGRE_dynamic(
    DG_unit_type;
    idx = 1:12,
    P_ref = t -> -2,
    Q_ref = t -> 0.0,
)
    busses_static, lines_static, T, elist, Zs, Yshs = CIGRE_static()
    busses = copy(busses_static)
    DG_locs = 2:12
    for i in DG_locs
        S_bkgrnd = zero(im)
        try
            S_bkgrnd = busses_static[i].S
        catch
            S_bkgrnd = complex(busses_static[i].P, busses_static[i].Q)
        end
        busses[i] = eval(Symbol(DG_unit_type))(;
            K_pll = 1632.993, #Hz/pu
            T_int = 2.,
            K_PT1 = 1.0,
            T_PT1 = 1e-12,
            K_FRT = 2.0,
            I_max = 1.0, # pu
            Vref = 1.0, #power_flow.V[i], #1.0, # pu
            Vdead = 0.1, #pu
            S_pq = S_bkgrnd, # TODO: check sign of PQ node
            Y_n = 0.0,
        )
    end

    lines = Array{AbstractLine,1}([])
    # # For ssp=0, replace it with PiModelLine
    # # assume t_km = 1, t_mk real
    # Y = PiModel(T) #./ admittance_m
    # ue = sqrt(-Y[2, 2] / Y[1, 1])
    # Z = inv(Y[1, 2] / ue)
    # TD = RLLine(
    #     from = 1,
    #     to = 2,
    #     R = real(Z),
    #     L = (imag(Z) / ω),
    #     ω0 = ω
    # )
    push!(lines, T)
    for (e, Z, Ysh) in zip(elist, Zs, Yshs)
        push!(
            lines,
            RLLine(
                from = first(e),
                to = last(e),
                R = real(Z) * base_admittance,
                L = (imag(Z) / ω)  * base_admittance,
                ω0 = ω
            ),
        )
    end
    pg = PowerGrid(busses[idx], lines[1:length(idx)-1])
    ADN(pg, DG_unit_type, P_ref, Q_ref)
end
