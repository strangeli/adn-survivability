using PowerDynamics: AbstractLine
###### knoten ######

function CIGRE_static()

    LoadP = [
        19839E3,
        0.0,
        501.7E3,
        431.65E3,
        727.5E3,
        548.05E3,
        76.5E3,
        586.85E3,
        573.75E3,
        543.3E3,
        329.8E3,
    ]


    cosϕ = [
        0.9737554,
        0.0,
        0.9231825,
        0.9700009,
        0.9699996,
        0.9700017,
        0.8500022,
        0.9699995,
        0.8499989,
        0.9586623,
        0.969997,
    ]

    LoadQ = LoadP .* sin.(acos.(cosϕ)) ./ cosϕ
    LoadQ[2] = 0.0

    begin
        busses_static = Array{PowerDynamics.AbstractNode,1}([])
        push!(busses_static, SlackAlgebraic(U = 110E3 / base_voltage_HV))

        for (P, Q) in zip(LoadP, LoadQ)
            try
                push!(
                    busses_static,
                    PQAlgebraic(S = - complex(P, Q) ./ base_power),
                )
            catch
                push!(
                    busses_static,
                    PQAlgebraic(P = - P / base_power, Q = - Q / base_power),
                )
            end
        end
    end



    ###### Kanten ######



    T = OLTC(
        from = 1,
        to = 2,
        Uos = 110E3 / base_voltage_HV, # Bemessungsspannung Oberspannungsseite in kV
        Uus = 20E3 / base_voltage,# Bemessungsspannung Unterspannungsseite in kV
        k = 0, # Kennzahl der Schaltgruppe, Bsp: YD5-> 5
        ssp = 10.0, # Stufenschalterposition
        stufe_us = 0.625, # Stufung pro Stufenschalterstellung in %
        Sr = 25E6 / base_power, #Bemessungsscheinleistung in MVA
        uk = 12.0, # Kurzschlussspannung in %
        Pvk = 25E3 / base_power, # Kupferverluste in MW
        Pvl = 0.0, # Eisenverluste in kW
        iLeer = 0.0, # Leerlaufstrom in %
    )

    #ldata = CSV.read("$dir/lines.csv"; header = true)
    ldata = [ #L_km(nicht ändern)
        2.82,
        4.42,
        0.61,
        1.3,
        0.56,
        1.54,
        1.67,
        0.32,
        0.77,
        0.33,
        0.24,
    ]

    elist = [
        (2, 3),
        (3, 4),
        (4, 5),
        (4, 9),
        (5, 6),
        (6, 7),
        (8, 9),
        (9, 10),
        (10, 11),
        (11, 12),
        #(7, 8), # Schalter
        #(12, 5), # Schalter
    ]

    R1 = 0.501 # Ω/km
    X1 = 0.716 # Ω/km
    C1 = 0.1511749E-6 # F/km

    Zs = complex(R1, X1) .* ldata #[!, Symbol("L_km(nicht ändern)")]
    Yshs = 1im .* ω .* C1 .* ldata #[!, Symbol("L_km(nicht ändern)")]

    begin
        lines = Array{AbstractLine,1}([])
        push!(lines, T)
        for (e, Z, Ysh) in zip(elist, Zs, Yshs)
            push!(
                lines,
                PiModelLine(
                    from = first(e),
                    to = last(e),
                    y = inv(Z) / base_admittance,
                    y_shunt_km = Ysh / 2.0 / base_admittance,
                    y_shunt_mk = Ysh / 2.0 / base_admittance,
                ),
            )
        end
    end




    busses_static, lines, T, elist, Zs, Yshs
end
