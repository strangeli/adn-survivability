const PERFORM_CHECKS = false

using MAT
using DataFrames
using CSV

using ForwardDiff: jacobian
using SparseArrays: sparse
using LinearAlgebra: eigvals

function map_complex_2D(array)
    n = length(array)
    out = zeros(2n)
    for i = 1:n
        out[2i-1] = real(array[i])
        out[2i] = imag(array[i])
    end
    return out
end

function map_2D_complex(array)
    return array[1:2:end] .+ 1im .* array[2:2:end]
end

function stable_power_flow(pg, ic_guess::AbstractArray; sparse_jac=false, project=false, eig=true, p=nothing, t0=0.0)
    pgr = rhs(pg)
    M = Array(pgr.mass_matrix)
    n = systemsize(pg)

    f!(dx, x) = pgr(dx, x, p, t0)
    f(x) = (dx = similar(x); f!(dx, x); dx)
    j!(J, x) = (J .= jacobian(f, x))
    j(x) = (dx = similar(x); jacobian(f!, dx, x))

    if project
        mpm = M * pinv(M)
        proj_rhs!(dx, x) = (f!(dx, x); dx .-= mpm * dx)
        proj_guess = nlsolve(proj_rhs!, ic_guess; xtol=1E-10)
    end

    if sparse_jac
        # we can specify the Jacobian to be sparse
        dx0 = similar(ic_guess)
        F0 = similar(ic_guess)
        J0 = sparse(zeros(n, n))
        df = OnceDifferentiable(f!, j!, dx0, F0, J0)
        if project
            fixed_point = nlsolve(df, proj_guess.zero; xtol = 1E-10)
        else
            fixed_point = nlsolve(df, ic_guess; xtol = 1E-10)
        end
    else
        if project
            fixed_point = nlsolve(f!, j!, proj_guess.zero; xtol = 1E-10)
        else
            fixed_point = nlsolve(f!, j!, ic_guess; xtol = 1E-10)
        end
    end

    if eig
        λ = eigvals(Array(M), j(fixed_point.zero)) .|> real |> extrema
        println("Jacobian spectrum \nmin : ", first(λ), "\nmax : ",last(λ), "\nstable : ", last(λ)==0)
    end

    return State(pg, fixed_point.zero)
end

stable_power_flow(
    pg,
    ic_guess::State;
    sparse_jac = false,
    project = false,
    eig = true,
    p = nothing,
    t0 = 0.0,
) = stable_power_flow(
    pg,
    ic_guess.vec;
    sparse_jac = sparse_jac,
    project = project,
    eig = eig,
    p = p,
    t0 = t0,
)




function pf_sol(pg, ic_guess, p)
    #### solve power flow
    root_rhs = (dx, x) -> rhs(pg)(dx, x, p, 0.)
    pf = nlsolve(root_rhs, ic_guess, autodiff = :forward, method = :newton, xtol = 1E-10).zero #
    return State(pg, pf)
end

function pf_sol(busses_static, lines; matfile="$dir/LFCigre.mat", outfile="$dir/results/power_flow_sol.csv")

    #### read from Matlab
    matlab_data = matread(matfile)
    uK = Array(matlab_data["uK"])[:]
    uK[1] *= 1e3/base_voltage_HV
    uK[2:end] .*= 1e3/base_voltage
    sK = Array(matlab_data["sK"])[:]

    #### solve power flow

    pg = PowerGrid(busses_static, lines)
    rpg! = rhs(pg)

    ic_guess = map_complex_2D(uK)
    root_rhs = RootRhs_ic(rhs(pg))
    pf = nlsolve(root_rhs, ic_guess, autodiff = :forward,  method = :newton, xtol = 1E-10).zero #

    line_currents = [PowerDynamics.get_current(State(pg, pf), k) for k = 1:length(busses_static)]
    sK_julia = map_2D_complex(pf) .* conj.(line_currents)

    power_flow_sol = DataFrame(
        V = abs.(map_2D_complex(pf)),
        V_matlab = abs.(uK),
        A = rad2deg.(angle.(map_2D_complex(pf))),
        A_matlab = rad2deg.(angle.(uK)),
        s_real = real.(sK_julia),
        s_real_matlab = real.(sK),
        s_imag = imag.(sK_julia),
        s_imag_matlab = imag.(sK),
        real = pf[1:2:end],
        imag = pf[2:2:end],)

    if outfile != nothing
        CSV.write(outfile, power_flow_sol)
    end

    return power_flow_sol
end

if PERFORM_CHECKS
    using GraphPlot
    using Plots
    using BenchmarkTools

    #### read from Matlab
    Q_marcel = matlab_data["LeistungsflussQuelle"]["Q0las"]
    P_marcel = matlab_data["LeistungsflussQuelle"]["P0las"]

    #### check output

    pf = map_complex_2D(complex.(power_flow_sol.real, power_flow_sol.imag))
    niceprint(pf)

    du = similar(pf)
    rpg!(du, ic_guess, 0, 0)
    println(sum(du))

    #### plot result

    # nodesize: Spannungsabfall [%]
    val = 100abs.(maximum(power_flow_sol.V[2:end]) .- power_flow_sol.V) ./ 20
    val[1] = val[2]
    gplot(pg.graph, nodesize = val, nodelabel = v_angle(pf))



    #### compare with Matlab

    @assert any(Q_marcel[2:end] .- LoadQ .== 0)
    @assert any(P_marcel[2:end] .- LoadP .== 0)

    # Trafo
    # need to multiply with .* [-1 -1; 1 1] to take the sign convention in PowerDynamics into account
    matlab_data["Y"]["TTVsoll"] .- PiModel(lines[1]) .* [-1 -1; 1 1] |> sum |> println

    # line admittances Vierpoldarstellung
    # need to multiply with .* [-1 -1; 1 1] to take the sign convention in PowerDynamics into account
    for (i, l) in enumerate(lines[2:end])
        PiModel(l) .* [-1 -1; 1 1] .- matlab_data["Y"]["LLsoll"][2i-1:2i, 2i-1:2i] |> sum |> println
    end

    @show abs.(map_2D_complex(pf)) ./ abs.(uK)
    @show angle.(map_2D_complex(pf)) - angle.(uK) .|> rad2deg

    currents = [PowerDynamics.get_current(State(pg, pf), k) for k = 1:length(busses)]

    @show sK_julia = map_2D_complex(pf) .* conj.(currents)

    @show real.(sK_julia) ./ real.(sK)
end
