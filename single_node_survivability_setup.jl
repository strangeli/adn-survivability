#using Pkg
#Pkg.activate(".")

using Distributed
using ClusterManagers

SLURM = false

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

@everywhere begin
    #using Pkg
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

@everywhere get_run(i, batch_size) =
    mod(i, batch_size) == 0 ? batch_size : mod(i, batch_size)
@everywhere get_batch(i, batch_size) = 1 + (i - 1) ÷ batch_size

@everywhere function FRT_curve(t; t_error = 0.1, Usoll = 1.0)
    FRT_lower_limit = 0.9 * Usoll
    FRT_upper_limit = 1.1 * Usoll

    if t_error > t
        FRT_lower_limit = 0.9 * Usoll
        FRT_upper_limit = 1.1 * Usoll
    end

    #TODO: more elegant way for this? (smooth step not needed since this function runs on a ready solution object)
    # ]0, 0.1]
    if t > t_error && t <= 0.1 + t_error
        FRT_lower_limit = 0.15 * Usoll
        FRT_upper_limit = 1.25 * Usoll
        # ]0.1, 0.15]
    elseif t > 0.1 + t_error && t <= 0.15 + t_error
        FRT_lower_limit = 0.15 * Usoll
        FRT_upper_limit = 1.2 * Usoll
        # ]0.15, 3]
    elseif t > 0.15 + t_error && t <= 3.0 + t_error
        FRT_lower_limit =
            ((0.85 - 0.15) * Usoll / (3.0 - 0.15)) * (t - t_error) +
            (0.15 * Usoll - ((0.85 - 0.15) * Usoll / (3.0 - 0.15)) * 0.15)  # 43/380
        FRT_upper_limit = 1.2 * Usoll
        # ]3, 5]
    elseif t > 3.0 + t_error && t <= 5.0 + t_error
        FRT_lower_limit = 0.85 * Usoll
        FRT_upper_limit = 1.2 * Usoll
        # ]5, 60]
    elseif t > 5.0 + t_error && t <= 60.0 + t_error
        FRT_lower_limit = 0.85 * Usoll
        FRT_upper_limit = 1.15 * Usoll
        # ]60, Inf]
    elseif t > 60.0 + t_error
        FRT_upper_limit = 1.1 * Usoll
        FRT_lower_limit = 0.9 * Usoll
    end

    (FRT_lower_limit, FRT_upper_limit)
end

@everywhere function CIGRE_prob(;
    P_ref = 0.0,
    Q_ref = 0.0,
    nsc_node = 9,
    t_duration = 0.15,
    resistance = 10.0,
    quad = 0.0,
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
            S_pq = V -> S_bkgrnd * (quad * V^2 + 1 - quad), #power_flow[i, :v]
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
        verbose ? println("found fault state\n", State(cpg, integrator.u)[:, :v]) : nothing
        #integrator.u = x2
        #u_modified!(integrator,true)
    end

    function regularState(integrator)
        sol2 = integrator.sol
        verbose ? println("searching post-fault state") : nothing
        #x3 = find_valid_initial_condition(powergrid, sol2[end]) # Jump the state to be valid for the new system.
        integrator.f = ode
        DiffEqBase.initialize_dae!(integrator, BrownFullBasicInit())
        verbose ? println("found post-fault state\n", State(cpg, integrator.u)[:, :v]) :
        nothing
        #integrator.u = x3
        #u_modified!(integrator,true)
    end

    # discrete callback requires tstop
    cb1 = DiscreteCallback(((u, t, integrator) -> t in nsc.tspan_fault[1]), errorState)
    cb2 = DiscreteCallback(((u, t, integrator) -> t in nsc.tspan_fault[2]), regularState)

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

## setup EnsembleProblem
@everywhere quad = 1

@everywhere cpg, prob = CIGRE_prob(;
    P_ref = 24.0,
    Q_ref = 5.0,
    nsc_node = 9,
    t_duration = 0.15,
    resistance = 3.5,
    quad = quad,
    get_cpg = true,
    verbose = true,
)

@everywhere begin
    # load flow sol

    S_total = complex(24.373141067954748, 6.115974820490637) * 1e6 / base_power
    nsc_nodes = 2:12
    n_nodes_LV = length(nsc_nodes)
    batch_size = 1000

    Tdist = Normal(150e-3, 10e-3)
    Rdist = Uniform(3.0, 10.0)


    NT = n_nodes_LV * batch_size

    # use the same perturbations for each node
    T = rand(Tdist, batch_size)
    R = rand(Rdist, batch_size)
end

## sampling, skip for plotting only

# naive binomial proportion
@everywhere function bp(sample)
    # add one failure + one success
    N = length(sample) + 2
    m = (1 + sum(sample)) / N
    m ± sqrt(m * (1 - m) / N)
end

@everywhere function prob_func(prob, i, repeat)
    return CIGRE_prob(;
        P_ref = real(S_total),
        Q_ref = imag(S_total),
        nsc_node = nsc_nodes[get_batch(i, batch_size)],
        t_duration = T[get_run(i, batch_size)],
        resistance = R[get_run(i, batch_size)],
        quad = quad,
        get_cpg = false,
        verbose = false,
    )
end

@everywhere begin
    r_idx = ones(Int, length(cpg.nodes))
    for n = 2:length(cpg.nodes)
        r_idx[n] = r_idx[n-1] + dimension(cpg.nodes[n-1])
    end
    i_idx = ones(Int, length(cpg.nodes)) .+ r_idx

    t_eval = range(prob.tspan...; step = 0.001)
    lower_lim = first.(FRT_curve.(t_eval))
    upper_lim = last.(FRT_curve.(t_eval))

    function observer(sol, i)
        v = sqrt.(sol(t_eval, idxs = r_idx) .^ 2 .+ sol(t_eval, idxs = i_idx) .^ 2)
        in_bounds = [all(lower_lim[i] .< v[:, i] .< upper_lim[i]) for i = 1:length(t_eval)]
        surv_idx = findfirst(.!in_bounds)
        survived = all(in_bounds)
        if isnothing(surv_idx)
            nodes_out = nothing
            surv_time = t_eval[end]
        else
            nodes_out = findall(
                (lower_lim[surv_idx] .≥ v[:, surv_idx]) .|
                (v[:, surv_idx] .≥ upper_lim[surv_idx]),
            )
            surv_time = t_eval[surv_idx]
        end
        if sol.retcode == :Success
            return ([survived, surv_time, nodes_out, sol.retcode], false)
        else
            # count solver error as indication of non-survival
            return ([false, sol.t[end], nodes_out, sol.retcode], false)
        end
    end
end

@everywhere function reduction(u, batch, I)
    successes = findall([last(b) == :Success for b in batch])
    survs = [first(b) for b in batch]
    surv = survs |> bp
    n_successes = length(successes)
    mean_surv_time = [b[2] for b in batch] |> mean

    out_of_bounds = vcat([b[3] for b in batch]...)
    nodes_out = findall(out_of_bounds .!== nothing)
    if isempty(nodes_out)
        fpn = (missing, missing)
        npn = (missing, missing)
    else
        unique_nodes = unique(out_of_bounds[nodes_out])
        if length(unique_nodes) == 1
            fpn = tuple(first(unique_nodes), missing) .- 1
            npn = tuple(1, missing)
        else
            score = [count(==(n), out_of_bounds) for n in unique_nodes]
            idx = sortperm(score)
            fpn = tuple(unique_nodes[idx[end]], unique_nodes[idx[end-1]]) .- 1
            npn = tuple(score[idx[end]], score[idx[end-1]])
        end
    end

    out = [surv.val surv.err n_successes mean_surv_time fpn npn]

    println(out)
    if isempty(u)
        return (out, false)
    else
        return ([u; out], false)
    end
end

mcprob = EnsembleProblem(
    prob;
    prob_func = prob_func,
    output_func = observer,
    reduction = reduction,
    u_init = [],
)

sim = solve(
    mcprob,
    Rodas4(),
    EnsembleDistributed(),
    trajectories = NT,
    batch_size = batch_size,
)
@show sim.elapsedTime

result = DataFrame(
    idx = [Symbol("MV_$a") for a = 1:length(nsc_nodes)],
    deg = length.(cpg.graph.fadjlist)[nsc_nodes],
    surv = sim.u[:, 1],
    err = sim.u[:, 2],
    ret_suc = sim.u[:, 3],
    time = sim.u[:, 4],
    first_location = first.(sim.u[:, 5]),
    first_freq = first.(sim.u[:, 6]),
    sec_location = last.(sim.u[:, 5]),
    sec_freq = last.(sim.u[:, 6]),
)


#now = unix2datetime(time())
CSV.write("$dir/results/single_node_mc_res_NT$(batch_size)_quad$(quad).csv", result) # _$(Dates.format(now, "e_d_u_yy-H_M_S"))

## plots

dir = @__DIR__
result_PQ = DataFrame(CSV.File("$dir/results/single_node_mc_res_NT1000_quad0.csv"))
result_ZIP = DataFrame(CSV.File("$dir/results/single_node_mc_res_NT1000_quad1.csv"))


using GraphPlot, Colors, Plots
using LaTeXStrings
#using StatsPlots

cmap = Colors.sequential_palette(200, 100; logscale = true);

#node_colours = result_PQ.surv .* length(cmap) .|> floor .|> Int;
#c=cmap[node_colours]

##
plot(
    legendfont = (18, "times"),
    guidefont = (22, "times"),
    tickfont = (18, "times"),
    framestyle = :box,
    dpi = 300,
    size = (800, 400),
)
scatter!(
    result_PQ.surv .± result_PQ.err,
    c = :blue,
    ms = 10,
    ylims = (0.85, 1.02),
    xlabel = "bus index i",
    ylabel = L"$\mu (i | R_{on}, \Delta t_{fault})$",
    xticks = nsc_nodes .- 1,  #[Symbol("MV-$a") for a = 1:length(nsc_nodes)]),
    legend = :bottomright,
    label = "PQ loads",
)
#node_colours = result_ZIP.surv .* length(cmap) .|> floor .|> Int;
scatter!(
    result_ZIP.surv .± result_ZIP.err,
    marker = :d,
    c = :orange,
    ms = 10,
    label = "Z loads",
)
savefig("$dir/plots/single_node_surv.pdf")
##

let res = result_PQ
    z = zeros(length(nsc_nodes), length(nsc_nodes))
    for i = 1:last(size(z))
        if ismissing(res.first_location[i])
            continue
        else
            z[i, res.first_location[i]] = res.first_freq[i]
        end
        if ismissing(res.sec_location[i])
            continue
        else
            z[i, res.sec_location[i]] = res.sec_freq[i]
        end
    end

    heatmap(
        z',
        xticks = (nsc_nodes .- 1, res.idx),
        yticks = (nsc_nodes .- 1, res.idx),
        xguide = "SC node",
        yguide = "likely perturbed node",
        colorbar_title = "counts (total=$(batch_size))",
    )
    savefig("$dir/plots/single_node_surv_most_affected_nodes.pdf")
end

node_colours = result_PQ.surv .* length(cmap) .|> floor .|> Int;
labels = ["ext. grid"; ["MV-$(lpad(i,2,"0"))" for i = 1:11]]
gplot(
    cpg.graph,
    nodelabeldist = 2,
    nodelabelangleoffset = π / 2,
    nodelabel = labels,
    nodelabelc = colorant"black",
    nodefillc = [RGB(0, 0, 0); cmap[node_colours]],
)

#using Cairo, Compose
# # save to pdf
#draw(PDF("karate.pdf", 16cm, 16cm), gplot(g))
