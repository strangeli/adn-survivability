begin
    using Dates
    using CSV
    using DataFrames
    using Query
    using StatsPlots
    using Measurements
    using Distributions
    using LaTeXStrings
end

##

function bp(sample)
    # add one failure + one success
    N = length(sample) + 2
    m = (1 + sum(sample)) / N
    m ± sqrt(m * (1 - m) / N)
end

function survPQ(P_range, Q_range, P_step, Q_step, result, threshold=100)
    surv = zeros(length(Q_range), length(P_range))
    serr = similar(surv)
    for (i, P) in enumerate(P_range)
        for (j, Q) in enumerate(Q_range)
            # symmetric window of width P_step x Q_step
            df = result[(P - P_step/2 .< result.P .< P + P_step/2) .& (Q - Q_step/2 .< result.Q .< Q + Q_step/2), :]
            if size(df, 1) > threshold
                s = bp(df.nsc_surv)
                surv[j, i] = s.val
                serr[j, i] = s.err
            else
                surv[j, i] = NaN
                serr[j, i] = NaN
            end
            # println("idx_P: ", idx_P, " idx_Q: ", idx_Q)
            # println(size(idx_P))
        end
    end
    surv, serr
end

##

const base_power = 1E6 # 1MW
S_total = complex(24.373141067954748, 6.115974820490637) * 1e6 / base_power
# S_total = complex(0.0, 0.0) * 1e6 / base_power


## load result from CSV
dir = @__DIR__

resultPQ = CSV.read("$dir/results/mc_res_Fri_21_Aug_20-16_43_13_step_nsc_node_4_quad0_error.csv")
resultPQ = sort(resultPQ, (:P, :Q))
resultZ = CSV.read("$dir/results/mc_res_Fri_21_Aug_20-7_17_29_step_nsc_node_4_quad1_error.csv")
resultZ = sort(resultZ, (:P, :Q))

##

P_min = -15.0 + real(S_total)
Q_min = -10.0 + imag(S_total)
P_step = .5 # window size
Q_step = .5 # window size
P_max = 15.0 + real(S_total)
Q_max = 10.0 + imag(S_total)
P_range = P_min:0.1:P_max
Q_range = Q_min:0.1:Q_max

##
surviPQ, err = survPQ(P_range, Q_range, P_step, Q_step, resultPQ);
@show err[.!isnan.(err)] |> extrema

##

surviZ, err = survPQ(P_range, Q_range, P_step, Q_step, resultZ);
@show err[.!isnan.(err)] |> extrema

##

let quad = 1
    default(
            legendfont = (12, "times"),
            titlefont = (14, "times"),
            guidefont = (16, "times"),
            tickfont = (14, "times"),
            framestyle = :box,
        )
    heatmap(
        P_range,
        Q_range,
        quad == 0 ? surviPQ : surviZ,
        xlabel = L"P_{ref} [MW]",
        ylabel = L"Q_{ref} [Mvar]",
        colorbar = :right,
        title = L"$\mu (3 | R_{on}, \Delta t_{fault})$",
        legend=:bottomleft,
        clims = (0, 1),
        size = (600, 400),
        xlims=(12, 30),
        ylims=(-6, 18),
        dpi=300,
        margin=Plots.mm,
    )
    scatter!(
        [real(S_total)],
        [imag(S_total)],
        c = :green,
        ms = 5,
        label = L"\{P_{PF}, Q_{PF}\}",
    )

    if quad == 1
        savefig("$dir/plots/survivability_Z_slide.png")
    else
        savefig("$dir/plots/survivability_PQ_slide.png")
    end
end

##




##







#result.step_success = result.step_success .|> Symbol |> categorical
result.nsc_success = result.nsc_success .|> Symbol |> categorical



# how many successful integrations?
result |>
@filter((_.nsc_success == "Success") | (_.nsc_success == :Success)) |>
DataFrame |>
size

result |>
@filter((_.nsc_success == "Success") | (_.nsc_success == :Success)) |>
@df scatter(:P, :Q, :nsc_surv, xlabel = "P_ref [pu]", ylabel = "Q_ref [pu]")



## minimum voltage during fault
result |>
@filter((_.nsc_success == "Success") | (_.nsc_success == :Success)) |>
@df scatter(
    :R,
    :nsc_min_v,
    group = :nsc_surv,
    ylabel = "avg min V [pu]",
    xlabel = "R [Ohm]",
    legend = :topleft,
    #xscale=:log10,
)
##

result |>
@filter((_.nsc_success == "Success") | (_.nsc_success == :Success)) |>
@df scatter(
    :T,
    :nsc_min_v,
    group = :nsc_surv,
    ylabel = "avg min V [pu]",
    xlabel = "t_d [s]",
    legend = :topleft,
    #xscale=:log10,
)

result |>
@filter((_.nsc_success == "Success") | (_.nsc_success == :Success)) |>
@df scatter(
    :P,
    :Q,
    group = :nsc_surv,
    xlabel = "P_ref [pu]",
    ylabel = "Q_ref [pu]",
    legend = :topright,
    ms = 3,
    #xscale=:log10,
)
scatter!([real(S_total)], [imag(S_total)], ms = 10, label = "power flow")
savefig("$dir/plots/surv.pdf")


k = 0.5
result[!, :catP] = round.(result.P .* k, digits = 1) ./ k
result[!, :catQ] = round.(result.Q .* k, digits = 1) ./ k
# assume that the problematic runs indicate non-survival
df = by(result, [:catP, :catQ], nsc_surv = :nsc_surv => mean)
dat = unstack(df, :catP, :catQ, :nsc_surv)
heatmap(
    sort(unique(df.catP)),
    sort(unique(df.catQ)),
    Array(dat)[:, 2:end]',
    xlabel = L"P_{ref} [MW]",
    ylabel = L"Q_{ref} [Mvar]",
    colorbar_title = "short-circuit survivability",
    clims = (0, 1),
)
scatter!([real(S_total)], [imag(S_total)], c = :green, ms = 5, label = "initial power flow")
savefig("$dir/plots/survivability_PQ.pdf")




df_err = by(result, [:catP, :catQ], err = :nsc_surv => x -> bp(x).err)
dat_err = unstack(df_err, :catP, :catQ, :err) |> Array |> transpose
heatmap(
    sort(unique(df_err.catP)),
    sort(unique(df_err.catQ)),
    dat_err[2:end, :],
    xlabel = "P_ref [pu]",
    ylabel = "Q_ref [pu]",
    colorbar_title = "NSC survivability std.err.",
)


# groups for binning
result[!, :catS] = round.(abs.(complex.(result.P, result.Q)), digits = 1)
result |>
@filter((_.nsc_success == "Success") | (_.nsc_success == :Success)) |>
@groupby(_.catS) |>
@map({S = key(_), Survivability = bp(_.nsc_surv)}) |>
DataFrame |>
@df scatter(
    :S,
    :Survivability,
    #err = :err,
    ylabel = "NSC survivability",
    xlabel = "|S_ref| [pu]",
    legend = false,
    ylims = (0.0, 1),
)
vline!([abs(S_total)], label = "power flow")

result[!, :catW] = round.(angle.(complex.(result.P, result.Q)), digits = 2)
result |>
@filter((_.nsc_success == "Success") | (_.nsc_success == :Success)) |>
@groupby(_.catW) |>
@map({S = key(_), Survivability = bp(_.nsc_surv)}) |>
DataFrame |>
@df scatter(
    :S,
    :Survivability,
    #err = :err,
    ylabel = "NSC survivability",
    xlabel = "angle(S_ref) [pu]",
    legend = false,
    ylims = (0.0, 1),
)
vline!([angle(S_total)], label = "power flow")

savefig("$dir/plots/surv_angle.pdf")



#cmap = map(s -> s == :Success ? :blue : :green, last.(sim.u))
#scatter(resistance, first.(sim.u), c=cmap, ylabel="norm(V) [pu]", xlabel="R [Ohm]", legend=false)
