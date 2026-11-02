# Plot results/timing.csv (written by timing.jl) into timing_vs_*.png next to it.
#
#     julia --project=test benchmark/verification/plot_timing.jl

using CairoMakie

csv = joinpath(@__DIR__, "results", "timing.csv")
lines = readlines(csv)
header = Symbol.(split(lines[1], ","))
rows = [Dict(zip(header, split(l, ","))) for l in lines[2:end]]
families = unique(r[:family] for r in rows)

for xcol in (:qubits, :t_count, :outcomes)
    fig = Figure(size=(640, 420))
    ax = Axis(fig[1, 1]; xlabel=string(xcol), ylabel="verification time (s)", yscale=log10,
              title="Verification time vs $(xcol)")
    for fam in families
        rs = sort([r for r in rows if r[:family] == fam]; by = r -> parse(Float64, r[xcol]))
        scatterlines!(ax, [parse(Float64, r[xcol]) for r in rs], [parse(Float64, r[:seconds]) for r in rs]; label=fam)
    end
    axislegend(ax; position=:lt)
    save(joinpath(@__DIR__, "results", "timing_vs_$(xcol).png"), fig)
end
println("plots written to ", joinpath(@__DIR__, "results"))
