# Verification timing versus qubit count, T count and outcome count.
#
#     julia -tauto --project=test benchmark/verification/timing.jl
#
# Families: approximate QFT-n, Cuccaro adders, seeded random Clifford+T circuits,
# and the repository fixture. Each program is compiled with PBCCompiler and then
# verified (normal-form path, oracles disabled). The best of three runs is
# recorded. Writes benchmark/verification/results/timing.csv. Plot it with
# plot_timing.jl.

using PBCCompiler
using PBCCompiler: parse_input
using PBCCompiler.Verification: random_clifford_t_circuit, _width, _opkind, OP_EIGHPI, OP_MEASUREMENT
using Random

include(joinpath(@__DIR__, "circuits.jl"))

count_kind(c, k) = count(op -> _opkind(op) == k, c)

function measure(family, param, input)
    compiled = PBCCompiler.preprocess_circuit(copy(input))
    verify_compilation(input, compiled; oracle=:none)           # warm-up / compile
    best = Inf; status = nothing
    for _ in 1:3
        t = @elapsed r = verify_compilation(input, compiled; oracle=:none)
        best = min(best, t); status = r.status
    end
    (family=family, param=param, qubits=_width(input), t_count=count_kind(input, OP_EIGHPI),
     outcomes=count_kind(compiled, OP_MEASUREMENT), input_ops=length(input), compiled_ops=length(compiled),
     status=string(status), seconds=best)
end

rows = []
for n in 2:2:16
    push!(rows, measure("approx_qft", n, approx_qft(n)))
end
for k in 1:5
    push!(rows, measure("cuccaro_adder", k, cuccaro_adder(k)))
end
rng = MersenneTwister(1234)
for n in (2, 4, 8, 12, 16), td in (0.1, 0.3)
    push!(rows, measure("random_t$(td)", n, random_clifford_t_circuit(rng, n; depth=8n, t_density=td, sorted_qubits=true)))
end
push!(rows, measure("fixture", 0, parse_input(joinpath(@__DIR__, "..", "..", "test", "fixtures", "sample_input.qasm"))))

mkpath(joinpath(@__DIR__, "results"))
out = joinpath(@__DIR__, "results", "timing.csv")
open(out, "w") do io
    println(io, join(keys(rows[1]), ","))
    for r in rows
        println(io, join(values(r), ","))
    end
end
for r in rows
    println(rpad(r.family, 14), " param=", rpad(r.param, 3), " qubits=", rpad(r.qubits, 3), " T=", rpad(r.t_count, 4),
            " outcomes=", rpad(r.outcomes, 4), " ", rpad(r.status, 11), " ", round(r.seconds * 1000; digits=2), " ms")
end
println("wrote ", out)
