# Mutation-testing report for the verifier.
#
#     julia -tauto --project=test benchmark/verification/mutation_report.jl [trials] [seed]
#
# For small random programs, compile with PBCCompiler, apply each mutation kind to a
# copy of the compiled output, decide with the dense oracle (identity outcome map,
# PBCCompiler's contract) whether the mutation changes the program's behaviour, and
# record whether the verifier reports NOT_EQUIVALENT:
#   * with the default oracle stack, and
#   * with the dense oracle disabled (refutation by the stabilizer oracle only).
# Writes benchmark/verification/results/mutation_report.csv and prints a markdown table.

using PBCCompiler
using PBCCompiler.Verification: random_clifford_t_circuit, dense_oracle, _width
using Random, Printf

include(joinpath(@__DIR__, "..", "..", "test", "verification", "mutations.jl"))

trials = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 300
seed = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2026
rng = MersenneTwister(seed)
fields = (:applied, :changing, :detected, :by_normal_form, :by_dense, :by_stabilizer, :detected_stab_only,
          :false_refutation, :not_certified)
stats = Dict(k => Dict(f => 0 for f in fields) for k in MUTATION_KINDS)
misses = String[]
t0 = time()
for trial in 1:trials
    n = rand(rng, 1:4)
    a = random_clifford_t_circuit(rng, n; depth=rand(rng, 2:3n+2), t_density=0.3, multiqubit_density=0.2,
                                  meas_density=0.1, cc_density=0.1, sorted_qubits=true)
    b = PBCCompiler.preprocess_circuit(copy(a))
    for kind in MUTATION_KINDS
        bm, ok = mutate(rng, b, kind; n_data=n)
        ok || continue
        max(n, _width(bm)) <= 8 || continue
        dv = dense_oracle(a, bm; n_data=n, outcome_map=:identity, max_qubits=8)
        dv.decided || continue
        s = stats[kind]; s[:applied] += 1
        r = verify_compilation(a, bm; outcome_map=:identity, minimize=false, max_dense_qubits=8)
        rs = verify_compilation(a, bm; outcome_map=:identity, minimize=false, max_dense_qubits=0, rng=rng)
        r.status == EQUIVALENT || (s[:not_certified] += 1)
        if !dv.agree
            s[:changing] += 1
            if r.status == NOT_EQUIVALENT
                s[:detected] += 1
                s[r.method == :normal_form ? :by_normal_form : r.method == :dense_oracle ? :by_dense : :by_stabilizer] += 1
            else
                push!(misses, "$kind, trial $trial: $(r.status) via $(r.method)")
            end
            rs.status == NOT_EQUIVALENT && (s[:detected_stab_only] += 1)
        else
            (r.status == NOT_EQUIVALENT || rs.status == NOT_EQUIVALENT) && (s[:false_refutation] += 1)
        end
    end
end
mkpath(joinpath(@__DIR__, "results"))
open(joinpath(@__DIR__, "results", "mutation_report.csv"), "w") do io
    println(io, join(("mutation", fields...), ","))
    for k in MUTATION_KINDS
        println(io, join((k, (stats[k][f] for f in fields)...), ","))
    end
end
pct(a, b) = b == 0 ? "n/a" : @sprintf("%.1f%%", 100a / b)
println("| mutation | applied | channel-changing | detected | rate | by dense | by stabilizer | stabilizer-only rate | false refutations |")
println("|---|---|---|---|---|---|---|---|---|")
for k in MUTATION_KINDS
    s = stats[k]
    println("| $k | $(s[:applied]) | $(s[:changing]) | $(s[:detected]) | $(pct(s[:detected], s[:changing])) | ",
            "$(s[:by_dense]) | $(s[:by_stabilizer]) | $(pct(s[:detected_stab_only], s[:changing])) | $(s[:false_refutation]) |")
end
println("\nmisses: ", isempty(misses) ? "none" : "")
foreach(m -> println("  - ", m), misses)
println("seed $seed, $trials programs, $(round(time() - t0; digits=1)) s")
