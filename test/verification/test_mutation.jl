@testitem "V9: mutation testing (100% detection of channel-changing mutants)" tags=[:verification, :v9] begin
using PBCCompiler
using PBCCompiler.Verification: random_clifford_t_circuit, dense_oracle, _width
using Random
include(joinpath(pkgdir(PBCCompiler), "test", "verification", "mutations.jl"))

rng = MersenneTwister(9090)
stats = Dict(k => Dict(:applied => 0, :changing => 0, :detected => 0, :false_refutation => 0,
                       :detected_stab => 0, :false_refutation_stab => 0) for k in MUTATION_KINDS)
misses = String[]
for trial in 1:60
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
        # independent refutation path: dense oracle disabled, stabilizer oracle only
        rs = verify_compilation(a, bm; outcome_map=:identity, minimize=false, max_dense_qubits=0, rng=rng)
        if !dv.agree
            rs.status == NOT_EQUIVALENT && (s[:detected_stab] += 1)
            s[:changing] += 1
            if r.status == NOT_EQUIVALENT
                s[:detected] += 1
            else
                push!(misses, "$kind (trial $trial): $(r.status) $(get(r.stats, :reasons, String[]))")
            end
        else
            r.status == NOT_EQUIVALENT && (s[:false_refutation] += 1)
            rs.status == NOT_EQUIVALENT && (s[:false_refutation_stab] += 1)
        end
    end
end
@test isempty(misses)
@test all(s -> s[:false_refutation] == 0 && s[:false_refutation_stab] == 0, values(stats))
# the stabilizer oracle alone (randomized) catches nearly all of them
@test sum(s[:detected_stab] for s in values(stats)) >= 0.9 * sum(s[:changing] for s in values(stats))
@test all(s -> s[:detected] == s[:changing], values(stats))
@test sum(s[:changing] for s in values(stats)) >= 100
# every mutation kind was actually exercised
@test all(k -> stats[k][:changing] >= 1, MUTATION_KINDS)
isempty(misses) || foreach(println, misses)
end
