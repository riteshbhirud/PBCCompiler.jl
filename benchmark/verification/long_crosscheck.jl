# Long-running cross-check: the normal-form verdict must never contradict the dense
# oracle. Usage (from the repository root):
#
#     julia -tauto --project=test benchmark/verification/long_crosscheck.jl [ncases] [seed]
#
# Prints a summary line. Exits with status 1 on any contradiction.

using PBCCompiler
using PBCCompiler.Verification: random_clifford_t_circuit, dense_oracle, _width
using Random

include(joinpath(@__DIR__, "..", "..", "test", "verification", "mutations.jl"))

ncases = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2000
seed = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
rng = MersenneTwister(seed)
checked = 0; nf_eq = 0; dense_neq = 0; contradictions = 0; t0 = time()
while checked < ncases
    n = rand(rng, 1:5)
    a = random_clifford_t_circuit(rng, n; depth=rand(rng, 1:4n+2), t_density=0.2, multiqubit_density=0.2,
                                  meas_density=0.1, cc_density=0.1, sorted_qubits=rand(rng, Bool))
    b = try PBCCompiler.preprocess_circuit(copy(a)) catch; continue end
    _width(b) <= 8 || continue
    rand(rng) < 0.5 && (b = first(mutate(rng, b, rand(rng, MUTATION_KINDS); n_data=n)))
    tc = time()
    dv = dense_oracle(a, b; n_data=n, outcome_map=:infer, max_qubits=8)
    time() - tc > 5 && (println("  slow dense case: n=$n, |a|=$(length(a)), |b|=$(length(b)), width=$(_width(b)), $(round(time()-tc;digits=1)) s"); flush(stdout))
    dv.decided || continue
    r = verify_compilation(a, b; oracle=:none)
    global checked += 1
    if checked % 100 == 0
        println("  ... $checked cases, $(round(time() - t0; digits=1)) s"); flush(stdout)
    end
    global nf_eq += r.status == EQUIVALENT
    global dense_neq += !dv.agree
    if r.status == EQUIVALENT && !dv.agree
        global contradictions += 1
        println("CONTRADICTION (seed $seed, case $checked)")
    end
end
println("long cross-check: $checked cases, $nf_eq proved EQUIVALENT by the normal form, ",
        "$dense_neq refuted by the dense oracle, $contradictions contradictions, ",
        round(time() - t0; digits=1), " s")
contradictions == 0 || exit(1)
