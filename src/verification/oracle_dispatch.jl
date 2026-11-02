# Oracle fallback used by `verify_compilation` when the normal form does not
# decide. The dense and stabilizer oracles are added in V6. Until then only
# `oracle = :none` behaviour exists.

function _oracle_fallback(input, compiled, n; oracle, max_dense_qubits, max_branches, rng, semantics,
                          magic_k, outcome_map, stabilizer_trials, divergence, stats, reasons)
    push!(reasons, "no oracle available to decide the remaining question")
    (UNKNOWN, :none, nothing)
end
