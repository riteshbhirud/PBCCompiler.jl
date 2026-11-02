function _oracle_fallback(input, compiled, n; oracle, max_dense_qubits, max_branches, rng, semantics,
                          magic_k, outcome_map, stabilizer_trials, divergence, stats, reasons, visible_bits = nothing)
    N = max(n, _width(compiled), _width(input))
    if oracle == :none
        push!(reasons, "oracles disabled (oracle = :none)")
        return (UNKNOWN, :none, nothing)
    end
    if oracle == :dense || (oracle == :auto && N <= max_dense_qubits)
        t = time_ns()
        dv = dense_oracle(input, compiled; n_data = n, semantics, outcome_map, magic_k,
                          max_qubits = max_dense_qubits, max_branches, visible_bits)
        stats[:time_dense_oracle] = (time_ns() - t) / 1e9
        stats[:dense_branches] = dv.leaves_input + dv.leaves_compiled
        if dv.decided && dv.agree
            stats[:dense_max_deviation] = dv.max_deviation
            push!(reasons, "dense oracle: the programs agree to within $(dv.max_deviation) under " *
                           "an admissible outcome map (numerical, so not a proof)")
            return (UNKNOWN, :dense_oracle, nothing)
        elseif dv.decided
            cex = Counterexample(:dense_oracle, divergence, dv.witness_state, nothing, dv.witness_record,
                                 dv.p_input, dv.p_compiled, nothing, nothing, true, dv.note)
            return (NOT_EQUIVALENT, :dense_oracle, cex)
        end
        push!(reasons, "dense oracle undecided: " * dv.note)
        oracle == :dense && return (UNKNOWN, :dense_oracle, nothing)
    end
    t = time_ns()
    sv = stabilizer_oracle(input, compiled; n_data = n, outcome_map = outcome_map === :infer ? :identity : outcome_map,
                           magic_k, trials = stabilizer_trials, rng, visible_bits)
    stats[:time_stabilizer_oracle] = (time_ns() - t) / 1e9
    stats[:stabilizer_trials] = sv.trials
    if !sv.applicable
        push!(reasons, "stabilizer oracle not applicable: " * sv.note)
        return (UNKNOWN, :stabilizer_oracle, nothing)
    elseif !sv.refuted
        push!(reasons, "stabilizer oracle: " * sv.note)
        return (UNKNOWN, :stabilizer_oracle, nothing)
    end
    cex = Counterexample(:stabilizer_oracle, divergence, ComplexF64[], sv.witness_state, sv.witness_record,
                         sv.p_input, sv.p_compiled, nothing, nothing, false, sv.note)
    if outcome_map === :infer
        stats[:identity_map_witness] = cex
        push!(reasons, "stabilizer oracle: the programs differ under the identity outcome map " *
                       "(p = $(sv.p_input) vs $(sv.p_compiled)); with outcome_map = :infer this does not rule out a parity map")
        return (UNKNOWN, :stabilizer_oracle, nothing)
    end
    (NOT_EQUIVALENT, :stabilizer_oracle, cex)
end
