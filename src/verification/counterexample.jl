# V7: counterexample minimization and confirmation.

# Generic delta debugging (ddmin): a 1-minimal sub-sequence of `xs` on which `test`
# still holds. `test` is assumed to hold on `xs`. Stops after `max_tests` calls.
function _ddmin(test, xs::Vector{T}; max_tests::Int = 500) where {T}
    ntests = 0
    k = 2
    while length(xs) >= 2 && ntests < max_tests
        chunk = cld(length(xs), k)
        reduced = false
        for s in 1:chunk:length(xs)
            cand = vcat(xs[1:s-1], xs[min(s + chunk, length(xs) + 1):end])
            ntests += 1
            if test(cand)
                xs = cand; k = max(k - 1, 2); reduced = true
                break
            end
            ntests >= max_tests && break
        end
        if !reduced
            k >= length(xs) && break
            k = min(2k, length(xs))
        end
    end
    xs
end

_has_measurement(c) = any(op -> _opkind(op) == OP_MEASUREMENT, c)

# Whether the pair is refuted by a cheap oracle run (dense if small, else stabilizer).
function _still_fails(input, compiled, n; outcome_map, semantics, magic_k, max_dense_qubits, rng, trials)
    (_has_measurement(input) && isempty(wellformedness_issues(input))) || return false
    N = max(n, _width(compiled))
    om = outcome_map === :infer ? :infer : outcome_map
    if N <= max_dense_qubits
        dv = dense_oracle(input, compiled; n_data = n, semantics, outcome_map = om, magic_k, max_qubits = max_dense_qubits)
        return dv.decided && !dv.agree
    end
    om === :infer && (om = :identity)
    sv = stabilizer_oracle(input, compiled; n_data = n, outcome_map = om, magic_k, trials, rng)
    sv.applicable && sv.refuted
end

"""
    minimize_counterexample(input; outcome_map = :identity, semantics = :outcomes, magic_k = 1,
                            max_dense_qubits = 10, max_tests = 300, rng) -> (input′, compiled′)

Delta-debug a failing *input* circuit: delete ops while, after recompiling the
remainder through PBCCompiler's public `preprocess_circuit`, the compiled program
is still refuted by an oracle. Returns the 1-minimal input found and its compiled
program. Requires the compiler to be deterministic, which it is.
"""
function minimize_counterexample(input::AbstractVector{CircuitOp.Type}; outcome_map = :identity, semantics::Symbol = :outcomes,
                                 magic_k::Int = 1, max_dense_qubits::Int = 10, max_tests::Int = 300,
                                 rng::AbstractRNG = Random.default_rng(), trials::Int = 10)
    n = _width(input)
    function test(sub)
        _width(sub) == n || return false          # keep the data register fixed
        comp = try _compile(sub) catch; return false end
        _still_fails(sub, comp, n; outcome_map, semantics, magic_k, max_dense_qubits, rng, trials)
    end
    m = _ddmin(test, collect(CircuitOp.Type, input); max_tests)
    (m, _compile(m))
end

"""
    minimize_pair(input, compiled; n_data, outcome_map = :identity, kwargs...) -> (input′, compiled′)

Minimize an inequivalent (input, compiled) pair that did not come from the
compiler: delta-debug the set of visible classical bits, dropping each removed
bit's measurements from *both* programs, while the pair stays refuted. A bit read
by a `BitConditional` is never removed.
"""
function minimize_pair(input::AbstractVector{CircuitOp.Type}, compiled::AbstractVector{CircuitOp.Type};
                       n_data::Int = _width(input), outcome_map = :identity, semantics::Symbol = :outcomes,
                       magic_k::Int = 1, max_dense_qubits::Int = 10, max_tests::Int = 300,
                       rng::AbstractRNG = Random.default_rng(), trials::Int = 10)
    read_bits = Set{Int}()
    for c in (input, compiled), op in c
        _opkind(op) == OP_BITCOND && push!(read_bits, _bitcond_parts(op)[2])
    end
    bits = sort!([b for b in keys(_final_writes(input)) if !(b in read_bits)])
    keep_only(c, kept) = [op for op in c if !(_opkind(op) == OP_MEASUREMENT && _measurement_bit(op) in bits &&
                                              !(_measurement_bit(op) in kept))]
    test(kept) = _still_fails(keep_only(input, kept), keep_only(compiled, kept), n_data;
                              outcome_map, semantics, magic_k, max_dense_qubits, rng, trials)
    kept = test(bits) ? _ddmin(test, bits; max_tests) : bits
    (keep_only(input, kept), keep_only(compiled, kept))
end

# Attach minimization and dense confirmation to a counterexample.
function _finish_counterexample(cex::Counterexample, input, compiled, n; compiled_by_us::Bool, minimize::Bool,
                                outcome_map, semantics, magic_k, max_dense_qubits, rng, stats)
    mi, mc = nothing, nothing
    if minimize
        t = time_ns()
        om = outcome_map === :infer ? :infer : outcome_map
        mi, mc = try
            compiled_by_us ? minimize_counterexample(input; outcome_map = om, semantics, magic_k, max_dense_qubits, rng) :
                             minimize_pair(input, compiled; n_data = n, outcome_map = om, semantics, magic_k, max_dense_qubits, rng)
        catch err
            err isa InterruptException && rethrow()
            (nothing, nothing)
        end
        stats[:time_minimization] = (time_ns() - t) / 1e9
    end
    confirmed = cex.confirmed
    note = cex.note
    # re-confirm with the dense oracle whenever the (minimized) pair is small enough
    ci, cc = mi === nothing ? (input, compiled) : (mi, mc)
    nn = _width(ci)
    if max(nn, _width(cc)) <= max_dense_qubits
        dv = dense_oracle(ci, cc; n_data = mi === nothing ? n : nn, semantics,
                          outcome_map = outcome_map === :infer ? :infer : outcome_map, magic_k, max_qubits = max_dense_qubits)
        if dv.decided
            confirmed = !dv.agree
            note *= confirmed ? " [confirmed by the dense oracle]" : " [NOT confirmed: the dense oracle finds agreement]"
        end
    elseif !confirmed
        note *= " [unconfirmed: too large for the dense oracle]"
    end
    Counterexample(cex.method, cex.divergence, cex.input_state, cex.input_stabilizer, cex.branch,
                   cex.prob_input, cex.prob_compiled, mi, mc, confirmed, note)
end
