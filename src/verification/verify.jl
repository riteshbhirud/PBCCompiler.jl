"""
    verify_compilation(input::Circuit, compiled::Circuit; kwargs...) -> VerificationResult
    verify_compilation(input::Circuit; kwargs...) -> VerificationResult

Check that `compiled` implements the same logical operation as `input`, and
return a [`VerificationResult`](@ref) whose status is `EQUIVALENT`,
`NOT_EQUIVALENT` or `UNKNOWN`. The one-argument form compiles a copy of `input`
with `preprocess_circuit` first, and defaults to `outcome_map = :identity`, since
the compiler is expected to keep every classical bit at its index.

`EQUIVALENT` is only returned when the two programs are proved equivalent.
`NOT_EQUIVALENT` comes with a [`Counterexample`](@ref). The numerical oracles can
only refute, never certify, so when they find no difference the result is `UNKNOWN`.

Keyword arguments:

- `semantics = :outcomes`: compare the distribution of classical outcomes for every input state; `:channel` also compares the final quantum state.
- `outcome_map = :infer`: how input bits relate to compiled outcomes: `:infer` (any parity map), `:identity` (same bit index), or a `Dict` from input bits to `AffineBit`s over compiled bits.
- `oracle = :auto`: which fallback oracle to use: `:auto`, `:dense`, `:stabilizer` or `:none`.
- `max_dense_qubits = 10`, `max_branches = 4096`: size limits for the dense oracle.
- `stabilizer_trials = 20`: number of random trials for the stabilizer oracle.
- `rng = Random.default_rng()`: random number generator used by the oracles.
- `magic_k = 1`: magic qubits are prepared in `(|0⟩ + exp(iπ magic_k/4)|1⟩)/√2`; 1 means `T|+⟩`.
- `n_data`: number of data qubits (default: the width of `input`); higher qubits are magic qubits.
- `visible_bits`: the input bits to compare (default: every bit the input writes).
- `minimize = true`: shrink a counterexample with [`minimize_counterexample`](@ref) or [`minimize_pair`](@ref) when one is found.
"""
function verify_compilation(input::AbstractVector{CircuitOp.Type}, compiled::AbstractVector{CircuitOp.Type};
                            outcome_map = :infer, oracle::Symbol = :auto, max_dense_qubits::Int = 10,
                            max_branches::Int = 4096, rng::AbstractRNG = Random.default_rng(),
                            semantics::Symbol = :outcomes, magic_k::Int = 1, n_data::Union{Nothing,Int} = nothing,
                            stabilizer_trials::Int = 20, minimize::Bool = true,
                            visible_bits::Union{Nothing,AbstractVector{Int}} = nothing, _compiled_by_us::Bool = false)
    semantics in (:outcomes, :channel) || throw(ArgumentError("semantics must be :outcomes or :channel"))
    oracle in (:auto, :dense, :stabilizer, :none) || throw(ArgumentError("oracle must be :auto, :dense, :stabilizer or :none"))
    t0 = time_ns()
    stats = Dict{Symbol,Any}()
    reasons = String[]
    n = n_data === nothing ? _width(input) : n_data
    stats[:n_data] = n
    stats[:n_qubits_compiled] = max(n, _width(compiled))
    stats[:input_ops] = length(input)
    stats[:compiled_ops] = length(compiled)

    wf = vcat(["input: " * s for s in check_wellformed(input; n_data = n)],
              ["compiled: " * s for s in check_wellformed(compiled; n_data = n)])

    # The input may contain gadgets too, e.g. an intermediate program in verify_passes.
    lr_in = lift_gadgets(input, n; magic_k)
    prog_in = lower_lifted(input, lr_in)
    for f in lr_in.failed
        push!(reasons, "input gadget at ops $(f.match.ops) not proved: $(f.reason)")
    end
    t1 = time_ns()
    lr = lift_gadgets(compiled, n; magic_k)
    prog_out = lower_lifted(compiled, lr)
    stats[:gadgets_lifted] = length(lr.lifted)
    stats[:gadgets_failed] = length(lr.failed)
    stats[:gadget_cache_hits] = lr.cache_hits
    stats[:time_lifting] = (time_ns() - t1) / 1e9
    unverified = Int[]
    for f in lr.failed
        append!(unverified, f.match.ops)
        push!(reasons, "gadget at ops $(f.match.ops) not proved: $(f.reason)")
    end
    for (i, r) in prog_in.unsupported
        push!(reasons, "input op $i unsupported: $r")
    end
    for (i, r) in prog_out.unsupported
        push!(unverified, i)
        push!(reasons, "compiled op $i unsupported: $r")
    end
    sort!(unique!(unverified))

    divergence = nothing
    if isempty(prog_in.unsupported) && isempty(prog_out.unsupported)
        t2 = time_ns()
        cf_in = canonicalize(normal_form(prog_in); semantics)
        cf_out = canonicalize(normal_form(prog_out); semantics)
        cmp = compare_forms(cf_in, cf_out; outcome_map, visible_bits)
        stats[:time_normal_form] = (time_ns() - t2) / 1e9
        stats[:items_input] = length(cf_in.nf.items)
        stats[:items_compiled] = length(cf_out.nf.items)
        stats[:merges] = cf_in.merges + cf_out.merges
        if cmp.equal
            if isempty(wf)
                ident = all(e == affine_var(get(cf_out.nf.final_var, b, -1)) for (b, e) in cmp.outcome_map)
                cert = Certificate(semantics, cf_in, cf_out, cmp.outcome_map, copy(cf_out.nf.var_bit), ident, lr.lifted)
                stats[:time_total] = (time_ns() - t0) / 1e9
                stats[:reasons] = reasons
                return VerificationResult(EQUIVALENT, :normal_form, cert, nothing, unverified, wf, stats)
            end
            push!(reasons, "normal forms match (under register-overwrite semantics), but the programs violate the IR contract; not certified")
            stats[:time_total] = (time_ns() - t0) / 1e9
            stats[:reasons] = reasons
            return VerificationResult(UNKNOWN, :wellformedness, nothing, nothing, unverified, wf, stats)
        end
        push!(reasons, "normal forms differ" * (cmp.position > 0 ? " at canonical position $(cmp.position)" : "") * ": " * cmp.reason)
        k = cmp.position
        if k > 0
            xi = k <= length(cf_in.nf.items) ? cf_in.nf.items[k] : nothing
            yi = k <= length(cf_out.nf.items) ? cf_out.nf.items[k] : nothing
            divergence = (k, xi, yi)
        end
    end

    res = _oracle_fallback(input, compiled, n; oracle, max_dense_qubits, max_branches, rng, semantics,
                           magic_k, outcome_map, stabilizer_trials, divergence, stats, reasons, visible_bits)
    status, method, cex = res
    if status == NOT_EQUIVALENT
        cex isa Counterexample || error("internal: NOT_EQUIVALENT without a counterexample")
        cex = _finish_counterexample(cex, input, compiled, n; compiled_by_us = _compiled_by_us, minimize,
                                     outcome_map, semantics, magic_k, max_dense_qubits, rng, stats)
    end
    stats[:time_total] = (time_ns() - t0) / 1e9
    stats[:reasons] = reasons
    VerificationResult(status, method, nothing, cex, unverified, wf, stats)
end

function verify_compilation(input::AbstractVector{CircuitOp.Type}; outcome_map = :identity, kwargs...)
    compiled = try
        _compile(input)
    catch err
        stats = Dict{Symbol,Any}(:reasons => ["the compiler threw an error: " * sprint(showerror, err)])
        return VerificationResult(UNKNOWN, :none, nothing, nothing, Int[], String[], stats)
    end
    verify_compilation(input, compiled; outcome_map, _compiled_by_us = true, kwargs...)
end

_compile(input) = PBCCompiler.preprocess_circuit(copy(collect(CircuitOp.Type, input)))

function Base.show(io::IO, ::MIME"text/plain", r::VerificationResult)
    println(io, "VerificationResult: ", r.status, " (method: ", r.method, ")")
    if r.certificate !== nothing
        c = r.certificate
        println(io, "  semantics: ", c.semantics, ", gadgets lifted: ", length(c.lifted_gadgets),
                ", identity outcome map: ", c.identity_map)
        c.identity_map || for b in sort!(collect(keys(c.outcome_map)))
            println(io, "    input bit ", b, " = ", c.outcome_map[b])
        end
    end
    if r.counterexample !== nothing
        cx = r.counterexample
        println(io, "  counterexample (", cx.method, ", ", cx.confirmed ? "confirmed" : "unconfirmed", "): ", cx.note)
    end
    isempty(r.unverified_components) || println(io, "  unverified compiled ops: ", r.unverified_components)
    for e in r.wellformedness_errors
        println(io, "  well-formedness: ", e)
    end
    for s in get(r.stats, :reasons, String[])
        println(io, "  reason: ", s)
    end
end
