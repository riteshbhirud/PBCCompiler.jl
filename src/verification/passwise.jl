# V8: per-pass translation validation.
#
# PBCCompiler's passes are individually callable functions. `preprocess_circuit`
# applies them in the order below, after returning early for circuits with fewer
# than two ops. The order is duplicated here because the compiler does not expose
# it as data (hook request H1 in docs/NOTES_FOR_COMPILER.md). A test checks that
# running these passes in sequence reproduces `preprocess_circuit` exactly, so a
# change in the compiler's pipeline is noticed.

"""The compiler passes, in the order `preprocess_circuit` runs them."""
const PREPROCESS_PASSES = (
    :remove_pauliconditional,
    :group_nonclifford,
    :merge_ops,
    :remove_clifford,
    :remove_nonclifford,
    :remove_post_measurement,
)

"""
    PassStep

One step of [`verify_passes`](@ref): the pass name, the program before and after
it, and the verification result of `after` against `before`.
"""
struct PassStep
    pass::Symbol
    before::Vector{CircuitOp.Type}
    after::Vector{CircuitOp.Type}
    result::VerificationResult
end

"""
    run_passes(input) -> Vector{Tuple{Symbol,Vector{CircuitOp.Type}}}

Run the compiler's passes one at a time on a copy of `input`, recording the
program after each one. Mirrors `preprocess_circuit`, including its early return
for circuits with fewer than two ops.
"""
function run_passes(input::AbstractVector{CircuitOp.Type})
    c = collect(CircuitOp.Type, input)
    out = Tuple{Symbol,Vector{CircuitOp.Type}}[]
    length(c) < 2 && return out
    for name in PREPROCESS_PASSES
        getfield(PBCCompiler, name)(c)
        push!(out, (name, copy(c)))
    end
    out
end

"""
    verify_passes(input; kwargs...) -> Vector{VerificationResult}

Translation validation of each compiler pass: run the passes one at a time and
verify the program after each pass against the program before it. The first
result that is not `EQUIVALENT` localizes a problem to that pass. Results come in
pass order, see [`PREPROCESS_PASSES`](@ref). Use [`verify_passes_detailed`](@ref)
to also get the intermediate programs.

The outcome record is always the original input's bits (gadget outcomes
introduced along the way are hidden), and the outcome map is the identity unless
`outcome_map` is given. Other keyword arguments are passed to
[`verify_compilation`](@ref) (`minimize` defaults to `false` here).
"""
verify_passes(input::AbstractVector{CircuitOp.Type}; kwargs...) = [s.result for s in verify_passes_detailed(input; kwargs...)]

"""
    verify_passes_detailed(input; kwargs...) -> Vector{PassStep}

Like [`verify_passes`](@ref), but also returns the program before and after each pass.
"""
function verify_passes_detailed(input::AbstractVector{CircuitOp.Type}; outcome_map = :identity, minimize::Bool = false, kwargs...)
    n = _width(input)
    visible = sort!(collect(keys(_final_writes(input))))
    steps = PassStep[]
    before = collect(CircuitOp.Type, input)
    trace = try
        run_passes(input)
    catch err
        stats = Dict{Symbol,Any}(:reasons => ["a compiler pass threw an error: " * sprint(showerror, err)])
        push!(steps, PassStep(:error, before, CircuitOp.Type[], VerificationResult(UNKNOWN, :none, nothing, nothing, Int[], String[], stats)))
        return steps
    end
    for (name, after) in trace
        r = verify_compilation(before, after; n_data = n, outcome_map, minimize, visible_bits = visible, kwargs...)
        r.stats[:pass] = name
        push!(steps, PassStep(name, before, after, r))
        before = after
    end
    steps
end
