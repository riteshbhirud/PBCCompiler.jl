"""
    VerificationStatus

`EQUIVALENT` (proved), `NOT_EQUIVALENT` (a concrete witness was found), or
`UNKNOWN` (neither; the reasons are in the result's `stats[:reasons]`).
"""
@enum VerificationStatus EQUIVALENT NOT_EQUIVALENT UNKNOWN

"""
    Certificate

Evidence for an `EQUIVALENT` verdict:

- `semantics` (`:outcomes` or `:channel`);
- the canonical forms of the input and compiled programs;
- the outcome map from each input bit to an affine function of compiled outcome variables. Each
  variable is the index of a compiled measurement op, and `variable_bits` gives
  its classical bit;
- whether that map is the identity on bits;
- the proved gadgets that were lifted.
"""
struct Certificate
    semantics::Symbol
    input_form::CanonicalForm
    compiled_form::CanonicalForm
    outcome_map::Dict{Int,AffineBit}
    variable_bits::Dict{Int,Int}
    identity_map::Bool
    lifted_gadgets::Vector{LiftedGadget}
end

"""
    Counterexample

A concrete witness for `NOT_EQUIVALENT`:

- `method`: the oracle that found it;
- `divergence`: the first diverging canonical position, with the input and compiled
  items there (when the normal forms were built);
- the input state, as a dense data-qubit vector (`input_state`) or as a stabilizer
  state (`input_stabilizer`);
- `branch`: the visible outcome record (input bit to value), and the probability of
  that record under the input (`prob_input`) and the compiled program
  (`prob_compiled`);
- the minimized failing input and compiled programs (when minimization ran);
- `confirmed`: whether the dense oracle re-confirmed the difference;
- `note`: a human-readable explanation.
"""
struct Counterexample
    method::Symbol
    divergence::Union{Nothing,Tuple{Int,Union{Nothing,NFItem},Union{Nothing,NFItem}}}
    input_state::Vector{ComplexF64}
    input_stabilizer::Union{Nothing,QuantumClifford.Stabilizer}
    branch::Dict{Int,Bool}
    prob_input::Float64
    prob_compiled::Float64
    minimized_input::Union{Nothing,Vector{CircuitOp.Type}}
    minimized_compiled::Union{Nothing,Vector{CircuitOp.Type}}
    confirmed::Bool
    note::String
end

"""
    VerificationResult

Result of [`verify_compilation`](@ref).

- `status`: a [`VerificationStatus`](@ref).
- `method`: what decided it (`:normal_form`, `:dense_oracle`, `:stabilizer_oracle`,
  `:wellformedness`, `:none`).
- `certificate` (`EQUIVALENT` only) and `counterexample` (`NOT_EQUIVALENT` only).
- `unverified_components`: compiled op indices that could not be proved (failed
  gadgets, unsupported ops).
- `wellformedness_errors`: messages from [`check_wellformed`](@ref) on both programs.
- `stats`: timings, sizes, branch counts, and `:reasons` (why the result is not
  `EQUIVALENT`, when it is not).
"""
struct VerificationResult
    status::VerificationStatus
    method::Symbol
    certificate::Union{Nothing,Certificate}
    counterexample::Union{Nothing,Counterexample}
    unverified_components::Vector{Int}
    wellformedness_errors::Vector{String}
    stats::Dict{Symbol,Any}
end
