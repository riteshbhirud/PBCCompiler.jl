"""
    WellformednessIssue

One violated rule found by [`wellformedness_issues`](@ref): the op index (in
program order, 1-based), a rule identifier, and a human-readable explanation.
"""
struct WellformednessIssue
    index::Int
    rule::Symbol
    message::String
end

Base.show(io::IO, w::WellformednessIssue) = print(io, "op ", w.index, " [", w.rule, "]: ", w.message)

"""
    wellformedness_issues(c::Circuit; nqubits=nothing, n_data=nothing) -> Vector{WellformednessIssue}

Run the static contract checks on a circuit and return every violation, in
program order. See [`check_wellformed`](@ref) for the rules.

- `nqubits`: if given, qubit indices must lie in `1:nqubits`.
- `n_data`: if given, qubits with index `> n_data` are magic qubits and the
  magic-qubit discipline is enforced.
"""
function wellformedness_issues(c::AbstractVector{CircuitOp.Type}; nqubits::Union{Nothing,Int}=nothing,
                               n_data::Union{Nothing,Int}=nothing)
    issues = WellformednessIssue[]
    written = Dict{Int,Int}()          # bit => op index of first write
    for (i, op) in enumerate(c)
        _check_op!(issues, i, op, nqubits, written; toplevel=true)
    end
    n_data === nothing || _check_magic!(issues, c, n_data)
    sort!(issues; by = w -> w.index)
    issues
end

"""
    check_wellformed(c::Circuit; nqubits=nothing, n_data=nothing) -> Vector{String}

Static contract checks on a circuit. Returns one message per violation (empty if
the circuit is well formed). Each message starts with the offending op index.

Rules:
- every qubit index is `≥ 1` (and `≤ nqubits` if given), each Pauli has one letter
  per listed qubit, and no op lists the same qubit twice;
- every measured or rotated Pauli is Hermitian (phase ±1, not ±i);
- `PauliConditional` control and target qubits are disjoint;
- every `BitConditional` reads a bit written by an earlier measurement, with a valid
  bit index (no read-before-write, no dangling bit);
- no classical bit is written by two measurements;
- with `n_data` given, each magic qubit (index `> n_data`) is used only by
  measurements: first by one joint measurement with data qubits, then by one
  single-qubit measurement that consumes it, and never again after that.
"""
check_wellformed(c::AbstractVector{CircuitOp.Type}; kwargs...) = map(string, wellformedness_issues(c; kwargs...))

function _check_qubits!(issues, i, qubits::Vector{Int}, nq, what)
    for q in qubits
        if q < 1 || (nq !== nothing && q > nq)
            push!(issues, WellformednessIssue(i, :qubit_range,
                "$what acts on qubit $q, outside the valid range 1:$(nq === nothing ? "∞" : nq)"))
        end
    end
    if length(unique(qubits)) != length(qubits)
        dup = [q for q in unique(qubits) if count(==(q), qubits) > 1]
        push!(issues, WellformednessIssue(i, :duplicate_qubit, "$what lists qubit(s) $dup more than once"))
    end
end

function _check_pauli!(issues, i, p::PauliOperator, qubits::Vector{Int}, what; hermitian::Bool=true)
    if nqubits(p) != length(qubits)
        push!(issues, WellformednessIssue(i, :pauli_length,
            "$what has a $(nqubits(p))-qubit Pauli but lists $(length(qubits)) qubits"))
    end
    if hermitian && !_ishermitian(p)
        push!(issues, WellformednessIssue(i, :non_hermitian,
            "$what uses the non-Hermitian Pauli $p (phase ±i); measurements and rotations need phase ±1"))
    end
end

function _check_op!(issues, i, op, nq, written; toplevel::Bool)
    k = _opkind(op)
    if k == OP_MEASUREMENT
        p, q = _pauli_and_qubits(op)
        b = _measurement_bit(op)
        _check_qubits!(issues, i, q, nq, "Measurement")
        _check_pauli!(issues, i, p, q, "Measurement")
        if b < 1
            push!(issues, WellformednessIssue(i, :bit_range, "Measurement writes invalid bit index $b"))
        elseif haskey(written, b)
            push!(issues, WellformednessIssue(i, :bit_rewrite,
                "Measurement writes bit $b, which was already written by op $(written[b])"))
        elseif toplevel
            written[b] = i
        end
        if !toplevel
            push!(issues, WellformednessIssue(i, :conditional_measurement,
                "a Measurement nested inside a BitConditional has no defined write semantics"))
        end
    elseif _isrotation(k)
        p, q = _pauli_and_qubits(op)
        _check_qubits!(issues, i, q, nq, _kindname(k))
        _check_pauli!(issues, i, p, q, _kindname(k))
    elseif k == OP_PAULICOND
        cp, cq, tp, tq = _paulicond_parts(op)
        _check_qubits!(issues, i, cq, nq, "PauliConditional control")
        _check_qubits!(issues, i, tq, nq, "PauliConditional target")
        _check_pauli!(issues, i, cp, cq, "PauliConditional control")
        _check_pauli!(issues, i, tp, tq, "PauliConditional target")
        ov = intersect(cq, tq)
        isempty(ov) || push!(issues, WellformednessIssue(i, :control_target_overlap,
            "PauliConditional control and target share qubit(s) $ov"))
        (isempty(cq) || isempty(tq)) && push!(issues, WellformednessIssue(i, :empty_pauli,
            "PauliConditional has an empty control or target"))
    elseif k == OP_BITCOND
        inner, b = _bitcond_parts(op)
        if b < 1
            push!(issues, WellformednessIssue(i, :bit_range, "BitConditional reads invalid bit index $b"))
        elseif !haskey(written, b)
            push!(issues, WellformednessIssue(i, :read_before_write,
                "BitConditional reads bit $b before any measurement writes it"))
        end
        _check_op!(issues, i, inner, nq, written; toplevel=false)
    else # OP_PREPMAGIC
        q, qs = _prepmagic_parts(op)
        _check_qubits!(issues, i, vcat([q], qs), nq, "PrepMagic")
    end
    nothing
end

function _check_magic!(issues, c, n_data::Int)
    # magic qubit => [(op index, role)] in program order
    uses = Dict{Int,Vector{Tuple{Int,Symbol}}}()
    for (i, op) in enumerate(c)
        qs = _support(op)
        magic = filter(>(n_data), qs)
        isempty(magic) && continue
        k = _opkind(op)
        if k != OP_MEASUREMENT
            for a in magic
                push!(issues, WellformednessIssue(i, :magic_nonmeasurement,
                    "magic qubit $a is used by a non-measurement op ($(_kindname(k))); magic qubits may only be consumed by measurements"))
                push!(get!(uses, a, Tuple{Int,Symbol}[]), (i, :other))
            end
            continue
        end
        if length(magic) > 1
            push!(issues, WellformednessIssue(i, :magic_multiple,
                "measurement touches several magic qubits $magic"))
        end
        role = length(qs) == 1 ? :single : :joint
        for a in magic
            push!(get!(uses, a, Tuple{Int,Symbol}[]), (i, role))
        end
    end
    for a in sort!(collect(keys(uses)))
        u = uses[a]
        first_i, first_role = u[1]
        if first_role != :joint
            push!(issues, WellformednessIssue(first_i, :magic_not_prepared,
                "magic qubit $a is first used by a $(first_role == :single ? "single-qubit measurement" : "non-measurement op") instead of a joint measurement with data qubits (it is consumed before being used as a resource)"))
        end
        consumed = findfirst(x -> x[2] == :single, u)
        if consumed === nothing
            push!(issues, WellformednessIssue(u[end][1], :magic_not_consumed,
                "magic qubit $a is never consumed by a single-qubit measurement before the end of the program"))
        else
            for (j, _) in u[consumed+1:end]
                push!(issues, WellformednessIssue(j, :magic_reuse,
                    "magic qubit $a is used again after being consumed at op $(u[consumed][1])"))
            end
            njoint = count(x -> x[2] == :joint, u[1:consumed])
            if njoint > 1
                push!(issues, WellformednessIssue(u[consumed][1], :magic_multiple_joint,
                    "magic qubit $a takes part in $njoint joint measurements before being consumed (expected 1)"))
            end
        end
    end
    nothing
end
