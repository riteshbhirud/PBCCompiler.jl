const PauliT = PauliOperator{Array{UInt8,0},Vector{UInt64}}

"""Kinds of `CircuitOp` as seen by the verifier."""
@enum OpKind OP_MEASUREMENT OP_PAULI OP_HALFPI OP_QUATPI OP_EIGHPI OP_PREPMAGIC OP_PAULICOND OP_BITCOND

function _opkind(op::CircuitOp.Type)::OpKind
    @match op begin
        CircuitOp.Measurement(_, _, _) => OP_MEASUREMENT
        CircuitOp.Pauli(_, _) => OP_PAULI
        CircuitOp.ExpHalfPiPauli(_, _) => OP_HALFPI
        CircuitOp.ExpQuatPiPauli(_, _) => OP_QUATPI
        CircuitOp.ExpEighPiPauli(_, _) => OP_EIGHPI
        CircuitOp.PrepMagic(_, _) => OP_PREPMAGIC
        CircuitOp.PauliConditional(_, _, _, _) => OP_PAULICOND
        CircuitOp.BitConditional(_, _) => OP_BITCOND
    end
end

const _KIND_NAMES = Dict(OP_MEASUREMENT => "Measurement", OP_PAULI => "Pauli", OP_HALFPI => "ExpHalfPiPauli",
    OP_QUATPI => "ExpQuatPiPauli", OP_EIGHPI => "ExpEighPiPauli", OP_PREPMAGIC => "PrepMagic",
    OP_PAULICOND => "PauliConditional", OP_BITCOND => "BitConditional")
_kindname(k::OpKind)::String = _KIND_NAMES[k]

_isrotation(k::OpKind) = k == OP_PAULI || k == OP_HALFPI || k == OP_QUATPI || k == OP_EIGHPI

# angle θ in units of π, for exp(-iθπ P)
function _rotation_angle(k::OpKind)::Rational{Int}
    k == OP_PAULI && return 1 // 2   # P ∝ exp(-iπ/2 P)
    k == OP_HALFPI && return 1 // 2
    k == OP_QUATPI && return 1 // 4
    k == OP_EIGHPI && return 1 // 8
    throw(ArgumentError("not a rotation kind: $k"))
end

function _pauli_and_qubits(op::CircuitOp.Type)::Tuple{PauliT,Vector{Int}}
    @match op begin
        CircuitOp.Measurement(p, _, q) => (p, q)
        CircuitOp.Pauli(p, q) => (p, q)
        CircuitOp.ExpHalfPiPauli(p, q) => (p, q)
        CircuitOp.ExpQuatPiPauli(p, q) => (p, q)
        CircuitOp.ExpEighPiPauli(p, q) => (p, q)
        _ => throw(ArgumentError("op has no single (pauli, qubits) pair: $op"))
    end
end

_measurement_bit(op::CircuitOp.Type)::Int = @match op begin
    CircuitOp.Measurement(_, b, _) => b
    _ => throw(ArgumentError("not a Measurement: $op"))
end

_bitcond_parts(op::CircuitOp.Type)::Tuple{CircuitOp.Type,Int} = @match op begin
    CircuitOp.BitConditional(inner, b) => (inner, b)
    _ => throw(ArgumentError("not a BitConditional: $op"))
end

_paulicond_parts(op::CircuitOp.Type)::Tuple{PauliT,Vector{Int},PauliT,Vector{Int}} = @match op begin
    CircuitOp.PauliConditional(cp, cq, tp, tq) => (cp, cq, tp, tq)
    _ => throw(ArgumentError("not a PauliConditional: $op"))
end

_prepmagic_parts(op::CircuitOp.Type)::Tuple{Int,Vector{Int}} = @match op begin
    CircuitOp.PrepMagic(q, qs) => (q, qs)
    _ => throw(ArgumentError("not a PrepMagic: $op"))
end

# may contain duplicates
function _raw_qubits(op::CircuitOp.Type)::Vector{Int}
    k = _opkind(op)
    if k == OP_MEASUREMENT || _isrotation(k)
        return copy(_pauli_and_qubits(op)[2])
    elseif k == OP_PAULICOND
        _, cq, _, tq = _paulicond_parts(op)
        return vcat(cq, tq)
    elseif k == OP_BITCOND
        return _raw_qubits(_bitcond_parts(op)[1])
    else # OP_PREPMAGIC
        q, qs = _prepmagic_parts(op)
        return vcat([q], qs)
    end
end

_support(op::CircuitOp.Type) = sort!(unique(_raw_qubits(op)))

function _width(c::AbstractVector{CircuitOp.Type})
    w = 0
    for op in c
        q = _raw_qubits(op)
        isempty(q) || (w = max(w, maximum(q)))
    end
    w
end

_ishermitian(p::PauliOperator) = iseven(p.phase[])

_isidentity(p::PauliOperator) = all(iszero, p.xz)
