# Mutation operators on compiled programs, shared by the mutation-testing test items
# and the benchmark scripts. They never edit compiler files: they copy a compiled
# output and change one thing in it. This file defines no test items.

using PBCCompiler: CircuitOp
import QuantumClifford
using QuantumClifford: PauliOperator, nqubits, embed
using Random

const _MV = PBCCompiler.Verification

"""Mutation kinds (spec V9)."""
const MUTATION_KINDS = (:flip_condition_bit, :drop_correction, :flip_sign, :swap_anticommuting,
                        :change_letter, :change_angle, :drop_gadget)

_rebuild_like(op, p, q) = begin
    k = _MV._opkind(op)
    k == _MV.OP_MEASUREMENT ? CircuitOp.Measurement(p, _MV._measurement_bit(op), q) :
    k == _MV.OP_EIGHPI ? CircuitOp.ExpEighPiPauli(p, q) :
    k == _MV.OP_QUATPI ? CircuitOp.ExpQuatPiPauli(p, q) :
    k == _MV.OP_HALFPI ? CircuitOp.ExpHalfPiPauli(p, q) : CircuitOp.Pauli(p, q)
end

_single_pauli_ops(c) = [i for (i, op) in enumerate(c) if (k = _MV._opkind(op); k == _MV.OP_MEASUREMENT || _MV._isrotation(k) ||
                        (k == _MV.OP_BITCOND && (ik = _MV._opkind(_MV._bitcond_parts(op)[1]); ik == _MV.OP_MEASUREMENT || _MV._isrotation(ik))))]

function _replace_pauli(op, f)
    if _MV._opkind(op) == _MV.OP_BITCOND
        inner, b = _MV._bitcond_parts(op)
        p, q = _MV._pauli_and_qubits(inner)
        return CircuitOp.BitConditional(_rebuild_like(inner, f(copy(p)), q), b)
    end
    p, q = _MV._pauli_and_qubits(op)
    _rebuild_like(op, f(copy(p)), q)
end

"""
    mutate(rng, c, kind; n_data = width of c) -> (mutated, applied::Bool)

Apply one mutation of the given kind to a copy of the compiled program `c`.
`applied` is false when the program has nothing to mutate for that kind.
`n_data` (qubits above it are magic qubits) is used by `:drop_gadget`.
"""
function mutate(rng::AbstractRNG, c, kind::Symbol; n_data::Int = _MV._width(c))
    c = copy(c)
    bcs = [i for (i, op) in enumerate(c) if _MV._opkind(op) == _MV.OP_BITCOND]
    if kind == :flip_condition_bit
        isempty(bcs) && return (c, false)
        i = rand(rng, bcs); inner, b = _MV._bitcond_parts(c[i])
        written = sort!(unique([_MV._measurement_bit(c[j]) for j in 1:i-1 if _MV._opkind(c[j]) == _MV.OP_MEASUREMENT]))
        others = filter(!=(b), written)
        isempty(others) && return (c, false)
        c[i] = CircuitOp.BitConditional(inner, rand(rng, others))
    elseif kind == :drop_correction
        isempty(bcs) && return (c, false)
        deleteat!(c, rand(rng, bcs))
    elseif kind == :flip_sign
        idx = _single_pauli_ops(c); isempty(idx) && return (c, false)
        i = rand(rng, idx)
        c[i] = _replace_pauli(c[i], p -> (p.phase[] = (p.phase[] + 0x02) & 0x03; p))
    elseif kind == :change_letter
        idx = _single_pauli_ops(c); isempty(idx) && return (c, false)
        i = rand(rng, idx)
        c[i] = _replace_pauli(c[i], p -> begin
            j = rand(rng, 1:nqubits(p)); x, z = p[j]
            new = rand(rng, setdiff([(true, false), (false, true), (true, true)], [(x, z)]))
            p[j] = new; p
        end)
    elseif kind == :change_angle
        idx = [i for (i, op) in enumerate(c) if (k = _MV._opkind(op); _MV._isrotation(k) ||
               (k == _MV.OP_BITCOND && _MV._isrotation(_MV._opkind(_MV._bitcond_parts(op)[1]))))]
        isempty(idx) && return (c, false)
        i = rand(rng, idx)
        op = c[i]
        cond = _MV._opkind(op) == _MV.OP_BITCOND
        inner, b = cond ? _MV._bitcond_parts(op) : (op, 0)
        p, q = _MV._pauli_and_qubits(inner)
        k = _MV._opkind(inner)
        choices = filter(!=(k), [_MV.OP_HALFPI, _MV.OP_QUATPI, _MV.OP_EIGHPI])
        nk = rand(rng, choices)
        ni = nk == _MV.OP_HALFPI ? CircuitOp.ExpHalfPiPauli(p, q) : nk == _MV.OP_QUATPI ? CircuitOp.ExpQuatPiPauli(p, q) :
             CircuitOp.ExpEighPiPauli(p, q)
        c[i] = cond ? CircuitOp.BitConditional(ni, b) : ni
    elseif kind == :swap_anticommuting
        meas = [i for (i, op) in enumerate(c) if _MV._opkind(op) == _MV.OP_MEASUREMENT]
        pairs = Tuple{Int,Int}[]
        n = _MV._width(c)
        for x in meas, y in meas
            x < y || continue
            px = embed(n, reverse(_MV._pauli_and_qubits(c[x]))...); py = embed(n, reverse(_MV._pauli_and_qubits(c[y]))...)
            QuantumClifford.comm(px, py) == 0x01 && push!(pairs, (x, y))
        end
        isempty(pairs) && return (c, false)
        x, y = rand(rng, pairs)
        c[x], c[y] = c[y], c[x]
    elseif kind == :drop_gadget
        gs = PBCCompiler.GadgetShapes.find_gadgets(c, n_data)
        isempty(gs) && return (c, false)
        deleteat!(c, collect(rand(rng, gs).ops))
    else
        throw(ArgumentError("unknown mutation $kind"))
    end
    (c, true)
end
