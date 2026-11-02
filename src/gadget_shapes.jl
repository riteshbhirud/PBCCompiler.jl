"""
    PBCCompiler.GadgetShapes

Find magic-state gadgets in a compiled PBC program by their shape. Used by the
verifier and the resource estimator. The recognized shapes are listed in
[`GADGET_TEMPLATES`](@ref).
"""
module GadgetShapes

using ..PBCCompiler: CircuitOp
using Moshi.Match: @match
using QuantumClifford: PauliOperator, nqubits

export GadgetTemplate, GADGET_TEMPLATES, GadgetMatch, find_gadgets, count_gadgets

"""
    GadgetTemplate(name, joint_magic_letter, magic_measure_letter, max_corrections, description)

The shape of a magic-state gadget, made of consecutive ops:

1. a joint `Measurement` on data qubits and exactly one magic qubit `a`, with letter `joint_magic_letter` (`'X'`, `'Y'` or `'Z'`) on `a`;
2. a single-qubit `Measurement` of `a` with letter `magic_measure_letter`;
3. up to `max_corrections` `BitConditional` rotations on data qubits that read one of the two outcome bits.

A template only describes the shape; it says nothing about whether the gadget is correct.
"""
struct GadgetTemplate
    name::Symbol
    joint_magic_letter::Char
    magic_measure_letter::Char
    max_corrections::Int
    description::String
end

"""
    GADGET_TEMPLATES

The gadget shapes recognized by default. It contains the shape produced by
PBCCompiler's `gadgetize`: measure `P⊗Z` with the magic qubit, measure the magic
qubit in `X`, then apply `P_{π/4}` if the first outcome is 1 and `P_{π/2}` if the
second is 1.
"""
const GADGET_TEMPLATES = GadgetTemplate[
    GadgetTemplate(:pi8_autocorrected, 'Z', 'X', 2,
                   "π/8 rotation by magic-state injection with π/4 and π/2 corrections (PBCCompiler gadgetize)"),
]

"""
    GadgetMatch

One recognized gadget: the template name, the op indices it spans (`ops`), the
magic qubit, the indices of the joint measurement, the magic measurement and the
corrections, and the two outcome bits (`joint_bit`, `magic_bit`).
"""
struct GadgetMatch
    template::Symbol
    ops::UnitRange{Int}
    magic_qubit::Int
    joint_index::Int
    magic_index::Int
    correction_indices::Vector{Int}
    joint_bit::Int
    magic_bit::Int
end

function _letter(p::PauliOperator, pos::Int)
    x, z = p[pos]
    x ? (z ? 'Y' : 'X') : (z ? 'Z' : 'I')
end

_measurement_fields(op) = @match op begin
    CircuitOp.Measurement(p, b, q) => (p, b, q)
    _ => nothing
end

_bitcond_fields(op) = @match op begin
    CircuitOp.BitConditional(inner, b) => (inner, b)
    _ => nothing
end

_rotation_qubits(op) = @match op begin
    CircuitOp.Pauli(_, q) => q
    CircuitOp.ExpHalfPiPauli(_, q) => q
    CircuitOp.ExpQuatPiPauli(_, q) => q
    CircuitOp.ExpEighPiPauli(_, q) => q
    _ => nothing
end

function _match_at(c, i::Int, n_data::Int, t::GadgetTemplate)
    m = _measurement_fields(c[i])
    m === nothing && return nothing
    p, b1, q = m
    nqubits(p) == length(q) || return nothing
    magic_pos = findall(>(n_data), q)
    length(magic_pos) == 1 || return nothing
    data_pos = findall(<=(n_data), q)
    any(j -> _letter(p, j) != 'I', data_pos) || return nothing
    _letter(p, magic_pos[1]) == t.joint_magic_letter || return nothing
    a = q[magic_pos[1]]
    i + 1 <= length(c) || return nothing
    m2 = _measurement_fields(c[i + 1])
    m2 === nothing && return nothing
    p2, b2, q2 = m2
    (q2 == [a] && nqubits(p2) == 1 && _letter(p2, 1) == t.magic_measure_letter) || return nothing
    corr = Int[]
    j = i + 2
    while j <= length(c) && length(corr) < t.max_corrections
        bc = _bitcond_fields(c[j])
        bc === nothing && break
        inner, b = bc
        (b == b1 || b == b2) || break
        rq = _rotation_qubits(inner)
        (rq === nothing || any(>(n_data), rq)) && break
        push!(corr, j)
        j += 1
    end
    GadgetMatch(t.name, i:j-1, a, i, i + 1, corr, b1, b2)
end

"""
    find_gadgets(c, n_data::Int; templates = GADGET_TEMPLATES) -> Vector{GadgetMatch}

Return the non-overlapping gadgets in a compiled program that match one of
`templates`, in program order. Qubits above `n_data` are magic qubits. When
several templates match at the same position, the first one is used.
"""
function find_gadgets(c::AbstractVector{CircuitOp.Type}, n_data::Int; templates = GADGET_TEMPLATES)
    out = GadgetMatch[]
    i = 1
    while i <= length(c)
        found = nothing
        for t in templates
            found = _match_at(c, i, n_data, t)
            found === nothing || break
        end
        if found === nothing
            i += 1
        else
            push!(out, found)
            i = last(found.ops) + 1
        end
    end
    out
end

"""
    count_gadgets(c, n_data::Int; templates = GADGET_TEMPLATES) -> Dict{Symbol,Int}

Count the gadgets of each template in a compiled program, which is the number of
magic states it consumes by gadget type.
"""
function count_gadgets(c::AbstractVector{CircuitOp.Type}, n_data::Int; templates = GADGET_TEMPLATES)
    d = Dict{Symbol,Int}()
    for g in find_gadgets(c, n_data; templates)
        d[g.template] = get(d, g.template, 0) + 1
    end
    d
end

end # module GadgetShapes
