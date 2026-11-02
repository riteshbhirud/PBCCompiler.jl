# V3: rotation normal form ("push Cliffords to the end").
#
# A program is rewritten as
#
#     U = PF · F · Iₖ ⋯ I₁
#
# where I₁…Iₖ are the recorded items (non-Clifford rotations, conditional
# rotations, and measurements, all expressed in *input* coordinates), F is an
# unconditional Clifford frame, and PF is an outcome-dependent Pauli frame (affine
# in the outcome variables, tracked up to phase). Processing a new operation G that
# acts after everything so far:
#
# * Pauli `Q^c` (θ ≡ π/2, possibly conditional): PF ← PF·Q^c.
# * Unconditional Clifford rotation exp(-iθ(-1)^s Q), θ ≡ ±π/4:
#   G·PF = PF·exp(-iθ(-1)^{s⊕c} Q) ∝ PF·Q^{s⊕c}·exp(-iθQ)  with c = [PF, Q],
#   so PF ← PF·Q^{s⊕c} and F ← exp(-iθQ)·F.
# * Any other rotation (non-Clifford angle, or conditional Clifford):
#   record exp(-iθ(-1)^{s⊕c} F†QF) as an item.
# * Measurement of Q: equivalent to measuring (-1)^c F†QF on the item state.
#
# Angles are exact rationals in units of π, normalized mod 1 (exp(-iπP) = −I is a
# global phase) into (0, 1/2] with the sign carried by `sign`.

@enum ItemKind ITEM_ROTATION ITEM_MEASURE

"""
    NFItem

One item of a [`NormalForm`](@ref), in input coordinates.

- `ITEM_ROTATION`: apply `exp(-iπ·angle·(-1)^sign · pauli)` in the branches where
  `cond == 1`. `angle ∈ (0, 1/2]`, `pauli` has phase `+1`.
- `ITEM_MEASURE`: measure `(-1)^sign · pauli`, writing the outcome to variable `var`
  (classical bit `bit`).

`sign` and `cond` are [`AffineBit`](@ref)s over earlier outcome variables. `src`
is the index of the originating op.
"""
struct NFItem
    kind::ItemKind
    pauli::PauliT
    angle::Rational{Int}
    sign::AffineBit
    cond::AffineBit
    var::Int
    bit::Int
    src::Int
end

"""
    NormalForm

A program in rotation normal form `U = PF · F · items[end] ⋯ items[1]`:

- `items`: the [`NFItem`](@ref)s in program order;
- `fwd`, `inv`: the Clifford frame `F` as the maps `P ↦ F P F†` and `P ↦ F† P F`;
- `pauli_frame`: the outcome-dependent Pauli frame `PF`. Its constant part is
  folded into `F`, so every bit is linear in the outcome variables.
- `final_var`: classical bit ⇒ variable of its final write.

Use [`frame_clifford`](@ref) to get `F` as a `CliffordOperator`.
"""
struct NormalForm
    n::Int
    items::Vector{NFItem}
    fwd::PauliMap
    inv::PauliMap
    pauli_frame::SymbolicPauli
    final_var::Dict{Int,Int}
end

"""
    frame_clifford(nf::NormalForm) -> CliffordOperator

The final unconditional Clifford frame `F` of the normal form.
"""
frame_clifford(nf::NormalForm) = to_clifford(nf.fwd)

"""Normalize an angle (units of π) modulo 1 into (−1/2, 1/2]."""
function _normalize_angle(θ::Rational{Int})
    t = mod(θ, 1)
    t > 1 // 2 ? t - 1 : t
end

# Mutable builder state.
mutable struct _NFBuilder
    n::Int
    items::Vector{NFItem}
    fwd::PauliMap
    inv::PauliMap
    pf::SymbolicPauli
end

_NFBuilder(n::Int) = _NFBuilder(n, NFItem[], identity_map(n), identity_map(n), SymbolicPauli(n))

function _push_op!(b::_NFBuilder, kind::VKind, q::PauliT, angle::Rational{Int}, sign::AffineBit,
                   cond::AffineBit, var::Int, bit::Int, src::Int)
    if kind == V_MEAS
        pp = apply_map(b.inv, q)
        s = commutation_parity(b.pf, q) ⊻ (pp.phase[] == 0x02) ⊻ sign
        push!(b.items, NFItem(ITEM_MEASURE, _set_phase(pp, 0), 0 // 1, s, AffineBit(true), var, bit, src))
        return
    end
    θ = _normalize_angle(angle)
    (iszero(θ) || _isidentity(q) || iszero(cond)) && return   # identity or global phase
    if θ == 1 // 2
        mul_conditional!(b.pf, q, cond)                       # Pauli (sign is a phase)
    elseif denominator(θ) == 4 && isone(cond)
        c = commutation_parity(b.pf, q) ⊻ sign
        mul_conditional!(b.pf, q, c)
        rotate!(b.fwd, b.inv, q, θ > 0 ? 1 : 3)
    else
        pp = apply_map(b.inv, q)
        s = sign ⊻ commutation_parity(b.pf, q) ⊻ (pp.phase[] == 0x02)
        t = θ
        if t < 0
            t = -t
            s = s ⊻ true
        end
        push!(b.items, NFItem(ITEM_ROTATION, _set_phase(pp, 0), t, s, cond, 0, 0, src))
    end
    nothing
end

# Fold the constant part of the Pauli frame into the Clifford frame.
function _finalize(b::_NFBuilder, final_var::Dict{Int,Int})
    n = b.n
    pc = PauliOperator(0x00, [b.pf.x[j].constant for j in 1:n], [b.pf.z[j].constant for j in 1:n])
    if !_isidentity(pc)
        pauli_sign_conjugate_fwd!(b.fwd, pc)
        pauli_sign_conjugate_inv!(b.inv, pc)
    end
    pf = SymbolicPauli([AffineBit(false, a.vars) for a in b.pf.x], [AffineBit(false, a.vars) for a in b.pf.z])
    NormalForm(n, b.items, b.fwd, b.inv, pf, copy(final_var))
end

"""
    normal_form(prog::LoweredProgram) -> NormalForm
    normal_form(c::Circuit; n = width of c) -> NormalForm

Rewrite a program into rotation normal form: all unconditional Clifford
operations are pushed to a final frame, and non-Clifford rotations, conditional
rotations and measurements are recorded in input coordinates. Exact: angles are
rationals, and Pauli signs and outcome dependences are tracked symbolically.

The `Circuit` method throws an `ArgumentError` if the circuit contains a
construct the verifier does not support (see [`lower`](@ref)).
"""
function normal_form(prog::LoweredProgram)
    b = _NFBuilder(prog.n)
    for v in prog.ops
        _push_op!(b, v.kind, v.pauli, v.angle, v.sign, v.cond, v.var, v.bit, v.src)
    end
    _finalize(b, prog.final_var)
end

function normal_form(c::AbstractVector{CircuitOp.Type}; n::Int = _width(c))
    prog = lower(c, n)
    isempty(prog.unsupported) || throw(ArgumentError("cannot build a normal form: " *
        join(("op $(i): $r" for (i, r) in prog.unsupported), "; ")))
    normal_form(prog)
end

function Base.show(io::IO, it::NFItem)
    if it.kind == ITEM_MEASURE
        print(io, "Measure((-1)^(", it.sign, ") ", it.pauli, " → m", it.var, " [bit ", it.bit, "])")
    else
        print(io, "Rot(", it.pauli, ", ", it.angle, "π, sign=", it.sign)
        isone(it.cond) || print(io, ", if ", it.cond)
        print(io, ")")
    end
end

function Base.show(io::IO, ::MIME"text/plain", nf::NormalForm)
    println(io, "NormalForm on ", nf.n, " qubits with ", length(nf.items), " items:")
    for (k, it) in enumerate(nf.items)
        println(io, "  [", k, "] ", it)
    end
    println(io, "  Clifford frame: ", isequal(nf.fwd, identity_map(nf.n)) ? "identity" : "non-trivial")
    print(io, "  Pauli frame: ", nf.pauli_frame)
end
