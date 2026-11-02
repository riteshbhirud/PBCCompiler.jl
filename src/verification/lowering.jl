@enum VKind V_ROT V_MEAS

"""
    VOp

One lowered operation.

- Rotation (`V_ROT`): apply `exp(-iπ·angle·(-1)^sign · pauli)` in the branches
  where `cond == 1`.
- Measurement (`V_MEAS`): measure `(-1)^sign · pauli` and write the outcome
  (1 means eigenvalue −1) to variable `var`, which stands for classical bit `bit`.

`src` is the index of the originating op in the input circuit.
"""
struct VOp
    kind::VKind
    pauli::PauliT
    angle::Rational{Int}
    sign::AffineBit
    cond::AffineBit
    var::Int
    bit::Int
    src::Int
end

_vrot(p, angle, cond, src) = VOp(V_ROT, p, angle, AffineBit(false), cond, 0, 0, src)

"""
    LoweredProgram

Result of lowering a circuit: the op list over `n` qubits, the outcome variable
of the final write to each classical bit (`final_var`), the classical bit of every
variable (`var_bit`), and every construct that could not be lowered
(`unsupported`, as `(op index, reason)` pairs).
"""
struct LoweredProgram
    n::Int
    ops::Vector{VOp}
    final_var::Dict{Int,Int}
    var_bit::Dict{Int,Int}
    unsupported::Vector{Tuple{Int,String}}
end

function _embed_checked(p::PauliT, q::Vector{Int}, n::Int, what::String)
    nqubits(p) == length(q) || return "$what: Pauli has $(nqubits(p)) letters for $(length(q)) qubits"
    length(unique(q)) == length(q) || return "$what: repeated qubit in $q"
    all(x -> 1 <= x <= n, q) || return "$what: qubit index outside 1:$n in $q"
    _ishermitian(p) || return "$what: non-Hermitian Pauli $p"
    embed(n, q, p)
end

"""
    lower(c, n::Int; replacements = Dict()) -> LoweredProgram

Lower the circuit `c` to rotations and measurements on qubits `1:n`. Unsupported
constructs are recorded in `unsupported`: `PrepMagic`, nested `BitConditional`,
conditional measurements, reads before writes, and malformed ops (including ops
on qubits outside `1:n`).

Each measurement gets a fresh outcome variable (its op index), and a
`BitConditional` reads the variable of the latest write to its bit.

`replacements[i] = (j, vops, writes)` replaces ops `i:j` by the given `vops`, and
records each `(bit, var)` in `writes` as a write of `bit` with variable `var`.
"""
function lower(c::AbstractVector{CircuitOp.Type}, n::Int;
               replacements::Dict{Int,Tuple{Int,Vector{VOp},Vector{Tuple{Int,Int}}}} =
                   Dict{Int,Tuple{Int,Vector{VOp},Vector{Tuple{Int,Int}}}}())
    prog = LoweredProgram(n, VOp[], Dict{Int,Int}(), Dict{Int,Int}(), Tuple{Int,String}[])
    i = 1
    while i <= length(c)
        if haskey(replacements, i)
            j, vops, writes = replacements[i]
            append!(prog.ops, vops)
            for (b, v) in writes
                prog.final_var[b] = v
                prog.var_bit[v] = b
            end
            i = j + 1
        else
            _lower_op!(prog, c[i], i, AffineBit(true))
            i += 1
        end
    end
    prog
end

function _lower_op!(prog::LoweredProgram, op, i::Int, cond::AffineBit)
    n = prog.n
    k = _opkind(op)
    if k == OP_MEASUREMENT
        if !isone(cond)
            push!(prog.unsupported, (i, "measurement inside a BitConditional"))
            return
        end
        p, q = _pauli_and_qubits(op)
        b = _measurement_bit(op)
        e = _embed_checked(p, q, n, "Measurement")
        e isa String && (push!(prog.unsupported, (i, e)); return)
        b >= 1 || (push!(prog.unsupported, (i, "Measurement writes invalid bit $b")); return)
        push!(prog.ops, VOp(V_MEAS, e, 0 // 1, AffineBit(false), AffineBit(true), i, b, i))
        prog.final_var[b] = i
        prog.var_bit[i] = b
    elseif _isrotation(k)
        p, q = _pauli_and_qubits(op)
        e = _embed_checked(p, q, n, _kindname(k))
        e isa String && (push!(prog.unsupported, (i, e)); return)
        push!(prog.ops, _vrot(e, _rotation_angle(k), cond, i))
    elseif k == OP_PAULICOND
        cp, cq, tp, tq = _paulicond_parts(op)
        isempty(intersect(cq, tq)) || (push!(prog.unsupported, (i, "PauliConditional control and target overlap")); return)
        (isempty(cq) || isempty(tq)) && (push!(prog.unsupported, (i, "PauliConditional with empty control or target")); return)
        p1 = _embed_checked(cp, cq, n, "PauliConditional control")
        p1 isa String && (push!(prog.unsupported, (i, p1)); return)
        p2 = _embed_checked(tp, tq, n, "PauliConditional target")
        p2 isa String && (push!(prog.unsupported, (i, p2)); return)
        # (I + P₁ + P₂ − P₁P₂)/2 ∝ exp(-iπ/4 P₁) exp(-iπ/4 P₂) exp(iπ/4 P₁P₂)
        push!(prog.ops, _vrot(p1, 1 // 4, cond, i))
        push!(prog.ops, _vrot(p2, 1 // 4, cond, i))
        push!(prog.ops, _vrot(p1 * p2, -1 // 4, cond, i))
    elseif k == OP_BITCOND
        inner, b = _bitcond_parts(op)
        if !isone(cond)
            push!(prog.unsupported, (i, "nested BitConditional (condition is a product of bits, not affine)"))
            return
        end
        if !haskey(prog.final_var, b)
            push!(prog.unsupported, (i, "BitConditional reads bit $b before it is written"))
            return
        end
        _lower_op!(prog, inner, i, affine_var(prog.final_var[b]))
    else # OP_PREPMAGIC
        push!(prog.unsupported, (i, "PrepMagic has no defined semantics in PBCCompiler"))
    end
    nothing
end
