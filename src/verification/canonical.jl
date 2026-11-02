# V5: canonicalization of normal forms.
#
# Two normal forms of equivalent programs can differ because commuting items were
# reordered, or because rotations about the same Pauli were merged or split. The
# canonical form:
#
# 1. (outcome semantics only) drops items that cannot affect any measurement: the
#    items that are not ancestors of a measurement in the dependency DAG. Their set
#    is closed under descendants, so they can all be moved after the last
#    measurement;
# 2. merges rotations about the same Pauli whenever every item between them commutes
#    with that Pauli (the earlier rotation is moved forward onto the later one). If
#    the merged angle becomes Clifford, it is re-absorbed into the frame and the
#    remaining items are rebuilt;
# 3. normalizes the angle of every unconditional, constant-sign non-Clifford rotation
#    modulo Cliffords: exp(-iπθP) = exp(-iπrP)·exp(-iπ(θ−r)P) with r ∈ (0, 1/4) and
#    θ − r a multiple of 1/4. The Clifford remainder is absorbed into the frame, so
#    for example T† = T·S† is always represented as T plus a frame change. Without
#    this, T†·T·T (input) and T†·S (compiler merged T·T) would have different normal
#    forms. Steps 1–3 are repeated to a fixed point;
# 4. orders items by Foata layers of the DAG (layer = 1 + max layer of the
#    predecessors), and sorts within a layer by a key that does not depend on
#    variable names.
#
# DAG edge i → j (i before j): their Paulis anticommute, or j's sign or condition
# reads the outcome variable written by i.

"""Whether two items must keep their relative order (`i` before `j`)."""
function _depends(a::NFItem, b::NFItem)
    comm(a.pauli, b.pauli) == 0x01 && return true
    if a.kind == ITEM_MEASURE
        (a.var in b.sign.vars || a.var in b.cond.vars) && return true
    end
    false
end

# Predecessor lists of the dependency DAG.
function _predecessors(items::Vector{NFItem})
    m = length(items)
    preds = [Int[] for _ in 1:m]
    for j in 1:m, i in 1:j-1
        _depends(items[i], items[j]) && push!(preds[j], i)
    end
    preds
end

"""Indices of the items that are measurements or ancestors of a measurement."""
function _measurement_ancestors(items::Vector{NFItem})
    preds = _predecessors(items)
    keep = falses(length(items))
    stack = [k for k in eachindex(items) if items[k].kind == ITEM_MEASURE]
    for k in stack
        keep[k] = true
    end
    while !isempty(stack)
        k = pop!(stack)
        for p in preds[k]
            keep[p] || (keep[p] = true; push!(stack, p))
        end
    end
    keep
end

# Rebuild a normal form from (possibly non-canonical) items in item coordinates,
# composing the new frames with the old ones: total = PF·F·(PF₂·F₂·items₂).
function _rebuild(nf::NormalForm, items::Vector{NFItem})
    b = _NFBuilder(nf.n)
    for it in items
        if it.kind == ITEM_MEASURE
            _push_op!(b, V_MEAS, it.pauli, 0 // 1, it.sign, AffineBit(true), it.var, it.bit, it.src)
        else
            _push_op!(b, V_ROT, it.pauli, it.angle, it.sign, it.cond, 0, 0, it.src)
        end
    end
    inner = _finalize(b, nf.final_var, nf.var_bit)
    fwd = compose(nf.fwd, inner.fwd)               # F·F₂
    inv = compose(inner.inv, nf.inv)               # (F·F₂)† · F·F₂ = F₂†(F†·F)F₂
    pf2 = copy(inner.pauli_frame)
    conjugate!(pf2, nf.fwd.img)                    # F·PF₂·F†
    pf = mul_frame!(copy(nf.pauli_frame), pf2)
    NormalForm(nf.n, inner.items, fwd, inv, pf, nf.final_var, nf.var_bit)
end

# Find the first mergeable pair (i, j): same Pauli and condition, sign difference
# constant, and every item strictly between commutes with the Pauli.
function _find_merge(items::Vector{NFItem})
    m = length(items)
    for j in 1:m
        b = items[j]
        b.kind == ITEM_ROTATION || continue
        for i in j-1:-1:1
            a = items[i]
            if comm(a.pauli, b.pauli) == 0x01
                break                                   # blocked
            end
            a.kind == ITEM_ROTATION || continue
            a.pauli.xz == b.pauli.xz || continue
            a.cond == b.cond || continue
            d = a.sign ⊻ b.sign
            isconstant(d) || continue
            return (i, j, d.constant)
        end
    end
    nothing
end

"""
    CanonicalForm

A [`NormalForm`](@ref) after pruning (outcome semantics), rotation merging to a
fixed point, and canonical Foata ordering. `layers[k]` is the Foata layer of
`nf.items[k]`. `merges` counts merged rotation pairs, `rebuilds` counts how
many merges produced a Clifford that was pushed back into the frame, `splits`
counts angle normalizations modulo Cliffords, and `pruned` counts items dropped
because they cannot affect any outcome.
"""
struct CanonicalForm
    nf::NormalForm
    layers::Vector{Int}
    merges::Int
    rebuilds::Int
    splits::Int
    pruned::Int
    semantics::Symbol
end

"""
    canonicalize(nf::NormalForm; semantics = :outcomes) -> CanonicalForm

Canonical form of a normal form, used for equivalence checking. With
`semantics = :outcomes`, items that cannot influence any measurement outcome are
dropped. With `:channel`, every item is kept.
"""
function canonicalize(nf::NormalForm; semantics::Symbol = :outcomes)
    semantics in (:outcomes, :channel) || throw(ArgumentError("semantics must be :outcomes or :channel"))
    merges = 0; rebuilds = 0; pruned = 0; splits = 0
    while true
        if semantics == :outcomes
            keep = _measurement_ancestors(nf.items)
            if !all(keep)
                pruned += count(!, keep)
                nf = NormalForm(nf.n, nf.items[keep], nf.fwd, nf.inv, nf.pauli_frame, nf.final_var, nf.var_bit)
            end
        end
        sp = _find_split(nf.items)
        if sp !== nothing
            k, r, c = sp
            it = nf.items[k]
            items = copy(nf.items)
            items[k] = NFItem(ITEM_ROTATION, it.pauli, r, AffineBit(false), it.cond, 0, 0, it.src)
            insert!(items, k + 1, NFItem(ITEM_ROTATION, it.pauli, c, AffineBit(false), it.cond, 0, 0, it.src))
            splits += 1
            nf = _rebuild(nf, items)
            continue
        end
        mg = _find_merge(nf.items)
        mg === nothing && break
        i, j, flip = mg
        a = nf.items[i]; b = nf.items[j]
        θ = _normalize_angle(b.angle + (flip ? -a.angle : a.angle))
        s = b.sign
        if θ < 0
            θ = -θ; s = s ⊻ true
        end
        items = copy(nf.items)
        merges += 1
        if iszero(θ)
            deleteat!(items, [i, j])
            nf = NormalForm(nf.n, items, nf.fwd, nf.inv, nf.pauli_frame, nf.final_var, nf.var_bit)
        else
            items[j] = NFItem(ITEM_ROTATION, b.pauli, θ, s, b.cond, 0, 0, b.src)
            deleteat!(items, i)
            if denominator(θ) <= 2 || (denominator(θ) == 4 && isone(b.cond))
                rebuilds += 1
                nf = _rebuild(nf, items)                # Clifford or Pauli: re-push into the frames
            else
                nf = NormalForm(nf.n, items, nf.fwd, nf.inv, nf.pauli_frame, nf.final_var, nf.var_bit)
            end
        end
    end
    # Foata layers and canonical order
    items = nf.items
    preds = _predecessors(items)
    layer = zeros(Int, length(items))
    for k in eachindex(items)
        layer[k] = 1 + (isempty(preds[k]) ? 0 : maximum(layer[p] for p in preds[k]))
    end
    order = sort(collect(eachindex(items)); by = k -> (layer[k], _item_key(items[k]), k))
    nf2 = NormalForm(nf.n, items[order], nf.fwd, nf.inv, nf.pauli_frame, nf.final_var, nf.var_bit)
    CanonicalForm(nf2, layer[order], merges, rebuilds, splits, pruned, semantics)
end

# First unconditional, constant-sign non-Clifford rotation whose signed angle is not
# already in (0, 1/4): returns (index, r, c) with signed angle = r + c, r ∈ (0, 1/4),
# c a nonzero multiple of 1/4.
function _find_split(items::Vector{NFItem})
    for (k, it) in enumerate(items)
        (it.kind == ITEM_ROTATION && isone(it.cond) && isconstant(it.sign)) || continue
        denominator(it.angle) > 4 || continue
        θ = it.sign.constant ? -it.angle : it.angle
        r = mod(θ, 1 // 4)
        (r == θ && !it.sign.constant) && continue
        return (k, r, θ - r)
    end
    nothing
end

# Ordering key that does not depend on the names of outcome variables.
_item_key(it::NFItem) = (Int(it.kind), it.pauli.xz, it.angle, it.sign.constant, length(it.sign.vars),
                         it.cond.constant, length(it.cond.vars))
