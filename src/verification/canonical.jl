function _depends(a::NFItem, b::NFItem)
    comm(a.pauli, b.pauli) == 0x01 && return true
    if a.kind == ITEM_MEASURE
        (a.var in b.sign.vars || a.var in b.cond.vars) && return true
    end
    false
end

function _predecessors(items::Vector{NFItem})
    m = length(items)
    preds = [Int[] for _ in 1:m]
    for j in 1:m, i in 1:j-1
        _depends(items[i], items[j]) && push!(preds[j], i)
    end
    preds
end

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

# The new frames are composed with the old ones: total = PF·F·(PF₂·F₂·items₂).
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
    fwd = compose(nf.fwd, inner.fwd)
    inv = compose(inner.inv, nf.inv)
    pf2 = copy(inner.pauli_frame)
    conjugate!(pf2, nf.fwd.img)
    pf = mul_frame!(copy(nf.pauli_frame), pf2)
    NormalForm(nf.n, inner.items, fwd, inv, pf, nf.final_var, nf.var_bit)
end

function _find_merge(items::Vector{NFItem})
    m = length(items)
    for j in 1:m
        b = items[j]
        b.kind == ITEM_ROTATION || continue
        for i in j-1:-1:1
            a = items[i]
            if comm(a.pauli, b.pauli) == 0x01
                break
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

The canonical form of a [`NormalForm`](@ref), returned by [`canonicalize`](@ref).
`nf` holds the canonical items and `layers[k]` is the dependency layer of
`nf.items[k]`. The counters `merges`, `rebuilds`, `splits` and `pruned` record
how many rotations were merged, how many merges were absorbed back into the
Clifford frame, how many angles were reduced modulo Cliffords, and how many items
were dropped because they cannot affect any outcome.
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

Return the canonical form of a normal form, for use with [`compare_forms`](@ref).
With `semantics = :outcomes` (the default), items that cannot influence any
measurement outcome are dropped. With `:channel`, every item is kept.
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
                nf = _rebuild(nf, items)
            else
                nf = NormalForm(nf.n, items, nf.fwd, nf.inv, nf.pauli_frame, nf.final_var, nf.var_bit)
            end
        end
    end
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

# Returns (k, r, c) with signed angle r + c, r in (0, 1/4) and c a multiple of 1/4.
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

_item_key(it::NFItem) = (Int(it.kind), it.pauli.xz, it.angle, it.sign.constant, length(it.sign.vars),
                         it.cond.constant, length(it.cond.vars))
