"""
    EstOp

One operation of the compiled program with its classification and cost:

- `index`: index of the op in the compiled program
- `kind`: `:measurement`, `:clifford`, `:t_rotation`, `:frame` or `:other`
- `class`: `:intra` (one block), `:cross` (several blocks), `:magic` (magic qubits only) or `:free`
- `pauli`: the full-width Pauli
- `blocks`, `weights`: data blocks touched and the logical weight on each
- `consumes_magic`: whether the op uses a magic state
- `reads`, `writes`: classical bit indices (`writes = 0` if none)
- `ancilla`, `rounds`, `round_time` (seconds), `d_joint`, `cost_source`, `note`
- `alt_ancilla`: the analytic model's ancilla count when the constructed model was used, else -1

Pi/4 rotations and `PauliConditional`s are costed as one Pauli product measurement;
Pauli and pi/2 rotations are Pauli-frame updates at no cost (Litinski 2019, Sec. 1).
"""
mutable struct EstOp
    index::Int
    kind::Symbol
    class::Symbol
    pauli::EPauli
    blocks::Vector{Int}
    weights::Vector{Int}
    consumes_magic::Bool
    reads::Vector{Int}
    writes::Int
    ancilla::Int
    rounds::Int
    round_time::Float64
    d_joint::Int
    cost_source::Symbol
    note::String
    alt_ancilla::Int
end

duration(op::EstOp) = op.rounds * op.round_time

function _blocks_and_weights(support::Vector{Int}, n_data::Int, layout::LogicalLayout)
    w = Dict{Int,Int}()
    for q in support
        q <= n_data || continue
        b = layout.qubit_to_block[q]
        w[b] = get(w, b, 0) + 1
    end
    bs = sort!(collect(keys(w)))
    (bs, [w[b] for b in bs])
end

"""
    classify_operations(compiled, layout, hw; n_data, adapter_model = nothing,
                        adapter_kwargs = (;)) -> Vector{EstOp}

Classify and cost every op of a compiled program, returning one [`EstOp`](@ref)
per op. Data qubits are `1:n_data`; higher indices are magic qubits. Throws if a
cross-block measurement exists but no `adapter_model` is given.
"""
function classify_operations(c::AbstractVector{CircuitOp.Type}, layout::LogicalLayout, hw::HardwareModel;
                             n_data::Int, adapter_model::Union{Nothing,AdapterCostModel} = nothing,
                             adapter_kwargs = (;))
    N = max(n_data, _est_width(c))
    ref_cycle = maximum(b -> cycle_time(b, hw), layout.code_blocks)
    ops = EstOp[]
    for (i, op) in enumerate(c)
        reads = Int[]; writes = 0
        kind = :other; p = nothing; q = nothing; magic = false
        m = _measurement_parts(op)
        inner = op
        if m !== nothing
            p, writes, q = m
            kind = :measurement
        else
            bc = _bitcond_parts_est(op)
            if bc !== nothing
                inner, b = bc
                push!(reads, b)
            end
            r = _rotation_parts(inner)
            if r !== nothing
                p, q, θ = r
                kind = θ == 1 // 2 ? :frame : θ == 1 // 4 ? :clifford : :t_rotation
            else
                pc = _paulicond_parts_est(inner)
                if pc !== nothing
                    cp, cq, tp, tq = pc
                    p = QuantumClifford.tensor(cp, tp); q = vcat(cq, tq)
                    kind = :clifford
                else
                    kind = :other
                end
            end
        end
        if p === nothing
            push!(ops, EstOp(i, kind, :free, zero(PauliOperator, N), Int[], Int[], false, reads, writes, 0, 0, 0.0, 0,
                             :none, "not costed", -1))
            continue
        end
        q === nothing && error("internal: op $i has a Pauli but no qubit list")
        P = embed(N, q, p)
        support = _pauli_support(p, q)
        magic = any(>(n_data), support)
        bs, ws = _blocks_and_weights(support, n_data, layout)
        consumes = (kind == :measurement && magic && !isempty(bs)) || kind == :t_rotation
        if kind == :frame
            push!(ops, EstOp(i, kind, :free, P, bs, ws, false, reads, writes, 0, 0, 0.0, 0, :none, "Pauli frame", -1))
        elseif isempty(bs)
            # measurement of an injected magic patch: one code cycle
            push!(ops, EstOp(i, kind, :magic, P, bs, ws, false, reads, writes, 0, 1, ref_cycle, 0, :magic, "", -1))
        elseif length(bs) == 1
            model = layout.code_blocks[bs[1]]
            anc, rounds = intra_measurement_cost(model, ws[1])
            push!(ops, EstOp(i, kind, :intra, P, bs, ws, consumes, reads, writes, anc, rounds, cycle_time(model, hw),
                             code_d(model), :intra, "", -1))
        else
            adapter_model === nothing &&
                throw(ArgumentError("op $i spans blocks $bs: cross-block measurements need an adapter_model " *
                                    "(e.g. AnalyticAdapterModel(cycle_time_factor = ...))"))
            ac = adapter_cost(adapter_model, layout.code_blocks[bs], ws, hw; layout, qubits = support, pauli = P,
                              n_data, adapter_kwargs...)
            push!(ops, EstOp(i, kind, :cross, P, bs, ws, consumes, reads, writes, ac.ancilla, ac.rounds, ac.round_time,
                             ac.d_joint, ac.source, ac.note, ac.alt_ancilla))
        end
    end
    ops
end

"""
    weight_histogram(ops) -> Dict{Int,Int}

Histogram of per-block logical weights over all costed measurements: maps each
weight to the number of (measurement, block) pairs with that weight.
"""
function weight_histogram(ops::Vector{EstOp})
    h = Dict{Int,Int}()
    for op in ops, w in op.weights
        op.class in (:intra, :cross) || continue
        h[w] = get(h, w, 0) + 1
    end
    h
end
