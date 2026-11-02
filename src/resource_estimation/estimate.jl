"""
    ResourceEstimate

Result of [`estimate_resources`](@ref).

- `physical_qubits`: the sum of `qubit_breakdown` (`:data`, `:syndrome`,
  `:surgery_ancilla_peak`, `:factories`, `:routing`)
- `time_seconds`; `time_cycles`, in units of `reference_cycle` (the slowest block's
  syndrome cycle); `logical_depth`
- `failure_probability` and `error_breakdown` (`:blocks`, `:measurements`, `:magic`);
  `meets_budget` compares the failure probability with `error_budget`
- `per_block` and `per_measurement` tables, `weight_histogram`, `magic_states` by
  gadget type, `n_factories` and the factory search table
- `schedule` (a [`ScheduleResult`](@ref)), the chosen `distances` by code family,
  `hardware_overhead`, and `notes` (recorded fallbacks and warnings)

[`to_dict`](@ref) gives a JSON-friendly `Dict`.
"""
struct ResourceEstimate
    physical_qubits::Int
    qubit_breakdown::Dict{Symbol,Int}
    time_seconds::Float64
    time_cycles::Float64
    reference_cycle::Float64
    logical_depth::Int
    failure_probability::Float64
    error_budget::Float64
    meets_budget::Bool
    error_breakdown::Dict{Symbol,Float64}
    per_block::Vector{NamedTuple}
    per_measurement::Vector{NamedTuple}
    weight_histogram::Dict{Int,Int}
    magic_states::Dict{Symbol,Int}
    n_factories::Int
    factory_table::Vector{Tuple{Int,Float64,Float64}}
    schedule::ScheduleResult
    distances::Dict{Symbol,Int}
    hardware_overhead::Vector{EmbeddingOverhead}
    notes::Vector{String}
end

"""
    estimate_resources(circuit, layout; hardware, adapter_model = nothing,
                       factory = FifteenToOneFactory(), error_budget = 1e-2,
                       optimize_distances = false, n_data = nothing,
                       factory_objective = :makespan, max_factories = 64,
                       routing = :bridge, embedding_iterations = 20_000, seed = 0,
                       adapter_kwargs = (;)) -> ResourceEstimate

Estimate the physical qubits, run time and failure probability of a compiled PBC
program (`circuit`, typically the output of `preprocess_circuit`) whose data qubits
`1:n_data` are placed by `layout`. Qubits above `n_data` are magic qubits; `n_data`
defaults to the largest qubit in the layout. Returns a [`ResourceEstimate`](@ref).

`hardware` is required, as hardware times and error rates have no published
defaults. `adapter_model` is needed only if some measurement spans several blocks.
On hardware that is not all-to-all, each block's Tanner graph is embedded and the
routing qubits and extra CNOT layers are included.

With `optimize_distances = true`, each code family starts at its smallest distance
and the family contributing most to the failure probability is increased until the
estimate is within `error_budget` or no distance can be increased. All surface-code
blocks share one distance; bivariate-bicycle blocks move between the presets.
"""
function estimate_resources(c::AbstractVector{CircuitOp.Type}, layout::LogicalLayout; hardware::HardwareModel,
                            adapter_model::Union{Nothing,AdapterCostModel} = nothing,
                            factory::FactoryModel = FifteenToOneFactory(), error_budget::Real = 1e-2,
                            optimize_distances::Bool = false, n_data::Union{Nothing,Int} = nothing,
                            factory_objective::Symbol = :makespan, max_factories::Int = 64, routing::Symbol = :bridge,
                            embedding_iterations::Int = 20_000, seed::Int = 0, adapter_kwargs = (;))
    nd = n_data === nothing ? maximum(keys(layout.qubit_to_block); init = 0) : n_data
    errs = validate_layout(layout, 1:nd)
    isempty(errs) || throw(ArgumentError("invalid layout: " * join(errs, "; ")))
    kw = (; hardware, adapter_model, factory, error_budget = Float64(error_budget), n_data = nd, factory_objective,
          max_factories, routing, embedding_iterations, seed, adapter_kwargs)
    optimize_distances || return _estimate(c, layout; kw...)
    fams = unique!([distance_family(b) for b in layout.code_blocks if distance_family(b) !== nothing])
    dist = Dict(f => min_distance(first(b for b in layout.code_blocks if distance_family(b) == f)) for f in fams)
    notes = String[]
    while true
        lay = _with_distances(layout, dist)
        est = _estimate(c, lay; kw...)
        est.failure_probability <= error_budget && return _add_notes(est, notes)
        contrib = Dict(f => 0.0 for f in fams)
        for r in est.per_block
            fam = distance_family(lay.code_blocks[r.block])
            fam === nothing || (contrib[fam] += r.failure)
        end
        for r in est.per_measurement
            for b in r.blocks
                fam = distance_family(lay.code_blocks[b])
                fam === nothing || (contrib[fam] += r.failure / length(r.blocks))
            end
        end
        order = sort(fams; by = f -> -contrib[f])
        bumped = false
        for f in order
            contrib[f] > 0 || continue
            model = first(b for b in lay.code_blocks if distance_family(b) == f)
            nd2 = try
                next_distance(model, dist[f])
            catch
                nothing
            end
            nd2 === nothing && continue
            push!(notes, "increased $f distance $(dist[f]) → $nd2 (contributed $(contrib[f]))")
            dist[f] = nd2; bumped = true
            break
        end
        if !bumped
            push!(notes, "error budget $error_budget not reachable: no code family can be increased further " *
                         "(remaining failure $(est.failure_probability); magic contributes $(est.error_breakdown[:magic]))")
            return _add_notes(est, notes)
        end
    end
end

_add_notes(e::ResourceEstimate, notes) = (append!(e.notes, notes); e)

function _with_distances(layout::LogicalLayout, dist::Dict)
    blocks = CodeModel[(f = distance_family(b); f === nothing ? b : with_distance(b, dist[f])) for b in layout.code_blocks]
    LogicalLayout(blocks, layout.qubit_to_block, layout.qubit_within_block)
end

function _estimate(c, layout::LogicalLayout; hardware, adapter_model, factory, error_budget, n_data, factory_objective,
                   max_factories, routing, embedding_iterations, seed, adapter_kwargs)
    hw = hardware
    notes = String[]
    blocks = layout.code_blocks
    ops = classify_operations(c, layout, hw; n_data, adapter_model, adapter_kwargs)
    for o in ops
        isempty(o.note) || o.class == :free || push!(notes, "op $(o.index): $(o.note)")
    end
    overheads = EmbeddingOverhead[]
    extra_time = zeros(length(blocks))
    routing_qubits = 0
    if hw.connectivity !== :all_to_all
        cache = Dict{Any,EmbeddingOverhead}()
        for (b, m) in enumerate(blocks)
            key = (typeof(m), code_n(m), code_d(m))
            oh = get!(cache, key) do
                embedding_overhead(m, hw; routing, iterations = embedding_iterations, rng = Random.MersenneTwister(seed))
            end
            push!(overheads, oh)
            routing_qubits += oh.extra_qubits
            extra_time[b] = oh.extra_2q_layers * hw.t_2q
        end
        for o in ops
            isempty(o.blocks) || (o.round_time += maximum(extra_time[o.blocks]))
        end
        push!(notes, "routing overhead for adapter/merged codes is not included (their Tanner graphs are not available to the analytic model)")
    end
    cyc = [cycle_time(m, hw) + extra_time[b] for (b, m) in enumerate(blocks)]
    magic = magic_state_counts(c, n_data)
    spec = factory_spec(factory, layout, hw)
    append!(notes, spec.notes)
    base = sum(code_n, blocks; init = 0) + sum(syndrome_ancillas, blocks; init = 0) + routing_qubits
    F, sch, table = choose_factory_count(ops, layout, hw, spec; objective = factory_objective, max_factories,
                                         base_qubits = base)
    makespan = sch.makespan
    ref = maximum(cyc)
    p = hw.p_phys
    per_block = NamedTuple[]
    fb = 0.0
    for (b, m) in enumerate(blocks)
        plc = logical_error_per_cycle(m, p)
        alive = makespan / cyc[b]
        f = plc * alive
        fb += f
        push!(per_block, (block = b, model = string(m), d = code_d(m), qubits = code_n(m) + syndrome_ancillas(m),
                          p_L_per_cycle = plc, cycles_alive = alive, failure = f, utilization = sch.utilization[b]))
    end
    per_meas = NamedTuple[]
    fm = 0.0
    pos = Dict(S => k for (k, S) in enumerate(sch.ops))
    for (k, o) in enumerate(ops)
        o.class == :free && continue
        f = sum((logical_error_per_cycle(blocks[b], p) * o.rounds for b in o.blocks); init = 0.0)
        fm += f
        j = pos[k]
        push!(per_meas, (index = o.index, kind = o.kind, class = o.class, blocks = o.blocks, weights = o.weights,
                         ancilla = o.ancilla, rounds = o.rounds, seconds = duration(o), start = sch.start[j],
                         failure = f, consumes_magic = o.consumes_magic, cost_source = o.cost_source,
                         analytic_ancilla = o.alt_ancilla))
    end
    nmagic = count(o -> o.consumes_magic, ops)
    fmag = nmagic * spec.output_error
    qb = Dict(:data => sum(code_n, blocks; init = 0), :syndrome => sum(syndrome_ancillas, blocks; init = 0),
              :surgery_ancilla_peak => sch.peak_ancilla, :factories => F * spec.qubits, :routing => routing_qubits)
    total_fail = fb + fm + fmag
    dists = Dict{Symbol,Int}()
    for m in blocks
        f = distance_family(m)
        f === nothing || (dists[f] = code_d(m))
    end
    factory isa FifteenToOneFactory && nmagic > 0 && (dists[:factory] = _factory_d(factory, layout))
    ResourceEstimate(sum(values(qb)), qb, makespan, makespan / ref, ref, sch.logical_depth, total_fail, error_budget,
                     total_fail <= error_budget, Dict(:blocks => fb, :measurements => fm, :magic => fmag),
                     per_block, per_meas, weight_histogram(ops), magic, F, table, sch, dists, overheads, notes)
end
