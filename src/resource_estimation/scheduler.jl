"""
    ScheduleResult

- `ops`: indices into the `EstOp` list of the scheduled operations
- `start`, `finish`: their start and finish times in seconds
- `makespan`: total time in seconds
- `logical_depth`: the longest dependency chain, counted in operations
- `critical_path`, `critical_path_time`: compiled op indices along the longest chain
  of durations and reaction delays, and its length in seconds; a makespan above it
  means blocks or magic states were the bottleneck
- `block_busy`, `utilization`: per block
- `peak_ancilla`: the most surgery and adapter ancilla qubits in use at once
- `magic_waits`: total time operations waited for a magic state
"""
struct ScheduleResult
    ops::Vector{Int}
    start::Vector{Float64}
    finish::Vector{Float64}
    makespan::Float64
    logical_depth::Int
    critical_path::Vector{Int}
    critical_path_time::Float64
    block_busy::Vector{Float64}
    utilization::Vector{Float64}
    peak_ancilla::Int
    magic_waits::Float64
end

"""
    dependency_dag(ops, t_reaction) -> (preds, delays, costed)

Dependencies between the costed operations (`costed`, the indices of ops whose class
is not `:free`). `preds[j]` lists the earlier operations that `j` must wait for, and
`delays[j]` the matching delay: 0 if their Paulis anticommute, `t_reaction` if `j`
reads a bit written by that operation. All indices in `preds` are positions in `costed`.
"""
function dependency_dag(ops::Vector{EstOp}, t_reaction::Float64)
    S = [k for k in eachindex(ops) if ops[k].class != :free]
    m = length(S)
    preds = [Int[] for _ in 1:m]; delays = [Float64[] for _ in 1:m]
    writer = Dict{Int,Int}()            # bit => position in S of its latest writer
    for j in 1:m
        oj = ops[S[j]]
        cls = Dict{Int,Float64}()
        for b in oj.reads
            haskey(writer, b) && (cls[writer[b]] = t_reaction)
        end
        for i in 1:j-1
            d = get(cls, i, -1.0)
            if comm(ops[S[i]].pauli, oj.pauli) == 0x01
                d = max(d, 0.0)
            end
            d >= 0 && (push!(preds[j], i); push!(delays[j], d))
        end
        oj.writes > 0 && (writer[oj.writes] = j)
    end
    (preds, delays, S)
end

"""
    schedule_operations(ops, layout, hw; n_factories = 1, time_per_state = 0.0) -> ScheduleResult

Schedule the costed operations under resource limits and return a deterministic
[`ScheduleResult`](@ref):
- each block runs at most `max_parallel_surgeries` operations at once, and a
  cross-block operation uses a slot on every block it touches;
- magic state `j` is ready at `ceil(j / n_factories) * time_per_state` seconds,
  with unlimited storage;
- a classically conditioned operation starts at least `t_reaction` after the
  measurement it reads.
"""
function schedule_operations(ops::Vector{EstOp}, layout::LogicalLayout, hw::HardwareModel;
                             n_factories::Int = 1, time_per_state::Float64 = 0.0)
    n_factories >= 1 || throw(ArgumentError("need at least one factory"))
    preds, delays, S = dependency_dag(ops, hw.t_reaction)
    m = length(S)
    dur = [duration(ops[S[j]]) for j in 1:m]
    succs = [Int[] for _ in 1:m]; sdel = [Float64[] for _ in 1:m]
    for j in 1:m, (i, d) in zip(preds[j], delays[j])
        push!(succs[i], j); push!(sdel[i], d)
    end
    # longest weighted path to a sink, used as the scheduling priority
    bl = zeros(m)
    for i in m:-1:1
        bl[i] = dur[i] + maximum((sdel[i][k] + bl[succs[i][k]] for k in eachindex(succs[i])); init = 0.0)
    end
    ef = zeros(m); depth = zeros(Int, m); via = zeros(Int, m)
    for j in 1:m
        best = 0.0; bd = 0
        for (i, d) in zip(preds[j], delays[j])
            if via[j] == 0 || ef[i] + d > best
                best = ef[i] + d; via[j] = i
            end
            bd = max(bd, depth[i])
        end
        ef[j] = best + dur[j]; depth[j] = bd + 1
    end
    nb = length(layout.code_blocks)
    slots = [zeros(max_parallel_surgeries(layout.code_blocks[b])) for b in 1:nb]
    start = zeros(m); finish = zeros(m)
    remaining = [length(preds[j]) for j in 1:m]
    ready = [j for j in 1:m if remaining[j] == 0]
    magic_used = 0; waits = 0.0
    busy = zeros(nb)
    scheduled = 0
    while !isempty(ready)
        k = argmax(i -> (bl[ready[i]], -ready[i]), eachindex(ready))
        j = ready[k]; deleteat!(ready, k)
        op = ops[S[j]]
        t = maximum((finish[i] + d for (i, d) in zip(preds[j], delays[j])); init = 0.0)
        for b in op.blocks
            t = max(t, minimum(slots[b]))
        end
        if op.consumes_magic
            magic_used += 1
            r = cld(magic_used, n_factories) * time_per_state
            r > t && (waits += r - t)
            t = max(t, r)
        end
        start[j] = t; finish[j] = t + dur[j]
        for b in op.blocks
            s = argmin(slots[b]); slots[b][s] = finish[j]
            busy[b] += dur[j]
        end
        scheduled += 1
        for i in succs[j]
            remaining[i] -= 1
            remaining[i] == 0 && push!(ready, i)
        end
    end
    scheduled == m || error("internal: dependency cycle")
    makespan = maximum(finish; init = 0.0)
    cp = Int[]
    if m > 0
        j = argmax(ef)
        while j != 0
            pushfirst!(cp, ops[S[j]].index)
            j = via[j]
        end
    end
    ev = sort!(vcat([(start[j], ops[S[j]].ancilla) for j in 1:m], [(finish[j], -ops[S[j]].ancilla) for j in 1:m]);
               by = e -> (e[1], e[2]))
    cur = 0; peak = 0
    for (_, a) in ev
        cur += a; peak = max(peak, cur)
    end
    util = makespan > 0 ? busy ./ (makespan .* [max_parallel_surgeries(b) for b in layout.code_blocks]) : zeros(nb)
    ScheduleResult(S, start, finish, makespan, maximum(depth; init = 0), cp, maximum(ef; init = 0.0), busy, util,
                   peak, waits)
end
