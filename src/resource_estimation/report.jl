_fmt_time(s) = s >= 3600 ? string(round(s / 3600; sigdigits = 3), " h") : s >= 60 ? string(round(s / 60; sigdigits = 3), " min") :
               s >= 1 ? string(round(s; sigdigits = 3), " s") : s >= 1e-3 ? string(round(s * 1e3; sigdigits = 3), " ms") :
               string(round(s * 1e6; sigdigits = 3), " µs")

function Base.show(io::IO, ::MIME"text/plain", e::ResourceEstimate)
    println(io, "ResourceEstimate")
    println(io, "  physical qubits: ", e.physical_qubits)
    for k in (:data, :syndrome, :surgery_ancilla_peak, :factories, :routing)
        println(io, "    ", rpad(string(k), 22), e.qubit_breakdown[k])
    end
    println(io, "  time: ", _fmt_time(e.time_seconds), " (", round(e.time_cycles; sigdigits = 4),
            " cycles of ", _fmt_time(e.reference_cycle), "), logical depth ", e.logical_depth)
    println(io, "  failure probability: ", round(e.failure_probability; sigdigits = 3), " (budget ", e.error_budget,
            e.meets_budget ? ", met)" : ", NOT met)")
    for k in (:blocks, :measurements, :magic)
        println(io, "    ", rpad(string(k), 22), round(e.error_breakdown[k]; sigdigits = 3))
    end
    println(io, "  distances: ", e.distances)
    println(io, "  magic states: ", e.magic_states, ", factories: ", e.n_factories)
    println(io, "  measurement weight histogram: ", sort!(collect(e.weight_histogram)))
    println(io, "  critical path: ", length(e.schedule.critical_path), " ops, ",
            _fmt_time(e.schedule.critical_path_time), " (makespan ", _fmt_time(e.time_seconds), ")")
    println(io, "  per block:")
    for r in e.per_block
        println(io, "    #", r.block, " ", r.model, ": ", r.qubits, " qubits, p_L/cycle ", round(r.p_L_per_cycle; sigdigits = 3),
                ", failure ", round(r.failure; sigdigits = 3), ", utilization ", round(r.utilization; sigdigits = 3))
    end
    for n in e.notes
        println(io, "  note: ", n)
    end
end

"""
    to_dict(e::ResourceEstimate) -> Dict{String,Any}

The estimate as nested `Dict`s, `Vector`s, numbers and strings, ready for JSON
serialization.
"""
function to_dict(e::ResourceEstimate)
    nt(x::NamedTuple) = Dict{String,Any}(string(k) => (v isa Symbol ? string(v) : v) for (k, v) in pairs(x))
    Dict{String,Any}(
        "physical_qubits" => e.physical_qubits,
        "qubit_breakdown" => Dict(string(k) => v for (k, v) in e.qubit_breakdown),
        "time_seconds" => e.time_seconds, "time_cycles" => e.time_cycles, "reference_cycle" => e.reference_cycle,
        "logical_depth" => e.logical_depth,
        "failure_probability" => e.failure_probability, "error_budget" => e.error_budget, "meets_budget" => e.meets_budget,
        "error_breakdown" => Dict(string(k) => v for (k, v) in e.error_breakdown),
        "per_block" => [nt(r) for r in e.per_block],
        "per_measurement" => [nt(r) for r in e.per_measurement],
        "weight_histogram" => Dict(string(k) => v for (k, v) in e.weight_histogram),
        "magic_states" => Dict(string(k) => v for (k, v) in e.magic_states),
        "n_factories" => e.n_factories,
        "factory_search" => [Dict("factories" => f, "makespan" => m, "objective" => o) for (f, m, o) in e.factory_table],
        "schedule" => Dict("makespan" => e.schedule.makespan, "critical_path" => e.schedule.critical_path,
                           "critical_path_time" => e.schedule.critical_path_time,
                           "peak_ancilla" => e.schedule.peak_ancilla, "magic_waits" => e.schedule.magic_waits,
                           "utilization" => e.schedule.utilization),
        "distances" => Dict(string(k) => v for (k, v) in e.distances),
        "routing" => [Dict("extra_qubits" => o.extra_qubits, "extra_2q_layers" => o.extra_2q_layers,
                           "non_native_edges" => o.non_native_edges, "model" => string(o.model)) for o in e.hardware_overhead],
        "notes" => e.notes,
    )
end
