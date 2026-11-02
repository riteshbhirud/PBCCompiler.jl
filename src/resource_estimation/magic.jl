# Litinski 2019: "A Game of Surface Codes", Quantum 3, 128, arXiv:1808.02892
# Litinski 2019b: "Magic State Distillation: Not as Costly as You Think", Quantum 3, 205, arXiv:1905.06903
# Fowler et al. 2012: "Surface codes: Towards practical large-scale quantum computation", arXiv:1208.0928

"""A magic-state factory model."""
abstract type FactoryModel end

"""
    FifteenToOneFactory(; d = nothing)

The 15-to-1 distillation block at code distance `d` (Litinski 2019, Secs. 3 and 4):
- qubits `= 11 * 2d^2`
- expected code cycles per output state `= 11d / (1 - p)^15`
- output error `= 35p^3`, with `p` the hardware's physical error rate
With `d = nothing`, the largest surface-code distance in the layout is used; `d` is
required if the layout has no surface-code block.
"""
Base.@kwdef struct FifteenToOneFactory <: FactoryModel
    d::Union{Nothing,Int} = nothing
end

"""
    TableFactory(name, qubits, cycles_per_state, output_error, p_phys)

A factory with fixed qubit count, surface-code cycles per state and output error,
valid at physical error rate `p_phys`. Using it at another `p_phys` is recorded as a
warning in the estimate. See [`litinski_factory`](@ref) for presets.
"""
struct TableFactory <: FactoryModel
    name::Symbol
    qubits::Int
    cycles_per_state::Float64
    output_error::Float64
    p_phys::Float64
end

const _MSD_TABLE1 = Dict(
    :f_7_3_3 => (810, 18.1, 4.4e-8, 1e-4),
    :f_9_3_3 => (1150, 18.1, 9.3e-10, 1e-4),
    :f_11_5_5 => (2070, 30.0, 1.9e-11, 1e-4),
    :f_17_7_7 => (4620, 42.6, 4.5e-8, 1e-3),
)

"""
    litinski_factory(name) -> TableFactory

Single-level 15-to-1 factories (Litinski 2019b, Table 1): `:f_7_3_3`, `:f_9_3_3`,
`:f_11_5_5` (at `p = 1e-4`) and `:f_17_7_7` (at `p = 1e-3`). The name gives the
distances `(d_X, d_Z, d_m)`.
"""
function litinski_factory(name::Symbol)
    haskey(_MSD_TABLE1, name) || throw(ArgumentError("unknown preset $name (have $(sort!(collect(keys(_MSD_TABLE1)))))"))
    q, c, e, p = _MSD_TABLE1[name]
    TableFactory(name, q, c, e, p)
end

_factory_d(f::FifteenToOneFactory, layout) = f.d !== nothing ? f.d :
    begin
        ds = [code_d(b) for b in layout.code_blocks if b isa SurfaceCodeModel]
        isempty(ds) && throw(ArgumentError("FifteenToOneFactory needs d when the layout has no surface-code block"))
        maximum(ds)
    end

_surface_cycle(hw) = 3hw.t_1q + 4hw.t_2q + hw.t_meas        # Fowler et al. 2012, as in SurfaceCodeModel

"""
    FactorySpec

A factory instantiated for an estimate: its qubits, seconds per state, output error
per state, and notes (for example a `p_phys` mismatch).
"""
struct FactorySpec
    qubits::Int
    time_per_state::Float64
    output_error::Float64
    notes::Vector{String}
end

function factory_spec(f::FifteenToOneFactory, layout, hw::HardwareModel)
    d = _factory_d(f, layout); p = hw.p_phys
    succ = (1 - p)^15
    FactorySpec(11 * 2d^2, 11d / succ * _surface_cycle(hw), 35p^3, String[])
end

function factory_spec(f::TableFactory, layout, hw::HardwareModel)
    notes = String[]
    isapprox(hw.p_phys, f.p_phys; rtol = 1e-9) ||
        push!(notes, "factory $(f.name) is tabulated at p = $(f.p_phys), hardware has p = $(hw.p_phys)")
    FactorySpec(f.qubits, f.cycles_per_state * _surface_cycle(hw), f.output_error, notes)
end

"""
    magic_state_counts(compiled, n_data) -> Dict{Symbol,Int}

Number of magic states consumed by a compiled program, by gadget type. Remaining
pi/8 rotations are counted as `:t_rotation`.
"""
function magic_state_counts(c::AbstractVector{CircuitOp.Type}, n_data::Int)
    d = count_gadgets(c, n_data)
    nt = count(op -> (r = _rotation_parts(op); r !== nothing && r[3] == 1 // 8), c)
    nt > 0 && (d[:t_rotation] = nt)
    d
end

"""
    choose_factory_count(ops, layout, hw, spec; objective = :makespan, max_factories = 64)
        -> (n_factories, schedule, table)

Choose the number of factories that minimizes `objective`:
- `:makespan`: shortest run time, preferring fewer factories on ties
- `:spacetime`: total qubits times run time
Returns the chosen count, its [`ScheduleResult`](@ref), and a table of
`(n, makespan, objective value)` for every count tried.
"""
function choose_factory_count(ops::Vector{EstOp}, layout::LogicalLayout, hw::HardwareModel, spec::FactorySpec;
                              objective::Symbol = :makespan, max_factories::Int = 64, base_qubits::Int = 0)
    objective in (:makespan, :spacetime) || throw(ArgumentError("objective must be :makespan or :spacetime"))
    max_factories >= 1 || throw(ArgumentError("max_factories must be ≥ 1"))
    nmagic = count(o -> o.consumes_magic, ops)
    nmagic == 0 && return (0, schedule_operations(ops, layout, hw; n_factories = 1, time_per_state = 0.0), Tuple{Int,Float64,Float64}[])
    best = nothing; table = Tuple{Int,Float64,Float64}[]
    for F in 1:min(max_factories, nmagic)
        sch = schedule_operations(ops, layout, hw; n_factories = F, time_per_state = spec.time_per_state)
        val = objective == :makespan ? sch.makespan : (base_qubits + sch.peak_ancilla + F * spec.qubits) * sch.makespan
        push!(table, (F, sch.makespan, val))
        if best === nothing || val < best[3] * (1 - 1e-12)
            best = (F, sch, val)
        end
        sch.magic_waits == 0 && break
    end
    best === nothing && error("internal: no factory count was evaluated")
    (best[1], best[2], table)
end
