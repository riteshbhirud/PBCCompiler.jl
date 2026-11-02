@testmodule EstUtils begin
    using PBCCompiler.ResourceEstimation
    const hw = HardwareModel(:all_to_all; t_1q=1e-8, t_2q=2e-8, t_meas=5e-7, t_reaction=1e-6, p_phys=1e-3)
    const cyc = 3e-8 + 4 * 2e-8 + 5e-7                 # surface cycle 6.1e-7 s
    const am = AnalyticAdapterModel(cycle_time_factor=2.0)
    two_surface(d=3; maxpar=1) = LogicalLayout([SurfaceCodeModel(d; max_parallel=maxpar), SurfaceCodeModel(d; max_parallel=maxpar)],
                                               Dict(1 => 1, 2 => 2), Dict(1 => 1, 2 => 1))
end

@testitem "Estimation E7: scheduler dependencies, capacity and reaction time" tags=[:estimation, :e7] setup=[EstUtils] begin
using PBCCompiler
using PBCCompiler: CircuitOp
using PBCCompiler.ResourceEstimation
using QuantumClifford: @P_str
using .EstUtils: hw, cyc, am, two_surface
const C = CircuitOp
sched(c, lay) = schedule_operations(classify_operations(c, lay, hw; n_data=2, adapter_model=am), lay, hw)
r3 = 3cyc                                              # one d = 3 lattice surgery

# commuting ops on different blocks run in parallel
s = sched([C.Measurement(P"Z", 1, [1]), C.Measurement(P"Z", 2, [2])], two_surface())
@test s.makespan ≈ r3 && s.logical_depth == 1
# anticommuting ops serialize
s = sched([C.Measurement(P"ZZ", 1, [1, 2]), C.Measurement(P"X", 2, [1])], two_surface())
@test s.makespan ≈ 2 * r3 + r3 && s.logical_depth == 2 && s.critical_path == [1, 2]   # cross op takes 2·r3 (factor 2)
# commuting ops on the same block serialize at capacity 1, run in parallel at capacity 2
c = [C.Measurement(P"Z", 1, [1]), C.Measurement(P"Z", 2, [1])]
@test sched(c, two_surface()).makespan ≈ 2r3
@test sched(c, two_surface(; maxpar=2)).makespan ≈ r3
# a classical dependency adds t_reaction and puts both ops on the critical path
c = [C.Measurement(P"Z", 1, [1]), C.BitConditional(C.ExpQuatPiPauli(P"X", [2]), 1)]
s = sched(c, two_surface())
@test s.makespan ≈ r3 + 1e-6 + r3
@test s.critical_path == [1, 2] && s.critical_path_time ≈ s.makespan
# a Pauli-frame correction costs nothing and creates no dependency
c = [C.Measurement(P"Z", 1, [1]), C.BitConditional(C.ExpHalfPiPauli(P"X", [2]), 1), C.Measurement(P"Z", 2, [2])]
@test sched(c, two_surface()).makespan ≈ r3
# deterministic
c = [C.Measurement(P"Z", k, [mod1(k, 2)]) for k in 1:6]
s1 = sched(c, two_surface()); s2 = sched(c, two_surface())
@test s1.start == s2.start && s1.makespan ≈ 3r3
# peak ancilla: two parallel intra ops of 17 qubits each
@test sched([C.Measurement(P"Z", 1, [1]), C.Measurement(P"Z", 2, [2])], two_surface()).peak_ancilla == 34
end

@testitem "Estimation E8: magic-state factories and throughput search" tags=[:estimation, :e8] setup=[EstUtils] begin
using PBCCompiler
using PBCCompiler: CircuitOp, preprocess_circuit
using PBCCompiler.ResourceEstimation
using QuantumClifford: @P_str
using .EstUtils: hw, cyc
const C = CircuitOp

# GoSC 15-to-1 block: 11 tiles of 2d², 11d cycles per state (over success probability), 35p³
sp = factory_spec(FifteenToOneFactory(d=5), nothing, hw)
@test sp.qubits == 11 * 2 * 25
@test sp.time_per_state ≈ 11 * 5 / (1 - 1e-3)^15 * cyc
@test sp.output_error ≈ 35e-9
# Litinski MSD Table 1 presets, and a p_phys mismatch note
f = litinski_factory(:f_17_7_7)
@test (f.qubits, f.cycles_per_state, f.output_error, f.p_phys) == (4620, 42.6, 4.5e-8, 1e-3)
@test isempty(factory_spec(f, nothing, hw).notes)
@test !isempty(factory_spec(litinski_factory(:f_11_5_5), nothing, hw).notes)
@test_throws ArgumentError litinski_factory(:nope)

# four independent T gates on four blocks: magic-bound, so more factories help until no one waits
n = 4
input = vcat([C.ExpEighPiPauli(P"Z", [q]) for q in 1:n], [C.Measurement(P"X", q, [q]) for q in 1:n])
out = preprocess_circuit(copy(input))
lay = LogicalLayout([SurfaceCodeModel(3) for _ in 1:n], Dict(q => q for q in 1:n), Dict(q => 1 for q in 1:n))
@test magic_state_counts(out, n) == Dict(:pi8_autocorrected => 4)
ops = classify_operations(out, lay, hw; n_data=n)
spec = factory_spec(FifteenToOneFactory(d=3), lay, hw)
F, sch, table = choose_factory_count(ops, lay, hw, spec)
ms = [t[2] for t in table]
@test issorted(ms; rev=true)                           # makespan never increases with more factories
@test sch.makespan == minimum(ms) && F == table[argmin(ms)][1]
@test sch.magic_waits == 0 || F == length(table)
@test table[1][2] >= 4 * spec.time_per_state            # one factory: four states in sequence
@test 1 <= F <= 4
# spacetime objective trades time for qubits: never more factories than the makespan optimum
F2, _, _ = choose_factory_count(ops, lay, hw, spec; objective=:spacetime)
@test F2 <= F
end

@testitem "Estimation E9: full estimate, error budget and distance selection" tags=[:estimation, :e9, :e10] setup=[EstUtils] begin
using PBCCompiler
using PBCCompiler: CircuitOp, preprocess_circuit
using PBCCompiler.ResourceEstimation
using QuantumClifford: @P_str
using .EstUtils: hw, cyc, am, two_surface
const C = CircuitOp

# hand-computed: one cross-block and one intra-block measurement (they anticommute)
c = [C.Measurement(P"ZZ", 1, [1, 2]), C.Measurement(P"X", 2, [1])]
e = estimate_resources(c, two_surface(); hardware=hw, adapter_model=am)
@test e.qubit_breakdown == Dict(:data => 18, :syndrome => 16, :surgery_ancilla_peak => 39, :factories => 0, :routing => 0)
@test e.physical_qubits == 73
@test e.time_seconds ≈ 2 * 3cyc + 3cyc
@test e.time_cycles ≈ 9
pL = 0.1 * 0.1^2                                       # d = 3, p = 1e-3
@test e.error_breakdown[:blocks] ≈ 2 * pL * 9
@test e.error_breakdown[:measurements] ≈ 2 * pL * 3 + pL * 3
@test e.error_breakdown[:magic] == 0
@test e.failure_probability ≈ 0.027 && e.meets_budget == false     # budget 1e-2
@test e.n_factories == 0 && e.logical_depth == 2
d = to_dict(e)
@test d["physical_qubits"] == 73 && d["qubit_breakdown"]["data"] == 18 && length(d["per_measurement"]) == 2

# single block, one measurement
e1 = estimate_resources([C.Measurement(P"Z", 1, [1])], LogicalLayout([SurfaceCodeModel(5)], Dict(1 => 1), Dict(1 => 1));
                        hardware=hw)
@test e1.physical_qubits == 25 + 24 + 49 && e1.time_seconds ≈ 5cyc
@test e1.failure_probability ≈ 0.1 * 0.1^3 * 5 + 0.1 * 0.1^3 * 5

# classical dependency: t_reaction lands on the critical path
e2 = estimate_resources([C.Measurement(P"Z", 1, [1]), C.BitConditional(C.ExpQuatPiPauli(P"X", [2]), 1)], two_surface();
                        hardware=hw)
@test e2.time_seconds ≈ 3cyc + 1e-6 + 3cyc && e2.schedule.critical_path == [1, 2]

# invariants: total = sum of the breakdown; distance never decreases qubits; failure decreases with d
input = vcat([C.ExpEighPiPauli(P"Z", [q]) for q in 1:3], [C.PauliConditional(P"Z", [1], P"X", [2])],
             [C.Measurement(P"Z", q, [q]) for q in 1:3])
out = preprocess_circuit(copy(input))
prevq = 0; prevf = Inf
for dd in (3, 5, 7, 9)
    lay = LogicalLayout([SurfaceCodeModel(dd) for _ in 1:3], Dict(q => q for q in 1:3), Dict(q => 1 for q in 1:3))
    est = estimate_resources(out, lay; hardware=hw, adapter_model=am, factory=FifteenToOneFactory(d=7))
    @test est.physical_qubits == sum(values(est.qubit_breakdown))
    @test est.physical_qubits >= prevq
    @test est.error_breakdown[:blocks] + est.error_breakdown[:measurements] < prevf
    global prevq = est.physical_qubits
    global prevf = est.error_breakdown[:blocks] + est.error_breakdown[:measurements]
end

# distance optimization: meets a reachable budget with the smallest distance that does
lay = LogicalLayout([SurfaceCodeModel(3) for _ in 1:3], Dict(q => q for q in 1:3), Dict(q => 1 for q in 1:3))
opt = estimate_resources(out, lay; hardware=hw, adapter_model=am, error_budget=1e-4, optimize_distances=true)
@test opt.meets_budget
dsel = opt.distances[:surface]
smaller = LogicalLayout([SurfaceCodeModel(dsel - 2) for _ in 1:3], Dict(q => q for q in 1:3), Dict(q => 1 for q in 1:3))
@test !estimate_resources(out, smaller; hardware=hw, adapter_model=am, error_budget=1e-4).meets_budget
# an unreachable budget (below the magic-state error, which no code distance changes) is reported
bad = estimate_resources(out, lay; hardware=hw, adapter_model=am, error_budget=1e-12, optimize_distances=true,
                         factory=FifteenToOneFactory(d=3))
@test !bad.meets_budget && any(n -> occursin("not reachable", n), bad.notes)

# p_L from a table of the user's own simulations
path = tempname() * ".csv"
write(path, "code,d,p,p_L,ci_low,ci_high\nsurf,3,0.0005,1e-4,0,0\nsurf,3,0.002,4e-3,0,0\n")
lay = LogicalLayout([SurfaceCodeModel(3; pl=TablePL(path; code="surf"))], Dict(1 => 1), Dict(1 => 1))
et = estimate_resources([C.Measurement(P"Z", 1, [1])], lay; hardware=hw)
@test 1e-4 < et.per_block[1].p_L_per_cycle < 4e-3
rm(path)

# all-to-all hardware has zero routing; grid hardware adds routing qubits and slows cycles
hwg = HardwareModel(:grid; t_1q=1e-8, t_2q=2e-8, t_meas=5e-7, t_reaction=1e-6, p_phys=1e-3)
eg = estimate_resources(c, two_surface(); hardware=hwg, adapter_model=am, embedding_iterations=3000)
@test eg.qubit_breakdown[:routing] >= 0 && eg.time_seconds >= e.time_seconds
@test e.qubit_breakdown[:routing] == 0

# adding a measurement never decreases the makespan
base = [C.Measurement(P"Z", 1, [1]), C.Measurement(P"X", 2, [2])]
for extra in (C.Measurement(P"Z", 3, [2]), C.Measurement(P"XX", 3, [1, 2]), C.Measurement(P"Y", 3, [1]))
    t0 = estimate_resources(base, two_surface(); hardware=hw, adapter_model=am).time_seconds
    t1 = estimate_resources(vcat(base, [extra]), two_surface(); hardware=hw, adapter_model=am).time_seconds
    @test t1 >= t0
end
end
