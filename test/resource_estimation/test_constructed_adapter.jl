@testitem "Estimation E5: constructed adapter model (fork-dependent)" tags=[:estimation, :e5] begin
using PBCCompiler
using PBCCompiler: CircuitOp
using PBCCompiler.ResourceEstimation
using QuantumClifford: @P_str, embed
const C = CircuitOp

hw = HardwareModel(:all_to_all; t_1q=1e-8, t_2q=2e-8, t_meas=5e-7, t_reaction=1e-6, p_phys=1e-3)
s5 = SurfaceCodeModel(5)
lay = LogicalLayout([s5, s5, s5], Dict(1 => 1, 2 => 2, 3 => 3), Dict(1 => 1, 2 => 1, 3 => 1))
cm = ConstructedAdapterModel(cycle_time_factor=2.0)
analytic = adapter_cost(AnalyticAdapterModel(cycle_time_factor=2.0), [s5, s5], [1, 1], hw).ancilla

# Always testable: three-block measurements fall back to the analytic model, and the
# fallback is recorded in the estimate's notes (never silent).
c3 = [C.Measurement(P"ZZZ", 1, [1, 2, 3])]
e = estimate_resources(c3, lay; hardware=hw, adapter_model=cm)
@test e.per_measurement[1].cost_source == :analytic
@test any(n -> occursin("constructed adapter fallback", n), e.notes)

if !adapters_available()
    @info "E5 construction tests skipped: QuantumClifford.ECC.build_adapter is not available. " *
          "Install the riteshbhirud/QuantumClifford.jl fork (branch master, PR #742) to run them."
    # without the fork every cross-block measurement falls back, with a recorded reason
    e = estimate_resources([C.Measurement(P"ZZ", 1, [1, 2])], lay; hardware=hw, adapter_model=cm)
    @test e.per_measurement[1].cost_source == :analytic
    @test any(n -> occursin("not loaded", n), e.notes)
    @test_skip false     # marks the skipped construction tests in the summary
else
    # surface d = 5 × surface d = 5, Z̄Z̄: +13 qubits (aux edges + adapter width 5) and
    # +14 checks (4 bridge X checks, 10 vertex Z checks): 27 ancilla
    ac = adapter_cost(cm, [s5, s5], [1, 1], hw; layout=lay, qubits=[1, 2], pauli=P"ZZ_", n_data=3)
    @test ac.source == :constructed && ac.ancilla == 27 && ac.rounds == 5
    @test ac.alt_ancilla == analytic                  # cross-check kept side by side
    @test occursin("+13 qubits, +14 checks", ac.note)
    @test length(cm.cache) == 1
    adapter_cost(cm, [s5, s5], [1, 1], hw; layout=lay, qubits=[1, 2], pauli=P"ZZ_", n_data=3)
    @test length(cm.cache) == 1                       # cached
    # X-type parts are handled via the Hadamard-conjugated code
    acx = adapter_cost(cm, [s5, s5], [1, 1], hw; layout=lay, qubits=[1, 2], pauli=P"XZ_", n_data=3)
    @test acx.source == :constructed && length(cm.cache) == 2
    # Y-type parts fall back, with the reason
    acy = adapter_cost(cm, [s5, s5], [1, 1], hw; layout=lay, qubits=[1, 2], pauli=P"YZ_", n_data=3)
    @test acy.source == :analytic && occursin("Y-type", acy.note)
    # two gross blocks: the merged code keeps k₁ + k₂ − 1 logicals (checked inside the construction)
    bb = BBCodeModel(:gross; intra_ancilla=40, intra_rounds=12, aux_graph_qubits=w -> 3w)
    layb = pack_sequential(24, [bb, bb])
    acb = adapter_cost(ConstructedAdapterModel(cycle_time_factor=2.0), [bb, bb], [1, 1], hw; layout=layb,
                       qubits=[1, 13], pauli=embed(24, [1, 13], P"ZZ"), n_data=24)
    @test acb.source == :constructed && acb.ancilla > 0 && acb.d_joint == 12
    # a construction that throws (the fork's inter-code path needs equal-length supports:
    # gross Z̄₁ has weight 26, Z̄₂ weight 38) falls back with the error recorded, and the failure is cached
    layb2 = LogicalLayout([bb, bb], Dict(1 => 1, 2 => 2), Dict(1 => 1, 2 => 2))
    cmb = ConstructedAdapterModel(cycle_time_factor=2.0)
    acf = adapter_cost(cmb, [bb, bb], [1, 1], hw; layout=layb2, qubits=[1, 2], pauli=embed(2, [1, 2], P"ZZ"), n_data=2)
    @test acf.source == :analytic && occursin("construction failed", acf.note) && occursin("unequal lengths", acf.note)
    @test length(cmb.cache) == 1
    # in a full estimate the per-measurement table reports both numbers
    e = estimate_resources([C.Measurement(P"ZZ", 1, [1, 2])], lay; hardware=hw, adapter_model=cm)
    r = e.per_measurement[1]
    @test r.cost_source == :constructed && r.ancilla == 27 && r.analytic_ancilla == analytic
end
end
