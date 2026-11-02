@testitem "Estimation E3/E4: classification, intra-block and analytic adapter costs (hand-computed)" tags=[:estimation, :e3, :e4] begin
using PBCCompiler
using PBCCompiler: CircuitOp, preprocess_circuit
using PBCCompiler.ResourceEstimation
using QuantumClifford: @P_str
const C = CircuitOp

hw = HardwareModel(:all_to_all; t_1q=1e-8, t_2q=2e-8, t_meas=5e-7, t_reaction=1e-6, p_phys=1e-3)
cyc = 3e-8 + 4 * 2e-8 + 5e-7                       # surface cycle = 6.1e-7 s
s3 = SurfaceCodeModel(3)                           # tile = 2·9 − 1 = 17 qubits
lay = LogicalLayout([s3, s3], Dict(1 => 1, 2 => 2), Dict(1 => 1, 2 => 1))
am = AnalyticAdapterModel(cycle_time_factor=2.0)

c = [C.Measurement(P"ZZ", 1, [1, 2]),              # cross-block
     C.Measurement(P"Z", 2, [1]),                  # intra
     C.BitConditional(C.ExpQuatPiPauli(P"Z", [1]), 1),  # conditional Clifford → PPM, reads bit 1
     C.BitConditional(C.ExpHalfPiPauli(P"X", [2]), 2)]  # Pauli frame
ops = classify_operations(c, lay, hw; n_data=2, adapter_model=am)
@test [o.class for o in ops] == [:cross, :intra, :intra, :free]
@test [o.kind for o in ops] == [:measurement, :measurement, :clifford, :frame]
x = ops[1]
@test x.blocks == [1, 2] && x.weights == [1, 1]
@test x.ancilla == 17 + 17 + (2 * 3 - 1)            # Σ aux tiles + (a−1)(2d−1) = 39
@test x.rounds == 3 && x.round_time ≈ 2 * cyc && x.d_joint == 3
@test (ops[2].ancilla, ops[2].rounds, ops[2].round_time) == (17, 3, cyc)
@test ops[3].reads == [1] && ops[2].writes == 2
@test weight_histogram(ops) == Dict(1 => 4)
@test_throws ArgumentError classify_operations(c, lay, hw; n_data=2)        # cross-block needs an adapter model

# a gadget: the joint measurement consumes a magic state, the magic measurement is magic-only
g = preprocess_circuit([C.ExpEighPiPauli(P"Z", [1]), C.Measurement(P"Z", 1, [1])])
ops = classify_operations(g, LogicalLayout([s3], Dict(1 => 1), Dict(1 => 1)), hw; n_data=1)
@test [o.class for o in ops] == [:intra, :magic, :intra, :free, :intra]
@test [o.consumes_magic for o in ops] == [true, false, false, false, false]
@test ops[2].rounds == 1 && ops[2].ancilla == 0

# three blocks: a − 1 = 2 joins, d_joint = min distance
s5 = SurfaceCodeModel(5)
lay3 = LogicalLayout([s3, s5, s5], Dict(1 => 1, 2 => 2, 3 => 3), Dict(1 => 1, 2 => 1, 3 => 1))
op = only(classify_operations([C.Measurement(P"XXX", 1, [1, 2, 3])], lay3, hw; n_data=3, adapter_model=am))
@test op.d_joint == 3 && op.ancilla == 17 + 49 + 49 + 2 * 5 && op.rounds == 3

# qLDPC block: in-block cost and auxiliary graph come from the model's required arguments
bb = BBCodeModel(:gross; intra_ancilla=40, intra_rounds=12, aux_graph_qubits=w -> 3w + 5)
layb = pack_sequential(24, [bb, bb])
op = only(classify_operations([C.Measurement(P"ZZZ", 1, [1, 2, 3])], layb, hw; n_data=24))
@test (op.class, op.ancilla, op.rounds) == (:intra, 40, 12)
op = only(classify_operations([C.Measurement(P"ZZZ", 1, [1, 2, 13])], layb, hw; n_data=24, adapter_model=am))
@test op.weights == [2, 1] && op.ancilla == (3 * 2 + 5) + (3 * 1 + 5) + (2 * 12 - 1)
@test op.round_time ≈ 2.0 * (7 * 2e-8 + 1e-8 + 5e-7)

# the sourced defaults can be overridden
am2 = AnalyticAdapterModel(cycle_time_factor=1.5, adapter_qubits_per_join=d -> 3d, rounds=d -> 2d)
c2 = adapter_cost(am2, [s3, s5], [1, 1], hw)
@test (c2.ancilla, c2.rounds, c2.round_time) == (17 + 49 + 9, 6, 1.5 * cyc)
end
