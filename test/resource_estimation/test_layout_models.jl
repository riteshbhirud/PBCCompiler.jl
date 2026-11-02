@testitem "Estimation E1: layout validation and packing" tags=[:estimation, :e1] begin
using PBCCompiler
using PBCCompiler: CircuitOp
using PBCCompiler.ResourceEstimation
using QuantumClifford: @P_str
const C = CircuitOp

s = SurfaceCodeModel(3)
bb = BBCodeModel(:gross; intra_ancilla=10, intra_rounds=12, aux_graph_qubits=w -> 3w)
# validation: each rule violated once
lay = LogicalLayout([s, s], Dict(1 => 1, 2 => 2), Dict(1 => 1, 2 => 1))
@test isempty(validate_layout(lay, 1:2))
@test only(validate_layout(lay, 1:3)) == "logical qubit 3 is not mapped to a block"
@test occursin("share slot", only(validate_layout(LogicalLayout([s, s], Dict(1 => 1, 2 => 1), Dict(1 => 1, 2 => 1)), 1:2)))
@test occursin("encodes only 1", only(validate_layout(LogicalLayout([s], Dict(1 => 1), Dict(1 => 2)), 1:1)))
@test occursin("there are 1 blocks", only(validate_layout(LogicalLayout([s], Dict(1 => 3), Dict(1 => 1)), 1:1)))

# sequential packing
lay = pack_sequential(14, [bb, bb])
@test isempty(validate_layout(lay, 1:14))
@test lay.qubit_to_block[12] == 1 && lay.qubit_to_block[13] == 2 && lay.qubit_within_block[13] == 1
@test_throws ArgumentError pack_sequential(3, [s, s])

# interaction packing: two tightly coupled pairs (1,3) and (2,4), weakly coupled across
c = vcat([C.Measurement(P"ZZ", i, [1, 3]) for i in 1:5], [C.Measurement(P"XX", 5 + i, [2, 4]) for i in 1:5],
         [C.Measurement(P"ZZ", 11, [1, 2])])
two = SurfaceCodeModel(3)
gen2 = [BBCodeModel(:gross; intra_ancilla=1, intra_rounds=1, aux_graph_qubits=identity) for _ in 1:2]
# blocks of capacity 2 via a small generic model is not needed: use two gross blocks but only 4 qubits
lay, w = pack_by_interaction(c, gen2; n_data=4)
@test isempty(validate_layout(lay, 1:4))
# capacity 12 each: everything fits in one block, so no cross-block weight remains
@test w == 0
# force capacity 2 per block with 2-qubit generic blocks
using QuantumClifford: Stabilizer
struct TwoSlot <: CodeModel end
PBCCompiler.ResourceEstimation.code_k(::TwoSlot) = 2
lay, w = pack_by_interaction(c, [TwoSlot(), TwoSlot()]; n_data=4)
@test w == 1                                   # only the single (1,2) measurement crosses
@test lay.qubit_to_block[1] == lay.qubit_to_block[3] && lay.qubit_to_block[2] == lay.qubit_to_block[4]
@test interaction_graph(c, 4) == Dict((1, 3) => 5, (2, 4) => 5, (1, 2) => 1)
end

@testitem "Estimation E2: code models" tags=[:estimation, :e2] begin
using PBCCompiler.ResourceEstimation
using PBCCompiler.ResourceEstimation: pl
hw = HardwareModel(:all_to_all; t_1q=1e-8, t_2q=2e-8, t_meas=5e-7, t_reaction=1e-6, p_phys=1e-3)

# surface code
s = SurfaceCodeModel(5)
@test (code_n(s), code_k(s), code_d(s), syndrome_ancillas(s)) == (25, 1, 5, 24)
@test cycle_time(s, hw) ≈ 3e-8 + 8e-8 + 5e-7                     # [Fowler] 8-step cycle
@test logical_error_per_cycle(s, 1e-3) ≈ 0.1 * 0.1^3             # [GoSC] Eq. (10)
@test intra_measurement_cost(s, 2) == (2 * 49, 5)
@test_throws ArgumentError SurfaceCodeModel(4)

# BB presets, with the published fit reproducing Table 1 of Bravyi et al.
g = BBCodeModel(:gross; intra_ancilla=30, intra_rounds=12, aux_graph_qubits=w -> 3w + 5)
tg = BBCodeModel(:two_gross; intra_ancilla=30, intra_rounds=18, aux_graph_qubits=w -> 3w)
@test (code_n(g), code_k(g), code_d(g), syndrome_ancillas(g)) == (144, 12, 12, 144)
@test (code_n(tg), code_k(tg), code_d(tg)) == (288, 12, 18)
@test 1.5e-7 < logical_error_per_cycle(g, 1e-3) < 2.5e-7          # Table 1: 2 × 10⁻⁷
@test 5e-13 < logical_error_per_cycle(g, 1e-4) < 1.2e-12          # Table 1: 8 × 10⁻¹³
@test 1.5e-12 < logical_error_per_cycle(tg, 1e-3) < 2.5e-12       # Table 1: 2 × 10⁻¹²
@test cycle_time(g, hw) ≈ 7 * 2e-8 + 1e-8 + 5e-7
@test aux_graph_qubits(g, 14) == 47
# code parameters follow from the explicit check matrices
for m in (g, tg, SurfaceCodeModel(3), SurfaceCodeModel(9))
    Hx, Hz = parity_check_matrices(m)
    n = size(Hx, 2)
    @test n == code_n(m)
    @test n - gf2_rank(Hx) - gf2_rank(Hz) == code_k(m)
    @test all(iszero, (Int.(Hx) * transpose(Int.(Hz))) .% 2)
end

# p_L models
@test pl(AnalyticPL(A=0.03, p_th=0.01), 1e-3, 3) ≈ 0.03 * 0.01
@test pl(FitPL((p, d) -> p * d), 0.5, 4) == 2.0
path = tempname() * ".csv"
write(path, "code,d,p,p_L,ci_low,ci_high\nsurf,3,0.001,1e-4,0,0\nsurf,3,0.01,1e-2,0,0\nother,3,0.001,1,0,0\n")
t = TablePL(path; code="surf")
@test pl(t, 0.001, 3) ≈ 1e-4 && pl(t, 0.01, 3) ≈ 1e-2
@test pl(t, sqrt(0.001 * 0.01), 3) ≈ 1e-3                         # log-log interpolation
@test_throws ArgumentError pl(t, 0.1, 3)
@test_throws ArgumentError pl(t, 0.001, 5)
rm(path)

# generic CSS wrapper reads n, k from a QuantumClifford code
using QuantumClifford.ECC: Steane7
m = GenericCSSModel(Steane7(); d=3, pl=AnalyticPL(), cycle_layers=(2, 4, 1), intra_ancilla=7, intra_rounds=3,
                    aux_graph_qubits=w -> 2w)
@test (code_n(m), code_k(m), code_d(m), syndrome_ancillas(m)) == (7, 1, 3, 6)
Hx, Hz = parity_check_matrices(m)
@test 7 - gf2_rank(Hx) - gf2_rank(Hz) == 1
@test_throws ArgumentError GenericCSSModel(Steane7(); pl=AnalyticPL(), cycle_layers=(1, 1, 1), intra_ancilla=1,
                                           intra_rounds=1, aux_graph_qubits=identity)
end
