@testitem "Verification V1: well-formedness" tags=[:verification, :v1] begin

using PBCCompiler
using PBCCompiler: Circuit, CircuitOp, preprocess_circuit
using PBCCompiler.Verification: wellformedness_issues, check_wellformed
using QuantumClifford: @P_str
const C = CircuitOp

# the single issue reported for `c`, asserting there is exactly one
function only_issue(c; kw...)
    iss = wellformedness_issues(Circuit(c); kw...)
    @test length(iss) == 1
    iss[1]
end

@testset "clean circuits" begin
    c = Circuit([C.ExpEighPiPauli(P"Z", [1]), C.PauliConditional(P"Z", [1], P"X", [2]),
                 C.Measurement(P"Z", 1, [1]), C.BitConditional(C.ExpHalfPiPauli(P"X", [2]), 1),
                 C.Measurement(P"-ZZ", 2, [2, 1])])
    @test isempty(check_wellformed(c))
    @test isempty(check_wellformed(c; nqubits=2, n_data=2))
end

@testset "each rule violated once, at the right index" begin
    m = C.Measurement(P"Z", 1, [1])
    w = only_issue([m, C.ExpQuatPiPauli(P"X", [0])])
    @test (w.index, w.rule) == (2, :qubit_range)
    w = only_issue([m, C.ExpQuatPiPauli(P"X", [3])]; nqubits=2)
    @test (w.index, w.rule) == (2, :qubit_range)
    w = only_issue([C.ExpQuatPiPauli(P"XX", [1, 1]), m])
    @test (w.index, w.rule) == (1, :duplicate_qubit)
    w = only_issue([m, C.ExpEighPiPauli(P"XZ", [1])])
    @test (w.index, w.rule) == (2, :pauli_length)
    w = only_issue([C.ExpHalfPiPauli(P"X", [1]), C.Measurement(P"iZ", 1, [1])])
    @test (w.index, w.rule) == (2, :non_hermitian)
    w = only_issue([m, C.PauliConditional(P"Z", [1], P"X", [1])])
    @test (w.index, w.rule) == (2, :control_target_overlap)
    w = only_issue([C.Measurement(P"Z", 0, [1])])
    @test (w.index, w.rule) == (1, :bit_range)
    w = only_issue([C.ExpQuatPiPauli(P"X", [1]), C.BitConditional(C.ExpHalfPiPauli(P"X", [1]), 1), m])
    @test (w.index, w.rule) == (2, :read_before_write)
    w = only_issue([m, C.ExpQuatPiPauli(P"X", [1]), C.Measurement(P"X", 1, [1])])
    @test (w.index, w.rule) == (3, :bit_rewrite)
    @test occursin("already written by op 1", w.message)
    w = only_issue([m, C.BitConditional(C.Measurement(P"X", 2, [1]), 1)])
    @test (w.index, w.rule) == (2, :conditional_measurement)
end

@testset "magic-qubit discipline (n_data = 1)" begin
    joint = C.Measurement(P"ZZ", 2, [1, 2])
    cons = C.Measurement(P"X", 3, [2])
    fin = C.Measurement(P"Z", 1, [1])
    @test isempty(wellformedness_issues(Circuit([joint, cons, fin]); n_data=1))
    w = only_issue([joint, C.ExpQuatPiPauli(P"Z", [2]), cons, fin]; n_data=1)
    @test (w.index, w.rule) == (2, :magic_nonmeasurement)
    w = only_issue([cons, fin]; n_data=1)
    @test (w.index, w.rule) == (1, :magic_not_prepared)
    w = only_issue([joint, fin]; n_data=1)
    @test (w.index, w.rule) == (1, :magic_not_consumed)
    w = only_issue([joint, cons, C.Measurement(P"Z", 4, [2]), fin]; n_data=1)
    @test (w.index, w.rule) == (3, :magic_reuse)
    w = only_issue([joint, C.Measurement(P"XZ", 4, [1, 2]), cons, fin]; n_data=1)
    @test (w.index, w.rule) == (3, :magic_multiple_joint)
    iss = wellformedness_issues(Circuit([C.Measurement(P"ZZZ", 2, [1, 2, 3]), C.Measurement(P"X", 3, [2]),
                                         C.Measurement(P"X", 4, [3]), fin]); n_data=1)
    @test [(w.index, w.rule) for w in iss] == [(1, :magic_multiple)]
end

@testset "compiler output" begin
    c = Circuit([C.PauliConditional(P"Z", [1], P"X", [2]), C.ExpEighPiPauli(P"Z", [2]),
                 C.ExpEighPiPauli(P"X", [1]), C.Measurement(P"Z", 1, [1]), C.Measurement(P"Z", 2, [2])])
    out = preprocess_circuit(copy(c))
    @test isempty(wellformedness_issues(out; n_data=2))
    # compiler bug B3: gadget bits collide with the program's own bits
    c = Circuit([C.Measurement(P"Z", 1, [1]), C.ExpEighPiPauli(P"X", [1]), C.Measurement(P"Z", 2, [1])])
    iss = wellformedness_issues(preprocess_circuit(copy(c)); n_data=1)
    @test any(w -> w.rule == :bit_rewrite, iss)
end

end
