# Reference comparison: surface-code-only estimate vs Litinski, "A Game of Surface
# Codes" (Quantum 3, 128, 2019), Sec. 4, minimal setup at p = 1e-4:
#   100 logical qubits, 10^8 T gates, compact data block (153 tiles) + one 15-to-1
#   block (11 tiles) = 164 tiles × 2d², d = 13 → ≈ 55,400 physical qubits,
#   13·11·10^8 code cycles ≈ 4 h at 1 µs per code cycle, final error ≈ 0.2%.
#
#     julia -tauto --project=test benchmark/estimation/gosc_comparison.jl
#
# A 10^8-T circuit cannot be built explicitly. We estimate programs with N random
# π/8 Pauli-product rotations on 100 qubits (compiled by PBCCompiler), check that time
# and error per T are independent of N, and extrapolate linearly to 10^8. Settings
# mirror GoSC's assumptions (each one is printed with the result):
#   * one surface-code patch per logical qubit (d = 13), p = 1e-4, 1 µs code cycle
#     (each of the 8 cycle steps is 1/8 µs), no classical reaction time;
#   * multi-patch Pauli products by lattice surgery: one routing tile (2d² − 1 qubits)
#     per involved patch and no adapter qubits (AnalyticAdapterModel with
#     adapter_qubits_per_join = 0, cycle_time_factor = 1);
#   * one GoSC 15-to-1 factory (11 tiles, 11d cycles per state, 35p³).
# Nothing is tuned to match the reference.

using PBCCompiler
using PBCCompiler: CircuitOp, preprocess_circuit
using PBCCompiler.ResourceEstimation
using QuantumClifford: PauliOperator
using Random, Printf

const NQ = 100
const D = 13
hw = HardwareModel(:all_to_all; t_1q=1e-6 / 8, t_2q=1e-6 / 8, t_meas=1e-6 / 8, t_reaction=0.0, p_phys=1e-4)
am = AnalyticAdapterModel(cycle_time_factor=1.0, adapter_qubits_per_join=d -> 0)

function program(NT, rng)
    c = CircuitOp.Type[]
    for _ in 1:NT
        x = rand(rng, Bool, NQ); z = rand(rng, Bool, NQ)
        any(x .| z) || (z[1] = true)
        push!(c, CircuitOp.ExpEighPiPauli(PauliOperator(0x00, x, z), collect(1:NQ)))
    end
    append!(c, [CircuitOp.Measurement(PauliOperator(0x00, [false], [true]), q, [q]) for q in 1:NQ])
    preprocess_circuit(c)
end

layout(d) = LogicalLayout([SurfaceCodeModel(d) for _ in 1:NQ], Dict(q => q for q in 1:NQ), Dict(q => 1 for q in 1:NQ))

rows = []
for NT in (100, 200, 400), d in (11, 13)
    out = program(NT, MersenneTwister(NT))
    e = estimate_resources(out, layout(d); hardware=hw, adapter_model=am, factory=FifteenToOneFactory(d=d), max_factories=1)
    nT = e.magic_states[:pi8_autocorrected]
    push!(rows, (NT=nT, d=d, qubits=e.physical_qubits, q=e.qubit_breakdown, cyc_per_T=e.time_cycles / nT,
                 fail_per_T=e.failure_probability / nT, e=e))
end

println("Scaling check (per-T quantities should not depend on N):")
for r in rows
    @printf("  N_T=%4d d=%d  qubits=%6d  cycles/T=%7.1f  failure/T=%.3e\n", r.NT, r.d, r.qubits, r.cyc_per_T, r.fail_per_T)
end
println("\nExtrapolation to 10^8 T gates (from the largest N):")
for d in (11, 13)
    r = last(filter(x -> x.d == d, rows))
    cycles = r.cyc_per_T * 1e8
    hours = cycles * 1e-6 / 3600
    fail = r.fail_per_T * 1e8
    @printf("  d=%d: %d physical qubits (data %d, syndrome %d, surgery ancilla peak %d, factory %d), %.2e cycles = %.2f h, failure %.3g\n",
            d, r.qubits, r.q[:data], r.q[:syndrome], r.q[:surgery_ancilla_peak], r.q[:factories], cycles, hours, fail)
end
println("\nReference (GoSC Sec. 4.3, d = 13): 55,400 qubits, 13·11·10^8 = 1.43e10 cycles = 3.97 h, failure 0.2%",
        " (d = 11: 19.8%; their failure counts memory errors on all 164 tiles; the T-gate error budget is separate)")
