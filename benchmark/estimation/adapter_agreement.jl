# E4 vs E5: analytic vs constructed adapter costs on several code pairs.
#
# Needs the universal-adapter construction (riteshbhirud/QuantumClifford.jl, master,
# PR #742) in the active environment:
#
#     julia -tauto --project=<env with the fork and PBCCompiler> benchmark/estimation/adapter_agreement.jl
#
# Prints one row per pair: the constructed ancilla (new qubits + new checks), the
# analytic join term alone (2d − 1, arXiv:2410.03628), and the analytic total (which
# adds each block's auxiliary-graph assumption: one routing tile for surface codes,
# and for BB codes the illustrative user choice `aux_graph_qubits = w -> 0`, so that
# only the join term is compared).

using PBCCompiler
using PBCCompiler.ResourceEstimation
using QuantumClifford: @P_str, embed, PauliOperator
using Printf

adapters_available() || error("the adapter construction is not loaded; see the header")
hw = HardwareModel(:all_to_all; t_1q=1e-8, t_2q=2e-8, t_meas=5e-7, t_reaction=1e-6, p_phys=1e-3)
bb(name) = BBCodeModel(name; intra_ancilla=0, intra_rounds=0, aux_graph_qubits=w -> 0)
pairs = Any[]
for d in (3, 5, 7, 9, 11)
    push!(pairs, ("surface d=$d × d=$d, ZZ", SurfaceCodeModel(d), SurfaceCodeModel(d), (1, 'Z'), (1, 'Z')))
    push!(pairs, ("surface d=$d × d=$d, XZ", SurfaceCodeModel(d), SurfaceCodeModel(d), (1, 'X'), (1, 'Z')))
end
push!(pairs, ("gross × gross, Z̄₁Z̄₁", bb(:gross), bb(:gross), (1, 'Z'), (1, 'Z')))
push!(pairs, ("gross × gross, Z̄₁Z̄₂", bb(:gross), bb(:gross), (1, 'Z'), (2, 'Z')))
push!(pairs, ("gross × two-gross, Z̄₁Z̄₁", bb(:gross), bb(:two_gross), (1, 'Z'), (1, 'Z')))
push!(pairs, ("two-gross × two-gross, Z̄₁Z̄₁", bb(:two_gross), bb(:two_gross), (1, 'Z'), (1, 'Z')))

cm = ConstructedAdapterModel(cycle_time_factor=1.0)
println("| pair | d_joint | constructed ancilla | analytic join (2d−1) | analytic total | construction note |")
println("|---|---|---|---|---|---|")
for (name, A, B, (ia, la), (ib, lb)) in pairs
    ka, kb = code_k(A), code_k(B)
    lay = LogicalLayout([A, B], Dict(1 => 1, 2 => 2), Dict(1 => ia, 2 => ib))
    P = embed(2, [1, 2], PauliOperator(0x00, [la == 'X', lb == 'X'], [la == 'Z', lb == 'Z']))
    t = @elapsed ac = adapter_cost(cm, [A, B], [1, 1], hw; layout=lay, qubits=[1, 2], pauli=P, n_data=2)
    dj = ac.d_joint
    @printf("| %s | %d | %d | %d | %d | %s (%.2f s) |\n", name, dj, ac.ancilla, 2dj - 1, ac.alt_ancilla,
            split(ac.note, ";")[1], t)
end
