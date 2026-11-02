@testitem "Estimation E6: hardware graphs and embedding overhead" tags=[:estimation, :e6] begin
using PBCCompiler.ResourceEstimation
using PBCCompiler: ResourceEstimation
using Graphs, Random

# generators
g = ResourceEstimation.grid(3, 4)
@test nv(g) == 12 && ne(g) == 3 * 3 + 2 * 4
h = heavy_hex(2)
@test is_connected(h) && maximum(degree(h)) <= 3
@test all(v -> degree(h, v) == 2, nv(h)-ne(h)÷2+1:nv(h))        # the heavy (edge) qubits have degree 2

# Tanner graphs
s3 = SurfaceCodeModel(3)
T = tanner_graph(s3)
@test nv(T) == 9 + 8 && ne(T) == 4 * 4 + 2 * 4                   # 4 weight-4 and 4 weight-2 checks

# all-to-all hardware has zero overhead
hw_a = HardwareModel(:all_to_all; t_1q=1e-8, t_2q=2e-8, t_meas=5e-7, t_reaction=1e-6, p_phys=1e-3)
o = embedding_overhead(s3, hw_a)
@test (o.extra_qubits, o.extra_2q_layers, o.non_native_edges) == (0, 0, 0)

# grid hardware: deterministic for a fixed seed, consistent bookkeeping
hw_g = HardwareModel(:grid; t_1q=1e-8, t_2q=2e-8, t_meas=5e-7, t_reaction=1e-6, p_phys=1e-3)
o1 = embedding_overhead(s3, hw_g; rng=MersenneTwister(5), iterations=5000)
o2 = embedding_overhead(s3, hw_g; rng=MersenneTwister(5), iterations=5000)
@test o1.embedding.placement == o2.embedding.placement && o1.extra_qubits == o2.extra_qubits
@test allunique(o1.embedding.placement)
@test o1.embedding.cost == sum(o1.embedding.edge_lengths) >= ne(T)
@test o1.non_native_edges == count(>(1), o1.embedding.edge_lengths)
@test o1.extra_qubits == sum(L - 1 for L in o1.embedding.edge_lengths if L > 1; init=0)
# annealing never does worse than the breadth-first start (0 iterations)
o0 = embedding_overhead(s3, hw_g; rng=MersenneTwister(5), iterations=0)
@test o1.embedding.cost <= o0.embedding.cost
# swap routing: no extra qubits, 3× the bridge layers
os = embedding_overhead(s3, hw_g; routing=:swap, rng=MersenneTwister(5), iterations=5000)
@test os.extra_qubits == 0 && os.extra_2q_layers == 3 * o1.extra_2q_layers

# a code that is natively a grid embeds with no overhead: a path graph on a line
line = path_graph(6)
emb = embed_tanner_graph(line, ResourceEstimation.grid(1, 6); iterations=2000, rng=MersenneTwister(1))
@test emb.cost == 5 && all(==(1), emb.edge_lengths)
# heavy-hex: BB gross code (288 nodes, degree 6) cannot be native on degree-3 hardware
bb = BBCodeModel(:gross; intra_ancilla=1, intra_rounds=1, aux_graph_qubits=identity)
hw_h = HardwareModel(:heavy_hex; t_1q=1e-8, t_2q=2e-8, t_meas=5e-7, t_reaction=1e-6, p_phys=1e-3)
ob = embedding_overhead(bb, hw_h; rng=MersenneTwister(2), iterations=2000)
@test ob.non_native_edges > 0 && ob.extra_qubits > 0
end
