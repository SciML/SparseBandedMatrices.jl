using SparseBandedMatrices, BenchmarkTools
using StableRNGs, LinearAlgebra

const SUITE = BenchmarkGroup()
const rng = StableRNG(123)

N = 1000
x = rand(rng, N)
y = zeros(N)

# =============================================================================
# Construction and diagonal writes
# =============================================================================

SUITE["construct"] = BenchmarkGroup()

SUITE["construct"]["undef"] = @benchmarkable SparseBandedMatrix{Float64}(
    undef, $N, $N
)

diagvals = rand(rng, N)
A = SparseBandedMatrix{Float64}(undef, N, N)

SUITE["construct"]["setdiagonal_main"] = @benchmarkable setdiagonal!(
    A, $diagvals, false
) setup = (A = SparseBandedMatrix{Float64}(undef, $N, $N))
SUITE["construct"]["setdiagonal_upper"] = @benchmarkable setdiagonal!(
    A, $(rand(rng, N - 1)), true
) setup = (A = SparseBandedMatrix{Float64}(undef, $N, $N))

# =============================================================================
# Element access and matvec
# =============================================================================

SUITE["ops"] = BenchmarkGroup()

SUITE["ops"]["setindex!"] = @benchmarkable $A[100, 200] = 3.0
SUITE["ops"]["getindex"] = @benchmarkable $A[100, 200]
SUITE["ops"]["matvec"] = @benchmarkable mul!($y, $A, $x)
SUITE["ops"]["matvec_oop"] = @benchmarkable $A * $x
