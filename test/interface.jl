using SparseBandedMatrices, LinearAlgebra, Test

@testset "BigFloat support" begin
    dim = 5
    diag_locs = [1, 3, 5, 7, 9]
    diag_vals = Vector{Vector{BigFloat}}(undef, length(diag_locs))
    for (j, loc) in enumerate(diag_locs)
        len = min(loc, 2 * dim - loc)
        diag_vals[j] = BigFloat.(rand(len))
    end

    A = SparseBandedMatrix{BigFloat}(diag_locs, diag_vals, dim, dim)
    @test eltype(A) == BigFloat
    @test size(A) == (dim, dim)

    # Test getindex returns BigFloat
    @test A[1, 1] isa BigFloat
    @test A[2, 3] isa BigFloat  # zero element

    # Test setindex! with BigFloat
    A[1, 1] = big"42.0"
    @test A[1, 1] == big"42.0"

    # Test mul! with BigFloat matrices
    B = Matrix{BigFloat}(ones(BigFloat, dim, 3))
    C = Matrix{BigFloat}(zeros(BigFloat, dim, 3))
    mul!(C, A, B, big"1.0", big"0.0")
    @test eltype(C) == BigFloat

    # Verify result
    A_dense = zeros(BigFloat, dim, dim)
    for i in 1:dim, j in 1:dim
        A_dense[i, j] = A[i, j]
    end
    expected = A_dense * B
    @test isapprox(C, expected)

    # Test mul! with Matrix * SparseBandedMatrix
    B2 = Matrix{BigFloat}(ones(BigFloat, 3, dim))
    C2 = Matrix{BigFloat}(zeros(BigFloat, 3, dim))
    mul!(C2, B2, A, big"1.0", big"0.0")
    expected2 = B2 * A_dense
    @test isapprox(C2, expected2)
end

@testset "ComplexF64 support" begin
    A = SparseBandedMatrix{ComplexF64}(undef, 5, 5)
    A[1, 1] = 2.0 + 1.0im
    A[2, 2] = 3.0 - 2.0im
    A[3, 3] = 4.0 + 0.5im

    @test eltype(A) == ComplexF64
    @test A[1, 1] == 2.0 + 1.0im
    @test A[4, 4] == 0.0 + 0.0im  # zero element

    B = Matrix{ComplexF64}(ones(ComplexF64, 5, 3))
    C = Matrix{ComplexF64}(zeros(ComplexF64, 5, 3))

    mul!(C, A, B, 1.0 + 0im, 0.0 + 0im)
    @test eltype(C) == ComplexF64

    # Verify result
    A_dense = zeros(ComplexF64, 5, 5)
    for i in 1:5, j in 1:5
        A_dense[i, j] = A[i, j]
    end
    expected = A_dense * B
    @test isapprox(C, expected)

    # Test mul! with Matrix * SparseBandedMatrix
    B2 = Matrix{ComplexF64}(ones(ComplexF64, 3, 5))
    C2 = Matrix{ComplexF64}(zeros(ComplexF64, 3, 5))
    mul!(C2, B2, A, 1.0 + 0im, 0.0 + 0im)
    expected2 = B2 * A_dense
    @test isapprox(C2, expected2)
end

@testset "Float32 support" begin
    A = SparseBandedMatrix{Float32}(undef, 5, 5)
    A[1, 1] = 2.0f0
    A[2, 2] = 3.0f0
    A[3, 3] = 4.0f0

    @test eltype(A) == Float32
    @test A[1, 1] == 2.0f0

    B = Matrix{Float32}(ones(Float32, 5, 3))
    C = Matrix{Float32}(zeros(Float32, 5, 3))

    mul!(C, A, B, 1.0f0, 0.0f0)
    @test eltype(C) == Float32

    # Verify result
    A_dense = zeros(Float32, 5, 5)
    for i in 1:5, j in 1:5
        A_dense[i, j] = A[i, j]
    end
    expected = A_dense * B
    @test isapprox(C, expected)
end

@testset "AbstractArray interface" begin
    A = SparseBandedMatrix{Float64}(undef, 5, 5)
    A[1, 1] = 2.0
    A[2, 2] = 3.0
    A[3, 3] = 4.0

    # Test size
    @test size(A) == (5, 5)
    @test size(A, 1) == 5
    @test size(A, 2) == 5

    # Test length
    @test length(A) == 25

    # Test eltype
    @test eltype(A) == Float64

    # Test firstindex/lastindex
    @test firstindex(A) == 1
    @test lastindex(A) == 25

    # Test IndexStyle
    @test Base.IndexStyle(typeof(A)) == Base.IndexCartesian()

    # Test iteration
    count = 0
    for _ in A
        count += 1
    end
    @test count == 25
end
