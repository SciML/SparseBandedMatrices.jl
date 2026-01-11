using Test
using SparseBandedMatrices
using LinearAlgebra

@testset "BigFloat Support" begin
    # Test 1: Basic constructor
    @testset "Constructor" begin
        A = SparseBandedMatrix{BigFloat}(undef, 5, 5)
        @test size(A) == (5, 5)
        @test eltype(A) == BigFloat
    end

    # Test 2: setindex!/getindex
    @testset "Indexing" begin
        A = SparseBandedMatrix{BigFloat}(undef, 5, 5)
        A[1, 1] = BigFloat(2.0)
        A[2, 3] = BigFloat(3.14)
        @test A[1, 1] == BigFloat(2.0)
        @test typeof(A[1, 1]) == BigFloat
        @test A[3, 3] == BigFloat(0.0)  # Unset element returns zero
        @test typeof(A[3, 3]) == BigFloat
    end

    # Test 3: setdiagonal!
    @testset "setdiagonal!" begin
        A = SparseBandedMatrix{BigFloat}(undef, 5, 5)
        diagvals = BigFloat[1.0, 2.0, 3.0]
        setdiagonal!(A, diagvals, true)
        setdiagonal!(A, BigFloat[4.0, 5.0], false)
        # Verify the diagonals were set (checking actual values)
        @test A[3, 1] == BigFloat(1.0)  # first element of lower diagonal
        @test A[4, 2] == BigFloat(2.0)  # second element of lower diagonal
        @test A[5, 3] == BigFloat(3.0)  # third element of lower diagonal
        @test A[1, 4] == BigFloat(4.0)  # first element of upper diagonal
        @test A[2, 5] == BigFloat(5.0)  # second element of upper diagonal
    end

    # Test 4: Constructor with diagonal values
    @testset "Constructor with diagonals" begin
        ind_vals = [2, 7]
        diag_vals = [BigFloat[1.0, 2.0], BigFloat[3.0, 4.0, 5.0]]
        B = SparseBandedMatrix{BigFloat}(ind_vals, diag_vals, 5, 5)
        @test eltype(B) == BigFloat
    end

    # Test 5: mul! with Matrix
    @testset "mul! Matrix operations" begin
        A = SparseBandedMatrix{BigFloat}(undef, 3, 3)
        A[1, 1] = BigFloat(1.0)
        A[2, 2] = BigFloat(2.0)
        A[3, 3] = BigFloat(3.0)

        B = ones(BigFloat, 3, 2)
        C = zeros(BigFloat, 3, 2)

        mul!(C, A, B, BigFloat(1.0), BigFloat(0.0))
        expected = Matrix(A) * B
        @test isapprox(C, expected)
    end

    # Test 6: * operator (this was the key bug fix)
    @testset "* operator" begin
        A = SparseBandedMatrix{BigFloat}(undef, 3, 3)
        A[1, 1] = BigFloat(1.0)
        A[2, 2] = BigFloat(2.0)
        A[3, 3] = BigFloat(3.0)

        B = ones(BigFloat, 3, 2)
        C = A * B
        expected = Matrix(A) * B

        @test eltype(C) == BigFloat
        @test isapprox(C, expected)
    end

    # Test 7: Matrix * SparseBandedMatrix
    @testset "Matrix * SparseBandedMatrix" begin
        A = SparseBandedMatrix{BigFloat}(undef, 3, 3)
        A[1, 1] = BigFloat(1.0)
        A[2, 2] = BigFloat(2.0)
        A[3, 3] = BigFloat(3.0)

        B = ones(BigFloat, 2, 3)
        C = B * A
        expected = B * Matrix(A)

        @test eltype(C) == BigFloat
        @test isapprox(C, expected)
    end

    # Test 8: SparseBandedMatrix * SparseBandedMatrix
    @testset "SparseBandedMatrix * SparseBandedMatrix" begin
        A = SparseBandedMatrix{BigFloat}(undef, 3, 3)
        A[1, 1] = BigFloat(1.0)
        A[2, 2] = BigFloat(2.0)

        B = SparseBandedMatrix{BigFloat}(undef, 3, 3)
        B[1, 1] = BigFloat(2.0)
        B[3, 3] = BigFloat(3.0)

        C = A * B
        expected = Matrix(A) * Matrix(B)
        @test isapprox(C, expected)
    end
end

# Note: ComplexF64 is not fully supported because fma() is not defined for complex numbers.
# This is a known limitation of the current implementation which uses fma for performance.

@testset "AbstractArray Interface" begin
    A = SparseBandedMatrix{Float64}(undef, 5, 5)
    A[1, 1] = 1.0
    A[2, 2] = 2.0
    A[3, 3] = 3.0
    A[1, 3] = 0.5

    @testset "size" begin
        @test size(A) == (5, 5)
        @test size(A, 1) == 5
        @test size(A, 2) == 5
    end

    @testset "length" begin
        @test length(A) == 25
    end

    @testset "eltype" begin
        @test eltype(A) == Float64
    end

    @testset "axes" begin
        @test axes(A) == (Base.OneTo(5), Base.OneTo(5))
    end

    @testset "firstindex/lastindex" begin
        @test firstindex(A) == 1
        @test lastindex(A) == 25
    end

    @testset "iteration" begin
        vals = collect(A)
        @test length(vals) == 25
        @test vals[1] == 1.0  # A[1,1]
    end

    @testset "Matrix conversion" begin
        B = Matrix(A)
        @test typeof(B) == Matrix{Float64}
        @test all(A[i, j] == B[i, j] for i in 1:5, j in 1:5)
    end

    @testset "IndexStyle" begin
        @test Base.IndexStyle(typeof(A)) == IndexCartesian()
    end
end
