using SafeTestsets, Test

@testset "SparseBandedMatrices" begin
    @safetestset "Quality Assurance" include("qa.jl")
    @safetestset "JET Static Analysis" include("jet.jl")

    @safetestset "Constructors" begin
        using SparseBandedMatrices

        A = SparseBandedMatrix{Float64}(undef, 5, 5)
        A[1, 1] = 2
        @test A[1, 1] == 2.0
        A[4, 1] = 0
        @test A[4, 1] == 0.0
        A[1, 3] = 5
        @test A[1, 3] == 5.0

        @test size(A) == (5, 5)
    end

    @safetestset "Multiplication" begin
        using SparseBandedMatrices, Random
        dim = 5000
        x = rand(10:75)
        diag_vals = Vector{Vector{Float64}}(undef, x)
        diag_locs = randperm(dim * 2 - 1)[1:x]
        for j in 1:x
            diag_vals[j] = rand(min(diag_locs[j], 2 * dim - diag_locs[j]))
        end

        x_butterfly = SparseBandedMatrix{Float64}(diag_locs, diag_vals, dim, dim)
        x_dense = copy(x_butterfly)

        y = rand(dim, dim)
        z = zeros(dim, dim)

        @test isapprox(x_dense * y, x_butterfly * y)
        @test isapprox(y * x_dense, y * x_butterfly)

        y = rand(10:75)
        diag_vals = Vector{Vector{Float64}}(undef, y)
        diag_locs = randperm(dim * 2 - 1)[1:y]
        for j in 1:y
            diag_vals[j] = rand(min(diag_locs[j], 2 * dim - diag_locs[j]))
        end

        y_butterfly = SparseBandedMatrix{Float64}(diag_locs, diag_vals, dim, dim)
        y_dense = copy(y_butterfly)

        @test isapprox(x_butterfly * y_butterfly, x_dense * y_dense)
    end

    @safetestset "Division" begin
        using SparseBandedMatrices, Random
        dim = 5000
        x = rand(10:75)
        diag_vals = Vector{Vector{Float64}}(undef, x)
        diag_locs = randperm(dim * 2 - 1)[1:x]
        for j in 1:x
            diag_vals[j] = rand(min(diag_locs[j], 2 * dim - diag_locs[j]))
        end

        x_butterfly = SparseBandedMatrix{Float64}(diag_locs, diag_vals, dim, dim)
        x_dense = copy(x_butterfly)

        y = rand(dim, dim)
        z = zeros(dim, dim)

        @test isapprox(x_dense / y, x_butterfly / y)
        @test isapprox(y / x_dense, y / x_butterfly)

        y = rand(10:75)
        diag_vals = Vector{Vector{Float64}}(undef, y)
        diag_locs = randperm(dim * 2 - 1)[1:y]
        for j in 1:y
            diag_vals[j] = rand(min(diag_locs[j], 2 * dim - diag_locs[j]))
        end

        y_butterfly = SparseBandedMatrix{Float64}(diag_locs, diag_vals, dim, dim)
        y_dense = copy(y_butterfly)

        @test isapprox(x_butterfly / y_butterfly, x_dense / y_dense)
    end
end
