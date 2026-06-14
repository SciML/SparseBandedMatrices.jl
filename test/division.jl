using SparseBandedMatrices, Random
using Test

@testset "Division" begin
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
