using SparseBandedMatrices, JET, LinearAlgebra

@testset "JET static analysis" begin
    # Test that there are no unresolved dispatch errors in key functions
    @testset "Error-free analysis" begin
        rep = JET.report_call(getindex, (SparseBandedMatrix{Float64}, Int, Int))
        @test length(JET.get_reports(rep)) == 0

        rep = JET.report_call(setindex!, (SparseBandedMatrix{Float64}, Float64, Int, Int))
        @test length(JET.get_reports(rep)) == 0

        rep = JET.report_call(size, (SparseBandedMatrix{Float64},))
        @test length(JET.get_reports(rep)) == 0

        rep = JET.report_call(setdiagonal!, (SparseBandedMatrix{Float64}, Vector{Float64}, Bool))
        @test length(JET.get_reports(rep)) == 0
    end

    @testset "Type stability for core operations" begin
        # Test getindex type stability
        rep = JET.report_opt(getindex, (SparseBandedMatrix{Float64}, Int, Int))
        @test length(JET.get_reports(rep)) == 0

        # Test size type stability
        rep = JET.report_opt(size, (SparseBandedMatrix{Float64},))
        @test length(JET.get_reports(rep)) == 0
    end
end
