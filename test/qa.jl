using SparseBandedMatrices, Aqua
@testset "Aqua" begin
    Aqua.test_all(SparseBandedMatrices; ambiguities = (recursive = false,))
end
