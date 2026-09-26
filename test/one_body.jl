@testitem "One-body weights" tags=[:one_body] begin
    weights = OrbitalRotationsTK.weights
    @test weights(NPL(), 3) == [1.0, 1.0, 1.0]
    @test weights(NPL(w=[1.0, 0.0, 2.0]), 3) == [1.0, 0.0, 2.0]
    @test weights(NPL(w=[1, 0, 2]), 3) isa Vector{Float64}
    @test_throws DimensionMismatch weights(NPL(w=[1.0, 0.0]), 3)
end
