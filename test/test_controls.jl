using Test
using QuantumPropagators.Controls

@testset "simple getcontrols" begin

    H = (nothing, (nothing, nothing))
    @test_throws TypeError getcontrols(H)

    ϵ₁ = t -> t
    ϵ₂ = t -> 0.0
    H = (nothing, (nothing, ϵ₁), (nothing, ϵ₂))
    @test getcontrols(H) == (ϵ₁, ϵ₂)

    u₁ = [0.1, 1.0]
    u₂ = [0.1, 2.0]
    H = (nothing, (nothing, u₁), (nothing, u₂))
    @test getcontrols(H) == (u₁, u₂)

end
