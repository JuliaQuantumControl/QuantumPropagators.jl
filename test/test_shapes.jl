using Test
using QuantumPropagators.Shapes


@testset "blackman" begin
    ≈(a, b) = isapprox(a, b; atol=1e-15)
    @test blackman(-1.0, 0.0, 10.0) == 0.0
    @test blackman(0.0, 0.0, 10.0) ≈ 0.0
    @test blackman(10.0, 0.0, 10.0) ≈ 0.0
    @test blackman(11.0, 0.0, 10.0) == 0.0
    @test blackman(5.0, 0.0, 10.0) ≈ 1.0
end


@testset "flattop-blackman" begin
    ≈(a, b) = isapprox(a, b; atol=1e-14)
    shape(v) = flattop(v, t₀=10, T=20, t_rise=2, func=:blackman)
    @test shape(9.9) == 0
    @test shape(10) ≈ 0.0
    @test shape(12) ≈ 1.0
    @test shape(18) ≈ 1.0
    @test shape(20) ≈ 0.0
    @test shape(20.1) == 0
    @test shape(15) == 1

    tlist = collect(range(0, 30, step=1.0))
    shape_vec = flattop.(tlist, t₀=10, T=20, t_rise=2, func=:blackman)
    @test shape_vec[11] == shape(10)
end


@testset "flattop-sinsq" begin
    ≈(a, b) = isapprox(a, b; atol=1e-14)
    shape(v) = flattop(v, t₀=10, T=20, t_rise=2, func=:sinsq)
    @test shape(9.9) == 0.0
    @test shape(10) ≈ 0.0
    @test shape(12) ≈ 1.0
    @test shape(18) ≈ 1.0
    @test shape(20) ≈ 0.0
    @test shape(20.1) == 0
    @test shape(15) == 1

    tlist = collect(range(0, 30, step=1.0))
    shape_vec = flattop.(tlist, t₀=10, T=20, t_rise=2, func=:sinsq)
    @test shape_vec[11] == shape(10)
end


@testset "flattop-box" begin
    # Test that a flattop with a zero rise time is equivalent to a box, and
    # specifically that it does not return NaN (#42)
    ≈(a, b) = isapprox(a, b; atol=1e-15)
    for func in [:blackman, :sinqu]
        for v ∈ [-1.0, 0.0, 0.5, 1.0, 2.0]
            f = flattop(v; T=1.0, t₀=0.0, t_rise=0.0)
            @test !isnan(f)
            @test f ≈ box(v, 0.0, 1.0)
        end
    end
end
