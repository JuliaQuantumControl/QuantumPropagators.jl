using Test
using QuantumPropagators.Shapes


@testset "flattop-blackman" begin
    shape(v) = flattop(v, t₀=10, T=20, t_rise=2, func=:blackman)
    @test shape(9.9) == 0
    @test shape(10) < 1e-14
    @test shape(20) < 1e-14
    @test shape(20.1) == 0
    @test shape(15) == 1

    tlist = collect(range(0, 30, step=1.0))
    shape_vec = flattop.(tlist, t₀=10, T=20, t_rise=2, func=:blackman)
    @test shape_vec[11] == shape(10)
end


@testset "flattop-sinsq" begin
    shape(v) = flattop(v, t₀=10, T=20, t_rise=2, func=:sinsq)
    @test shape(9.9) == 0
    @test shape(10) < 1e-14
    @test shape(20) < 1e-14
    @test shape(20.1) == 0
    @test shape(15) == 1

    tlist = collect(range(0, 30, step=1.0))
    shape_vec = flattop.(tlist, t₀=10, T=20, t_rise=2, func=:sinsq)
    @test shape_vec[11] == shape(10)
end
