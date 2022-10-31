using Test
using LinearAlgebra
using Logging
using QuantumPropagators
using QuantumPropagators.Generators
using QuantumPropagators.Controls
using QuantumPropagators: Generator, Operator
using QuantumControlBase.TestUtils: random_hermitian_matrix


@testset "Simple getcontrols" begin

    H₀ = random_hermitian_matrix(5, 1.0)
    H₁ = random_hermitian_matrix(5, 1.0)
    H₂ = random_hermitian_matrix(5, 1.0)

    @test length(getcontrols(H₀)) == 0

    H = (nothing, (nothing, nothing))
    @test_throws MethodError getcontrols(H)

    ϵ₁ = t -> t
    ϵ₂ = t -> 0.0
    H = (H₀, (H₁, ϵ₁), (H₂, ϵ₂))
    @test getcontrols(H) == (ϵ₁, ϵ₂)

    H = (H₀, (H₁, ϵ₁), (H₂, ϵ₁))
    @test getcontrols(H) == (ϵ₁,)

    u₁ = [0.1, 1.0]
    u₂ = [0.1, 2.0]
    H = (H₁, (H₂, u₁), (H₂, u₂))
    @test getcontrols(H) == (u₁, u₂)

end

@testset "Tuple evalcontrols" begin
    H₀ = random_hermitian_matrix(5, 1.0)
    H₁ = random_hermitian_matrix(5, 1.0)
    H₂ = random_hermitian_matrix(5, 1.0)
    ϵ₁(t) = 1.0
    ϵ₂(t) = 1.0
    H = (H₀, (H₁, ϵ₁), (H₂, ϵ₂))

    Op = evalcontrols(H, IdDict(ϵ₁ => 1.2, ϵ₂ => 1.3))
    @test Op isa Matrix{ComplexF64}
    @test norm(Op - (H₀ + 1.2 * H₁ + 1.3 * H₂)) < 1e-15

    evalcontrols!(Op, H, IdDict(ϵ₁ => 2.2, ϵ₂ => 2.3))
    @test norm(Op - (H₀ + 2.2 * H₁ + 2.3 * H₂)) < 1e-15

    H = ((H₁, ϵ₁), (H₂, ϵ₂))

    Op = evalcontrols(H, IdDict(ϵ₁ => 1.2, ϵ₂ => 1.3))
    @test Op isa Matrix{ComplexF64}
    @test norm(Op - (1.2 * H₁ + 1.3 * H₂)) < 1e-15

    evalcontrols!(Op, H, IdDict(ϵ₁ => 2.2, ϵ₂ => 2.3))
    @test norm(Op - (2.2 * H₁ + 2.3 * H₂)) < 1e-15

    H = ((H₁, ϵ₁), H₂)

    Op = evalcontrols(H, IdDict(ϵ₁ => 1.2, ϵ₂ => 1.3))
    @test Op isa Matrix{ComplexF64}
    @test norm(Op - (1.2 * H₁ + H₂)) < 1e-15

    evalcontrols!(Op, H, IdDict(ϵ₁ => 2.2, ϵ₂ => 2.3))
    @test norm(Op - (2.2 * H₁ + H₂)) < 1e-15

end

@testset "Generator evalcontrols" begin

    H₀ = random_hermitian_matrix(5, 1.0)
    H₁ = random_hermitian_matrix(5, 1.0)
    H₂ = random_hermitian_matrix(5, 1.0)
    ϵ₁ = t -> 1.0
    ϵ₂ = t -> 1.0

    H = hamiltonian(H₀, (H₁, ϵ₁), (H₂, ϵ₂))

    Op = evalcontrols(H, IdDict(ϵ₁ => 1.2, ϵ₂ => 1.3))
    @test Op isa Operator{Matrix{ComplexF64},Float64}
    @test norm(Array(Op) - (H₀ + 1.2 * H₁ + 1.3 * H₂)) < 1e-15

    evalcontrols!(Op, H, IdDict(ϵ₁ => 2.2, ϵ₂ => 2.3))
    @test norm(Array(Op) - (H₀ + 2.2 * H₁ + 2.3 * H₂)) < 1e-15

    H = hamiltonian((H₁, ϵ₁), (H₂, ϵ₂))

    Op = evalcontrols(H, IdDict(ϵ₁ => 1.2, ϵ₂ => 1.3))
    @test Op isa Operator{Matrix{ComplexF64},Float64}
    @test norm(Array(Op) - (1.2 * H₁ + 1.3 * H₂)) < 1e-15

    evalcontrols!(Op, H, IdDict(ϵ₁ => 2.2, ϵ₂ => 2.3))
    @test norm(Array(Op) - (2.2 * H₁ + 2.3 * H₂)) < 1e-15

end

@testset "Tuple substitute_controls" begin
    H₀ = random_hermitian_matrix(5, 1.0)
    H₁ = random_hermitian_matrix(5, 1.0)
    H₂ = random_hermitian_matrix(5, 1.0)
    ϵ₁ = t -> 1.0
    ϵ₂ = t -> 1.0
    H = (H₀, (H₁, ϵ₁), (H₂, ϵ₂))

    ϵ₁_new = t -> 2.0
    ϵ₂_new = t -> 2.0

    H_new = substitute_controls(H, IdDict(ϵ₁ => ϵ₁_new, ϵ₂ => ϵ₂_new))

    @test H_new isa Tuple

    H_new[1] ≡ H₀
    H_new[2][1] ≡ H₁
    H_new[2][2] ≡ ϵ₁_new
    H_new[3][1] ≡ H₂
    H_new[3][2] ≡ ϵ₂_new

end


@testset "Generator substitute_controls" begin

    H₀ = random_hermitian_matrix(5, 1.0)
    H₁ = random_hermitian_matrix(5, 1.0)
    H₂ = random_hermitian_matrix(5, 1.0)
    ϵ₁ = t -> 1.0
    ϵ₂ = t -> 1.0
    H = hamiltonian(H₀, (H₁, ϵ₁), (H₂, ϵ₂))

    ϵ₁_new = t -> 2.0
    ϵ₂_new = t -> 2.0

    @test substitute_controls(H₀, IdDict(ϵ₁ => ϵ₁_new, ϵ₂ => ϵ₂_new)) ≡ H₀

    @test H isa Generator
    H_new = substitute_controls(H, IdDict(ϵ₁ => ϵ₁_new, ϵ₂ => ϵ₂_new))

    @test H_new isa Generator
    @test H_new.ops[1] ≡ H₀
    @test H_new.ops[2] ≡ H₁
    @test H_new.ops[3] ≡ H₂
    @test H_new.amplitudes[1] ≡ ϵ₁_new
    @test H_new.amplitudes[2] ≡ ϵ₂_new

end

struct MySquareAmpl
    control::Function
end

struct MyScaledAmpl
    c::Number
    control::Function
end


function QuantumPropagators.Generators.substitute_controls(a::MySquareAmpl, controls_map)
    return MySquareAmpl(get(controls_map, a.control, a.control))
end

function QuantumPropagators.Generators.evalcontrols(a::MyScaledAmpl, vals_dict, args...)
    return a.c * evalcontrols(a.control, vals_dict, args...)
end


@testset "Nonlinear substitute_controls" begin

    H₀ = random_hermitian_matrix(5, 1.0)
    H₁ = random_hermitian_matrix(5, 1.0)
    H₂ = random_hermitian_matrix(5, 1.0)
    ϵ₁ = t -> 1.0
    ϵ₂ = t -> 1.0

    with_logger(NullLogger()) do

        H = hamiltonian(H₀, (H₁, MySquareAmpl(ϵ₁)), (H₂, MySquareAmpl(ϵ₂)))

        ϵ₁_new = t -> 2.0
        ϵ₂_new = t -> 2.0

        @test H isa Generator
        H_new = substitute_controls(H, IdDict(ϵ₁ => ϵ₁_new, ϵ₂ => ϵ₂_new))

        @test H_new isa Generator
        @test H_new.ops[1] ≡ H₀
        @test H_new.ops[2] ≡ H₁
        @test H_new.ops[3] ≡ H₂

        @test H_new.amplitudes[1] isa MySquareAmpl
        @test H_new.amplitudes[2] isa MySquareAmpl
        @test H_new.amplitudes[1].control ≡ ϵ₁_new
        @test H_new.amplitudes[2].control ≡ ϵ₂_new

    end

end
