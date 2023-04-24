using Test
using LinearAlgebra
using QuantumControlTestUtils.RandomObjects: random_matrix, random_state_vector
using QuantumControlTestUtils.Interfaces: check_operator

using QuantumPropagators: Generator, Operator, ScaledOperator


@testset "Operator mul!" begin

    H₀ = random_matrix(5; hermitian=true)
    H₁ = random_matrix(5; hermitian=true)
    H₂ = random_matrix(5; hermitian=true)

    Op = Operator([H₀, H₁, H₂], [2.1, 1.1])
    H = H₀ + 2.1 * H₁ + 1.1 * H₂

    @test norm(convert(Matrix{ComplexF64}, Op) - H) < 1e-12

    Ψ = random_state_vector(5)
    ϕ0 = random_state_vector(5)

    ϕ = copy(ϕ0)
    mul!(ϕ, Op, Ψ, true, false)
    ϕ_expected = H * Ψ
    @test norm(ϕ - ϕ_expected) < 1e-12

    ϕ = copy(ϕ0)
    mul!(ϕ, Op, Ψ, true, true)
    ϕ_expected = ϕ0 + H * Ψ
    @test norm(ϕ - ϕ_expected) < 1e-12

    ϕ = copy(ϕ0)
    mul!(ϕ, Op, Ψ, 2.0, true)
    ϕ_expected = ϕ0 + (2 * H) * Ψ
    @test norm(ϕ - ϕ_expected) < 1e-12

    ϕ = copy(ϕ0)
    mul!(ϕ, Op, Ψ, 2.0, 2.0)
    ϕ_expected = 2 * ϕ0 + (2 * H) * Ψ
    @test norm(ϕ - ϕ_expected) < 1e-12


end


@testset "Operator interface" begin
    H₀ = random_matrix(5; hermitian=true)
    H₁ = random_matrix(5; hermitian=true)
    H₂ = random_matrix(5; hermitian=true)
    Op = Operator([H₀, H₁, H₂], [2.1, 1.1])
    Ψ = random_state_vector(5)
    @test check_operator(Op; state=Ψ)
end


@testset "Operator copy" begin

    H₀ = random_matrix(5; hermitian=true)
    H₁ = random_matrix(5; hermitian=true)
    H₂ = random_matrix(5; hermitian=true)

    Op = Operator([H₀, H₁, H₂], [2.1, 1.1])
    Op2 = copy(Op)

    @test length(Op.ops) == length(Op2.ops)
    @test length(Op.coeffs) == length(Op2.coeffs)

    @test norm(Op.coeffs - Op2.coeffs) < 1e-14

    @test norm(Op2.ops[1] - H₀) < 1e-14
    @test norm(Op2.ops[2] - H₁) < 1e-14
    @test norm(Op2.ops[3] - H₂) < 1e-14

    @test Op.ops[1] ≡ H₀
    @test Op2.ops[1] ≢ H₀

    Op3 = Operator([H₁, H₂, H₀], [1.1, 2.1])
    copyto!(Op3, Op)
    @test norm(Op.coeffs - Op3.coeffs) < 1e-14
    @test norm(Op3.ops[1] - H₀) < 1e-14
    @test norm(Op3.ops[2] - H₁) < 1e-14
    @test norm(Op3.ops[3] - H₂) < 1e-14

end


@testset "Operator array conversion" begin

    H₀ = random_matrix(5; hermitian=true)
    H₁ = random_matrix(5; hermitian=true)
    H₂ = random_matrix(5; hermitian=true)

    Op = Operator([H₀, H₁, H₂], [2.1, 1.1])
    H = H₀ + 2.1 * H₁ + 1.1 * H₂
    @test ishermitian(Op) == ishermitian(H) == true
    @test eltype(Array(Op)) ≡ ComplexF64
    @test size(Op) == size(H)
    @test norm(Array{ComplexF64}(Op) - H) < 1e-12

    Op = Operator([H₁, H₂], [2.1, 1.1])
    H = 2.1 * H₁ + 1.1 * H₂
    @test eltype(Array(Op)) ≡ ComplexF64
    @test size(Op) == size(H)
    @test norm(Array{ComplexF64}(Op) - H) < 1e-12

    Op = Operator([H₀, H₁, H₂], [1.1,])
    H = H₀ + H₁ + 1.1 * H₂
    @test eltype(Array(Op)) ≡ ComplexF64
    @test size(Op) == size(H)
    @test norm(Array{ComplexF64}(Op) - H) < 1e-12

    Op = ScaledOperator(2, Operator([H₀, H₁, H₂], [2.1, 1.1]))
    @test repr(Op) == "ScaledOperator{Int64,Operator}(2, …)"
    H = 2 * (H₀ + 2.1 * H₁ + 1.1 * H₂)
    @test eltype(Array(Op)) ≡ ComplexF64
    @test size(Op) == size(H)
    @test norm(Array{ComplexF64}(Op) - H) < 1e-12

    Op = ScaledOperator(2, Operator([H₁, H₂], [2.1, 1.1]))
    H = 2 * (2.1 * H₁ + 1.1 * H₂)
    @test eltype(Array(Op)) ≡ ComplexF64
    @test ishermitian(Op) == ishermitian(H) == true
    @test size(Op) == size(H)
    @test norm(Array{ComplexF64}(Op) - H) < 1e-12

    Op = ScaledOperator(2, Operator([H₀, H₁, H₂], [1.1,]))
    H = 2 * (H₀ + H₁ + 1.1 * H₂)
    @test eltype(Array(Op)) ≡ ComplexF64
    @test size(Op) == size(H)
    @test norm(Array{ComplexF64}(Op) - H) < 1e-12

    Op = ScaledOperator(2, H₀)
    H = 2 * H₀
    @test eltype(Array(Op)) ≡ ComplexF64
    @test size(Op) == size(H)
    @test norm(Array{ComplexF64}(Op) - H) < 1e-12

end


@testset "Operator dot" begin

    H₀ = random_matrix(5; hermitian=true)
    H₁ = random_matrix(5; hermitian=true)
    H₂ = random_matrix(5; hermitian=true)

    Op = Operator([H₀, H₁, H₂], [2.1, 1.1])
    H = H₀ + 2.1 * H₁ + 1.1 * H₂

    Ψ = random_state_vector(5)
    ϕ = random_state_vector(5)

    c1 = dot(ϕ, Op, Ψ)
    c2 = dot(ϕ, H, Ψ)
    @test abs(c1 - c2) < 1e-12

end


@testset "ScaledOperator mul!" begin

    H₀ = random_matrix(5; hermitian=true)
    H₁ = random_matrix(5; hermitian=true)
    H₂ = random_matrix(5; hermitian=true)

    @test norm(Array(ScaledOperator(0.5, H₀)) - 0.5 * H₀) < 1e-12

    Op = 0.5 * Operator([H₀, H₁, H₂], [2.1, 1.1])
    @test Op isa ScaledOperator
    H = 0.5 * (H₀ + 2.1 * H₁ + 1.1 * H₂)
    @test norm(convert(Matrix{ComplexF64}, Op) - H) < 1e-12

    Op2 = 2 * Op
    @test norm(Array(Op2) - Array(Operator([H₀, H₁, H₂], [2.1, 1.1]))) < 1e-12
    Op2 = Op * 2
    @test norm(Array(Op2) - Array(Operator([H₀, H₁, H₂], [2.1, 1.1]))) < 1e-12

    Ψ = random_state_vector(5)
    ϕ0 = random_state_vector(5)

    ϕ = copy(ϕ0)
    mul!(ϕ, Op, Ψ, true, false)
    ϕ_expected = H * Ψ
    @test norm(ϕ - ϕ_expected) < 1e-12

    ϕ = copy(ϕ0)
    mul!(ϕ, Op, Ψ, true, true)
    ϕ_expected = ϕ0 + H * Ψ
    @test norm(ϕ - ϕ_expected) < 1e-12

    ϕ = copy(ϕ0)
    mul!(ϕ, Op, Ψ, 2.0, true)
    ϕ_expected = ϕ0 + (2 * H) * Ψ
    @test norm(ϕ - ϕ_expected) < 1e-12

    ϕ = copy(ϕ0)
    mul!(ϕ, Op, Ψ, 2.0, 2.0)
    ϕ_expected = 2 * ϕ0 + (2 * H) * Ψ
    @test norm(ϕ - ϕ_expected) < 1e-12

end


@testset "ScaledOperator dot" begin

    H₀ = random_matrix(5; hermitian=true)
    H₁ = random_matrix(5; hermitian=true)
    H₂ = random_matrix(5; hermitian=true)

    Op = Operator([H₀, H₁, H₂], [2.1, 1.1]) * 0.5
    @test Op isa ScaledOperator
    H = 0.5 * (H₀ + 2.1 * H₁ + 1.1 * H₂)

    Ψ = random_state_vector(5)
    ϕ = random_state_vector(5)

    c1 = dot(ϕ, Op, Ψ)
    c2 = dot(ϕ, H, Ψ)
    @test abs(c1 - c2) < 1e-12

end


@testset "ScaledOperator interface" begin
    H₀ = random_matrix(5; hermitian=true)
    H₁ = random_matrix(5; hermitian=true)
    H₂ = random_matrix(5; hermitian=true)
    Op = Operator([H₀, H₁, H₂], [2.1, 1.1]) * 0.5
    @test Op isa ScaledOperator
    Ψ = random_state_vector(5)
    @test check_operator(Op; state=Ψ)
end
