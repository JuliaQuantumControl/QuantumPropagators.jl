using Test
using LinearAlgebra
using QuantumControlTestUtils.RandomObjects: random_matrix, random_state_vector
using QuantumPropagators.Interfaces:
    check_operator, check_generator, supports_matrix_interface
import QuantumPropagators.Interfaces: supports_inplace
import QuantumPropagators.Controls: get_controls, evaluate
import ArrayInterface
using StaticArrays: SMatrix, SVector

using QuantumPropagators: Generator, Operator, ScaledOperator

# A minimal operator type that does not subtype AbstractMatrix, to test that
# supports_matrix_interface propagates correctly through Operator/ScaledOperator.
struct MatrixFreeOp
    mat::Matrix{ComplexF64}
end

Base.size(op::MatrixFreeOp) = size(op.mat)
Base.size(op::MatrixFreeOp, dim::Integer) = size(op.mat, dim)
Base.:*(α::Number, op::MatrixFreeOp) = MatrixFreeOp(α * op.mat)
Base.:*(op::MatrixFreeOp, Ψ::AbstractVector) = op.mat * Ψ
LinearAlgebra.mul!(ϕ, op::MatrixFreeOp, Ψ, α::Number, β::Number) = mul!(ϕ, op.mat, Ψ, α, β)
LinearAlgebra.dot(ϕ, op::MatrixFreeOp, Ψ) = dot(ϕ, op.mat, Ψ)
supports_inplace(::Type{MatrixFreeOp}) = true
get_controls(::MatrixFreeOp) = ()
evaluate(op::MatrixFreeOp, args...; kwargs...) = op


@testset "Operator mul!" begin

    H₀ = random_matrix(5; hermitian = true)
    H₁ = random_matrix(5; hermitian = true)
    H₂ = random_matrix(5; hermitian = true)

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


@testset "Operator * Ψ" begin

    N = 5
    H₀ = SMatrix{N,N}(random_matrix(N; hermitian = true))
    H₁ = SMatrix{N,N}(random_matrix(N; hermitian = true))
    H₂ = SMatrix{N,N}(random_matrix(N; hermitian = true))

    Op = Operator([H₀, H₁, H₂], [2.1, 1.1])
    H = H₀ + 2.1 * H₁ + 1.1 * H₂

    @test norm(convert(Matrix{ComplexF64}, Op) - H) < 1e-12

    Ψ = SVector{N}(random_state_vector(N))

    ϕ = Op * Ψ
    ϕ_expected = H * Ψ
    @test norm(ϕ - ϕ_expected) < 1e-12

end


@testset "Operator interface" begin
    H₀ = random_matrix(5; hermitian = true)
    H₁ = random_matrix(5; hermitian = true)
    H₂ = random_matrix(5; hermitian = true)
    Op = Operator([H₀, H₁, H₂], [2.1, 1.1])
    Ψ = random_state_vector(5)
    @test supports_matrix_interface(typeof(Op))
    @test check_operator(Op; state = Ψ)
end


@testset "Operator copy" begin

    H₀ = random_matrix(5; hermitian = true)
    H₁ = random_matrix(5; hermitian = true)
    H₂ = random_matrix(5; hermitian = true)

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

    H₀ = random_matrix(5; hermitian = true)
    H₁ = random_matrix(5; hermitian = true)
    H₂ = random_matrix(5; hermitian = true)

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
    @test startswith(repr(Op), "ScaledOperator(2, Operator(Matrix{ComplexF64}")
    @test summary(Op) == "ScaledOperator with coeff=2"
    repl_repr = repr("text/plain", Op; context = (:limit => true))
    @test contains(repl_repr, "operator.ops::Vector{Matrix{ComplexF64}}")
    @test contains(repl_repr, "operator.coeffs: [2.1, 1.1]")
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

    H₀ = random_matrix(5; hermitian = true)
    H₁ = random_matrix(5; hermitian = true)
    H₂ = random_matrix(5; hermitian = true)

    Op = Operator([H₀, H₁, H₂], [2.1, 1.1])
    H = H₀ + 2.1 * H₁ + 1.1 * H₂

    Ψ = random_state_vector(5)
    ϕ = random_state_vector(5)

    c1 = dot(ϕ, Op, Ψ)
    c2 = dot(ϕ, H, Ψ)
    @test abs(c1 - c2) < 1e-12

end


@testset "ScaledOperator mul!" begin

    H₀ = random_matrix(5; hermitian = true)
    H₁ = random_matrix(5; hermitian = true)
    H₂ = random_matrix(5; hermitian = true)

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


@testset "ScaledOperator * Ψ" begin

    N = 5
    H₀ = SMatrix{N,N}(random_matrix(N; hermitian = true))
    H₁ = SMatrix{N,N}(random_matrix(N; hermitian = true))
    H₂ = SMatrix{N,N}(random_matrix(N; hermitian = true))

    Op = 0.5 * Operator([H₀, H₁, H₂], [2.1, 1.1])
    @test Op isa ScaledOperator

    H = 0.5 * (H₀ + 2.1 * H₁ + 1.1 * H₂)

    @test norm(convert(Matrix{ComplexF64}, Op) - H) < 1e-12

    Ψ = SVector{N}(random_state_vector(N))

    ϕ = Op * Ψ
    ϕ_expected = H * Ψ
    @test norm(ϕ - ϕ_expected) < 1e-12

end


@testset "ScaledOperator dot" begin

    H₀ = random_matrix(5; hermitian = true)
    H₁ = random_matrix(5; hermitian = true)
    H₂ = random_matrix(5; hermitian = true)

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
    H₀ = random_matrix(5; hermitian = true)
    H₁ = random_matrix(5; hermitian = true)
    H₂ = random_matrix(5; hermitian = true)
    Op = Operator([H₀, H₁, H₂], [2.1, 1.1]) * 0.5
    @test Op isa ScaledOperator
    Ψ = random_state_vector(5)
    @test supports_matrix_interface(typeof(Op))
    @test check_operator(Op; state = Ψ)
end


@testset "supports_matrix_interface" begin

    H₀ = random_matrix(5; hermitian = true)
    H₁ = random_matrix(5; hermitian = true)
    H₂ = random_matrix(5; hermitian = true)
    Ψ = random_state_vector(5)

    # Operator wrapping standard Matrix → true
    Op = Operator([H₀, H₁, H₂], [2.1, 1.1])
    @test supports_matrix_interface(typeof(Op))

    # ScaledOperator wrapping Operator{Matrix} → true
    SOp = Op * 0.5
    @test SOp isa ScaledOperator
    @test supports_matrix_interface(typeof(SOp))

    # Operator wrapping MatrixFreeOp (not AbstractMatrix) → false, but valid operator
    h₀ = MatrixFreeOp(H₀)
    h₁ = MatrixFreeOp(H₁)
    h₂ = MatrixFreeOp(H₂)
    FreeOp = Operator([h₀, h₁, h₂], [2.1, 1.1])
    @test !supports_matrix_interface(typeof(FreeOp))
    @test check_operator(FreeOp; state = Ψ)

    # ScaledOperator wrapping Operator{MatrixFreeOp} → false, but valid operator
    SFreeOp = FreeOp * 0.5
    @test SFreeOp isa ScaledOperator
    @test !supports_matrix_interface(typeof(SFreeOp))
    @test check_operator(SFreeOp; state = Ψ)

end


@testset "Hermitian matrix supports in-place operations" begin

    # Test the resolution of
    # https://github.com/JuliaQuantumControl/QuantumPropagators.jl/issues/102

    Ψ0 = ComplexF64[1, 0]
    Ĥ = ComplexF64[
         0   0.5
        0.5   0
    ]
    tlist = collect(range(0.0, 1.0, length = 101))
    generator = (Hermitian(Ĥ),)
    op = evaluate(generator, tlist, 1)
    T = Hermitian{ComplexF64,Matrix{ComplexF64}}
    @test op isa T
    @test supports_inplace(T)
    op2 = similar(op)
    @test op2 isa T
    @test ArrayInterface.ismutable(T)
    @test check_generator(generator; state = Ψ0, tlist)

end
