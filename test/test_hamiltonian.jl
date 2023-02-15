using Test
using Logging
using LinearAlgebra

using QuantumPropagators
using QuantumPropagators: Generator, Operator
using QuantumControlBase
using QuantumControlTestUtils.RandomObjects: random_matrix, random_state_vector
using QuantumControlTestUtils: QuantumTestLogger

_OT(::Generator{OT,AT}) where {OT,AT} = OT
_AT(::Generator{OT,AT}) where {OT,AT} = AT
_OT(::Operator{OT,CT}) where {OT,CT}  = OT
_CT(::Operator{OT,CT}) where {OT,CT}  = CT

#=

@testset "Operator interface" begin

    H₀ = random_matrix(5; hermitian=true)
    H₁ = random_matrix(5; hermitian=true)
    H₂ = random_matrix(5; hermitian=true)
    Ψ = random_state_vector(5)

    H = Operator([H₀, H₁, H₂], [1.1, 2.1])

    test_logger = NullLogger()  # hide errors
    #test_logger = global_logger()  # show errors

    @test check_operator(H₀, Ψ)
    @test check_operator(H, Ψ)

end


@testset "Amplitude interface" begin

    ϵ(t) = 1.0

    test_logger = NullLogger()  # hide errors
    #test_logger = global_logger()  # show errors

    @test check_amplitude(ϵ)

    struct BrokenAmplitude end
    @test_throws ErrorException("Invalid control amplitude") begin
        with_logger(test_logger) do
            @info "BrokenAmplitude"
            check_amplitude(BrokenAmplitude())
        end
    end

    struct BrokenAmplitude2 end
    QuantumPropagators.Controls.get_controls(::BrokenAmplitude2) = (ϵ,)
    @test_throws ErrorException("Invalid control amplitude") begin
        with_logger(test_logger) do
            @info "BrokenAmplitude2"
            check_amplitude(BrokenAmplitude2())
        end
    end

    struct BrokenAmplitude3 end
    QuantumPropagators.Controls.get_controls(::BrokenAmplitude3) = (ϵ,)
    QuantumPropagators.Controls.evalcontrols(::BrokenAmplitude3, args...; kwargs...) = nothing
    QuantumControlBase.get_control_deriv(::BrokenAmplitude3, _) = nothing
    @test_throws ErrorException("Invalid control amplitude") begin
        with_logger(test_logger) do
            @info "BrokenAmplitude3"
            check_amplitude(BrokenAmplitude3())
        end
    end

    struct BrokenAmplitude4 end
    QuantumPropagators.Controls.get_controls(::BrokenAmplitude4) = error("Not implemented")
    @test_throws ErrorException("Not implemented") begin
        with_logger(test_logger) do
            @info "BrokenAmplitude4"
            check_amplitude(BrokenAmplitude4())
        end
    end

    struct BrokenAmplitude5 end
    QuantumPropagators.Controls.get_controls(::BrokenAmplitude5) = (ϵ,)
    QuantumPropagators.Controls.evalcontrols(::BrokenAmplitude5, args...; kwargs..) =
        error("Not implemented")
    @test_throws ErrorException("Not implemented") begin
        with_logger(test_logger) do
            @info "BrokenAmplitude5"
            check_amplitude(BrokenAmplitude5())
        end
    end

    struct BrokenAmplitude6 end
    QuantumPropagators.Controls.get_controls(::BrokenAmplitude6) = (ϵ,)
    QuantumPropagators.Controls.evalcontrols(::BrokenAmplitude6, args...) = 1.0
    QuantumControlBase.get_control_deriv(::BrokenAmplitude6, _) = error("Not implemented")
    @test_throws ErrorException("Invalid control amplitude") begin
        with_logger(test_logger) do
            @info "BrokenAmplitude6"
            check_amplitude(BrokenAmplitude6())
        end
    end

    struct BrokenAmplitude7 end
    QuantumPropagators.Controls.get_controls(::BrokenAmplitude7) = (ϵ,)
    QuantumPropagators.Controls.evalcontrols(::BrokenAmplitude7, args...) = 1.0
    @test_throws ErrorException("Invalid control amplitude") begin
        with_logger(test_logger) do
            @info "BrokenAmplitude7"
            check_amplitude(BrokenAmplitude7())
        end
    end

end


@testset "Generator interface" begin

    H₀ = random_matrix(5; hermitian=true)
    H₁ = random_matrix(5; hermitian=true)
    ϵ₁ = t -> 1.0
    H₂ = random_matrix(5; hermitian=true)
    ϵ₂ = t -> 1.0

    H = hamiltonian(H₀)
    Ψ = random_state_vector(5)

    test_logger = NullLogger()  # hide errors
    #test_logger = global_logger()  # show errors

    @test check_generator(H₀, Ψ)

    @test check_generator(H, Ψ)

    @test_throws ErrorException("Invalid generator") begin
        with_logger(test_logger) do
            @info "nothing-Generator"
            check_generator(nothing, Ψ)
        end
    end

    struct BrokenGenerator end
    QuantumPropagators.Controls.get_controls(::BrokenGenerator) = (ϵ₁, ϵ₂)
    @test_throws ErrorException("Invalid generator") begin
        with_logger(test_logger) do
            @info "BrokenGenerator"
            check_generator(BrokenGenerator(), Ψ)
        end
    end

    struct BrokenGenerator2 end
    QuantumPropagators.Controls.get_controls(::BrokenGenerator2) = (ϵ₁, ϵ₂)
    QuantumPropagators.Controls.evalcontrols(::BrokenGenerator2, _...) = nothing
    @test_throws ErrorException("Invalid generator") begin
        with_logger(test_logger) do
            @info "BrokenGenerator2"
            check_generator(BrokenGenerator2(), Ψ)
        end
    end

    struct BrokenGenerator3 end
    QuantumPropagators.Controls.get_controls(::BrokenGenerator3) = (ϵ₁, ϵ₂)
    QuantumPropagators.Controls.evalcontrols(::BrokenGenerator3, _...) = H₀
    QuantumControlBase.get_control_deriv(::BrokenGenerator3, _) =
        random_matrix(5; hermitian=true)
    @test_throws ErrorException("Invalid generator") begin
        with_logger(test_logger) do
            @info "BrokenGenerator3"
            check_generator(BrokenGenerator3(), Ψ)
        end
    end

    struct BrokenGenerator4 end
    QuantumPropagators.Controls.get_controls(::BrokenGenerator4) = (ϵ₁, ϵ₂)
    QuantumPropagators.Controls.evalcontrols(::BrokenGenerator4, _...) = H₀
    QuantumControlBase.get_control_deriv(::BrokenGenerator4, _) = nothing
    @test_throws ErrorException("Invalid generator") begin
        with_logger(test_logger) do
            @info "BrokenGenerator4"
            check_generator(BrokenGenerator4(), Ψ)
        end
    end

    struct BrokenGenerator5 end
    QuantumPropagators.Controls.get_controls(::BrokenGenerator5) = (ϵ₁, ϵ₂)
    QuantumPropagators.Controls.evalcontrols(::BrokenGenerator5, _...) =
        error("Not Implemented")
    @test_throws ErrorException("Not Implemented") begin
        with_logger(test_logger) do
            @info "BrokenGenerator5"
            check_generator(BrokenGenerator5(), Ψ)
        end
    end

    struct BrokenGenerator6 end
    QuantumPropagators.Controls.get_controls(::BrokenGenerator6) = error("Not Implemented")
    @test_throws ErrorException("Not Implemented") begin
        with_logger(test_logger) do
            @info "BrokenGenerator6"
            check_generator(BrokenGenerator6(), Ψ)
        end
    end

    struct BrokenGenerator7 end
    QuantumPropagators.Controls.get_controls(::BrokenGenerator7) = (ϵ₁, ϵ₂)
    QuantumPropagators.Controls.evalcontrols(::BrokenGenerator7, _...) = H₀
    QuantumControlBase.get_control_deriv(::BrokenGenerator7, control) =
        error("Not Implemented")
    @test_throws ErrorException("Invalid generator") begin
        with_logger(test_logger) do
            @info "BrokenGenerator7"
            check_generator(BrokenGenerator7(), Ψ)
        end
    end

end

=#


@testset "standard Hamiltonians" begin

    H₀ = random_matrix(5; hermitian=true)
    H₁ = random_matrix(5; hermitian=true)
    ϵ₁ = t -> 1.0
    H₂ = random_matrix(5; hermitian=true)
    ϵ₂ = t -> 1.0

    H = hamiltonian(H₀)
    @test H ≡ H₀

    H = hamiltonian(H₀, H₁)
    @test H isa Matrix
    @test norm(H - H₀ - H₁) < 1e-14

    H = hamiltonian(H₀, (H₁, ϵ₁))
    @test H isa Generator
    @test _OT(H) == Matrix{ComplexF64}
    @test _AT(H) <: Function
    @test _AT(H) == typeof(ϵ₁)
    @test length(H.ops) == 2
    @test length(H.amplitudes) == 1
    @test (length(H.ops) - length(H.amplitudes)) == 1
    @test H.ops[1] ≡ H₀
    @test H.ops[2] ≡ H₁
    @test H.amplitudes[1] ≡ ϵ₁

    H = hamiltonian(H₀, (H₁, ϵ₁), (H₂, ϵ₂))
    @test H isa Generator
    @test _OT(H) == Matrix{ComplexF64}
    @test _AT(H) == Function
    @test length(H.ops) == 3
    @test length(H.amplitudes) == 2
    @test (length(H.ops) - length(H.amplitudes)) == 1
    @test H.ops[1] ≡ H₀
    @test H.ops[2] ≡ H₁
    @test H.ops[3] ≡ H₂
    @test H.amplitudes[1] ≡ ϵ₁
    @test H.amplitudes[2] ≡ ϵ₂
    @test repr(H) == "Generator{Matrix{ComplexF64}, Function}(<3 ops>, <2 amplitudes>)"

    H = hamiltonian(H₀, (H₁, ϵ₁), (H₂, ϵ₁))
    @test H isa Generator
    @test length(H.ops) == 2
    @test length(H.amplitudes) == 1
    @test (length(H.ops) - length(H.amplitudes)) == 1
    @test H.ops[1] ≡ H₀
    @test norm(H.ops[2] - H₁ - H₂) < 1e-14
    @test H.amplitudes[1] ≡ ϵ₁

    H = hamiltonian(H₀, (H₁, ϵ₁), (H₁, ϵ₂))
    @test H isa Generator
    @test length(H.ops) == 3
    @test length(H.amplitudes) == 2
    @test (length(H.ops) - length(H.amplitudes)) == 1
    @test H.ops[1] ≡ H₀
    @test H.ops[2] ≡ H₁
    @test H.ops[3] ≡ H₁
    @test H.amplitudes[1] ≡ ϵ₁
    @test H.amplitudes[2] ≡ ϵ₂

    H = hamiltonian(H₀, (H₁, ϵ₁), H₂)
    @test H isa Generator
    @test length(H.ops) == 2
    @test length(H.amplitudes) == 1
    @test (length(H.ops) - length(H.amplitudes)) == 1
    @test norm(H.ops[1] - H₀ - H₂) < 1e-14
    @test H.amplitudes[1] ≡ ϵ₁

end

@testset "pathologial Hamiltonians" begin

    H₀ = random_matrix(5; hermitian=true)

    @test_throws ErrorException("Generator has no terms") begin
        H = hamiltonian()
    end

    @test_throws ErrorException("A Generator requires at least one amplitude") begin
        G = Generator([H₀], [])
    end

end


@testset "vector control Hamiltonians" begin

    H₀ = random_matrix(5; hermitian=true)
    H₁ = random_matrix(5; hermitian=true)
    ϵ₁ = rand(10)
    H₂ = random_matrix(5; hermitian=true)
    ϵ₂ = rand(10)

    H = hamiltonian(H₀, (H₁, ϵ₁), (H₂, ϵ₂))
    @test H isa Generator
    @test _OT(H) == Matrix{ComplexF64}
    @test _AT(H) == Vector{Float64}
    @test repr(H) ==
          "Generator{Matrix{ComplexF64}, Vector{Float64}}(<3 ops>, <2 amplitudes>)"

end


@testset "pathological Hamiltonians" begin

    H₀_r = random_matrix(5; complex=false)
    H₀ = random_matrix(5; hermitian=true)
    H₁ = random_matrix(5; hermitian=true)
    ϵ₁ = t -> 1.0
    H₂ = random_matrix(5; hermitian=true)
    ϵ₂ = t -> 1.0
    ϵ₂_v = rand(10)

    @test_logs (:warn, "Collected amplitudes are of disparate types") begin
        H = hamiltonian(H₀, (H₁, ϵ₁), (H₂, ϵ₂_v))
    end
    @test H isa Generator
    @test _OT(H) == Matrix{ComplexF64}
    @test _AT(H) == Any
    @test repr(H) == "Generator{Matrix{ComplexF64}, Any}(<3 ops>, <2 amplitudes>)"

    H = hamiltonian(nothing, (nothing, ϵ₁), (nothing, ϵ₂))
    @test repr(H) == "Generator{Nothing, Function}(<3 ops>, <2 amplitudes>)"

    logger = QuantumTestLogger()
    with_logger(logger) do
        try
            H = hamiltonian((H₀, (H₁, ϵ₁), (H₂, ϵ₂)))
        catch
        end
    end
    @test "Warn: Generator terms may not have been properly expanded" ∈ logger

    logger = QuantumTestLogger()
    with_logger(logger) do
        try
            H = hamiltonian([H₀, (H₁, ϵ₁), (H₂, ϵ₂)])
        catch
        end
    end
    @test "Warn: Generator terms may not have been properly expanded" ∈ logger

    logger = QuantumTestLogger()
    with_logger(logger) do
        H = hamiltonian((H₀, (H₁, ϵ₁)))
    end
    @test "Warn: Generator terms may not have been properly expanded" ∈ logger

    logger = QuantumTestLogger()
    with_logger(logger) do
        H = hamiltonian([H₀, (H₁, ϵ₁)])
    end
    @test "Warn: Generator terms may not have been properly expanded" ∈ logger

    logger = QuantumTestLogger()
    with_logger(logger) do
        try
            H = hamiltonian(H₀_r, (H₁, ϵ₁, ϵ₂))
            @test "raised Exception" == ""
        catch exc
            @test exc == ErrorException("time-dependent term must be 2-tuple")
        end
    end

    logger = QuantumTestLogger()
    with_logger(logger) do
        try
            H = hamiltonian(H₀_r, ϵ₁, (H₂, ϵ₂))
        catch exc
        end
    end
    @test "Error: Collected drift operators are of a disparate type" ∈ logger

    logger = QuantumTestLogger()
    with_logger(logger) do
        H = hamiltonian((ϵ₁, H₁), (ϵ₂, H₂))
    end
    @test "Warn: It looks like (op, ampl) in term are reversed" ∈ logger
    @test "Warn: Collected operators are not of a concrete type" ∈ logger

    logger = QuantumTestLogger()
    with_logger(logger) do
        H = hamiltonian(H₀, (ϵ₁, H₁), (ϵ₂, H₂))
    end
    @test "Warn: It looks like (op, ampl) in term are reversed" ∈ logger
    @test "Warn: Collected operators are not of a concrete type: Any" ∈ logger
    #@test r"Error: evaluate.* for amplitude does not return a number" ∈ logger
    #@test "Warn: Collected amplitude #1 is invalid" ∈ logger
    #@test "Warn: Collected amplitude #2 is invalid" ∈ logger
    @test repr(H) == "Generator{Any, Matrix{ComplexF64}}(<3 ops>, <2 amplitudes>)"

end
