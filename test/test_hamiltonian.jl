using Test
using Logging
using LinearAlgebra

using QuantumPropagators
using QuantumPropagators: Generator, Operator
using QuantumPropagators.Generators: check_generator, check_operator, check_amplitude
using QuantumControlBase.TestUtils:
    random_hermitian_matrix, random_real_matrix, random_state_vector

_OT(::Generator{OT,AT}) where {OT,AT} = OT
_AT(::Generator{OT,AT}) where {OT,AT} = AT
_OT(::Operator{OT,CT}) where {OT,CT}  = OT
_CT(::Operator{OT,CT}) where {OT,CT}  = CT


@testset "Operator interface" begin

    H₀ = random_hermitian_matrix(5, 1.0)
    H₁ = random_hermitian_matrix(5, 1.0)
    H₂ = random_hermitian_matrix(5, 1.0)
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
    QuantumPropagators.Generators.getcontrols(::BrokenAmplitude2) = (ϵ,)
    @test_throws ErrorException("Invalid control amplitude") begin
        with_logger(test_logger) do
            @info "BrokenAmplitude2"
            check_amplitude(BrokenAmplitude2())
        end
    end

    struct BrokenAmplitude3 end
    QuantumPropagators.Generators.getcontrols(::BrokenAmplitude3) = (ϵ,)
    QuantumPropagators.Generators.evalcontrols(::BrokenAmplitude3, _...) = nothing
    QuantumPropagators.Generators.getcontrolderiv(::BrokenAmplitude3, _) = nothing
    @test_throws ErrorException("Invalid control amplitude") begin
        with_logger(test_logger) do
            @info "BrokenAmplitude3"
            check_amplitude(BrokenAmplitude3())
        end
    end

    struct BrokenAmplitude4 end
    QuantumPropagators.Generators.getcontrols(::BrokenAmplitude4) = error("Not implemented")
    @test_throws ErrorException("Not implemented") begin
        with_logger(test_logger) do
            @info "BrokenAmplitude4"
            check_amplitude(BrokenAmplitude4())
        end
    end

    struct BrokenAmplitude5 end
    QuantumPropagators.Generators.getcontrols(::BrokenAmplitude5) = (ϵ,)
    QuantumPropagators.Generators.evalcontrols(::BrokenAmplitude5, _...) =
        error("Not implemented")
    @test_throws ErrorException("Not implemented") begin
        with_logger(test_logger) do
            @info "BrokenAmplitude5"
            check_amplitude(BrokenAmplitude5())
        end
    end

    struct BrokenAmplitude6 end
    QuantumPropagators.Generators.getcontrols(::BrokenAmplitude6) = (ϵ,)
    QuantumPropagators.Generators.evalcontrols(::BrokenAmplitude6, _...) = 1.0
    QuantumPropagators.Generators.getcontrolderiv(::BrokenAmplitude6, _) =
        error("Not implemented")
    @test_throws ErrorException("Invalid control amplitude") begin
        with_logger(test_logger) do
            @info "BrokenAmplitude6"
            check_amplitude(BrokenAmplitude6())
        end
    end

    struct BrokenAmplitude7 end
    QuantumPropagators.Generators.getcontrols(::BrokenAmplitude7) = (ϵ,)
    QuantumPropagators.Generators.evalcontrols(::BrokenAmplitude7, _...) = 1.0
    @test_throws ErrorException("Invalid control amplitude") begin
        with_logger(test_logger) do
            @info "BrokenAmplitude7"
            check_amplitude(BrokenAmplitude7())
        end
    end

end


@testset "Generator interface" begin

    H₀ = random_hermitian_matrix(5, 1.0)
    H₁ = random_hermitian_matrix(5, 1.0)
    ϵ₁ = t -> 1.0
    H₂ = random_hermitian_matrix(5, 1.0)
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
    QuantumPropagators.Generators.getcontrols(::BrokenGenerator) = (ϵ₁, ϵ₂)
    @test_throws ErrorException("Invalid generator") begin
        with_logger(test_logger) do
            @info "BrokenGenerator"
            check_generator(BrokenGenerator(), Ψ)
        end
    end

    struct BrokenGenerator2 end
    QuantumPropagators.Generators.getcontrols(::BrokenGenerator2) = (ϵ₁, ϵ₂)
    QuantumPropagators.Generators.evalcontrols(::BrokenGenerator2, _...) = nothing
    @test_throws ErrorException("Invalid generator") begin
        with_logger(test_logger) do
            @info "BrokenGenerator2"
            check_generator(BrokenGenerator2(), Ψ)
        end
    end

    struct BrokenGenerator3 end
    QuantumPropagators.Generators.getcontrols(::BrokenGenerator3) = (ϵ₁, ϵ₂)
    QuantumPropagators.Generators.evalcontrols(::BrokenGenerator3, _...) = H₀
    QuantumPropagators.Generators.getcontrolderiv(::BrokenGenerator3, _) =
        random_hermitian_matrix(5, 1.0)
    @test_throws ErrorException("Invalid generator") begin
        with_logger(test_logger) do
            @info "BrokenGenerator3"
            check_generator(BrokenGenerator3(), Ψ)
        end
    end

    struct BrokenGenerator4 end
    QuantumPropagators.Generators.getcontrols(::BrokenGenerator4) = (ϵ₁, ϵ₂)
    QuantumPropagators.Generators.evalcontrols(::BrokenGenerator4, _...) = H₀
    QuantumPropagators.Generators.getcontrolderiv(::BrokenGenerator4, _) = nothing
    @test_throws ErrorException("Invalid generator") begin
        with_logger(test_logger) do
            @info "BrokenGenerator4"
            check_generator(BrokenGenerator4(), Ψ)
        end
    end

    struct BrokenGenerator5 end
    QuantumPropagators.Generators.getcontrols(::BrokenGenerator5) = (ϵ₁, ϵ₂)
    QuantumPropagators.Generators.evalcontrols(::BrokenGenerator5, _...) =
        error("Not Implemented")
    @test_throws ErrorException("Not Implemented") begin
        with_logger(test_logger) do
            @info "BrokenGenerator5"
            check_generator(BrokenGenerator5(), Ψ)
        end
    end

    struct BrokenGenerator6 end
    QuantumPropagators.Generators.getcontrols(::BrokenGenerator6) = error("Not Implemented")
    @test_throws ErrorException("Not Implemented") begin
        with_logger(test_logger) do
            @info "BrokenGenerator6"
            check_generator(BrokenGenerator6(), Ψ)
        end
    end

    struct BrokenGenerator7 end
    QuantumPropagators.Generators.getcontrols(::BrokenGenerator7) = (ϵ₁, ϵ₂)
    QuantumPropagators.Generators.evalcontrols(::BrokenGenerator7, _...) = H₀
    QuantumPropagators.Generators.getcontrolderiv(::BrokenGenerator7, control) =
        error("Not Implemented")
    @test_throws ErrorException("Invalid generator") begin
        with_logger(test_logger) do
            @info "BrokenGenerator7"
            check_generator(BrokenGenerator7(), Ψ)
        end
    end

end


@testset "standard Hamiltonians" begin

    H₀ = random_hermitian_matrix(5, 1.0)
    H₁ = random_hermitian_matrix(5, 1.0)
    ϵ₁ = t -> 1.0
    H₂ = random_hermitian_matrix(5, 1.0)
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

    H₀ = random_hermitian_matrix(5, 1.0)

    @test_throws ErrorException("Generator has no terms") begin
        H = hamiltonian()
    end

    @test_throws ErrorException("A Generator requires at least one amplitude") begin
        G = Generator([H₀], [])
    end

end


@testset "vector control Hamiltonians" begin

    H₀ = random_hermitian_matrix(5, 1.0)
    H₁ = random_hermitian_matrix(5, 1.0)
    ϵ₁ = rand(10)
    H₂ = random_hermitian_matrix(5, 1.0)
    ϵ₂ = rand(10)

    H = hamiltonian(H₀, (H₁, ϵ₁), (H₂, ϵ₂))
    @test H isa Generator
    @test _OT(H) == Matrix{ComplexF64}
    @test _AT(H) == Vector{Float64}
    @test repr(H) ==
          "Generator{Matrix{ComplexF64}, Vector{Float64}}(<3 ops>, <2 amplitudes>)"

end


@testset "pathological Hamiltonians" begin

    H₀_r = random_real_matrix(5, 1.0)
    H₀ = random_hermitian_matrix(5, 1.0)
    H₁ = random_hermitian_matrix(5, 1.0)
    ϵ₁ = t -> 1.0
    H₂ = random_hermitian_matrix(5, 1.0)
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

    @test_logs (:warn, "Generator terms may not have been properly expanded") try
        H = hamiltonian((H₀, (H₁, ϵ₁), (H₂, ϵ₂)))
    catch
    end

    @test_logs (:warn, "Generator terms may not have been properly expanded") try
        H = hamiltonian([H₀, (H₁, ϵ₁), (H₂, ϵ₂)])
    catch
    end

    @test_warn "Generator terms may not have been properly expanded" begin
        H = hamiltonian((H₀, (H₁, ϵ₁)))
    end

    @test_warn "Generator terms may not have been properly expanded" begin
        H = hamiltonian([H₀, (H₁, ϵ₁)])
    end

    @test_throws ErrorException("time-dependent term must be 2-tuple") begin
        H = hamiltonian(H₀_r, (H₁, ϵ₁, ϵ₂))
    end

    @test_logs (:error, r"Collected drift operators are of a disparate type.*") try
        H = hamiltonian(H₀_r, ϵ₁, (H₂, ϵ₂))
    catch
    end

    @test_warn "Collected operators are not of a concrete type: Function" begin
        H = hamiltonian((ϵ₁, H₁), (ϵ₂, H₂))
    end

    #! format: off
    _w1 = "It looks like (op, ampl) in term are reversed"
    _w2 = "It looks like (op, ampl) in term are reversed"
    _w3 = "Collected operators are not of a concrete type: Any"
    _e1 = "evalcontrols(::Matrix{ComplexF64}, …) for amplitude does not return a number"
    _e2 = "getcontrolderiv(::Matrix{ComplexF64}, …) for amplitude does not return `0.0` for a control that the amplitude does not depend on"
    _w4 = "Collected amplitude #1 is invalid"
    _e3 = "evalcontrols(::Matrix{ComplexF64}, …) for amplitude does not return a number"
    _e4 = "getcontrolderiv(::Matrix{ComplexF64}, …) for amplitude does not return `0.0` for a control that the amplitude does not depend on"
    _w5 = "Collected amplitude #2 is invalid"
    @test_logs (:warn, _w1) (:warn, _w2) (:warn, _w3) (:error, _e1) (:error, _e2) (:warn, _w4) (:error, _e3) (:error, _e4) (:warn, _w5) begin
        H = hamiltonian(H₀, (ϵ₁, H₁), (ϵ₂, H₂))
    end
    @test repr(H) == "Generator{Any, Matrix{ComplexF64}}(<3 ops>, <2 amplitudes>)"
    #! format: on

end
