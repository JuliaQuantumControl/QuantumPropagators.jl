# SPDX-FileCopyrightText: © 2021 Michael Goerz <mail@michaelgoerz.net>
#
# SPDX-License-Identifier: MIT

using Test
using Logging: with_logger, NullLogger
using LinearAlgebra: LinearAlgebra, I, norm, dot, lmul!
using StableRNGs: StableRNG
using IOCapture: IOCapture
using StaticArrays: SVector, SMatrix
using QuantumControlTestUtils.RandomObjects: random_state_vector, random_matrix

using QuantumPropagators: QuantumPropagators, Cheby, init_prop, prop_step!, propagate
using QuantumPropagators.Storage:
    init_storage, write_to_storage!, get_from_storage, get_from_storage!, map_observables
using QuantumPropagators.Interfaces: check_storage
import QuantumPropagators.Interfaces: supports_inplace
import QuantumPropagators.Storage


# A state that wraps a `Matrix`, to exercise the generic `Vector`-of-states
# storage path with a type that is not an `AbstractArray`.
struct GateState
    U::Matrix{ComplexF64}
end

supports_inplace(::Type{GateState}) = true
Base.copy(Ψ::GateState) = GateState(copy(Ψ.U))
Base.similar(Ψ::GateState) = GateState(similar(Ψ.U))
Base.copyto!(Ψ::GateState, ϕ::GateState) = (copyto!(Ψ.U, ϕ.U); Ψ)
Base.:*(c::Number, Ψ::GateState) = GateState(c * Ψ.U)
Base.:-(Ψ::GateState, ϕ::GateState) = GateState(Ψ.U - ϕ.U)
LinearAlgebra.norm(Ψ::GateState) = norm(Ψ.U)


# An immutable state that is not a `StaticArray`: `supports_inplace` is
# `false`, but the object is not `isbits` either, so the storage must still
# take a copy.
struct FrozenState
    v::Vector{ComplexF64}
end

supports_inplace(::Type{FrozenState}) = false
Base.copy(Ψ::FrozenState) = FrozenState(copy(Ψ.v))
Base.:*(c::Number, Ψ::FrozenState) = FrozenState(c * Ψ.v)
Base.:-(Ψ::FrozenState, ϕ::FrozenState) = FrozenState(Ψ.v - ϕ.v)
LinearAlgebra.norm(Ψ::FrozenState) = norm(Ψ.v)


# A deliberately broken storage: it keeps a *reference* to the written data.
# This is the bug that `check_storage` must detect.
struct AliasingStorage{T}
    data::Vector{T}
end

Storage.init_storage(Ψ::GateState, nt::Integer) =
    AliasingStorage(Vector{GateState}(undef, nt))
Storage.write_to_storage!(storage::AliasingStorage, i::Integer, data) =
    (storage.data[i] = data; storage)  # BUG: no copy
Storage.get_from_storage(storage::AliasingStorage, i) = storage.data[i]
Storage.get_from_storage!(data, storage::AliasingStorage, i) =
    copyto!(data, storage.data[i])


quiet_check(args...; kwargs...) =
    with_logger(NullLogger()) do
        check_storage(args...; kwargs..., quiet = true)
    end


@testset "Ownership of stored data" begin

    # The contract violation from #119, in isolation
    Û = ComplexF64[1 0; 0 1]
    storage = init_storage(Û, 3)
    write_to_storage!(storage, 1, Û)
    lmul!(2.0, Û)
    @test get_from_storage(storage, 1) == ComplexF64[1 0; 0 1]
    @test get_from_storage(storage, 1) ≢ Û

    # A `Vector` state goes into a `Matrix` storage, which always copies
    Ψ = ComplexF64[1, 0]
    storage = init_storage(Ψ, 3)
    write_to_storage!(storage, 1, Ψ)
    lmul!(2.0, Ψ)
    @test get_from_storage(storage, 1) == ComplexF64[1, 0]

    # Immutable `isbits` data is stored by reference (no copy)
    Ψ = SVector{2}(ComplexF64[1, 0])
    storage = init_storage(Ψ, 3)
    write_to_storage!(storage, 1, Ψ)
    @test get_from_storage(storage, 1) === Ψ

    # A non-`isbits` immutable wrapper is *not* provably safe, so it is copied
    Ψ = FrozenState(ComplexF64[1, 0])
    storage = init_storage(Ψ, 3)
    write_to_storage!(storage, 1, Ψ)
    @test get_from_storage(storage, 1) ≢ Ψ
    @test get_from_storage(storage, 1).v == Ψ.v

    # Writing to one slot must not disturb another
    Û = ComplexF64[1 0; 0 1]
    storage = init_storage(Û, 3)
    write_to_storage!(storage, 1, Û)
    lmul!(2.0, Û)
    write_to_storage!(storage, 2, Û)
    @test get_from_storage(storage, 1) == ComplexF64[1 0; 0 1]
    @test get_from_storage(storage, 2) == ComplexF64[2 0; 0 2]

    # `write_to_storage!` returns the storage
    @test write_to_storage!(storage, 3, Û) ≡ storage
    storage = init_storage(ComplexF64[1, 0], 3)
    @test write_to_storage!(storage, 1, ComplexF64[1, 0]) ≡ storage

end


@testset "check_storage for states" begin

    tlist = collect(range(0, 10; length = 10))

    Ψ = random_state_vector(4; rng = StableRNG(2821563937))
    @test check_storage(Ψ, tlist; rng = StableRNG(610341865))

    Û = random_matrix(4; rng = StableRNG(1671185571))
    @test check_storage(Û, tlist; rng = StableRNG(1321808964))

    @test check_storage(SVector{4}(Ψ), tlist; rng = StableRNG(3722552871))
    @test check_storage(SMatrix{4,4}(Û), tlist; rng = StableRNG(2428674146))

    @test check_storage(FrozenState(Ψ), tlist; rng = StableRNG(3139269722))

    # `nt` must not matter: with a long time grid, no slot may drift below
    # `atol` and thus become effectively unchecked
    @test check_storage(
        Û,
        collect(range(0, 10; length = 1000));
        rng = StableRNG(3688239531)
    )

end


@testset "check_storage detects a reference-storing implementation" begin

    tlist = collect(range(0, 10; length = 10))
    Ψ = GateState(random_matrix(4; rng = StableRNG(2821563937)))

    # `AliasingStorage` stores a reference, which the ownership loop detects
    @test !quiet_check(Ψ, tlist; rng = StableRNG(610341865))

    # The failure is reported
    c = IOCapture.capture(rethrow = Union{}) do
        check_storage(Ψ, tlist; rng = StableRNG(610341865))
    end
    @test c.value ≡ false
    @test contains(c.output, "the storage must own its data")

end


@testset "check_storage for observables" begin

    tlist = collect(range(0, 10; length = 10))
    Ψ = random_state_vector(4; rng = StableRNG(2821563937))
    Ô₁ = random_matrix(4; hermitian = true, rng = StableRNG(610341865))
    Ô₂ = random_matrix(4; hermitian = true, rng = StableRNG(1671185571))

    # a number
    @test check_storage(Ψ, tlist, (Ψ -> dot(Ψ, Ô₁, Ψ),); rng = StableRNG(1321808964))
    # a matrix observable (expectation value via three-argument `dot`)
    @test check_storage(Ψ, tlist, (Ô₁,); rng = StableRNG(3722552871))
    # a vector (e.g., populations)
    @test check_storage(Ψ, tlist, (Ψ -> abs.(Ψ) .^ 2,); rng = StableRNG(2428674146))
    # uniform multiple observables (stored as the columns of a matrix)
    @test check_storage(Ψ, tlist, (Ô₁, Ô₂); rng = StableRNG(3139269722))
    # non-uniform multiple observables (stored as a vector of tuples)
    observables = (Ψ -> real(dot(Ψ, Ô₁, Ψ)), Ψ -> count(abs.(Ψ) .^ 2 .> 0.1))
    @test check_storage(Ψ, tlist, observables; rng = StableRNG(3688239531))
    # a tuple containing mutable data
    observables = (Ψ -> abs.(Ψ) .^ 2, Ψ -> count(abs.(Ψ) .^ 2 .> 0.1))
    @test check_storage(Ψ, tlist, observables; rng = StableRNG(2821563937))

    # the state itself: the `_StoreState` path used by `propagate`
    @test check_storage(
        Ψ,
        tlist,
        QuantumPropagators._StoreState();
        rng = StableRNG(610341865)
    )
    Û = random_matrix(4; rng = StableRNG(1671185571))
    @test check_storage(
        Û,
        tlist,
        QuantumPropagators._StoreState();
        rng = StableRNG(1321808964)
    )

end


@testset "check_storage for plain data" begin

    @test check_storage(1.0, 10; transform = (x -> 0.9x))
    @test check_storage(1.0 + 0.5im, 10; transform = (x -> 0.9x))
    @test check_storage(ComplexF64[1, 2, 3], 10; transform = (v -> 0.9v))
    @test check_storage((1.0, 2), 10; transform = (t -> (0.9t[1], t[2] + 1)))

    # An in-place `transform` on mutable non-`Vector` data exercises ownership
    # in the generic `Vector`-of-data storage path
    M = ComplexF64[1 0; 0 1]
    @test check_storage(M, 10; transform = (M -> (lmul!(0.9, M); M)))

end


@testset "Storing the states of an in-place propagation (#119)" begin

    # The MWE from https://github.com/JuliaQuantumControl/QuantumPropagators.jl/issues/119
    GHz = 1.0
    MHz = 0.001GHz
    ns = 1.0

    Û₀ = Matrix{ComplexF64}(I(4)) / 2

    Ĥ = let δ₁ = 100MHz, δ₂ = -100MHz, J = 3MHz, Ω₀ = 35MHz
        Matrix{ComplexF64}([
                0 Ω₀/2 Ω₀/2 0
                Ω₀/2 δ₂ J Ω₀/2
                Ω₀/2 J δ₁ Ω₀/2
                0 Ω₀/2 Ω₀/2 δ₁+δ₂
            ])
    end

    T = 400ns
    tlist = collect(range(0, T, step = 0.1ns))
    nt = length(tlist)

    propagator = init_prop(Û₀, Ĥ, tlist; method = Cheby)
    storage = init_storage(Û₀, tlist)
    expected = [copy(Û₀)]
    write_to_storage!(storage, 1, Û₀)
    for n = 1:(nt-1)
        state = prop_step!(propagator)
        write_to_storage!(storage, n + 1, state)
        push!(expected, copy(state))
    end

    # Every slot must hold the state as of the time it was written
    @test all(get_from_storage(storage, n) ≈ expected[n] for n = 1:nt)
    # Consecutive stored states must be distinct objects with different values
    @test get_from_storage(storage, nt) ≢ get_from_storage(storage, nt - 1)
    @test norm(get_from_storage(storage, nt) - get_from_storage(storage, nt - 1)) > 1e-3
    @test norm(get_from_storage(storage, nt) - get_from_storage(storage, 1)) > 0.1
    @test get_from_storage(storage, 1) ≈ Û₀

    # The same via `propagate`
    storage = propagate(Û₀, Ĥ, tlist; method = Cheby, storage = true)
    @test size(storage) == (nt,)
    @test all(storage[n] ≈ expected[n] for n = 1:nt)
    @test storage[nt] ≢ storage[nt-1]
    @test norm(storage[nt] - storage[nt-1]) > 1e-3

end
