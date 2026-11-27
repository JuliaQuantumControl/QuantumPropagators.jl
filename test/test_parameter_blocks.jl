# SPDX-FileCopyrightText: © 2026 Michael Goerz <mail@michaelgoerz.net>
#
# SPDX-License-Identifier: MIT

using Test
using ComponentArrays: ComponentVector, Axis
using RecursiveArrayTools: ArrayPartition
using OffsetArrays: OffsetVector
using LinearAlgebra: norm
using StableRNGs: StableRNG
using IOCapture: IOCapture
using QuantumControlTestUtils.RandomObjects: random_matrix, random_state_vector
using QuantumPropagators: QuantumPropagators, hamiltonian, propagate, ExpProp
using QuantumPropagators.Controls:
    ParameterizedFunction,
    StaticParameter,
    get_controls,
    get_parameters,
    combine_parameters,
    evaluate
using QuantumPropagators.Interfaces:
    check_parameterized_function, check_parameterized, check_generator


# A linear control `a + b * t`, with `a` and `b` read from
# `parameters[offset+1]` and `parameters[offset+2]`.
struct LinearControl{P} <: ParameterizedFunction
    parameters::P
    offset::Int
end
LinearControl(parameters) = LinearControl(parameters, 0)

function (control::LinearControl)(t)
    p = control.parameters
    k = control.offset
    return p[begin+k] + p[begin+k+1] * t
end


# A custom array type that is a window into a vector, without a `Base.dataids`
# method of its own. To `Base.mightalias`, two windows into the same vector
# look like independent memory.
struct Window <: AbstractVector{Float64}
    data::Vector{Float64}
    offset::Int
    length::Int
end
Base.size(w::Window) = (w.length,)
Base.getindex(w::Window, i::Int) = w.data[w.offset+i]
Base.setindex!(w::Window, v, i::Int) = (w.data[w.offset+i] = v)


# A custom generator type that holds a parameter block for its drift and
# contains controls.
struct ScaledDriftGenerator{OT,CT}
    H₀::OT
    parameters::Vector{Float64}  # [scale]
    H₁::OT
    control::CT
end

QuantumPropagators.Controls.get_controls(G::ScaledDriftGenerator) = (G.control,)

function QuantumPropagators.Controls.get_parameters(G::ScaledDriftGenerator)
    return combine_parameters(
        G.parameters,
        (get_parameters(control) for control in get_controls(G))...
    )
end

function QuantumPropagators.Controls.evaluate(
    G::ScaledDriftGenerator,
    args...;
    vals_dict = IdDict()
)
    return G.parameters[1] * G.H₀ + evaluate(G.control, args...; vals_dict) * G.H₁
end


@testset "combine_parameters" begin
    a = [1.0, 2.0]
    b = [3.0]
    c = ComponentVector(x = 4.0, y = 5.0)
    empty = combine_parameters()
    @test empty isa ArrayPartition{Float64,Tuple{}}
    @test empty == Float64[]
    @test combine_parameters(Float64[], Float64[]) === empty
    @test combine_parameters(a) isa ArrayPartition
    @test combine_parameters(a).x === (a,)
    # A block is never mistaken for parameters (`c` has a component named `x`)
    @test c.x == 4.0
    @test combine_parameters(c).x === (c,)
    @test combine_parameters(a, Float64[], a) === combine_parameters(a)
    p = combine_parameters(a, b, a, c)
    @test p isa ArrayPartition
    @test p.x === (a, b, c)
    @test collect(p) == [1.0, 2.0, 3.0, 4.0, 5.0]
    # The result is flat
    p2 = combine_parameters(combine_parameters(a, b), combine_parameters(b, c))
    @test all(x -> !(x isa ArrayPartition), p2.x)
    @test p2.x === (a, b, c)
    @test p2 === p  # ArrayPartition is immutable
    # Round trip
    @test combine_parameters(p.x...) === p
    @test combine_parameters(empty.x...) === empty
    # A nested ArrayPartition is not a block
    nested = ArrayPartition(ArrayPartition(a, b), c)
    @test_throws ArgumentError combine_parameters(nested)
    # Writing to the combined parameters writes into the blocks
    p[3] = 30.0
    @test b[1] == 30.0
end


@testset "combine_parameters with aliasing blocks" begin
    x = collect(1.0:6.0)
    # overlapping views
    @test_throws ArgumentError combine_parameters(view(x, 1:4), view(x, 3:6))
    # non-overlapping views are fine
    p = combine_parameters(view(x, 1:3), view(x, 4:6))
    @test length(p) == 6
    @test length(p.x) == 2
    # two ComponentVectors over the same data
    data = [1.0, 2.0, 3.0]
    c1 = ComponentVector(data, Axis(A = 1, t₀ = 2, σ = 3))
    c2 = ComponentVector(data, Axis(ω = 1, ϕ = 2, E₀ = 3))
    @test_throws ArgumentError combine_parameters(c1, c2)
    # an array and a view into it
    @test_throws ArgumentError combine_parameters(x, view(x, 5:6))
end


@testset "issue #148 example" begin
    # Two generators with the controls `(ϵ1, ϵ2)` and `(ϵ2, ϵ3)`, each with
    # three parameters.
    ϵ1 = StaticParameter([1.0, 2.0, 3.0])
    ϵ2 = StaticParameter([4.0, 5.0, 6.0])
    ϵ3 = StaticParameter([7.0, 8.0, 9.0])
    H₀ = [1.0 0.0; 0.0 -1.0]
    H₁ = [0.0 1.0; 1.0 0.0]
    H₂ = [0.0 -1im; 1im 0.0]
    G1 = hamiltonian(H₀, (H₁, ϵ1), (H₂, ϵ2))
    G2 = hamiltonian(H₀, (H₁, ϵ2), (H₂, ϵ3))
    @test length(get_parameters(G1)) == 6
    @test length(get_parameters(G2)) == 6
    p = combine_parameters(get_parameters(G1), get_parameters(G2))
    @test length(p) == 9
    @test collect(p) == collect(1.0:9.0)
    @test p.x === (ϵ1.parameters, ϵ2.parameters, ϵ3.parameters)
end


@testset "overlapping blocks without dataids" begin
    data = [1.0, 2.0, 3.0]
    ϵ1 = LinearControl(Window(data, 0, 2))  # data[1:2]
    ϵ2 = LinearControl(Window(data, 1, 2))  # data[2:3]
    @test check_parameterized_function(ϵ1; tlist = [0.0, 1.0])
    @test check_parameterized_function(ϵ2; tlist = [0.0, 1.0])
    H₀ = [1.0 0.0; 0.0 -1.0]
    H₁ = [0.0 1.0; 1.0 0.0]
    H = hamiltonian(H₀, (H₀, ϵ1), (H₁, ϵ2))
    # Base.mightalias does not catch the overlap…
    p = get_parameters(H)
    @test length(p) == 4
    # … but check_parameterized does
    captured = IOCapture.capture(passthrough = false) do
        check_parameterized(H)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "Parameter blocks must not share memory")
    @test data == [1.0, 2.0, 3.0]  # restored
    # With independent windows, the check passes
    data2 = [1.0, 2.0, 3.0, 4.0]
    H2 = hamiltonian(
        H₀,
        (H₀, LinearControl(Window(data2, 0, 2))),
        (H₁, LinearControl(Window(data2, 2, 2)))
    )
    @test check_parameterized(H2)
    @test data2 == [1.0, 2.0, 3.0, 4.0]
end


@testset "check_parameterized_function for invalid blocks" begin
    tlist = [0.0, 1.0]
    @test check_parameterized_function(LinearControl([1.0, 2.0]); tlist)
    ϵ_ap = LinearControl(ArrayPartition([1.0], [2.0]))
    captured = IOCapture.capture(passthrough = false) do
        check_parameterized_function(ϵ_ap; tlist)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "must be a single parameter block")
    ϵ_offset = LinearControl(OffsetVector([1.0, 2.0], 0:1))
    @test ϵ_offset(1.0) == 3.0
    captured = IOCapture.capture(passthrough = false) do
        check_parameterized_function(ϵ_offset; tlist)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "must use one-based indexing")
    ϵ_int = LinearControl([1, 2])
    captured = IOCapture.capture(passthrough = false) do
        check_parameterized_function(ϵ_int; tlist)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "must be a vector of Float64")
end


@testset "custom generator holding a block" begin
    rng = StableRNG(1266426139)
    N = 4
    H₀ = random_matrix(N; hermitian = true, rng)
    H₁ = random_matrix(N; hermitian = true, rng)
    ϵ = LinearControl([0.5, 0.1])
    G = ScaledDriftGenerator(H₀, [2.0], H₁, ϵ)
    p = get_parameters(G)
    @test p.x === (G.parameters, ϵ.parameters)
    @test check_parameterized(G)
    p .= [3.0, 0.2, 0.0]
    @test norm(evaluate(G, 1.0) - (3.0 * H₀ + 0.2 * H₁)) < 1e-12
    # a control holding the same block as the generator
    G2 = ScaledDriftGenerator(H₀, ϵ.parameters, H₁, ϵ)
    @test get_parameters(G2).x === (ϵ.parameters,)
end


@testset "StaticParameter" begin
    w = [0.3, -0.1]
    ν = StaticParameter(w, 1)
    α = StaticParameter(w, 2)
    @test StaticParameter(w).index == 1
    @test ν(0.0) == 0.3
    @test α(10.0) == -0.1
    @test get_parameters(ν).x === (w,)
    @test get_parameters(ν) === get_parameters(α)  # same block
    @test get_controls(ν) === (ν,)
    @test repr(α) == "StaticParameter([0.3, -0.1], 2)"
    @test_throws BoundsError StaticParameter(w, 3)
    tlist = collect(range(0, 1; length = 11))
    @test check_parameterized_function(ν; tlist)
    w[1] = 0.4
    @test ν(0.5) == 0.4
    @test evaluate(ν, tlist, 3) == 0.4
    @test evaluate(ν, 0.5) == 0.4
    w[1] = 0.3

    rng = StableRNG(2683712180)
    N = 4
    H₀ = random_matrix(N; hermitian = true, rng)
    H₁ = random_matrix(N; hermitian = true, rng)
    Ψ₀ = random_state_vector(N; rng)
    H = hamiltonian((H₀, ν), (H₁, α))
    @test get_parameters(H) === get_parameters(ν)
    @test check_generator(
        H;
        state = Ψ₀,
        tlist,
        for_time_continuous = true,
        for_parameterization = true
    )
    T = tlist[end]
    Ψ = propagate(Ψ₀, H, tlist; method = ExpProp)
    Ψ_exact = exp(-1im * (w[1] * H₀ + w[2] * H₁) * T) * Ψ₀
    @test norm(Ψ - Ψ_exact) < 1e-10
    # Mutating the parameters changes the dynamics
    get_parameters(H) .= [0.6, 0.2]
    Ψ = propagate(Ψ₀, H, tlist; method = ExpProp)
    Ψ_exact = exp(-1im * (0.6 * H₀ + 0.2 * H₁) * T) * Ψ₀
    @test norm(Ψ - Ψ_exact) < 1e-10
end


# A mutable wrapper around a vector: a new wrapper is not identical (`===`)
# to an earlier one, even if it wraps the same data
mutable struct MutableWindow <: AbstractVector{Float64}
    data::Vector{Float64}
end
Base.size(w::MutableWindow) = size(w.data)
Base.getindex(w::MutableWindow, i::Int) = w.data[i]
Base.setindex!(w::MutableWindow, v, i::Int) = (w.data[i] = v)

struct FreshParametersGenerator
    data::Vector{Float64}
    wrap::Bool
end
QuantumPropagators.Controls.get_controls(::FreshParametersGenerator) = ()
function QuantumPropagators.Controls.get_parameters(G::FreshParametersGenerator)
    return G.wrap ? combine_parameters(MutableWindow(G.data)) : G.data
end


@testset "get_parameters contract" begin
    # repeated calls must return identical results
    G = FreshParametersGenerator([1.0, 2.0], true)
    @test get_parameters(G) == [1.0, 2.0]
    @test get_parameters(G) !== get_parameters(G)
    captured = IOCapture.capture(passthrough = false) do
        check_parameterized(G)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "must return identical (`===`) results")
    # the result must be an ArrayPartition
    G = FreshParametersGenerator([1.0, 2.0], false)
    captured = IOCapture.capture(passthrough = false) do
        check_parameterized(G)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "must return an ArrayPartition")
    # objects without parameters
    @test get_parameters(t -> 1.0) === get_parameters(zeros(10))
    @test isempty(get_parameters(hamiltonian([1.0 0.0; 0.0 -1.0])))
end


@testset "check_parameterized edge cases" begin
    H₀ = [1.0 0.0; 0.0 -1.0]
    H₁ = [0.0 1.0; 1.0 0.0]
    # The original values are restored even if writing fails partway
    a = [5.0, 6.0]
    H = hamiltonian(H₀, (H₀, LinearControl(a)), (H₁, LinearControl(7.0:8.0)))
    captured = IOCapture.capture(passthrough = false) do
        check_parameterized(H)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "Cannot mutate")
    @test a == [5.0, 6.0]
    # quiet=true suppresses all messages, also from nested checks
    captured = IOCapture.capture(passthrough = false) do
        check_generator(
            H;
            state = ComplexF64[1, 0],
            tlist = [0.0, 1.0],
            for_parameterization = true,
            quiet = true
        )
    end
    @test captured.value ≡ false
    @test captured.output == ""
    # Offset indexing results in `false`, not an exception
    ϵ_offset = LinearControl(OffsetVector([1.0, 2.0], 0:1))
    captured = IOCapture.capture(passthrough = false) do
        check_parameterized(ϵ_offset)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "must use one-based indexing")
    # An empty ArrayPartition is not a block
    ϵ_empty = LinearControl(ArrayPartition(Float64[]))
    captured = IOCapture.capture(passthrough = false) do
        check_parameterized_function(ϵ_empty; tlist = [0.0, 1.0])
    end
    @test captured.value ≡ false
    @test contains(captured.output, "must be a single parameter block")
    # Every block is checked, not only the combined parameters (whose `eltype`
    # is promoted to `Float64`)
    H = hamiltonian(H₀, (H₀, LinearControl([1.0, 2.0])), (H₁, LinearControl(Float32[3, 4])))
    @test eltype(get_parameters(H)) == Float64
    captured = IOCapture.capture(passthrough = false) do
        check_parameterized(H)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "parameter block 2")
    @test contains(captured.output, "must be a vector of Float64")
    # A control in a generator that holds an `ArrayPartition`
    ϵ_ap = LinearControl(ArrayPartition([3.0], [4.0]))
    H = hamiltonian(H₀, (H₀, LinearControl([1.0, 2.0])), (H₁, ϵ_ap))
    captured = IOCapture.capture(passthrough = false) do
        check_parameterized(H)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "must be a single parameter block")
end
