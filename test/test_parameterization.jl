# SPDX-FileCopyrightText: © 2021 Michael Goerz <mail@michaelgoerz.net>
#
# SPDX-License-Identifier: MIT

using Test
using ComponentArrays
using RecursiveArrayTools  # to ensure extension is loaded
using UnPack: @unpack
using QuantumPropagators: hamiltonian
using QuantumPropagators.Controls: ParameterizedFunction, get_parameters, get_controls
using QuantumPropagators.Interfaces: check_parameterized_function
using QuantumPropagators.Interfaces: check_parameterized
using QuantumPropagators.Interfaces: check_generator

struct GaussianControl <: ParameterizedFunction
    parameters::ComponentVector{Float64,Vector{Float64},Tuple{Axis{(A = 1, t₀ = 2, σ = 3)}}}
    GaussianControl(; kwargs...) = new(ComponentVector(; kwargs...))
end

function (control::GaussianControl)(t)
    @unpack A, t₀, σ = control.parameters
    return A * exp(-(t - t₀)^2 / (2 * σ^2))
end


const 𝕚 = 1im


function total_enantiomer_ham(
    parameters;
    sign,
    a,
    independent_parameters = false,
    kwargs...
)

    μ = (sign == "-" ? -1 : 1)
    H₁Re = μ * ComplexF64[0 1 0; 1 0 0; 0  0 0]
    H₁Im = μ * ComplexF64[0 𝕚 0; -𝕚 0 0; 0  0 0]
    H₂Re = μ * ComplexF64[0 0 0; 0 0 1; 0  1 0]
    H₂Im = μ * ComplexF64[0 0 0; 0 0 𝕚; 0 -𝕚 0]
    H₃Re = μ * ComplexF64[0 0 1; 0 0 0; 1  0 0]
    H₃Im = μ * ComplexF64[0 0 𝕚; 0 0 0; -𝕚  0 0]

    if independent_parameters
        # This doesn't make sense physically, but it's a good way to test
        # collecting multiple parameter arrays
        return hamiltonian(
            (H₁Re, TEH_field1Re(copy(parameters), a)),
            (H₁Im, TEH_field1Im(copy(parameters), a)),
            (H₂Re, TEH_field2Re(copy(parameters), a)),
            (H₂Im, TEH_field2Im(copy(parameters), a)),
            (H₃Re, TEH_field3Re(copy(parameters), a)),
            (H₃Im, TEH_field3Im(copy(parameters), a));
            kwargs...
        )
    else
        return hamiltonian(
            (H₁Re, TEH_field1Re(parameters, a)),
            (H₁Im, TEH_field1Im(parameters, a)),
            (H₂Re, TEH_field2Re(parameters, a)),
            (H₂Im, TEH_field2Im(parameters, a)),
            (H₃Re, TEH_field3Re(parameters, a)),
            (H₃Im, TEH_field3Im(parameters, a));
            kwargs...
        )
    end

end

struct TEH_field1Re <: ParameterizedFunction
    parameters::ComponentVector{
        Float64,
        Vector{Float64},
        Tuple{
            Axis{(
                ΔT₁ = 1,
                ΔT₂ = 2,
                ΔT₃ = 3,
                ϕ₁ = 4,
                ϕ₂ = 5,
                ϕ₃ = 6,
                E₀₁ = 7,
                E₀₂ = 8,
                E₀₃ = 9
            )}
        }
    }
    a::Float64
end

function (E::TEH_field1Re)(t)
    @unpack E₀₁, ΔT₁, ϕ₁ = E.parameters
    _tanhfield(t; E₀ = E₀₁, t₁ = 0.0, t₂ = ΔT₁, a = E.a) * cos(ϕ₁)
end

struct TEH_field2Re <: ParameterizedFunction
    parameters::ComponentVector{
        Float64,
        Vector{Float64},
        Tuple{
            Axis{(
                ΔT₁ = 1,
                ΔT₂ = 2,
                ΔT₃ = 3,
                ϕ₁ = 4,
                ϕ₂ = 5,
                ϕ₃ = 6,
                E₀₁ = 7,
                E₀₂ = 8,
                E₀₃ = 9
            )}
        }
    }
    a::Float64
end

function (E::TEH_field2Re)(t)
    @unpack E₀₂, ΔT₁, ΔT₂, ϕ₂ = E.parameters
    _tanhfield(t; E₀ = E₀₂, t₁ = ΔT₁, t₂ = (ΔT₁ + ΔT₂), a = E.a) * cos(ϕ₂)
end

struct TEH_field3Re <: ParameterizedFunction
    parameters::ComponentVector{
        Float64,
        Vector{Float64},
        Tuple{
            Axis{(
                ΔT₁ = 1,
                ΔT₂ = 2,
                ΔT₃ = 3,
                ϕ₁ = 4,
                ϕ₂ = 5,
                ϕ₃ = 6,
                E₀₁ = 7,
                E₀₂ = 8,
                E₀₃ = 9
            )}
        }
    }
    a::Float64
end

function (E::TEH_field3Re)(t)
    @unpack E₀₃, ΔT₁, ΔT₂, ΔT₃, ϕ₃ = E.parameters
    _tanhfield(t; E₀ = E₀₃, t₁ = (ΔT₁ + ΔT₂), t₂ = (ΔT₁ + ΔT₂ + ΔT₃), a = E.a) * cos(ϕ₃)
end

struct TEH_field1Im <: ParameterizedFunction
    parameters::ComponentVector{
        Float64,
        Vector{Float64},
        Tuple{
            Axis{(
                ΔT₁ = 1,
                ΔT₂ = 2,
                ΔT₃ = 3,
                ϕ₁ = 4,
                ϕ₂ = 5,
                ϕ₃ = 6,
                E₀₁ = 7,
                E₀₂ = 8,
                E₀₃ = 9
            )}
        }
    }
    a::Float64
end

function (E::TEH_field1Im)(t)
    @unpack E₀₁, ΔT₁, ϕ₁ = E.parameters
    _tanhfield(t; E₀ = E₀₁, t₁ = 0.0, t₂ = ΔT₁, a = E.a) * sin(ϕ₁)
end

struct TEH_field2Im <: ParameterizedFunction
    parameters::ComponentVector{
        Float64,
        Vector{Float64},
        Tuple{
            Axis{(
                ΔT₁ = 1,
                ΔT₂ = 2,
                ΔT₃ = 3,
                ϕ₁ = 4,
                ϕ₂ = 5,
                ϕ₃ = 6,
                E₀₁ = 7,
                E₀₂ = 8,
                E₀₃ = 9
            )}
        }
    }
    a::Float64
end

function (E::TEH_field2Im)(t)
    @unpack E₀₂, ΔT₁, ΔT₂, ϕ₂ = E.parameters
    _tanhfield(t; E₀ = E₀₂, t₁ = ΔT₁, t₂ = (ΔT₁ + ΔT₂), a = E.a) * sin(ϕ₂)
end

struct TEH_field3Im <: ParameterizedFunction
    parameters::ComponentVector{
        Float64,
        Vector{Float64},
        Tuple{
            Axis{(
                ΔT₁ = 1,
                ΔT₂ = 2,
                ΔT₃ = 3,
                ϕ₁ = 4,
                ϕ₂ = 5,
                ϕ₃ = 6,
                E₀₁ = 7,
                E₀₂ = 8,
                E₀₃ = 9
            )}
        }
    }
    a::Float64
end

function (E::TEH_field3Im)(t)
    @unpack E₀₃, ΔT₁, ΔT₂, ΔT₃, ϕ₃ = E.parameters
    _tanhfield(t; E₀ = E₀₃, t₁ = (ΔT₁ + ΔT₂), t₂ = (ΔT₁ + ΔT₂ + ΔT₃), a = E.a) * sin(ϕ₃)
end

_tanhfield(t; E₀, t₁, t₂, a) = (E₀ / 2) * (tanh(a * (t - t₁)) - tanh(a * (t - t₂)));


@testset "get_parameters of controls" begin

    ϵ(t) = 1.0
    p = get_parameters(ϵ)
    @test isempty(p)
    @test check_parameterized(ϵ)

    ϵ = [0.0, 1.0, 0.0]
    p = get_parameters(ϵ)
    @test isempty(p)
    @test check_parameterized(ϵ)

    gaussian = GaussianControl(A = 1.0, t₀ = 0.5, σ = 0.2)
    p = get_parameters(gaussian)
    @test p isa AbstractVector
    @test eltype(p) == Float64
    @test length(p) == 3
    @test Array(p) == [1.0, 0.5, 0.2]
    @test check_parameterized_function(gaussian; tlist = [0.0, 1.0])

end


@testset "enantiomer_ham" begin

    parameters = ComponentVector(
        ΔT₁ = 0.3,
        ΔT₂ = 0.4,
        ΔT₃ = 0.3,
        ϕ₁ = 0.0,
        ϕ₂ = 0.0,
        ϕ₃ = 0.0,
        E₀₁ = 4.5,
        E₀₂ = 4.0,
        E₀₃ = 5.0
    )
    ϵ = TEH_field1Re(parameters, 100.0)
    tlist = [0.0, 0.5, 1.0]
    @test check_parameterized_function(ϵ; tlist)
    H = total_enantiomer_ham(parameters; sign = "+", a = 100)
    @test length(get_controls(H)) == 6
    Ψ₀ = ComplexF64[1, 0, 0]
    @test check_generator(H; state = Ψ₀, tlist, for_parameterization = true)
    @test check_parameterized(H)
    @test get_parameters(H) === parameters

end


@testset "enantiomer_ham - collect independent parameters" begin

    parameters = ComponentVector(
        ΔT₁ = 0.3,
        ΔT₂ = 0.4,
        ΔT₃ = 0.3,
        ϕ₁ = 0.0,
        ϕ₂ = 0.0,
        ϕ₃ = 0.0,
        E₀₁ = 4.5,
        E₀₂ = 4.0,
        E₀₃ = 5.0
    )
    H = total_enantiomer_ham(parameters; independent_parameters = true, sign = "+", a = 100)
    @test length(get_controls(H)) == 6
    Ψ₀ = ComplexF64[1, 0, 0]
    tlist = [0.0, 0.5, 1.0]
    @test check_generator(H; state = Ψ₀, tlist, for_parameterization = true)
    @test check_parameterized(H)
    p = get_parameters(H)
    @test length(p) == length(get_controls(H)) * length(parameters)

end
