# SPDX-FileCopyrightText: © 2021 Michael Goerz <mail@michaelgoerz.net>
#
# SPDX-License-Identifier: MIT

using Test
using LinearAlgebra: Diagonal, norm
using QuantumPropagators.Shapes: flattop
using QuantumPropagators: ExpProp, hamiltonian, init_prop
using QuantumPropagators: propagate_sequence, Propagation
using IOCapture: IOCapture

const 𝕚 = 1im;

function ramsey_scheme(;
    T = 37,                 # duration of each pulse (μs)
    nt = 1001,              # number of time steps for each pulse
    τ = 1428,               # the free evolution time (μs) = 1.428ms
    μ = 1.0,                # extra prefactor for pulse amplitudes
    signs = (+1, +1),       # the signs of the π-pulses
    storage = nothing,      # optional: 4 storage objects
)

    tlist = collect(range(0, T; length = 1001))
    S1, S2 = signs
    @assert S1 in [-1, 1]
    @assert S2 in [-1, 1]

    propagations = [
        Propagation(
            ramsey_rwa_hamiltonian(; Ωₚ = (t -> μ * pi_shape(t; T)), Ωₛ = (t -> 0.0)),
            collect(range(0, T; length = nt));
            pre_propagation = apply_global_phase,
            post_propagation = nothing
        ),
        Propagation(
            ramsey_rwa_hamiltonian(;
                Ωₚ = (t -> (√2 / 2) * μ * pi_shape(t - T; T)),
                Ωₛ = (t -> (√2 / 2) * μ * pi_shape(t - T; T))
            ),
            collect(range(T, 2T; length = nt));
            post_propagation = ramsey_state_to_lab
        ),
        Propagation(
            ramsey_lab_free_hamiltonian(),
            [2T, 2T + τ];
            post_propagation = ramsey_state_to_rwa
        ),
        Propagation(
            ramsey_rwa_hamiltonian(;
                Ωₚ = (t -> S1 * (√2 / 2) * μ * pi_shape(t - (2T + τ); T)),
                Ωₛ = (t -> S2 * (√2 / 2) * μ * pi_shape(t - (2T + τ); T))
            ),
            collect(range(2T + τ, 3T + τ; length = nt));
            post_propagation = nothing
        )
    ]
    if !isnothing(storage)
        for (i, item) in enumerate(storage)
            propagations[i].kwargs[:storage] = item
        end
    end
    return propagations
end


function ramsey_scheme_pre_initialized(;
    T = 37,                 # duration of each pulse (μs)
    nt = 1001,              # number of time steps for each pulse
    τ = 1428,               # the free evolution time (μs) = 1.428ms
    μ = 1.0,                # extra prefactor for pulse amplitudes
    signs = (+1, +1),       # the signs of the π-pulses
)

    tlist = collect(range(0, T; length = 1001))
    S1, S2 = signs
    @assert S1 in [-1, 1]
    @assert S2 in [-1, 1]

    Ψ = zeros(ComplexF64, 3)  # storage for state in propagator

    propagators = [
        init_prop(
            copy(Ψ),
            ramsey_rwa_hamiltonian(; Ωₚ = (t -> μ * pi_shape(t; T)), Ωₛ = (t -> 0.0)),
            collect(range(0, T; length = nt));
            method = ExpProp,
            check = false,
        ),
        init_prop(
            copy(Ψ),
            ramsey_rwa_hamiltonian(;
                Ωₚ = (t -> (√2 / 2) * μ * pi_shape(t - T; T)),
                Ωₛ = (t -> (√2 / 2) * μ * pi_shape(t - T; T))
            ),
            collect(range(T, 2T; length = nt));
            method = ExpProp,
            check = false,
        ),
        init_prop(
            copy(Ψ),
            ramsey_lab_free_hamiltonian(),
            [2T, 2T + τ];
            method = ExpProp,
            check = false,
        ),
        init_prop(
            copy(Ψ),
            ramsey_rwa_hamiltonian(;
                Ωₚ = (t -> S1 * (√2 / 2) * μ * pi_shape(t - (2T + τ); T)),
                Ωₛ = (t -> S2 * (√2 / 2) * μ * pi_shape(t - (2T + τ); T))
            ),
            collect(range(2T + τ, 3T + τ; length = nt));
            method = ExpProp,
            check = false,
        ),
    ]

    propagations = [
        Propagation(propagators[1]; post_propagation = nothing),
        Propagation(propagators[2]; post_propagation = ramsey_state_to_lab),
        Propagation(propagators[3]; post_propagation = ramsey_state_to_rwa),
        Propagation(propagators[4]; post_propagation = nothing)
    ]

    return propagations

end


function ramsey_U_RWA(
    t;
    f_DQ = 0.293 * 2π,
    f_RF = 4.94 * 2π,
    f1 = (f_RF + f_DQ / 2),
    f2 = (f_RF - f_DQ / 2)
)
    Diagonal([1.0, exp(𝕚 * f1 * t), exp(𝕚 * (f1 - f2) * t)])
end


function apply_global_phase(Ψ, generator, tlist; _...)
    if get(ENV, "DEBUG_RAMSEY_STATES", "") == "true"
        println("Applying phase factor")
    end
    return 𝕚 * Ψ
end

function ramsey_state_to_lab(Ψ, generator, tlist; _...)
    U = ramsey_U_RWA(tlist[end])
    if get(ENV, "DEBUG_RAMSEY_STATES", "") == "true"
        println("RWA frame $(round.(Ψ; digits=2)) → lab frame $(round.(U'*Ψ; digits=2))")
    end
    return U' * Ψ
end

function ramsey_state_to_lab(Ψ, propagator; _...)
    U = ramsey_U_RWA(propagator.tlist[end])
    return U' * Ψ
end

function ramsey_state_to_rwa(Ψ, generator, tlist; _...)
    U = ramsey_U_RWA(tlist[end])
    if get(ENV, "DEBUG_RAMSEY_STATES", "") == "true"
        println("lab frame $(round.(Ψ; digits=2)) → RWA frame $(round.(U*Ψ; digits=2))")
    end
    return U * Ψ
end

function ramsey_state_to_rwa(Ψ, propagator; _...)
    U = ramsey_U_RWA(propagator.tlist[end])
    return U * Ψ
end


pi_shape(t; T) = (2π / T) * flattop(t; T, t_rise = 0.5 * T, func = :sinsq)


function ramsey_rwa_hamiltonian(; Δ = 0.0, δ = 0.0, ν = 0.0, Ωₚ, Ωₛ)
    H₀ = [
        0 0 0
        0 Δ+ν 0
        0 0  δ+2*ν
    ]
    H₁ = [
        0    0.5 0
        0.5  0   0
        0    0   0
    ]
    H₂ = [
        0   0    0
        0   0    0.5
        0   0.5  0
    ]
    return hamiltonian(H₀, (H₁, Ωₚ), (H₂, Ωₛ))
end


function ramsey_lab_free_hamiltonian(;
    f_DQ = 0.293 * 2π,
    f_RF = 4.94 * 2π,
    ν = 0.0,
    Δ = 0.0,
    δ = 0.0
)
    Diagonal([0.0, f_DQ / 2 + f_RF + Δ + ν, f_DQ + δ + 2ν])
end


@testset "Ramsey Sequence Propagation" begin

    Ψ₀ = ComplexF64[1, 0, 0]
    ENV["DEBUG_RAMSEY_STATES"] = "true"
    captured = IOCapture.capture() do
        propagate_sequence(
            Ψ₀,
            ramsey_scheme(signs = (+1, +1), τ = 1380),
            method = ExpProp,
            check = false,
        )
    end
    ENV["DEBUG_RAMSEY_STATES"] = ""
    result = captured.value
    @test captured.output == """
    Applying phase factor
    RWA frame ComplexF64[0.0 - 0.71im, 0.0 + 0.0im, 0.0 - 0.71im] → lab frame ComplexF64[0.0 - 0.71im, -0.0 - 0.0im, 0.64 + 0.29im]
    lab frame ComplexF64[0.0 - 0.71im, 0.0 + 0.0im, -0.1 - 0.7im] → RWA frame ComplexF64[0.0 - 0.71im, 0.0 + 0.0im, 0.0 - 0.71im]
    """
    @test result isa Vector
    @test length(result) == 4
    Ψout = result[end]
    @test Ψout isa Vector{ComplexF64}
    #=
    println("[1] After initial π pulse:")
    @show round.(result[1]; digits=2)
    println("[2] After π/2 pulse (transformed to lab):")
    @show round.(result[2]; digits=2)
    println("[3] After free time evolution (transformed to RWA):")
    @show round.(result[3]; digits=2)
    println("[4] After π/2 pulse:")
    @show round.(result[4]; digits=2)
    =#

    result = propagate_sequence(
        Ψ₀,
        ramsey_scheme(signs = (+1, +1), τ = 1380),
        method = ExpProp,
        check = false,
        storage = true,
    )
    @test result isa Vector
    @test length(result) == 4
    @test result[end] isa Matrix{ComplexF64}
    @test size(result[end]) == (3, 1001)
    @test norm(result[end][:, end] - Ψout) < 1e-14

    P₁ = Float64[1 0 0; 0 0 0; 0 0 0]
    P₂ = Float64[0 0 0; 0 1 0; 0 0 0]
    P₃ = Float64[0 0 0; 0 0 0; 0 0 1]
    result = propagate_sequence(
        Ψ₀,
        ramsey_scheme(signs = (+1, +1), τ = 1380),
        method = ExpProp,
        check = false,
        storage = true,
        observables = [P₁, P₂, P₃]
    )
    @test result isa Vector
    @test length(result) == 4
    @test result[end] isa Matrix{ComplexF64}
    @test size(result[end]) == (3, 1001)
    @test norm(result[end][:, end] - abs2.(Ψout)) < 1e-14

    storage = [
        zeros(Float64, 3, 1001),
        zeros(Float64, 3, 1001),
        zeros(Float64, 3, 1001),
        zeros(Float64, 3, 1001)
    ]
    result = propagate_sequence(
        Ψ₀,
        ramsey_scheme(signs = (+1, +1), τ = 1380, storage = storage),
        method = ExpProp,
        check = false,
        observables = [Ψ -> abs2.(Ψ)]
    )
    @test result isa Vector
    @test length(result) == 4
    @test result[end] isa Vector{ComplexF64}
    @test norm(result[end] - Ψout) < 1e-12
    @test norm(storage[end][:, end] - abs2.(Ψout)) < 1e-12

    result = propagate_sequence(
        𝕚 * Ψ₀,
        ramsey_scheme_pre_initialized(signs = (+1, +1), τ = 1380),
    )
    @test result isa Vector
    @test length(result) == 4
    @test norm(result[end] - Ψout) < 1e-12

end
