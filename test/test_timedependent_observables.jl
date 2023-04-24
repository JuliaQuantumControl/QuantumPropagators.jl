using Test
using LinearAlgebra
using StaticArrays
using QuantumPropagators: propagate
#using QuantumPropagators.Storage: init_storage, get_from_storage!

@testset "TLS Lab Frame" begin

    𝕚 = 1im
    ≈(a, b) = isapprox(a, b; atol=1e-12)
    ket(α, β) = SVector{2}(ComplexF64[α, β])

    Ψ₀ = ket(1, 0)
    Ω = 1
    Ĥ = SMatrix{2,2}([0 -Ω/2; -Ω/2 0])
    ω_l = 0.5  # rotating frame frequency

    σ̂_y = [
        0  -𝕚
        𝕚   0
    ]

    tlist = collect(range(0, 1.5π, length=101)) # 3π/2 pulse

    # states in rotating frame

    states = propagate(Ψ₀, Ĥ, tlist; storage=true, method=:expprop, inplace=false)
    g(t) = cos(Ω * t / 2)
    e(t) = 𝕚 * sin(Ω * t / 2)
    for i in eachindex(tlist)
        Ψ = states[i]
        t = tlist[i]
        @test norm(Ψ - ket(g(t), e(t))) ≈ 0.0
    end


    # expectation values in rotating frame

    expvals = propagate(
        Ψ₀,
        Ĥ,
        tlist;
        observables=(σ̂_y,),
        storage=true,
        method=:expprop,
        inplace=false
    )
    @test norm(expvals .- [sin(Ω * t) for t ∈ tlist]) ≈ 0.0


    # states in lab frame

    Û(t) = SMatrix{2,2}([   # lab frame to rotating frame
        1  0
        0  exp(-𝕚 * ω_l * t)
    ])

    Û⁺(t) = SMatrix{2,2}([  # rotating frame to lab frame
        1  0
        0  exp(𝕚 * ω_l * t)
    ])

    to_lab(Ψ, tlist, i) = Û⁺(tlist[i]) * Ψ
    states_lab = propagate(
        Ψ₀,
        Ĥ,
        tlist;
        observables=(to_lab,),
        storage=true,
        method=:expprop,
        inplace=false
    )
    for (Ψ_lab, Ψ_rot, t) in zip(states_lab, states, tlist)
        @test norm(Ψ_lab - ket(cos(Ω * t / 2), 𝕚 * exp(𝕚 * ω_l * t) * sin(Ω * t / 2))) ≈ 0.0
        @test norm(Û(t) * Ψ_lab - Ψ_rot) ≈ 0.0
    end


    # expectation values in lab frame

    function σ̄_y_lab(Ψ_rot, tlist, i)
        t = tlist[i]
        Ψ_lab = Û⁺(t) * Ψ_rot
        return dot(Ψ_lab, σ̂_y, Ψ_lab)
    end

    expvals_lab = propagate(
        Ψ₀,
        Ĥ,
        tlist;
        observables=(σ̄_y_lab,),
        storage=true,
        method=:expprop,
        inplace=false
    )

    @test norm(expvals_lab .- [sin(Ω * t) * cos(ω_l * t) for t ∈ tlist]) ≈ 0.0

end
