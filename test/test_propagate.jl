using Test
using QuantumPropagators
using QuantumPropagators: Cheby
using QuantumPropagators.Storage
using ExponentialUtilities
using LinearAlgebra
using UnicodePlots
using StaticArrays: @SMatrix, SVector, @SVector

@testset "TLS Rabi Cycling" begin

    # We use Rabi cycling in a TLS as a test case, which allows to to compare
    # the propagation with the known analytic solution.

    SHOWPLOT = false

    Ψ0 = ComplexF64[1, 0]
    Ĥ = ComplexF64[
         0   0.5
        0.5   0
    ]
    tlist = collect(range(0, 1.5π, length = 101)) # 3π/2 pulse

    generator = (Ĥ,)

    storage = init_storage(Ψ0, tlist)
    @test isa(storage, Matrix)
    @test eltype(storage) == ComplexF64

    Ψ_out = propagate(Ψ0, generator, tlist; method = :expprop, storage = storage)
    Ψ_out_expv = propagate(Ψ0, generator, tlist; method = ExponentialUtilities)
    Ψ_out_expv_ni =
        propagate(Ψ0, generator, tlist; method = ExponentialUtilities, inplace = false)
    generator_herm = (Hermitian(Ĥ),)
    Ψ_out_expv_herm =
        propagate(Ψ0, generator_herm, tlist; method = ExponentialUtilities, inplace = true)
    Ψ_expected = ComplexF64[-1/√2, -1im/√2]  # note the phases

    pop0 = abs.(storage[1, :]) .^ 2
    SHOWPLOT && println(lineplot(tlist ./ π, pop0, ylim = [0, 1], title = "fw prop"))

    @test norm(Ψ_out - Ψ_expected) < 1e-12
    @test norm(Ψ_out_expv - Ψ_expected) < 1e-9
    @test norm(Ψ_out_expv_ni - Ψ_expected) < 1e-9
    @test norm(Ψ_out_expv_herm - Ψ_expected) < 1e-9
    @test pop0[end] ≈ 0.5

    # Propagating backward in time should exactly reverse the dynamics (since
    # they are unitary). Thus, we should end up back at the initial state, and
    # the stored states (being filled back-to-front) should exactly match the
    # stored states from the forward propagation.

    storage_bw = init_storage(Ψ0, tlist)
    Ψ_out_bw = propagate(
        Ψ_out,
        generator,
        tlist;
        method = :expprop,
        backward = true,
        storage = storage_bw
    )

    pop0_bw = abs.(storage_bw[1, :]) .^ 2
    SHOWPLOT && println(lineplot(tlist ./ π, pop0_bw, ylim = [0, 1], title = "bw prop"))

    @test norm(Ψ_out_bw - Ψ0) < 1e-12
    @test pop0_bw[1] ≈ 1.0

    @test norm(storage - storage_bw) < 1e-12

end


@testset "Immutable TLS" begin

    # We use Rabi cycling in a TLS as a test case, which allows to to compare
    # the propagation with the known analytic solution.

    SHOWPLOT = false

    Ψ0 = @SVector ComplexF64[1, 0]
    Ĥ = @SMatrix ComplexF64[
         0   0.5
        0.5   0
    ]
    tlist = collect(range(0, 1.5π, length = 101)) # 3π/2 pulse

    generator = (Ĥ,)

    storage = init_storage(Ψ0, tlist)
    @test isa(storage, Vector)
    @test eltype(storage) == SVector{2,ComplexF64}

    Ψ_out = propagate(
        Ψ0,
        generator,
        tlist;
        inplace = false,
        method = :expprop,
        storage = storage
    )
    Ψ_expected = @SVector ComplexF64[-1/√2, -1im/√2]  # note the phases

    pop0 = abs2.(hcat((Array.(storage))...))
    SHOWPLOT && println(lineplot(tlist ./ π, pop0', ylim = [0, 1], title = "fw prop"))

    @test norm(Ψ_out - Ψ_expected) < 1e-12
    @test pop0[end] ≈ 0.5

    Ψ_out_cheby =
        propagate(Ψ0, generator, tlist; inplace = false, method = Cheby, storage = storage,)
    @test norm(Ψ_out_cheby - Ψ_expected) < 1e-12

    # Propagating backward in time should exactly reverse the dynamics (since
    # they are unitary). Thus, we should end up back at the initial state, and
    # the stored states (being filled back-to-front) should exactly match the
    # stored states from the forward propagation.

    storage_bw = init_storage(Ψ0, tlist)
    Ψ_out_bw = propagate(
        Ψ_out,
        generator,
        tlist;
        inplace = false,
        method = :expprop,
        backward = true,
        storage = storage_bw
    )

    pop0_bw = abs2.(hcat((Array.(storage))...))
    SHOWPLOT && println(lineplot(tlist ./ π, pop0_bw', ylim = [0, 1], title = "bw prop"))

    @test norm(Ψ_out_bw - Ψ0) < 1e-12
    @test pop0_bw[1] ≈ 1.0

    @test norm(storage - storage_bw) < 1e-12

    Ψ_out_bw_cheby = propagate(
        Ψ_out,
        generator,
        tlist;
        inplace = false,
        method = Cheby,
        backward = true,
        storage = storage_bw
    )
    @test norm(Ψ_out_bw_cheby - Ψ0) < 1e-12


end


@testset "Optomech propagation" begin
    include("optomech.jl")
    Ψ0 = ket(0, N_cav) ⊗ ket(2, N_mech)
    H = H_cav + H_mech + H_int
    generator = (H,)
    tlist = collect(range(0, 50, step = 0.2))
    Ψ1 = propagate(Ψ0, generator, tlist, method = :newton)
    @test (norm(Ψ1) - 1.0) < 1e-12
    Ψ2 = propagate(Ψ0, generator, tlist, method = :cheby)
    @test norm(Ψ1 - Ψ2) < 1e-10
end
