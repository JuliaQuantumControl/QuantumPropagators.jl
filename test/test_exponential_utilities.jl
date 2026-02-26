using Test

using QuantumPropagators.Controls: evaluate
using QuantumPropagators.Generators: ScaledOperator
using StableRNGs: StableRNG
using QuantumControlTestUtils.RandomObjects:
    random_state_vector, random_dynamic_generator, random_matrix
using LinearAlgebra: norm, ishermitian
import ExponentialUtilities
import StaticArrays: SMatrix, SVector
import OffsetArrays: OffsetVector
using BenchmarkTools: @btimed

RUN_BENCHMARKS = false


@testset "expv with dense complex non-Hermitian matrix" begin

    N = 10
    rng = StableRNG(677968056)
    A = random_matrix(N; rng, complex = true, hermitian = false)
    Ψ₀ = random_state_vector(N; rng)
    dt = 0.1
    t = -1im * dt
    Ψ = ExponentialUtilities.expv(t, A, Ψ₀)
    Ψ_expected = exp(t * A) * Ψ₀
    @test norm(Ψ - Ψ_expected) < 1e-12

end


@testset "expv with dense real non-Hermitian matrix" begin

    N = 10
    rng = StableRNG(2945201349)
    A = random_matrix(N; rng, complex = false, hermitian = false)
    Ψ₀ = random_state_vector(N; rng)
    dt = 0.1
    t = -1im * dt
    Ψ = ExponentialUtilities.expv(t, A, Ψ₀)
    Ψ_expected = exp(t * A) * Ψ₀
    @test norm(Ψ - Ψ_expected) < 1e-12

end


@testset "expv with dense complex Hermitian matrix and complex/real time" begin

    # Hermitian A triggers Lanczos instead of Arnoldi
    N = 10
    rng = StableRNG(697691674)
    A = random_matrix(N; rng, hermitian = true)
    @test ishermitian(A)
    Ψ₀ = random_state_vector(N; rng)
    dt = 0.1
    t = -1im * dt
    Ψ = ExponentialUtilities.expv(t, A, Ψ₀)
    Ψ_expected = exp(t * A) * Ψ₀
    @test norm(Ψ - Ψ_expected) < 1e-12
    if RUN_BENCHMARKS
        b1 = @btimed ExponentialUtilities.expv($t, $A, $Ψ₀)
    end

    # We can also check against the code path where we explicitly disable
    # `ishermitian` detection
    Ψ = ExponentialUtilities.expv(t, A, Ψ₀; ishermitian = false)
    @test norm(Ψ - Ψ_expected) < 1e-12
    if RUN_BENCHMARKS
        b2 = @btimed ExponentialUtilities.expv($t, $A, $Ψ₀; ishermitian = false)
        @test b2.time > b1.time
    end

    # Instead of using complex time, we can also use real time and absorb the
    # factor `-1im` in `A`
    O = ScaledOperator(-1im, A)
    @test !ishermitian(O)
    Ψ = ExponentialUtilities.expv(dt, O, Ψ₀)
    @test norm(Ψ - Ψ_expected) < 1e-12
    if RUN_BENCHMARKS
        b3 = @btimed ExponentialUtilities.expv($dt, $O, $Ψ₀)
        @test b3.time > b1.time
    end

end


@testset "expv with dense real symmetric matrix" begin

    # Symmetric (Hermitian) `A` triggers Lanczos instead of Arnoldi
    N = 10
    rng = StableRNG(1056344208)
    A = random_matrix(N; rng, complex = false, hermitian = true)
    @test ishermitian(A)
    Ψ₀ = random_state_vector(N; rng)
    dt = 0.1
    t = -1im * dt
    Ψ = ExponentialUtilities.expv(t, A, Ψ₀)
    Ψ_expected = exp(t * A) * Ψ₀
    @test norm(Ψ - Ψ_expected) < 1e-12
    if RUN_BENCHMARKS
        b1 = @btimed ExponentialUtilities.expv($t, $A, $Ψ₀)
    end
    # We can also check against the code path where we explicitly disable
    # `ishermitian` detection
    Ψ = ExponentialUtilities.expv(t, A, Ψ₀; ishermitian = false)
    @test norm(Ψ - Ψ_expected) < 1e-12
    if RUN_BENCHMARKS
        b2 = @btimed ExponentialUtilities.expv($t, $A, $Ψ₀; ishermitian = false)
        @test b2.time > b1.time
    end

end


@testset "expv with sparse complex non-Hermitian matrix" begin

    N = 50
    rng = StableRNG(1788065498)
    A = random_matrix(N; rng, density = 0.3, hermitian = false)
    @test !ishermitian(A)
    Ψ₀ = random_state_vector(N; rng)
    dt = 0.1
    t = -1im * dt
    Ψ = ExponentialUtilities.expv(t, A, Ψ₀)
    Ψ_expected = exp(t * Matrix(A)) * Ψ₀
    @test norm(Ψ - Ψ_expected) < 1e-12

end


@testset "expv with sparse complex Hermitian matrix" begin

    # Sparse Hermitian A triggers Lanczos
    N = 50
    rng = StableRNG(2648750134)
    A = random_matrix(N; rng, density = 0.3, hermitian = true)
    @test ishermitian(A)
    Ψ₀ = random_state_vector(N; rng)
    dt = 0.1
    t = -1im * dt
    Ψ = ExponentialUtilities.expv(t, A, Ψ₀)
    Ψ_expected = exp(t * Matrix(A)) * Ψ₀
    @test norm(Ψ - Ψ_expected) < 1e-12

end


@testset "expv with Operator" begin

    N = 10
    rng = StableRNG(677968056)
    tlist = collect(range(0, 10, length = 101))
    H = random_dynamic_generator(N, tlist; rng)
    Ψ₀ = random_state_vector(N; rng)
    A = evaluate(H, tlist, 1)
    dt = tlist[2] - tlist[1]
    t = -1im * dt
    Ψ = ExponentialUtilities.expv(t, A, Ψ₀)
    Ψ_expected = exp(t * Array(A)) * Ψ₀
    @test norm(Ψ - Ψ_expected) < 1e-12

end


@testset "expv with StaticArrays SMatrix and SVector" begin

    N = 4
    rng = StableRNG(814242294)
    A = SMatrix{N,N}(random_matrix(N; rng, hermitian = true))
    Ψ₀ = SVector{N}(random_state_vector(N; rng))
    dt = 0.1
    t = -1im * dt
    Ψ = ExponentialUtilities.expv(t, A, Ψ₀)
    Ψ_expected = exp(t * Matrix(A)) * Vector(Ψ₀)
    @test Ψ isa SVector{N}
    @test norm(Ψ - Ψ_expected) < 1e-12

end


@testset "expv with pre-built KrylovSubspace" begin

    # arnoldi builds the Krylov subspace once; expv(t, Ks) reuses it,
    # which is efficient when applying the same operator at many time steps
    N = 10
    rng = StableRNG(3750746111)
    A = random_matrix(N; rng, hermitian = true)
    Ψ₀ = random_state_vector(N; rng)
    dt = 0.1
    t = -1im * dt
    Ks = ExponentialUtilities.arnoldi(A, Ψ₀)
    Ψ = ExponentialUtilities.expv(t, Ks)
    Ψ_expected = exp(t * A) * Ψ₀
    @test norm(Ψ - Ψ_expected) < 1e-12

end


@testset "expv! in-place with pre-built KrylovSubspace" begin

    N = 10
    rng = StableRNG(1187394246)
    A = random_matrix(N; rng, hermitian = true)
    Ψ₀ = random_state_vector(N; rng)
    dt = 0.1
    t = -1im * dt
    Ks = ExponentialUtilities.arnoldi(A, Ψ₀)
    Ψ = similar(Ψ₀, promote_type(typeof(t), eltype(A), eltype(Ψ₀)))
    ExponentialUtilities.expv!(Ψ, t, Ks)
    Ψ_expected = exp(t * A) * Ψ₀
    @test norm(Ψ - Ψ_expected) < 1e-12

end


@testset "expv! in-place with dense complex Hermitian matrix and complex/real time" begin

    N = 10
    rng = StableRNG(3503257612)
    A = random_matrix(N; rng, hermitian = true)
    @test ishermitian(A)
    Ψ₀ = random_state_vector(N; rng)
    dt = 0.1
    t = -1im * dt
    Ψ_expected = exp(t * A) * Ψ₀

    # Baseline: expv
    Ψ = ExponentialUtilities.expv(t, A, Ψ₀)
    @test norm(Ψ - Ψ_expected) < 1e-12
    if RUN_BENCHMARKS
        b0 = @btimed ExponentialUtilities.expv($t, $A, $Ψ₀)
    end

    # Lanczos path: arnoldi detects ishermitian(A) and builds a symmetric
    # tridiagonal H; expv! then uses the fast SymTridiagonal eigen path.
    # For Hermitian A the H matrix is real: ExpvCache element type = real(eltype(A)).
    Ks1 = ExponentialUtilities.arnoldi(A, Ψ₀)
    w = similar(Ψ₀)
    ExponentialUtilities.expv!(w, t, Ks1)
    @test norm(w - Ψ_expected) < 1e-12
    if RUN_BENCHMARKS
        cache1 = ExponentialUtilities.ExpvCache{real(eltype(A))}(Ks1.m)
        b1 = @btimed ExponentialUtilities.expv!($w, $t, $Ks1; cache = $cache1)
        @test b1.time < b0.time
    end

    # Arnoldi path: force full Arnoldi; H is upper Hessenberg, so expv! falls
    # back to the generic Higham matrix exponential, which is slower.
    # Arnoldi (including when forced on a Hermitian A) builds complex H:
    # ExpvCache element type = eltype(A).
    Ks2 = ExponentialUtilities.arnoldi(A, Ψ₀; ishermitian = false)
    ExponentialUtilities.expv!(w, t, Ks2)
    @test norm(w - Ψ_expected) < 1e-12
    if RUN_BENCHMARKS
        cache2 = ExponentialUtilities.ExpvCache{eltype(A)}(Ks2.m)
        b2 = @btimed ExponentialUtilities.expv!($w, $t, $Ks2; cache = $cache2)
        # The timings here are a lot flakier, so we don't
        # @test b2.time > b1.time
    end

    # Real-time path: absorb `-1im` into the operator via ScaledOperator;
    # the scaled operator is not Hermitian so arnoldi uses Arnoldi iteration,
    # giving an upper Hessenberg H (ComplexF64), same as Ks2.
    O = ScaledOperator(-1im, A)
    @test !ishermitian(O)
    Ks3 = ExponentialUtilities.arnoldi(O, Ψ₀)
    ExponentialUtilities.expv!(w, dt, Ks3)
    @test norm(w - Ψ_expected) < 1e-12
    if RUN_BENCHMARKS
        cache3 = ExponentialUtilities.ExpvCache{eltype(O)}(Ks3.m)
        b3 = @btimed ExponentialUtilities.expv!($w, $dt, $Ks3; cache = $cache3)
        # The timings here are a lot flakier, so we don't
        # @test b3.time > b1.time
    end

end


@testset "expv propagation loop over time steps" begin

    # Propagate Ψ through the first few time steps of a random dynamic
    # generator, reusing a single KrylovSubspace.
    # At each step: arnoldi! updates Ks in-place for the current (A, Ψ),
    # then expv(t, Ks) returns the new Ψ without mutating the input.
    N = 10
    rng = StableRNG(1373376837)
    tlist = collect(range(0, 10, length = 101))
    H = random_dynamic_generator(N, tlist; rng)
    Ψ₀ = random_state_vector(N; rng)
    dt = tlist[2] - tlist[1]
    t = -1im * dt
    n_steps = 5

    # Reference: sequential dense matrix exponentials
    Ψ_expected = copy(Ψ₀)
    for n = 1:n_steps
        A = evaluate(H, tlist, n)
        Ψ_expected = exp(t * Array(A)) * Ψ_expected
    end

    # Loop with arnoldi! + expv(t, Ks), reusing Ks across all steps
    A₀ = evaluate(H, tlist, 1)
    Ks = ExponentialUtilities.arnoldi(A₀, Ψ₀)
    Ψ = copy(Ψ₀)
    for n = 1:n_steps
        A = evaluate(H, tlist, n)
        ExponentialUtilities.arnoldi!(Ks, A, Ψ)
        Ψ = ExponentialUtilities.expv(t, Ks)
    end
    @test norm(Ψ - Ψ_expected) < 1e-12

end


@testset "expv! propagation loop over time steps" begin

    # Propagate Ψ through the first few time steps of a random dynamic
    # generator, reusing a single KrylovSubspace and ExpvCache.
    # At each step n: arnoldi! updates Ks in-place for the current (A, Ψ),
    # then expv! overwrites Ψ. arnoldi! reads b before expv! writes w, so
    # passing the same vector for both is safe.
    # Ks and cache must remain type-consistent across the loop: the H element
    # type is fixed at construction (real for Lanczos, complex for Arnoldi),
    # so all arnoldi! calls must use the same ishermitian setting as the
    # initial arnoldi call.
    N = 10
    rng = StableRNG(2266238742)
    tlist = collect(range(0, 10, length = 101))
    H = random_dynamic_generator(N, tlist; rng)
    Ψ₀ = random_state_vector(N; rng)
    dt = tlist[2] - tlist[1]
    t = -1im * dt
    n_steps = 5

    # Reference: sequential dense matrix exponentials
    Ψ_expected = copy(Ψ₀)
    for n = 1:n_steps
        A = evaluate(H, tlist, n)
        Ψ_expected = exp(t * Array(A)) * Ψ_expected
    end

    # Loop with arnoldi! + expv!, reusing Ks and cache across all steps.
    # The H element type is real(T) for Lanczos and T for Arnoldi, where
    # T = promote_type(eltype(A), eltype(Ψ)) is the Krylov basis vector type.
    # Hermitian A₀ → Lanczos → real H → ExpvCache{real(T)}
    A₀ = evaluate(H, tlist, 1)
    @test ishermitian(A₀)
    Ks = ExponentialUtilities.arnoldi(A₀, Ψ₀)
    T = promote_type(eltype(A₀), eltype(Ψ₀))
    cache = ExponentialUtilities.ExpvCache{real(T)}(Ks.m)
    Ψ = copy(Ψ₀)
    if RUN_BENCHMARKS
        # One might think that these in-place versions are completely
        # non-allocating. Unfortunately, that's not the case. The allocations
        # don't disappear if this is wrapped into a function.
        b1 = @btimed ExponentialUtilities.arnoldi!($Ks, $A₀, $Ψ)
        @test_broken b1.alloc == 0
        b2 = @btimed ExponentialUtilities.expv!($Ψ, $t, $Ks; cache = $cache)
        @test_broken b2.alloc == 0
    end
    Ψ = copy(Ψ₀)
    for n = 1:n_steps
        A = evaluate(H, tlist, n)
        ExponentialUtilities.arnoldi!(Ks, A, Ψ)
        ExponentialUtilities.expv!(Ψ, t, Ks; cache = cache)
    end
    @test norm(Ψ - Ψ_expected) < 1e-12

end


@testset "expv! propagation loop over time steps with non-Hermitian generator" begin

    # Non-Hermitian generators require Arnoldi iteration (not Lanczos), giving
    # an upper Hessenberg H matrix (ComplexF64). The ExpvCache element type
    # must be ComplexF64 accordingly. Otherwise the pattern is the same as for
    # Hermitian generators.
    N = 10
    rng = StableRNG(1573710086)
    tlist = collect(range(0, 10, length = 101))
    H = random_dynamic_generator(N, tlist; rng, hermitian = false)
    Ψ₀ = random_state_vector(N; rng)
    dt = tlist[2] - tlist[1]
    t = -1im * dt
    n_steps = 5

    # Reference: sequential dense matrix exponentials
    Ψ_expected = copy(Ψ₀)
    for n = 1:n_steps
        A = evaluate(H, tlist, n)
        Ψ_expected = exp(t * Array(A)) * Ψ_expected
    end

    # Loop with arnoldi! + expv!, reusing Ks and cache across all steps.
    # The H element type is T for Arnoldi, where
    # T = promote_type(eltype(A), eltype(Ψ)) is the Krylov basis vector type.
    # Non-Hermitian A₀ → Arnoldi → H element type T → ExpvCache{T}
    A₀ = evaluate(H, tlist, 1)
    @test !ishermitian(A₀)
    Ks = ExponentialUtilities.arnoldi(A₀, Ψ₀)
    T = promote_type(eltype(A₀), eltype(Ψ₀))
    cache = ExponentialUtilities.ExpvCache{T}(Ks.m)
    Ψ = copy(Ψ₀)
    for n = 1:n_steps
        A = evaluate(H, tlist, n)
        @test !ishermitian(A)
        ExponentialUtilities.arnoldi!(Ks, A, Ψ)
        ExponentialUtilities.expv!(Ψ, t, Ks; cache = cache)
    end
    @test norm(Ψ - Ψ_expected) < 1e-12

end


@testset "expv! propagation loop using error-estimate interface" begin

    # Variation of the propagation loop using expv!(w, t, A, b, Ks, cache):
    # this interface builds the Krylov subspace internally via Lanczos and
    # terminates adaptively when Saad's error estimate drops below the
    # requested tolerance (default atol=1e-8, rtol=1e-4).
    # Only Hermitian operators are supported: get_subspace_cache requires a
    # Lanczos (real H) subspace, which arnoldi builds automatically when
    # ishermitian(A) is true.
    # Passing w == b (same Ψ vector) is safe: b is fully consumed into
    # V[:,1] before w is written at the end of the call.
    N = 10
    rng = StableRNG(106614078)
    tlist = collect(range(0, 10, length = 101))
    H = random_dynamic_generator(N, tlist; rng)
    Ψ₀ = random_state_vector(N; rng)
    dt = tlist[2] - tlist[1]
    t = -1im * dt
    n_steps = 5

    # Reference: sequential dense matrix exponentials
    Ψ_expected = copy(Ψ₀)
    for n = 1:n_steps
        A = evaluate(H, tlist, n)
        Ψ_expected = exp(t * Array(A)) * Ψ_expected
    end

    A₀ = evaluate(H, tlist, 1)
    @test ishermitian(A₀)
    Ks = ExponentialUtilities.arnoldi(A₀, Ψ₀)  # auto-detects Hermitian → Lanczos
    cache = ExponentialUtilities.get_subspace_cache(Ks)
    Ψ = copy(Ψ₀)
    for n = 1:n_steps
        A = evaluate(H, tlist, n)
        @test ishermitian(A)
        ExponentialUtilities.expv!(Ψ, t, A, Ψ, Ks, cache)
    end
    # accuracy is governed by the per-step rtol=1e-4 default
    @test norm(Ψ - Ψ_expected) < 1e-3

end


@testset "expv! propagation loop using error-estimate interface with tight tolerance" begin

    # Same error-estimate interface, but with atol=rtol=1e-14. The termination
    # condition ε = atol + rtol*‖b‖ is now so tight that the algorithm uses
    # all m = min(Ks.maxiter, N) Lanczos steps, spanning the full N-dimensional
    # Krylov space and giving near-machine-precision results.
    N = 10
    rng = StableRNG(1653744503)
    tlist = collect(range(0, 10, length = 101))
    H = random_dynamic_generator(N, tlist; rng)
    Ψ₀ = random_state_vector(N; rng)
    dt = tlist[2] - tlist[1]
    t = -1im * dt
    n_steps = 5

    # Reference: sequential dense matrix exponentials
    Ψ_expected = copy(Ψ₀)
    for n = 1:n_steps
        A = evaluate(H, tlist, n)
        Ψ_expected = exp(t * Array(A)) * Ψ_expected
    end

    A₀ = evaluate(H, tlist, 1)
    Ks = ExponentialUtilities.arnoldi(A₀, Ψ₀)  # auto-detects Hermitian → Lanczos
    cache = ExponentialUtilities.get_subspace_cache(Ks)
    Ψ = copy(Ψ₀)
    for n = 1:n_steps
        A = evaluate(H, tlist, n)
        ExponentialUtilities.expv!(Ψ, t, A, Ψ, Ks, cache; atol = 1e-14, rtol = 1e-14)
    end
    @test norm(Ψ - Ψ_expected) < 1e-12

end


@testset "expv! propagation loop using error-estimate interface with non-Hermitian generator" begin

    # The error-estimate interface expv!(w, t, A, b, Ks, cache) only supports
    # Hermitian operators. For non-Hermitian generators, get_subspace_cache
    # throws an ErrorException because HermitianSubspaceCache requires a
    # Lanczos subspace with real H. Additionally, the method signature
    # requires KrylovSubspace{T,B,B} (second and third type params equal),
    # which only holds for Lanczos (Float64, Float64) not Arnoldi
    # (ComplexF64, Float64).
    N = 10
    rng = StableRNG(208995058)
    tlist = collect(range(0, 10, length = 101))
    H = random_dynamic_generator(N, tlist; rng, hermitian = false)
    Ψ₀ = random_state_vector(N; rng)
    dt = tlist[2] - tlist[1]
    t = -1im * dt
    n_steps = 5

    # Reference: sequential dense matrix exponentials
    Ψ_expected = copy(Ψ₀)
    for n = 1:n_steps
        A = evaluate(H, tlist, n)
        Ψ_expected = exp(t * Array(A)) * Ψ_expected
    end

    # get_subspace_cache throws ErrorException for non-Hermitian Ks;
    # a MethodError on expv! dispatch would also occur since the method
    # requires KrylovSubspace{T,B,B} but non-Hermitian Ks has B=ComplexF64 ≠ Float64
    A₀ = evaluate(H, tlist, 1)
    @test !ishermitian(A₀)
    Ks = ExponentialUtilities.arnoldi(A₀, Ψ₀)
    @test_broken begin
        cache = ExponentialUtilities.get_subspace_cache(Ks)
        Ψ = copy(Ψ₀)
        for n = 1:n_steps
            A = evaluate(H, tlist, n)
            ExponentialUtilities.expv!(Ψ, t, A, Ψ, Ks, cache)
        end
        norm(Ψ - Ψ_expected) < 1e-3
    end

end


@testset "expv with mode=:error_estimate" begin

    # mode=:error_estimate uses adaptive Krylov iteration with an on-the-fly
    # error estimate; default rtol ≈ 3e-4 (= sqrt(tol) with tol=1e-7)
    N = 10
    rng = StableRNG(3050745481)
    A = random_matrix(N; rng, hermitian = true)
    Ψ₀ = random_state_vector(N; rng)
    dt = 0.1
    t = -1im * dt
    Ψ = ExponentialUtilities.expv(t, A, Ψ₀; mode = :error_estimate)
    Ψ_expected = exp(t * A) * Ψ₀
    @test norm(Ψ - Ψ_expected) < 1e-4

end


@testset "expv with custom Krylov subspace size m" begin

    N = 20
    rng = StableRNG(1375059820)
    A = random_matrix(N; rng, hermitian = true)
    Ψ₀ = random_state_vector(N; rng)
    dt = 0.1
    t = -1im * dt
    Ψ = ExponentialUtilities.expv(t, A, Ψ₀; m = 8)
    Ψ_expected = exp(t * A) * Ψ₀
    @test norm(Ψ - Ψ_expected) < 1e-12

end


@testset "expv with iop (incomplete orthogonalization)" begin

    # iop=k uses only the last k vectors for re-orthogonalization,
    # reducing cost for large operators at the price of some accuracy
    N = 20
    rng = StableRNG(1918914874)
    A = random_matrix(N; rng, hermitian = false)
    Ψ₀ = random_state_vector(N; rng)
    dt = 0.1
    t = -1im * dt
    Ψ = ExponentialUtilities.expv(t, A, Ψ₀; iop = 2)
    Ψ_expected = exp(t * A) * Ψ₀
    @test norm(Ψ - Ψ_expected) < 1e-10

end


@testset "expv with custom opnorm" begin

    # opnorm can be provided to avoid recomputing it internally,
    # useful when the norm is expensive or already known
    N = 10
    rng = StableRNG(1355318934)
    A = random_matrix(N; rng, hermitian = false)
    Ψ₀ = random_state_vector(N; rng)
    dt = 0.1
    t = -1im * dt
    Ψ = ExponentialUtilities.expv(t, A, Ψ₀; opnorm = norm(A, Inf))
    Ψ_expected = exp(t * A) * Ψ₀
    @test norm(Ψ - Ψ_expected) < 1e-12

end


@testset "expv with OffsetArray state vector (broken)" begin

    # expv does not support OffsetVector as state vector b: arnoldi builds
    # V = similar(b, T, (n, m+1)) which yields an OffsetMatrix, and the
    # subsequent broadcast in expv! throws a DimensionMismatch
    N = 9
    rng = StableRNG(2841259413)
    A = random_matrix(N; rng, hermitian = true)
    Ψ₀ = OffsetVector(random_state_vector(N; rng), (-(N÷2)):(N÷2))
    dt = 0.1
    t = -1im * dt
    @test_broken begin
        Ψ = ExponentialUtilities.expv(t, A, Ψ₀)
        Ψ_expected = exp(t * A) * collect(Ψ₀)
        norm(Ψ - Ψ_expected) < 1e-12
    end

end
