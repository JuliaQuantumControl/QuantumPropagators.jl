module Cheby

export cheby_coeffs, cheby_coeffs!, ChebyWrk, cheby!, cheby

using SpecialFunctions
using LinearAlgebra
using TimerOutputs: @timeit_debug, TimerOutput


"""Calculate Chebychev coefficients.

```julia
a::Vector{Float64} = cheby_coeffs(Δ, dt; limit=1e-12)
```

return an array of coefficients larger than `limit`.

# Arguments

* `Δ`: the spectral radius of the underlying operator
* `dt`: the time step

See also [`cheby_coeffs!`](@ref) for an in-place version.
"""
function cheby_coeffs(Δ, dt; limit=1e-12)
    α = abs(0.5 * Δ * dt)
    coeffs = Float64[]
    a = besselj(0, α)
    append!(coeffs, a)
    ϵ = abs(a)
    i = 1
    while ϵ > limit
        a = 2 * besselj(i, α)
        append!(coeffs, a)
        ϵ = abs(a)
        i += 1
    end
    return coeffs
end


"""Calculate Chebychev coefficients in-place.

```julia
n::Int = cheby_coeffs!(coeffs, Δ, dt, limit=1e-12)
```

overwrites the first `n` values in `coeffs` with new coefficients larger than
`limit` for the given new spectral radius `Δ` and time step `dt`. The `coeffs`
array will be resized if necessary, and may have a length > `n` on exit.

See also [`cheby_coeffs`](@ref) for an non-in-place version.
"""
function cheby_coeffs!(coeffs, Δ, dt, limit=1e-12)
    α = abs(0.5 * Δ * dt)
    a = besselj(0, α)
    coeffs[1] = a
    N = length(coeffs)
    ϵ = abs(a)
    n = 1
    while ϵ > limit
        n += 1
        a = 2 * besselj(n - 1, α)
        coeffs[n] = a
        ϵ = abs(a)
        if (n >= N)
            N *= 2
            resize!(coeffs, N)
        end
    end
    return n
end


"""
Workspace for the Chebychev propagation routine.

```julia
ChebyWrk(Ψ, Δ, E_min, dt; limit=1e-12)
```

initializes the workspace for the propagation of a state similar to `Ψ` under a
Hamiltonian with eigenvalues between `E_min` and `E_min + Δ`,
and a time step `dt`. Chebychev coefficients smaller than the given `limit` are
discarded.
"""
mutable struct ChebyWrk{ST,CFS,FT<:AbstractFloat}
    v0::ST
    v1::ST
    v2::ST
    coeffs::CFS
    n_coeffs::Int64
    Δ::FT
    E_min::FT
    dt::FT
    limit::FT
    timing_data::TimerOutput
    function ChebyWrk(
        Ψ::ST,
        Δ::FT,
        E_min::FT,
        dt::FT;
        limit::FT=1e-12,
        _timing_data=TimerOutput()
    ) where {ST,FT}
        v0::ST = similar(Ψ)
        v1::ST = similar(Ψ)
        v2::ST = similar(Ψ)
        coeffs = cheby_coeffs(Δ, dt; limit=limit)
        n_coeffs = length(coeffs)
        new{ST,typeof(coeffs),FT}(
            v0,
            v1,
            v2,
            coeffs,
            n_coeffs,
            Δ,
            E_min,
            dt,
            limit,
            _timing_data
        )
    end
end


"""Evaluate `Ψ = exp(-𝕚 * H * dt) Ψ` in-place.

```julia
cheby!(Ψ, H, dt, wrk; E_min=nothing, check_normalization=false)
```

# Arguments

* `Ψ`: on input, initial vector. Will be overwritten with result.
* `H`: Hermitian operator
* `dt`: time step
* `wrk`: internal workspace
* `E_min`: minimum eigenvalue of H, to be used instead of the `E_min` from the
   initialization of `wrk`. The same `wrk` may be used for different values
   `E_min`, as long as the spectra radius `Δ` and the time step `dt` are the
   same as those used for the initialization of `wrk`.
* `check_normalizataion`: perform checks that the H does not exceed the
  spectral radius for which the workspace was initialized.

The routine will not allocate any internal storage. This implementation
requires `copyto!` `lmul!`, and `axpy!` to be implemented for `Ψ`, and the
three-argument `mul!` for `Ψ` and `H`.
"""
function cheby!(Ψ, H, dt, wrk; kwargs...)

    E_min = get(kwargs, :E_min, wrk.E_min)
    check_normalization = get(kwargs, :check_normalization, false)

    Δ = wrk.Δ
    β::Float64 = (Δ / 2) + E_min  # "normfactor"
    @assert abs(dt) ≈ abs(wrk.dt) "wrk was initialized for dt=$(wrk.dt), not dt=abs($dt)"
    if dt > 0
        c = -2im / Δ
    else
        c = 2im / Δ
    end
    a = wrk.coeffs
    ϵ = wrk.limit
    @assert length(a) > 1 "Need at least 2 Chebychev coefficients"
    v0 = wrk.v0
    v1 = wrk.v1
    v2 = wrk.v2

    # v0 ↔ Ψ; Ψ = a[1] * v0
    copyto!(v0, Ψ)
    lmul!(a[1], Ψ)

    # v1 = -i * H_norm * v0 = c * (H * v0 - β * v0)
    @timeit_debug wrk.timing_data "matrix-vector product" begin
        mul!(v1, H, v0)
    end
    axpy!(-β, v0, v1)
    lmul!(c, v1)

    # Ψ += a[2] * v1
    axpy!(a[2], v1, Ψ)

    c *= 2

    for i = 3:wrk.n_coeffs

        # v2 = -2i * H_norm * v1 + v0 = c * (H * v1 - β * v1) + v0
        @timeit_debug wrk.timing_data "matrix-vector product" begin
            mul!(v2, H, v1)
        end
        axpy!(-β, v1, v2)
        lmul!(c, v2)
        if check_normalization
            map_norm = abs(dot(v1, v2)) / (2 * norm(v1)^2)
            @assert(
                map_norm <= (1.0 + ϵ),
                "Incorrect normalization (E_min=$(E_min), Δ=$(Δ))"
            )
        end
        # v2 += v0
        axpy!(true, v0, v2)

        # Ψ += a[i] * v2
        axpy!(a[i], v2, Ψ)

        v0, v1, v2 = v1, v2, v0  # switch w/o copying

    end

    lmul!(exp(-im * β * dt), Ψ)

end


"""Evaluate `Ψ = exp(-𝕚 * H * dt) Ψ`.

```julia
Ψ_out = cheby(Ψ, H, dt, wrk; E_min=nothing, check_normalization=false)
```

acts like [`cheby!`](@ref) but does not modify `Ψ` in-place.
"""
function cheby(Ψ, H, dt, wrk; kwargs...)

    E_min = get(kwargs, :E_min, wrk.E_min)
    check_normalization = get(kwargs, :check_normalization, false)

    Δ = wrk.Δ
    β::Float64 = (Δ / 2) + E_min  # "normfactor"
    @assert abs(dt) ≈ wrk.dt "wrk was initialized for dt=$(wrk.dt), not dt=abs($dt)"
    if dt > 0
        c = -2im / Δ
    else
        c = 2im / Δ
    end
    a = wrk.coeffs
    ϵ = wrk.limit
    @assert length(a) > 1 "Need at least 2 Chebychev coefficients"

    v0 = Ψ
    Ψ = a[1] * v0

    @timeit_debug wrk.timing_data "matrix-vector product" begin
        v1 = c * (H * v0 - β * v0)
    end
    Ψ += a[2] * v1

    c *= 2

    for i = 3:wrk.n_coeffs

        @timeit_debug wrk.timing_data "matrix-vector product" begin
            v2 = H * v1
        end
        if check_normalization
            v2 = c * (v2 - v1 * β)
            map_norm = abs(dot(v1, v2)) / (2 * norm(v1)^2)
            @assert(
                map_norm <= (1.0 + ϵ),
                "Incorrect normalization (E_min=$(E_min), Δ=$(Δ))"
            )
            v2 += v0
        else
            v2 = c * (v2 - β * v1) + v0
        end

        Ψ += a[i] * v2

        v0, v1, v2 = v1, v2, v0  # switch w/o copying

    end

    return exp(-im * β * dt) * Ψ

end

end
