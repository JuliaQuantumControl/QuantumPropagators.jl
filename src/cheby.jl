"""Implementation of the Chebychev propagator."""
module ChebychevPropagator

export cheby_coeffs, cheby_coeffs!, ChebyWrk, cheby!

using SpecialFunctions
using LinearAlgebra


"""Calculate Chebychev coefficients.

Return an array of coefficiencts larger than `limit`.

Args:

* `Δ`: the spectral radius of the underlying operator
* `dt`: the time step
"""
function cheby_coeffs(Δ, dt; limit=1e-12)
    α = abs(0.5 * Δ * dt)
    coeffs = Float64[]
    a = besselj(0, α)
    append!(coeffs, a)
    ϵ = abs(a)
    i = 1
    while ϵ > limit
        a = 2 * besselj(i,α)
        append!(coeffs, a)
        ϵ = abs(a)
        i += 1
    end
    return coeffs
end


"""Calculate Chebychev coefficients in-place.

Writes the coefficients into `coeffs`. Returns the number `n` of coefficients
required for convergence, or the size of `coeffs`, whichever is less.

* `coeffs`: An array to hold the Chebychev cofficients. On output, the elements
   1 through `n` will hold the calculated coefficients, while the remaining
   elements will be unchanged
* `dt`: the time step
"""
function cheby_coeffs!(coeffs, Δ, dt; limit=1e-12)
    α = abs(0.5 * Δ * dt)
    coeffs[1] = besselj0(α)
    n = length(coeffs)
    for i = 2:n
        coeffs[i] = 2 * besselj(i-1,α)
        if abs(coeffs[i]) < limit
            n = i
            break
        end
    end
    return n
end


"""
Workspace for the Chebychev propagation routine.

```julia
    ChebyWrk(Ψ, Δ, E_min, dt; limit=1e-12)
```

initializes the workspace for the propagation of a state similar to Ψ under a
Hamiltonian with eigenvalues between `E_min` and `E_min + Δ`, and a time step
dt. Chebychev coefficients smaller than the given `limit` are discarded.
"""
struct ChebyWrk{T}
    v0 :: T
    v1 :: T
    v2 :: T
    coeffs :: Vector{Float64}
    Δ :: Float64
    E_min :: Float64
    dt :: Float64
    limit :: Float64
    function ChebyWrk(Ψ::T, Δ::Float64, E_min::Float64, dt::Float64;
                      limit::Float64=1e-12) where T
        v0::T = similar(Ψ)
        v1::T = similar(Ψ)
        v2::T = similar(Ψ)
        coeffs = cheby_coeffs(Δ, dt; limit=limit)
        new{T}(v0, v1, v2, coeffs, Δ, E_min, dt, limit)
    end
end


"""Evaluate `Ψ = exp(-i H dt) Ψ` in-place.

Args:

* `Ψ`: on input, initial vector. Will be overwritten with result.
* `H`: Hermitian operator
* `dt`: time step
* `E_min`: minimum eigenvalue of H, to be used instead of the `E_min` from the
   initialization of `wrk`. The same `wrk` may be used for different values
   `E_min`, as long as the spectra radius `Δ` and the time step `dt` are the
   same as those used for the initialization of `wrk`.

The routine will not allocate any internal storage. This implementation
requires `copyto!` `lmul!`, and `axpy!` to be implemented for `Ψ`, and the
three-argument `mul!` for `Ψ` and `H`.
"""
function cheby!(Ψ, H, dt, wrk; E_min::Union{Float64, Nothing}=nothing,
                check_normalization=false)

    Δ = wrk.Δ
    β :: Float64 = (Δ / 2) + wrk.E_min  # "normfactor"
    if E_min ≠ nothing
        β = (Δ / 2) + E_min
    end
    @assert dt ≈ wrk.dt "wrk was initialized for dt=$(wrk.dt), not dt=$dt"
    if dt > 0
        c = -2im / Δ
    else
        c = 2im / Δ
    end
    a = wrk.coeffs
    @assert length(a) > 1 "Need at least 2 Chebychev coefficients"
    v0 = wrk.v0
    v1 = wrk.v1
    v2 = wrk.v2

    # v0 ↔ Ψ; Ψ = a[1] * v0
    copyto!(v0, Ψ)
    lmul!(a[1], Ψ)

    # v1 = -i * H_norm * v0 = c * (H * v0 - β * v0)
    mul!(v1, H, v0)
    axpy!(-β, v0, v1)
    lmul!(c, v1)

    # Ψ += a[2] * v1
    axpy!(a[2], v1, Ψ)

    c *= 2

    for i = 3 : length(a)

        # v2 = -2i * H_norm * v1 + v0 = c * (H * v1 - β * v1) + v0
        mul!(v2, H, v1)
        axpy!(-β, v1, v2)
        lmul!(c, v2)
        if check_normalization
            map_norm = abs(dot(v1, v2)) / (2 * norm(v1)^2)
            @assert(
                map_norm <= (1.0 + wrk.limit),
                "Incorrect normalization (E_min, wrk.Δ)"
            )
        end
        v2 .+= v0

        # Ψ += a[i] * v2
        axpy!(a[i], v2, Ψ)

        v0, v1, v2 = v1, v2, v0  # switch w/o copying

    end

    lmul!(exp(-im * β * dt), Ψ)

end


end
