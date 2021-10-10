"""Implementation of propagation via direct matrix exponentiation."""

using LinearAlgebra
import StaticArrays


"""
```julia
    ExpPropWrk(v0)
```

Workspace for propagation via direct matrix exponentiation.

Initializes the workspace for the propagation of a vector `v0`
"""
struct ExpPropWrk{T}
    v :: T
    function ExpPropWrk(v0::T) where T
        new{T}(similar(v0))
    end
end


"""
```julia
expprop!(Ψ, H, dt, wrk; func=(H_dt -> exp(-1im * H_dt)))
```
Evaluate `Ψ = func(H*dt) Ψ` by directly evaluating `U = func(H*dt)`, i.e. by
matrix exponentiation for the default `func`, and then multiplying `U` and
`Ψ` in-place with `mul!`.

The workspace `wrk` must be initialized with [`ExpPropWrk`](@ref) to provide
storage for a temporary state.
"""
function expprop!(Ψ, H, dt, wrk; kwargs...)
    func = get(kwargs, :func, H_dt -> exp(-1im * H_dt))
    copyto!(wrk.v, Ψ)
    U = func(H*dt)
    mul!(Ψ, U, wrk.v)
end


function expprop(Ψ::StaticArrays.SVector, H, dt, wrk; kwargs...)
    func = get(kwargs, :func, H_dt -> exp(-1im * H_dt))
    return func(H*dt) * Ψ
end


function expprop(Ψ, H, dt, wrk; kwargs...)
    func = get(kwargs, :func, H_dt -> exp(-1im * H_dt))
    return func(H*dt) * Ψ
end
