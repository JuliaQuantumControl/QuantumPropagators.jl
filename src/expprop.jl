# Implementation of propagation via direct matrix exponentiation
module ExpProp

export ExpPropWrk, expprop!

using LinearAlgebra
using TimerOutputs: @timeit_debug, TimerOutput


"""
```julia
ExpPropWrk(v0)
```

Workspace for propagation via direct matrix exponentiation.

Initializes the workspace for the propagation of a vector `v0`
"""
struct ExpPropWrk{T}
    v::T
    timing_data::TimerOutput
    function ExpPropWrk(v0::T; _timing_data=TimerOutput()) where {T}
        new{T}(similar(v0), _timing_data)
    end
end


"""
```julia
expprop!(Ψ, H, dt, wrk; func=(H_dt -> exp(-1im * H_dt)), _...)
```
Evaluate `Ψ = func(H*dt) Ψ` by directly evaluating `U = func(H*dt)`, i.e. by
matrix exponentiation for the default `func`, and then multiplying `U` and
`Ψ` in-place with `mul!`.

The workspace `wrk` must be initialized with [`ExpPropWrk`](@ref) to provide
storage for a temporary state.

Keyword arguments besides `func` are ignored.
"""
function expprop!(Ψ, H, dt, wrk; func=(H_dt -> exp(-1im * H_dt)), _...)
    copyto!(wrk.v, Ψ)
    @timeit_debug wrk.timing_data "matrix exponentiation" begin
        U = func(H * dt)
    end
    @timeit_debug wrk.timing_data "matrix-vector product" begin
        mul!(Ψ, U, wrk.v)
    end
end


"""
```julia
Ψ_out = expprop(Ψ, H, dt, wrk; func=(H_dt -> exp(-1im * H_dt)), _...)
```

evaluates `Ψ_out = func(H*dt) Ψ` as in [`expprop!`](@ref), but not acting
in-place.
"""
function expprop(Ψ, H, dt, wrk; func=(H_dt -> exp(-1im * H_dt)), _...)
    @timeit_debug wrk.timing_data "matrix exponentiation" begin
        U = func(H * dt)
    end
    @timeit_debug wrk.timing_data "matrix-vector product" begin
        return U * Ψ
    end
end

end
