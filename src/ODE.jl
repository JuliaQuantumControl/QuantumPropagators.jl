"""Implementation of the ODE propagator."""

using LinearAlgebra
using OrdinaryDiffEq

"""
Workspace for the ODE propagation routine.

```julia
    ODEWrk(Ψ, dt; alg=DP5())
```

initializes the workspace for the propagation of a state similar to Ψ, under a
Hamiltonian similar to H, and for a time step similar to dt. 
"""
mutable struct ODEWrk{T, FT<:AbstractFloat}
    integrator
    function ODEWrk(Ψ::T, dt::FT, H; alg=DP5()) where {T, FT}
        tspan = (0.0,dt)
        H = zeros(eltype(T), length(Ψ), length(Ψ))
        function f!(dΨ::T, Ψ::T, H, t::FT)
            mul!(dΨ, H, Ψ)
            lmul!(-1im, dΨ)
        end
        prob = ODEProblem{true}(f!, Ψ, tspan, H)
        integrator = init(prob, alg; dt=dt, save_everystep=false)
        new{T, FT}(integrator)
    end
end


"""Evaluate `Ψ = exp(-i H dt) Ψ` in-place.

```julia
    ODE!(Ψ, H, dt, wrk)
```

Args:

* `Ψ`: on input, initial vector. Will be overwritten with result.
* `H`: Hermitian operator
* `dt`: time step
* `wrk`: internal workspace

The routine will not allocate any internal storage. This
implementation requires `copyto!` to be implemented for `Ψ`, and the
three-argument `mul!` for `Ψ` and `H`.  
"""
function ODE!(Ψ, H, dt, wrk; kwargs...)
    wrk.integrator.dt = dt
    copyto!(wrk.integrator.p, H)
    
    # this is probably not a good practice (but otherwise it will interpolate the controls)
    reinit!(wrk.integrator, Ψ) 
    
    step!(wrk.integrator, dt)
    copyto!(Ψ, wrk.integrator.u)
end
