using TimerOutputs: @timeit_debug, TimerOutput

@doc raw"""
Wrap around a [`Generator`](@ref), for use as an ODE function.

```julia
f = ode_function(generator, tlist; c=-1im)
```

creates a function suitable to be passed to
[`ODEProblem`](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEProblem).

```math
\gdef\op#1{\hat{#1}}
\gdef\ket#1{\vert{#1}\rangle}
```

With `generator` corresponding to ``\op{H}(t)``, this implicitly encodes the
ODE

```math
\frac{\partial}{\partial t} \ket{\Psi(t)} = c \op{H}(t) \ket{\Psi(t)}
```

for the state ``\ket{\Psi(t)}``. With the default ``c = -i``, this corresponds
to the Schrödinger equation, or the Liouville equation with
[`convention=:LvN`](@ref liouvillian).

The resulting `f` works both in-place and not-in-place, as

```julia
f(ϕ, Ψ, vals_dict, t)   # in-place `f(du, u, p, t)`
ϕ = f(Ψ, vals_dict, t)  # not-in-place `f(u, p, t)`
```

Calling `f` as above is functionally equivalent to calling [`evaluate`](@ref)
to obtain an operator `H` from the original time-dependent `generator`, and
then applying `H` to the current quantum state `Ψ`:

```julia
H = evaluate(f.generator, t; vals_dict=vals_dict)
ϕ = c * H * Ψ
```

where `vals_dict` may be a dictionary mapping controls to values (set as the
parameters `p` of the underlying ODE solver).

If [`QuantumPropagators.enable_timings()`](@ref
QuantumPropagators.enable_timings) has been called,
profiling data is collected in `f.timing_data`.
"""
function ode_function(generator::GT, tlist; c=-1im, _timing_data=TimerOutput()) where {GT}
    H = evaluate(generator, tlist, 1)
    OT = typeof(H)
    return QuantumODEFunction{GT,OT}(generator, H, c, _timing_data)
end


struct QuantumODEFunction{GT,OT}
    generator::GT
    operator::OT
    c::ComplexF64
    timing_data::TimerOutput
end

function (f::QuantumODEFunction)(du, u, p, t)
    @timeit_debug f.timing_data "operator evaluation" begin
        evaluate!(f.operator, f.generator, t; vals_dict=p)
    end
    @timeit_debug f.timing_data "matrix-vector product" begin
        return mul!(du, f.operator, u, f.c, false)
    end
end


function (f::QuantumODEFunction)(u, p, t)
    @timeit_debug f.timing_data "operator evaluation" begin
        H = evaluate(f.generator, t; vals_dict=p)
    end
    @timeit_debug f.timing_data "matrix-vector product" begin
        return f.c * H * u
    end
end
