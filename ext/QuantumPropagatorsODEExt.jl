module QuantumPropagatorsODEExt

using LinearAlgebra
using TimerOutputs: @timeit_debug, reset_timer!, TimerOutput
using OrdinaryDiffEq: OrdinaryDiffEq as ODE
using OrdinaryDiffEq.SciMLBase: FullSpecialize
using QuantumPropagators:
    QuantumPropagators,
    AbstractPropagator,
    PWCPropagator,
    _pwc_process_parameters,
    ode_function
using QuantumPropagators.Controls:
    get_controls, get_parameters, evaluate, evaluate!, discretize_on_midpoints
import QuantumPropagators: init_prop, reinit_prop!, prop_step!, set_state!, set_t!


"""
```julia
using OrdinaryDiffEq  # or: `using DifferentialEquations`

ode_propagator = init_prop(
    state,
    generator,
    tlist;
    method=OrdinaryDiffEq,  # or: `method=DifferentialEquations`
    inplace=true,
    backward=false,
    verbose=false,
    parameters=nothing,
    piecewise=false,
    pwc=false,
    alg=OrdinaryDiffEq.Tsit5(),
    solver_options...
)
```

initializes a [propagator](@ref AbstractPropagator) that uses an ODE solver
from the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl)
package as a backend.

By default, the resulting propagator is for time-continuous controls that can
be evaluated with [`evaluate(control, t)`](@ref evaluate) for any `t` in the
range of `tlist[begin]` to `tlist[end]`. The controls may be parametrized, see
[`get_parameters`](@ref). Any parameters will be available in the `parameters`
attribute of the resulting `ode_propagator`, as a dictionary mapping controls
to a vector of parameter values. Mutating `ode_propagator.parameters[control]`
will be reflected in any subsequent call to [`prop_step!`](@ref).

If `pwc=true` (or, equivalently `piecewise=true`), all controls will be
discretized with [`discretize_on_midpoints`](@ref) and the propagation will be
for piecewise constant dynamics. The resulting `ode_propagator` will be an
instance of [`PWCPropagator`](@ref), with the corresponding semantics. In
particular, the `ode_propatator.parameters` will be a mapping of controls to
discretized pulse values, *not* the analytical parameters obtained with
`get_parameters(control)` as in the default case.

Internally, the `generator` will be wrapped with
[`QuantumPropagators.ode_function`](@ref). The resulting function `f` will be
called internally as `f(ϕ, Ψ, vals_dict, t)` or `f(Ψ, vals_dict, t)` depending
on the `inplace` keyword argument.


# Method-specific keyword arguments

* `pwc`: Whether to propagate for piecewise-constant controls or, with the
  default `pwc=false`, for time-continuous controls.
* `piecewise`: Currently equivalent to `pwc`, but future version may change
   this to allow for other piecewise (e.g., piecewise-linear) controls.
* `parameters`: If given, a mapping of controls to parameter values
  (`pwc=false`) or pulse values on the intervals of the time grid (`pwc=true`).
  By default, the `parameters` are determined automatically using
  [`get_parameters`](@ref), respectively [`discretize_on_midpoints`](@ref) if
  `pwc=true`. If they are given manually, they must follow the exact same
  semantics. In particular, for `pwc=false`, any parameters must alias the
  parameters in the controls, such that mutating `parameters` is automatically
  reflected in [`evaluate`](@ref). The `parameters` will be available as an
  attribute of the `ode_propagator`.
* `alg`: The algorithm to use for the ODE Solver, see the list of solvers in
  the [DifferentialEquations manual](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
  The default `Tsit5()` method is the recommended choice for non-stiff
  problems.
* `solver_options`: All other keyword arguments are passed to the ODE solver,
  see the list of [Solve Keyword Arguments](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#solver_options)
  in the DifferentialEquations manual. Note that the options for "Default
  Algorithm Hinting" do not apply, since `alg` must be specified manually.
  Also, the "Output Control" is managed by the `ode_propagator`, so these
  options should not be used.
"""
function init_prop(
    state,
    generator,
    tlist,
    method::Val{:OrdinaryDiffEq};
    inplace=true,
    backward=false,
    verbose=false,
    parameters=nothing,
    piecewise=false,
    pwc=false,
    alg=ODE.Tsit5(),
    _ode_function_extra=Dict{Symbol,Any}(),  # undocumented
    solver_options...
)
    pwc = (pwc ≡ true) || (piecewise ≡ true)
    # In theory, there could be options (`piecewise_kind`?) for piecewise but
    # not pwc (e.g., piecewise-linear). We don't have anything like that
    # implemented at the moment, to piecewise and pwc are equivalent.
    timing_data = TimerOutput()
    f = ode_function(generator, tlist; _timing_data=timing_data)
    u0 = state
    tspan = (tlist[begin], tlist[end])
    if backward
        tspan = (tlist[end], tlist[begin])
    end
    controls = get_controls(generator)
    n = backward ? (length(tlist) - 1) : 1
    if pwc
        parameters = _pwc_process_parameters(parameters, controls, tlist)
        vals_dict = IdDict(control => parameters[control][n] for control in controls)
        # vals_dict will be available in integrator.p
    else  # time-continuous
        if isnothing(parameters)
            parameters = IdDict(control => get_parameters(control) for control in controls)
        else
            for control in controls
                @assert haskey(parameters, control)
            end
        end
        vals_dict = IdDict()
    end
    f_ode = ODE.ODEFunction{inplace,FullSpecialize}(f; _ode_function_extra...)
    prob = ODE.ODEProblem(f_ode, u0, tspan, vals_dict)
    integrator = ODE.init(prob, alg; save_everystep=false, solver_options...)
    IT = typeof(integrator)
    if pwc
        return ODEPWCPropagator{IT}(
            integrator,
            tlist,
            parameters,
            backward,
            inplace,
            timing_data,
            n
        )
    else
        return ODEContinuousPropagator{IT}(
            integrator,
            tlist,
            parameters,
            backward,
            inplace,
            timing_data,
            n
        )
    end
end

# Aliases

# Any `using DifferentialEquations` also loads `OrdinaryDiffEq`
init_prop(state, generator, tlist, method::Val{:DifferentialEquations}; kwargs...) =
    init_prop(state, generator, tlist, Val(:OrdinaryDiffEq); kwargs...)

###############################################################################


mutable struct ODEContinuousPropagator{IT} <: AbstractPropagator
    const integrator::IT
    const tlist::Vector{Float64}
    const parameters::AbstractDict
    const backward::Bool
    const inplace::Bool
    const timing_data::TimerOutput
    n::Int64 # index of next interval to propagate
end


mutable struct ODEPWCPropagator{IT} <: PWCPropagator
    # This is identical to `ODEContinuousPropagator`, but needs to be written
    # out because it's in a different place in the type hierarchy
    const integrator::IT
    const tlist::Vector{Float64}
    const parameters::AbstractDict
    const backward::Bool
    const inplace::Bool
    const timing_data::TimerOutput
    n::Int64 # index of next interval to propagate
end


const ODEPropagator = Union{ODEContinuousPropagator,ODEPWCPropagator}


function Base.getproperty(propagator::ODEPropagator, name::Symbol)
    if name ≡ :state
        integrator = getfield(propagator, :integrator)
        return integrator.u
    elseif name ≡ :t
        integrator = getfield(propagator, :integrator)
        return integrator.t
    end
    return getfield(propagator, name)
end


function reinit_prop!(propagator::ODEPropagator, state)
    ODE.reinit!(propagator.integrator, state)
    tlist = getfield(propagator, :tlist)
    n = propagator.backward ? (length(tlist) - 1) : 1
    setfield!(propagator, :n, n)
    reset_timer!(propagator.timing_data)
end


function prop_step!(propagator::ODEPropagator)
    @timeit_debug propagator.timing_data "prop_step!" begin
        tlist = getfield(propagator, :tlist)
        n = propagator.n  # index of interval we're going to propagate
        (0 < n < length(tlist)) || return nothing
        dt = tlist[n+1] - tlist[n]
        if propagator.backward
            dt = -dt
        end
        integrator = propagator.integrator
        if propagator isa ODEPWCPropagator
            controls = keys(integrator.p)
            for control in controls
                # `integrator.p` is the `vals_dict` in the call to
                # `QuantumODEFunction`
                integrator.p[control] = propagator.parameters[control][n]
            end
        end
        ODE.step!(integrator, dt, true) # stop_at_tdt
        n += propagator.backward ? -1 : +1
        setfield!(propagator, :n, n)
        return integrator.u
    end
end


function set_state!(propagator::ODEPropagator, state)
    if state ≢ propagator.state
        if propagator.inplace
            # ODE.set_u! does not work in-place
            copyto!(propagator.integrator.u, state)
            ODE.u_modified!(propagator.integrator, true)
        else
            ODE.set_u!(propagator.integrator, state)
        end
    end
end


function set_t!(propagator::ODEPropagator, t)
    tlist = getfield(propagator, :tlist)
    if t <= tlist[1]
        n = 1
    else
        N = length(tlist)
        if t >= tlist[end]
            n = N
        else
            n = min(searchsortedfirst(tlist, t), N)
        end
    end
    (t ≈ tlist[n]) || (@warn ("Snapping t=$t to time grid value $(tlist[n])"))
    setfield!(propagator, :n, propagator.backward ? n - 1 : n)
    ODE.set_t!(propagator.integrator, t)
end


end
