using .Controls: get_controls
using TimerOutputs: reset_timer!, @timeit_debug

"""Propagator for Newton propagation (`method=QuantumPropagators.Newton`).

This is a [`PWCPropagator`](@ref).
"""
mutable struct NewtonPropagator{GT,OT,ST} <: PWCPropagator
    generator::GT
    state::ST
    t::Float64  # time at which current `state` is defined
    n::Int64 # index of next interval to propagate
    tlist::Vector{Float64}
    parameters::AbstractDict
    controls
    genop::OT
    wrk::Newton.NewtonWrk{ST}
    backward::Bool
    inplace::Bool
    func::Function
    norm_min::Float64
    relerr::Float64
    max_restarts::Int64
end

set_t!(propagator::NewtonPropagator, t) = _pwc_set_t!(propagator, t)


"""
```julia
using QuantumPropagators: Newton

newton_propagator = init_prop(
    state,
    generator,
    tlist;
    method=Newton,
    inplace=true,
    backward=false,
    verbose=false,
    parameters=nothing,
    m_max=10,
    func=(z -> exp(-1im * z)),
    norm_min=1e-14,
    relerr=1e-12,
    max_restarts=50,
    _...
)
```

initializes a [`NewtonPropagator`](@ref).

# Method-specific keyword arguments

* `m_max`: maximum Krylov dimension, cf. [`NewtonWrk`](@ref
  QuantumPropagators.Newton.NewtonWrk)
* `func`, `norm_min`, `relerr`, `max_restarts`: parameter to pass to
  [`newton!`](@ref QuantumPropagators.Newton.newton!)
"""
function init_prop(
    state,
    generator,
    tlist,
    method::Val{:Newton};
    inplace=true,
    backward=false,
    verbose=false,
    parameters=nothing,
    m_max=10,
    func=(z -> exp(-1im * z)),
    norm_min=1e-14,
    relerr=1e-12,
    max_restarts=50,
    _...
)
    tlist = convert(Vector{Float64}, tlist)
    controls = get_controls(generator)
    G = _pwc_get_max_genop(generator, controls, tlist)

    parameters = _pwc_process_parameters(parameters, controls, tlist)
    wrk = Newton.NewtonWrk(state; m_max=m_max)
    n = 1
    t = tlist[1]
    if backward
        n = length(tlist) - 1
        t = float(tlist[n+1])
    end
    GT = typeof(generator)
    OT = typeof(G)
    ST = typeof(state)
    inplace || error("The Newton propagator is only implemented in-place")
    return NewtonPropagator{GT,OT,ST}(
        generator,
        inplace ? copy(state) : state,
        t,
        n,
        tlist,
        parameters,
        controls,
        G,
        wrk,
        backward,
        inplace,
        func,
        norm_min,
        relerr,
        max_restarts
    )
end

# Aliases (for backwards compatibility)
init_prop(state, generator, tlist, method::Val{:newton}; kwargs...) =
    init_prop(state, generator, tlist, Val(:Newton); kwargs...)


function reinit_prop!(propagator::NewtonPropagator, state; _...)
    set_state!(propagator, state)
    if propagator.backward
        set_t!(propagator, propagator.tlist[end])
    else
        set_t!(propagator, propagator.tlist[begin])
    end
    reset_timer!(propagator.wrk.timing_data)
end


function prop_step!(propagator::NewtonPropagator)
    @timeit_debug propagator.wrk.timing_data "prop_step!" begin
        Ψ = propagator.state
        H = propagator.genop
        n = propagator.n  # index of interval we're going to propagate
        tlist = getfield(propagator, :tlist)
        (0 < n < length(tlist)) || return nothing
        dt = tlist[n+1] - tlist[n]
        if propagator.backward
            dt = -dt
        end
        _pwc_set_genop!(propagator, n)
        if propagator.inplace
            Newton.newton!(
                Ψ,
                H,
                dt,
                propagator.wrk;
                func=propagator.func,
                norm_min=propagator.norm_min,
                relerr=propagator.relerr,
                max_restarts=propagator.max_restarts
            )
        else
            error("The Newton propagator is only implemented in-place")
        end
        _pwc_advance_time!(propagator)
        return propagator.state
    end
end
