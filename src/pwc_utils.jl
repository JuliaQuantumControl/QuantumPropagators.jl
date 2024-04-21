# Internal utility functions for PiecewisePropagators
#
# These allow to eliminate code duplication between the various sub-types of
# PiecewisePropagators, e.g., ChebyPropagator, NewtonPropagator, ExpPropagator.
# All of these types share the following conventions:
#
# * There is a `parameters` field that is a mapping `control` => `amplitude`,
#   where `amplitude` is the values of the control at the midpoints of the time
#   interval
# * There is a field `n` that stores the index of the next interval to
#   propagate
# * There is a field `t` that stores the time at which the current state
#   is defined
# * There is a field `generator` that contains the time-dependent generator,
#   e.g. a nested tuple `(H0, (H1, eps))`
# * There is a field `controls` that contains a list of all the controls in the
#   `generator` (usually functions or amplitude vectors)
# * There is a field `genop` that is pre-allocated for a static operator that
#   is the `generator` evaluated for a specific interval
#
# These conventions are not universally enforced. That is, a user may define a
# new PiecewisePropagator that does not follow these conventions. Thus, the
# functions defined in this file are always called explicitly, not implicitly
# via dispatch on PiecewisePropagator.

using .Controls: discretize, discretize_on_midpoints, evaluate, evaluate!


function _pwc_process_parameters(parameters, controls, tlist)
    if isnothing(parameters)
        parameters = IdDict(
            control => discretize_on_midpoints(control, tlist) for
            (i, control) in enumerate(controls)
        )
    else  # check that user-supplied parameters are for pulse parametrization
        for control ∈ controls
            amplitude = parameters[control]
            @assert amplitude isa AbstractVector
            @assert firstindex(amplitude) == 1
            @assert lastindex(amplitude) == length(amplitude)
            @assert length(amplitude) == length(tlist) - 1
        end
    end
    return parameters
end


function _pwc_set_t!(propagator::PiecewisePropagator, t)
    #
    #   | 1 | 2 | 3 |   ← intervals
    #   1   2   3   4   ← time grid points
    #
    # First, find `n` as the index on the time grid (not the intervals)
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
    # Set `propagator.n` as the index on the next interval to be propagated.
    # Note that we allow the invalid `propagator.n == N` if `t == tlist[end]`
    # (the end of the propagation).
    setfield!(propagator, :n, propagator.backward ? n - 1 : n)
    setfield!(propagator, :t, tlist[n])
end


function _pwc_get_max_genop(generator, controls, tlist)
    controlvals = [discretize(control, tlist) for control in controls]
    # TODO: this does not take into account explicit time dependencies (control
    # amplitude ≠ control function). For now, we just take any explicit time
    # dependency in the middle of the time grid.
    n = length(tlist) ÷ 2
    max_vals =
        IdDict(control => maximum(controlvals[i]) for (i, control) in enumerate(controls))
    return evaluate(generator, tlist, n; vals_dict=max_vals)
end


function _pwc_set_genop!(propagator::PiecewisePropagator, n)
    vals_dict = IdDict(c => propagator.parameters[c][n] for c in propagator.controls)
    generator = getfield(propagator, :generator)
    tlist = propagator.tlist
    evaluate!(propagator.genop, generator, tlist, n; vals_dict=vals_dict)
    return propagator.genop
end

function _pwc_get_genop(propagator::PiecewisePropagator, n)
    vals_dict = IdDict(c => propagator.parameters[c][n] for c in propagator.controls)
    generator = getfield(propagator, :generator)
    tlist = propagator.tlist
    return evaluate(generator, tlist, n; vals_dict=vals_dict)
end


function _pwc_advance_time!(propagator::PiecewisePropagator)
    n = getfield(propagator, :n)
    tlist = getfield(propagator, :tlist)
    if propagator.backward
        setfield!(propagator, :t, tlist[n])
        setfield!(propagator, :n, n - 1)
    else
        setfield!(propagator, :t, tlist[n+1])
        setfield!(propagator, :n, n + 1)
    end
end
