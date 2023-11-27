using .Controls: get_controls, evaluate, discretize
using TimerOutputs: reset_timer!, @timeit_debug

"""Propagator for Chebychev propagation (`method=QuantumPropagators.Cheby`).

This is a [`PWCPropagator`](@ref).
"""
mutable struct ChebyPropagator{GT,OT,ST} <: PWCPropagator
    generator::GT
    state::ST
    t::Float64  # time at which current `state` is defined
    n::Int64 # index of next interval to propagate
    tlist::Vector{Float64}
    parameters::AbstractDict
    controls
    control_ranges::AbstractDict
    genop::OT
    wrk::Cheby.ChebyWrk{ST,Vector{Float64},Float64}
    backward::Bool
    inplace::Bool
    specrange_method::Symbol
    specrange_buffer::Float64
    check_normalization::Bool
    specrange_options::Dict{Symbol,Any}
end


set_t!(propagator::ChebyPropagator, t) = _pwc_set_t!(propagator, t)


"""
```julia
using QuantumPropagators: Cheby

cheby_propagator = init_prop(
    state,
    generator,
    tlist;
    method=Cheby,
    inplace=true,
    backward=false,
    verbose=false,
    parameters=nothing,
    control_ranges=nothing,
    specrange_method=:auto,
    specrange_buffer=0.01,
    cheby_coeffs_limit=1e-12,
    check_normalization=false,
    specrange_kwargs...
)
```

initializes a [`ChebyPropagator`](@ref).

# Method-specific keyword arguments

* `control_ranges`: a dict the maps the controls in `generator` (see
  [`get_controls`](@ref QuantumPropagators.Controls.get_controls)) to a tuple
  of min/max values. The Chebychev coefficients will be calculated based on a
  spectral envelope that assumes that each control can take arbitrary values
  within the min/max range. If not given, the ranges are determined
  automatically. Specifying manual control ranges can be useful when the the
  control amplitudes (`parameters`) may change during the propagation, e.g. in
  a sequential-update control scheme.
* `specrange_method`: Method to pass to the [`specrange`](@ref
  QuantumPropagators.SpectralRange.specrange) function
* `specrange_buffer`: An additional factor by which to enlarge the estimated
   spectral range returned by [`specrange`](@ref
  QuantumPropagators.SpectralRange.specrange), in order to ensure that
  Chebychev coefficients are based on an overestimation of the spectral range.
* `cheby_coeffs_limit`: The maximum magnitude of Chebychev coefficients that
  should be treated as non-zero
* `check_normalization`: Check whether the Hamiltonian has been properly
  normalized, i.e., that the spectral range of `generator` has not been
  underestimated. This slowes down the propagation, but is advisable for novel
  `generators`.
* `uniform_dt_tolerance=1e-12`: How much the intervals of `tlist` are allowed
  to vary while still being considered constant.
* `specrange_kwargs`: All further keyword arguments are passed to the
  [`specrange`](@ref QuantumPropagators.SpectralRange.specrange) function
"""
function init_prop(
    state,
    generator,
    tlist,
    method::Val{:Cheby};
    inplace=true,
    backward=false,
    verbose=false,
    parameters=nothing,
    control_ranges=nothing,
    specrange_method=:auto,
    specrange_buffer=0.01,
    cheby_coeffs_limit=1e-12,
    check_normalization=false,
    uniform_dt_tolerance=1e-12,
    specrange_kwargs...
)
    tlist = convert(Vector{Float64}, tlist)
    controls = get_controls(generator)
    controlvals = [discretize(control, tlist) for control in controls]

    G = _pwc_get_max_genop(generator, controls, tlist)

    parameters = _pwc_process_parameters(parameters, controls, tlist)
    if isnothing(control_ranges)
        control_ranges = IdDict(
            control => (minimum(controlvals[i]), maximum(controlvals[i])) for
            (i, control) in enumerate(controls)
        )
    else  # check that control_ranges is a dict `control => (E_min, E_max)`
        for control ∈ controls
            @assert control ∈ keys(control_ranges)
            @assert control_ranges[control][1] ≤ control_ranges[control][2]
        end
    end
    E_min, E_max = cheby_get_spectral_envelope(
        generator,
        tlist,
        control_ranges,
        specrange_method;
        specrange_kwargs...
    )
    Δ = E_max - E_min
    @assert Δ > 0.0
    δ = specrange_buffer * Δ
    E_min = E_min - δ / 2
    Δ = Δ + δ
    dt = _get_uniform_dt(tlist; tol=uniform_dt_tolerance, warn=true)
    if isnothing(dt)
        error("Chebychev propagation only works on a uniform time grid")
    end
    wrk = Cheby.ChebyWrk(state, Δ, E_min, dt; limit=cheby_coeffs_limit)
    n = 1
    t = tlist[1]
    if backward
        n = length(tlist) - 1
        t = float(tlist[n+1])
    end
    GT = typeof(generator)
    OT = typeof(G)
    ST = typeof(state)
    return ChebyPropagator{GT,OT,ST}(
        generator,
        inplace ? copy(state) : state,
        t,
        n,
        tlist,
        parameters,
        controls,
        control_ranges,
        G,
        wrk,
        backward,
        inplace,
        specrange_method,
        specrange_buffer,
        check_normalization,
        specrange_kwargs,
    )
end

# Aliases (for backwards compatibility)
init_prop(state, generator, tlist, method::Val{:cheby}; kwargs...) =
    init_prop(state, generator, tlist, Val(:Cheby); kwargs...)


_transform_control_ranges(c, ϵ_min, ϵ_max, check) = (ϵ_min, ϵ_max)


"""
```julia
reinit_prop!(
    propagator::ChebyPropagator,
    state;
    transform_control_ranges=((c, ϵ_min, ϵ_max, check) => (ϵ_min, ϵ_max)),
    kwargs...
)
```

re-initializes an existing [`ChebyPropagator`](@ref). This may or may not
involve recalculating the Chebychev coefficients based on the current control
amplitudes in `propagator.parameters`.

# Method-specific keyword arguments

* `transform_control_ranges`: a function
  `(c, ϵ_min, ϵ_max, check) => (ϵ_min′, ϵ_max′)`.
  For each control `c`, the function is called with `check=true` and `ϵ_min`
  (`ϵ_max`) the current minimum (maximum) values for the control from
  `propagator.parameters`). The Chebychev coefficients will be recalculated if
  the existing coefficients were obtained assuming a range for `c` outside the
  *returned* `ϵ_min′, ϵ_max′`.

  If the coefficients do need to be recalculated, `transform_control_ranges` is
  called a second time with `check=false`, and the returned `(ϵ_min′, ϵ_max′)`
  are used for estimating the new spectral range.

  For example,

  ```
  function transform_control_ranges(c, ϵ_min, ϵ_max, check)
      if check
          return (min(ϵ_min, 2 * ϵ_min), max(ϵ_max, 2 * ϵ_max))
      else
          return (min(ϵ_min, 5 * ϵ_min), max(ϵ_max, 5 * ϵ_max))
      end
  end
  ```

  will re-calculate the Chebychev coefficients only if the current amplitudes
  differ by more than a factor of two from the ranges that were used when
  initializing the propagator (`control_ranges` parameter in
  [`init_prop`](@ref), which would have had to overestimate the actual
  amplitudes by at least a factor of two).  When re-calculating, the
  `control_ranges` will overestimate the amplitudes by a factor of five. With
  this `transform_control_ranges`, the propagation will be stable as long as
  the amplitudes do not change dynamically by more than a factor of 2.5 from
  their original range, while also not re-calculating coefficients
  unnecessarily in each pass because of modest changes in the amplitudes.

  The `transform_control_ranges` argument is only relevant in the context of
  optimal control, where the same `propagator` will be used for many iterations
  with changing control field amplitudes.


All other keyword arguments are ignored.
"""
function reinit_prop!(
    propagator::ChebyPropagator,
    state;
    transform_control_ranges=_transform_control_ranges,
    _...
)
    set_state!(propagator, state)

    wrk = propagator.wrk
    need_to_recalculate_cheby_coeffs = false
    control_ranges = IdDict(
        control => (
            minimum(propagator.parameters[control]),
            maximum(propagator.parameters[control])
        ) for control in propagator.controls
    )
    for control in propagator.controls
        ϵ_min = control_ranges[control][1]
        ϵ_max = control_ranges[control][2]
        ϵ_min_check, ϵ_max_check = transform_control_ranges(control, ϵ_min, ϵ_max, true)
        if (
            (ϵ_min_check < propagator.control_ranges[control][1]) ||
            (ϵ_max_check > propagator.control_ranges[control][2])
        )
            need_to_recalculate_cheby_coeffs = true
            break
        end
    end
    tlist = propagator.tlist
    if need_to_recalculate_cheby_coeffs
        for control in propagator.controls
            ϵ_min = control_ranges[control][1]
            ϵ_max = control_ranges[control][2]
            control_ranges[control] = transform_control_ranges(control, ϵ_min, ϵ_max, false)
        end
        E_min, E_max = cheby_get_spectral_envelope(
            getfield(propagator, :generator),
            tlist,
            control_ranges,
            propagator.specrange_method;
            propagator.specrange_options...
        )
        Δ = E_max - E_min
        @assert Δ > 0.0
        δ = propagator.specrange_buffer * Δ
        E_min = E_min - δ / 2
        Δ = Δ + δ
        dt = float(tlist[2] - tlist[1])
        propagator.control_ranges = control_ranges
        wrk = Cheby.ChebyWrk(state, Δ, E_min, dt; limit=wrk.limit)
    else
        reset_timer!(propagator.wrk.timing_data)
    end
    t = float(propagator.backward ? tlist[end] : tlist[1])
    _pwc_set_t!(propagator, t)
    propagator.wrk = wrk
end


"""Determine the spectral envelope of a `generator`.

```julia
E_min, E_max = cheby_get_spectral_envelope(
    generator, tlist, control_ranges, method; kwargs...
)
```

estimates a lower bound `E_min` the lowest eigenvalue of the generator for any
values of the controls specified by `control_ranges`, and an upper bound
`E_max` for the highest eigenvalue.

This is done by constructing operators from the extremal values for the
controls as specified in `control_ranges` and taking the smallest/largest
return values from [`specrange`](@ref
QuantumPropagators.SpectralRange.specrange) for those operators.

# Arguments

* `generator`: dynamical generator, e.g. a time-dependent
* `tlist`: The time grid for the propagation
* `control_ranges`: a dict that maps controls that occur in `generator` (cf.
  [`get_controls`](@ref) to a tuple of minimum and maximum amplitude for that
  control
* `method`: method name to pass to  [`specrange`](@ref
  QuantumPropagators.SpectralRange.specrange)
* `kwargs`: Any remaining keyword arguments are passed to [`specrange`](@ref
  QuantumPropagators.SpectralRange.specrange)
"""
function cheby_get_spectral_envelope(generator, tlist, control_ranges, method; kwargs...)
    min_vals = IdDict(control => r[1] for (control, r) ∈ control_ranges)
    # TODO: this does not take into account explicit time dependencies (control
    # amplitude ≠ control function). For now, we just take any explicit time
    # dependency in the middle of the time grid.
    n = length(tlist) ÷ 2
    G_min = evaluate(generator, tlist, n; vals_dict=min_vals)
    max_vals = IdDict(control => r[2] for (control, r) ∈ control_ranges)
    G_max = evaluate(generator, tlist, n; vals_dict=max_vals)
    E_min, E_max = SpectralRange.specrange(G_max, method; kwargs...)
    _E_min, _E_max = SpectralRange.specrange(G_min, method; kwargs...)
    E_min = (_E_min < E_min) ? _E_min : E_min
    E_max = (_E_max > E_max) ? _E_max : E_max
    return E_min, E_max
end


function prop_step!(propagator::ChebyPropagator)
    @timeit_debug propagator.wrk.timing_data "prop_step!" begin
        Ψ = propagator.state
        H = propagator.genop
        n = propagator.n
        dt = propagator.wrk.dt
        if propagator.backward
            dt = -dt
        end
        tlist = getfield(propagator, :tlist)
        (0 < n < length(tlist)) || return nothing
        _pwc_set_genop!(propagator, n)
        if propagator.inplace
            Cheby.cheby!(
                Ψ,
                H,
                dt,
                propagator.wrk;
                check_normalization=propagator.check_normalization
            )
        else
            Ψ = Cheby.cheby(
                Ψ,
                H,
                dt,
                propagator.wrk;
                check_normalization=propagator.check_normalization
            )
            setfield!(propagator, :state, Ψ)
        end
        _pwc_advance_time!(propagator)
        return propagator.state
    end
end
