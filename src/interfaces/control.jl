using Test

using ..Controls: discretize, discretize_on_midpoints, get_parameters, ParameterizedFunction


"""Check that `control` can be evaluated on a time grid.

```julia
@test check_control(
    control;
    tlist,
    for_parameterization=true,
    for_time_continuous=(control isa Function),
    quiet=false
)
```

verifies the given `control` (one of the elements of the tuple returned by
[`get_controls`](@ref)):

* [`evaluate(control, tlist, n)`](@ref evaluate) must be defined and return a
  `Float64`
* [`evaluate(control, tlist, n; vals_dict=IdDict(control => v))`](@ref evaluate)
  must be defined and return `v`
* [`discretize(control, tlist)`](@ref discretize) must be defined and return a
  vector of floats of the same size as `tlist`. Only if `length(tlist) > 2`.
* all values in [`discretize(control, tlist)`](@ref discretize) must be finite
  (`isfinite`).
* [`discretize_on_midpoints(control, tlist)`](@ref discretize_on_midpoints)
  must be defined and return a vector of floats with one element less than
  `tlist`. Only if `length(tlist) > 2`.
* all values in [`discretize_on_midpoints(control, tlist)`](@ref discretize)
  must be finite (`isfinite`)

If `for_time_continuous`:

* [`evaluate(control, t)`](@ref evaluate) must be defined and return a
  `Float64`
* [`evaluate(control, t; vals_dict=IdDict(control => v))`](@ref evaluate)
  must be defined and return `v`

If `for_parameterization`:

* [`get_parameters(control)`](@ref get_parameters) must be defined and return a
  vector of floats. Mutating that vector must mutate the control.

The function returns `true` for a valid control and `false` for an invalid
control. Unless `quiet=true`, it will log an error to indicate which of the
conditions failed.
"""
function check_control(
    control;
    tlist,
    for_parameterization=true,
    for_time_continuous=(control isa Function),
    quiet=false,
    _message_prefix=""  # for recursive calling
)

    @assert tlist isa Vector{Float64}

    px = _message_prefix
    success = true

    if for_parameterization
        if control isa ParameterizedFunction
            success &=
                check_parameterized_function(control; tlist, quiet, _message_prefix=px)
        else
            success &= check_parameterized(control; quiet, _message_prefix=px)
        end
    end

    try
        v = evaluate(control, tlist, 1)
        if !(v isa Float64)
            quiet ||
                @error "$(px)`evaluate(control, tlist, n)` must return a Float64, not $(typeof(v))"
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`evaluate(control, tlist, n)` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end
    try
        v₀ = rand()
        v = evaluate(control, tlist, 1; vals_dict=IdDict(control => v₀))
        if v != v₀
            msg = "$(px)`evaluate(control, tlist, n; vals_dict=IdDict(control => v))` must return `v`"
            quiet || @error msg v = v₀ result = v
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`evaluate(control, tlist, n; vals_dict=IdDict(control => v))`  must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    if for_time_continuous
        try
            v = evaluate(control, tlist[1])
            if !(v isa Float64)
                quiet ||
                    @error "$(px)`evaluate(control, t)` must return a Float64, not $(typeof(v))"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`evaluate(control, t)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end
        try
            v₀ = rand()
            v = evaluate(control, tlist[1]; vals_dict=IdDict(control => v₀))
            if v != v₀
                msg = "$(px)`evaluate(control, t; vals_dict=IdDict(control => v))` must return `v`"
                quiet || @error msg v = v₀ result = v
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`evaluate(control, t; vals_dict=IdDict(control => v))`  must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end
    end

    if length(tlist) > 2

        try
            a = discretize(control, tlist)
            if !(a isa Vector{Float64})
                quiet ||
                    @error "$(px)`discretize(control, tlist)` must return a Vector{Float64}, not $(typeof(a))"
                success = false
            end
            if !(length(a) == length(tlist))
                quiet ||
                    @error "$(px)`discretize(control, tlist)` must return a vector of the same size as `tlist` ($(length(a)) ≠ $(length(tlist)))"
                success = false
            end
            a = discretize(control, tlist)
            if !all(isfinite.(a))
                quiet ||
                    @error "$(px) all values in `discretize(control, tlist)` must be finite"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`discretize(control, tlist)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            a = discretize_on_midpoints(control, tlist)
            if !(a isa Vector{Float64})
                quiet ||
                    @error "$(px)`discretize_on_midpoints(control, tlist)` must return a Vector{Float64}, not $(typeof(a))"
                success = false
            end
            if length(a) ≠ length(tlist) - 1
                quiet ||
                    @error "$(px)`discretize_on_midpoints(control, tlist)` must return a vector of length one less than `tlist` ($(length(a)) ≠ $(length(tlist))-1)"
                success = false
            end
            if !all(isfinite.(a))
                quiet ||
                    @error "$(px)all values in `discretize_on_midpoints(control, tlist)` must be finite"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`discretize_on_midpoints(control, tlist)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

    end

    return success

end
