using Test

using ..Controls: discretize, discretize_on_midpoints


"""Check that `control` can be evaluated on a time grid.

```julia
@test check_control(control; tlist, quiet=false)
```

verifies the given `control` (one of the elements of the tuple returned by
[`get_controls`](@ref)):

* [`discretize(control, tlist)`](@ref discretize) must be defined and return a
  vector of floats of the same size as `tlist`.
* all values in [`discretize(control, tlist)`](@ref discretize) must be finite
  (`isfinite`)
* [`discretize_on_midpoints(control, tlist)`](@ref discretize_on_midpoints)
  must be defined and return a vector of floats with one element less than
  `tlist`.
* all values in [`discretize_on_midpoints(control, tlist)`](@ref discretize)
  must be finite (`isfinite`)

The function returns `true` for a valid control and `false` for an invalid
control. Unless `quiet=true`, it will log an error to indicate which of the
conditions failed.
"""
function check_control(
    control;
    tlist,
    quiet=false,
    _message_prefix=""  # for recursive calling
)

    @assert tlist isa Vector{Float64}

    px = _message_prefix
    success = true

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

    return success

end
