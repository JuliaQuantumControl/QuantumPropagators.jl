using Test

using ..Controls: discretize, discretize_on_midpoints


"""Check that `control` can be evaluated on a time grid.

```julia
@test check_control(control; tlist)
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
"""
function check_control(control; tlist)

    @assert tlist isa Vector{Float64}

    success = true

    try
        a = discretize(control, tlist)
        if !(a isa Vector{Float64})
            @error "`discretize(control, tlist)` must return a Vector{Float64}, not $(typeof(a))"
            success = false
        end
        if !(length(a) == length(tlist))
            @error "`discretize(control, tlist)` must return a vector of the same size as `tlist` ($(length(a)) ≠ $(length(tlist)))"
            success = false
        end
        a = discretize(control, tlist)
        if !all(isfinite.(a))
            @error "all values in `discretize(control, tlist)` must be finite"
            success = false
        end
    catch exc
        @error "`discretize(control, tlist)` must be defined: $exc"
        success = false
    end

    try
        a = discretize_on_midpoints(control, tlist)
        if !(a isa Vector{Float64})
            @error "`discretize_on_midpoints(control, tlist)` must return a Vector{Float64}, not $(typeof(a))"
            success = false
        end
        if length(a) ≠ length(tlist) - 1
            @error "`discretize_on_midpoints(control, tlist)` must return a vector of length one less than `tlist` ($(length(a)) ≠ $(length(tlist))-1)"
            success = false
        end
        if !all(isfinite.(a))
            @error "all values in `discretize_on_midpoints(control, tlist)` must be finite"
            success = false
        end
    catch exc
        @error "`discretize_on_midpoints(control, tlist)` must be defined: $exc"
        success = false
    end

    return success

end
