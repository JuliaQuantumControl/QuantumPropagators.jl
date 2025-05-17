"""Check that the given `tlist` is valid.

```julia
@test check_tlist(tlist; quiet=false)
```

verifies the given time grid. A valid time grid must

* be a `Vector{Float64}`,
* contain at least two points (beginning and end),
* be monotonically increasing

The function returns `true` for a valid time grid and `false` for an invalid
time grid. Unless `quiet=true`, it will log an error to indicated which of the
conditions failed.
"""
function check_tlist(tlist; quiet=false, _message_prefix="")

    px = _message_prefix
    success = true

    if tlist isa Vector{Float64}
        if length(tlist) >= 2
            try
                _t = tlist[begin]
                for t in tlist[(begin+1):end]
                    if t <= _t
                        quiet || @error "$(px)`tlist` must be monotonically increasing"
                        success = false
                    end
                    _t = t
                end
            catch exc
                quiet || @error(
                    "$(px)Cannot verify `tlist`.",
                    exception = (exc, catch_abbreviated_backtrace())
                )
                success = false
            end
        else
            quiet || @error "$(px)`tlist` must contain at least two points"
            success = false
        end
    else
        quiet || @error "$(px)`tlist` must be a Vector{Float64}, not $(typeof(tlist))"
        success = false
    end

    return success

end
