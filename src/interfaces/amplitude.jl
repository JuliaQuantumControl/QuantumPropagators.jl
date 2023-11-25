using Test

using ..Generators: Generator
using ..Controls: get_controls, evaluate, substitute

"""Check amplitude appearing in [`Generator`](@ref).

```julia
@test check_amplitude(ampl; tlist, quiet=false)
```

verifies that the given `ampl` is a valid element in the list of `amplitudes`
of a [`Generator`](@ref) object. Specifically:

* [`get_controls(ampl)`](@ref get_controls) must be defined and return a tuple
* all controls in `ampl` must pass [`check_control`](@ref)
* [`substitute(ampl, controls_replacements)`](@ref substitute) must be defined
* [`evaluate(ampl, tlist, n)`](@ref evaluate) must be defined and return a
  Number
* [`evaluate(ampl, tlist, n; vals_dict)`](@ref evaluate) must be defined and
  return a Number

The function returns `true` for a valid amplitude and `false` for an invalid
amplitude. Unless `quiet=true`, it will log an error to indicate which of the
conditions failed.
"""
function check_amplitude(
    ampl;
    tlist,
    quiet=false,
    _message_prefix=""  # for recursive calling
)

    px = _message_prefix
    success = true

    try
        controls = get_controls(ampl)
        if !(controls isa Tuple)
            quiet ||
                @error "$(px)get_controls(ampl) must return a tuple, not $(typeof(controls))"
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`get_controls(ampl)` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        controls = get_controls(ampl)
        for (i, control) ∈ enumerate(controls)
            valid_control = check_control(
                control;
                tlist,
                quiet,
                _message_prefix="On control $i ($(typeof(control))) in `ampl`: "
            )
            if !valid_control
                quiet ||
                    @error "$(px)control $i ($(typeof(control)) in `ampl` must pass `check_control`"
                success = false
            end
        end
    catch exc
        quiet || @error(
            "$(px)all controls in `ampl` must pass `check_control`.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        controls = get_controls(ampl)
        replacements = IdDict(ϵ => ϵ for ϵ ∈ controls)
        ampl_copy = substitute(ampl, replacements)
        if !(get_controls(ampl_copy) == get_controls(ampl))
            quiet || @error "$(px)`substitute(ampl, replacements)` must replace controls"
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`substitute(ampl, replacements)` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        val = evaluate(ampl, tlist, 1)
        if !(val isa Number)
            quiet ||
                @error "$(px)`evaluate(ampl, tlist, 1)` must return a Number, not $(typeof(val))"
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`evaluate(ampl, tlist, n)` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        controls = get_controls(ampl)
        vals_dict = IdDict(ϵ => 1.0 for ϵ ∈ controls)
        val = evaluate(ampl, tlist, 1; vals_dict)
        if !(val isa Number)
            quiet ||
                @error "$(px)evaluate(ampl, tlist, 1; vals_dict) must return a Number, not $(typeof(val))"
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`evaluate(ampl, tlist, n; vals_dict)` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    return success

end
