using Test

using ..Generators: Generator
using ..Controls: get_controls, evaluate, substitute

"""Check amplitude appearing in [`Generator`](@ref).

```julia
@test check_amplitude(ampl; tlist)
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
"""
function check_amplitude(ampl; tlist)

    success = true

    try
        controls = get_controls(ampl)
        if !(controls isa Tuple)
            @error "get_controls(ampl) must return a tuple, not $(typeof(controls))"
            success = false
        end
    catch exc
        @error "`get_controls(ampl)` must be defined: $exc"
        success = false
    end

    try
        controls = get_controls(ampl)
        for (i, control) ∈ enumerate(controls)
            if !check_control(control; tlist)
                @error "control $i in `ampl` must pass `check_control`"
                success = false
            end
        end
    catch exc
        @error "all controls in `ampl` must pass `check_control`: $exc"
        success = false
    end

    try
        controls = get_controls(ampl)
        replacements = IdDict(ϵ => ϵ for ϵ ∈ controls)
        ampl_copy = substitute(ampl, replacements)
        if !(get_controls(ampl_copy) == get_controls(ampl))
            @error "`substitute(ampl, replacments)` must replace controls"
            success = false
        end
    catch exc
        @error "`substitute(ampl, replacments)` must be defined: $exc"
        success = false
    end

    try
        val = evaluate(ampl, tlist, 1)
        if !(val isa Number)
            @error "`evaluate(ampl, tlist, 1)` must return a Number, not $(typeof(val))"
            success = false
        end
    catch exc
        @error "`evaluate(ampl, tlist, n)` must be defined: $exc"
        success = false
    end

    try
        controls = get_controls(ampl)
        vals_dict = IdDict(ϵ => 1.0 for ϵ ∈ controls)
        val = evaluate(ampl, tlist, 1; vals_dict)
        if !(val isa Number)
            @error "evaluate(ampl, tlist, 1; vals_dict) must return a Number, not $(typeof(val))"
            success = false
        end
    catch exc
        @error "`evaluate(ampl, tlist, n; vals_dict)` must be defined: $exc"
        success = false
    end

    return success

end
