using Test

using ..Controls: get_controls, evaluate, evaluate!
using ..Generators: Generator

"""Check the dynamical `generator` for propagating `state` over `tlist`.

```
@test check_generator(generator; state, tlist,
                     for_mutable_state=true, for_immutable_state=true,
                     for_expval=true, atol=1e-15)
```

verifies the given `generator`:

* [`get_controls(generator)`](@ref get_controls) must be defined and return a
  tuple
* all controls returned by [`get_controls(generator)`](@ref get_controls) must
  pass [`check_control`](@ref)
* [`evaluate(generator, tlist, n)`](@ref evaluate) must return a valid
  operator ([`check_operator`](@ref)).
* [`evaluate!(op, generator, tlist, n)`](@ref evaluate!) must be defined
* [`substitute(generator, replacements)`](@ref substitute) must be defined
* If `generator` is a [`Generator`](@ref) instance, all elements of
  `generator.amplitudes` must be valid, i.e., pass [`check_amplitude`](@ref).

"""
function check_generator(
    generator;
    state,
    tlist,
    for_mutable_state=true,
    for_immutable_state=true,
    for_expval=true,
    atol=1e-15
)

    @assert check_state(state; for_mutable_state, for_immutable_state, atol)
    @assert tlist isa Vector{Float64}
    @assert length(tlist) >= 2

    success = true

    try
        controls = get_controls(generator)
        if !(controls isa Tuple)
            @error "`get_controls(generator)` must return a tuple, not $(typeof(controls))"
            success = false
        end
    catch exc
        @error "`get_controls(generator)` must be defined: $exc"
        success = false
    end

    try
        op = evaluate(generator, tlist, 1)
        if !check_operator(
            op;
            state,
            tlist,
            for_mutable_state,
            for_immutable_state,
            for_expval,
            atol
        )
            @error "`evaluate(generator, tlist, n)` must return an operator that passes `check_operator`"
            success = false
        end
    catch exc
        @error "`evaluate(generator, tlist, n)` must return a valid operator: $exc"
        success = false
    end

    try
        op = evaluate(generator, tlist, 1)
        evaluate!(op, generator, tlist, length(tlist) - 1)
    catch exc
        @error "`evaluate!(op, generator, tlist, n)` must be defined: $exc"
        success = false
    end

    try
        controls = get_controls(generator)
        for (i, control) ∈ enumerate(controls)
            if !check_control(control; tlist)
                @error "control $i in `generator` must pass check_control"
                success = false
            end
        end
    catch exc
        @error "all controls in `generator` must pass `check_control`: $exc"
        success = false
    end

    try
        controls = get_controls(generator)
        replacements = IdDict(ϵ => ϵ for ϵ ∈ controls)
        generator2 = substitute(generator, replacements)
        if get_controls(generator) ≠ get_controls(generator2)
            @error "`substitute(generator, replacements)` must replace the controls in `generator`"
            success = false
        end
    catch exc
        @error "`substitute(generator, replacements)` must be defined: $exc"
        success = false
    end

    if generator isa Generator
        try
            for (i, ampl) in enumerate(generator.amplitudes)
                if !check_amplitude(ampl; tlist)
                    @error "amplitude $i in `generator` does not pass `check_amplitude`"
                    success = false
                end
            end
        catch exc
            @error "all elements of `generator.amplitudes` must be valid: $exc"
            success = false
        end
    end

    return success

end
