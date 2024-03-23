using ..Controls: get_controls, evaluate, evaluate!
using ..Generators: Generator

"""Check the dynamical `generator` for propagating `state` over `tlist`.

```julia
@test check_generator(generator; state, tlist,
                     for_mutable_state=true, for_immutable_state=true,
                     for_expval=true, for_parameterization=false,
                     atol=1e-14, quiet=false)
```

verifies the given `generator`:

* [`get_controls(generator)`](@ref get_controls) must be defined and return a
  tuple
* all controls returned by [`get_controls(generator)`](@ref get_controls) must
  pass [`check_control`](@ref)
* [`evaluate(generator, tlist, n)`](@ref evaluate) must return a valid
  operator ([`check_operator`](@ref)), with forwarded keyword arguments
  (including `for_expval`)
* [`evaluate!(op, generator, tlist, n)`](@ref evaluate!) must be defined
* [`substitute(generator, replacements)`](@ref substitute) must be defined
* If `generator` is a [`Generator`](@ref) instance, all elements of
  `generator.amplitudes` must pass [`check_amplitude`](@ref) with
  `for_parameterization`.

If `for_parameterization` (may require the `RecursiveArrayTools` package to be
loaded):

* [`get_parameters(generator)`](@ref get_parameters) must be defined and return
  a vector of floats. Mutating that vector must mutate the controls inside the
  `generator`.

The function returns `true` for a valid generator and `false` for an invalid
generator. Unless `quiet=true`, it will log an error to indicate which of the
conditions failed.
"""
function check_generator(
    generator;
    state,
    tlist,
    for_mutable_state=true,
    for_immutable_state=true,
    for_expval=true,
    for_parameterization=false,
    atol=1e-14,
    quiet=false,
    _message_prefix="",  # for recursive calling
    _check_amplitudes=true  # undocumented (internal use)
)

    @assert check_state(state; for_mutable_state, for_immutable_state, atol, quiet=true)
    @assert tlist isa Vector{Float64}
    @assert length(tlist) >= 2

    px = _message_prefix
    success = true

    try
        controls = get_controls(generator)
        if !(controls isa Tuple)
            quiet ||
                @error "$(px)`get_controls(generator)` must return a tuple, not $(typeof(controls))"
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`get_controls(generator)` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    if for_parameterization
        success &= check_parameterized(generator; _message_prefix=px)
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
            atol,
            quiet,
            _message_prefix="On `op = evaluate(generator, tlist, 1)`: "
        )
            quiet ||
                @error "$(px)`evaluate(generator, tlist, n)` must return an operator that passes `check_operator`"
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`evaluate(generator, tlist, n)` must return a valid operator.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        op = evaluate(generator, tlist, 1)
        evaluate!(op, generator, tlist, length(tlist) - 1)
    catch exc
        quiet || @error(
            "$(px)`evaluate!(op, generator, tlist, n)` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        controls = get_controls(generator)
        for (i, control) ∈ enumerate(controls)
            valid_control = check_control(
                control;
                tlist,
                for_parameterization=false,
                # We already check the parametrization for the entire
                # generator, so it would be redundant to do it again for
                # the controls
                quiet,
                _message_prefix="On control $i ($(typeof(control))) in `generator`: "
            )
            if !valid_control
                quiet ||
                    @error "$(px)control $i ($(typeof(control))) in `generator` must pass check_control"
                success = false
            end
        end
    catch exc
        quiet || @error(
            "$(px)all controls in `generator` must pass `check_control`.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        controls = get_controls(generator)
        replacements = IdDict(ϵ => ϵ for ϵ ∈ controls)
        generator2 = substitute(generator, replacements)
        if get_controls(generator) ≠ get_controls(generator2)
            quiet ||
                @error "$(px)`substitute(generator, replacements)` must replace the controls in `generator`"
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`substitute(generator, replacements)` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    if (generator isa Generator) && _check_amplitudes
        try
            for (i, ampl) in enumerate(generator.amplitudes)
                valid_ampl = check_amplitude(
                    ampl;
                    tlist,
                    for_parameterization=false,
                    # We already check the parametrization for the entire
                    # generator, so it would be redundant to do it again for
                    # the amplitudes
                    quiet,
                    _message_prefix="On ampl $i ($(typeof(ampl))) in `generator`: "
                )
                if !valid_ampl
                    quiet ||
                        @error "$(px)amplitude $i in `generator` does not pass `check_amplitude`"
                    success = false
                end
            end
        catch exc
            quiet || @error(
                "$(px)all elements of `generator.amplitudes` must be valid.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end
    end

    return success

end
