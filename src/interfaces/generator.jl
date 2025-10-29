using ..Controls: get_controls, evaluate, evaluate!
using ..Generators: Generator

"""Check the dynamical `generator` for propagating `state` over `tlist`.

```julia
@test check_generator(
    generator; state, tlist,
    for_pwc=true, for_time_continuous=false,
    for_expval=true, for_parameterization=false,
    atol=1e-14, quiet=false)
```

verifies the given `generator`:

* [`get_controls(generator)`](@ref get_controls) must be defined and return a
  tuple
* all controls returned by [`get_controls(generator)`](@ref get_controls) must
  pass [`check_control`](@ref)
* [`substitute(generator, replacements)`](@ref substitute) must be defined
* If `generator` is a [`Generator`](@ref) instance, all elements of
  `generator.amplitudes` must pass [`check_amplitude`](@ref) with
  `for_parameterization`.

If `for_pwc` (default):

* [`op = evaluate(generator, tlist, n)`](@ref evaluate) must return a valid
  operator ([`check_operator`](@ref)), with forwarded keyword arguments
  (including `for_expval`)
* If `QuantumPropagators.Interfaces.supports_inplace(op)` is `true`,
  [`evaluate!(op, generator, tlist, n)`](@ref evaluate!) must be defined

If `for_time_continuous`:

* [`evaluate(generator, t)`](@ref evaluate) must return a valid
  operator ([`check_operator`](@ref)), with forwarded keyword arguments
  (including `for_expval`)

* If `QuantumPropagators.Interfaces.supports_inplace(op)` is `true`,
  [`evaluate!(op, generator, t)`](@ref evaluate!) must be defined

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
    for_expval = true,
    for_pwc = true,
    for_time_continuous = false,
    for_parameterization = false,
    atol = 1e-14,
    quiet = false,
    _message_prefix = "",  # for recursive calling
    _check_amplitudes = true  # undocumented (internal use)
)

    @assert check_state(state; atol, quiet = true)
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
        success &= check_parameterized(generator; _message_prefix = px)
    end

    try
        controls = get_controls(generator)
        for (i, control) ∈ enumerate(controls)
            valid_control = check_control(
                control;
                tlist,
                for_parameterization = false,
                # We already check the parametrization for the entire
                # generator, so it would be redundant to do it again for
                # the controls
                quiet,
                _message_prefix = "On control $i ($(typeof(control))) in `generator`: "
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

    vals_dict = IdDict()
    try
        if success
            if for_pwc
                vals_dict = IdDict(
                    control => evaluate(control, tlist, 1) for
                    control in get_controls(generator)
                )
            elseif for_time_continuous
                vals_dict = IdDict(
                    control => evaluate(control, tlist[1]) for
                    control in get_controls(generator)
                )
            end
        end
    catch exc
        quiet || @error(
            "$(px)`evaluate(control, …)` must be defined for all controls in generator.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    if for_pwc

        try
            op = evaluate(generator, tlist, 1)
            if !check_operator(
                op;
                state,
                tlist,
                for_expval,
                atol,
                quiet,
                _message_prefix = "On `op = evaluate(generator, tlist, 1)`: "
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

        op = evaluate(generator, tlist, 1)

        try
            if !check_operator(
                op;
                state,
                tlist,
                for_expval,
                atol,
                quiet,
                _message_prefix = "On `op = evaluate(generator, tlist, 1; vals_dict)`: "
            )
                quiet ||
                    @error "$(px)`evaluate(generator, tlist, n; vals_dict)` must return an operator that passes `check_operator`"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`evaluate(generator, tlist, n; vals_dict)` must return a valid operator.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        if success && supports_inplace(op)
            try
                evaluate!(op, generator, tlist, length(tlist) - 1)
            catch exc
                quiet || @error(
                    "$(px)`evaluate!(op, generator, tlist, n)` must be defined.",
                    exception = (exc, catch_abbreviated_backtrace())
                )
                success = false
            end
            try
                op = evaluate(generator, tlist, 1)
                evaluate!(op, generator, tlist, length(tlist) - 1; vals_dict)
            catch exc
                quiet || @error(
                    "$(px)`evaluate!(op, generator, tlist, n; vals_dict)` must be defined.",
                    exception = (exc, catch_abbreviated_backtrace())
                )
                success = false
            end
        end

    end

    if for_time_continuous

        try
            op = evaluate(generator, tlist[1])
            if !check_operator(
                op;
                state,
                tlist,
                for_expval,
                atol,
                quiet,
                _message_prefix = "On `op = evaluate(generator, tlist[1])`: "
            )
                quiet ||
                    @error "$(px)`evaluate(generator, t)` must return an operator that passes `check_operator`"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`evaluate(generator, t)` must return a valid operator.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            op = evaluate(generator, tlist[1]; vals_dict)
            if !check_operator(
                op;
                state,
                tlist,
                for_expval,
                atol,
                quiet,
                _message_prefix = "On `op = evaluate(generator, tlist[1]; vals_dict)`: "
            )
                quiet ||
                    @error "$(px)`evaluate(generator, t; vals_dict)` must return an operator that passes `check_operator`"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`evaluate(generator, t; vals_dict)` must return a valid operator.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        op = evaluate(generator, tlist[begin])
        if success && supports_inplace(op)
            try
                evaluate!(op, generator, tlist[end])
            catch exc
                quiet || @error(
                    "$(px)`evaluate!(op, generator, t)` must be defined.",
                    exception = (exc, catch_abbreviated_backtrace())
                )
                success = false
            end
            try
                op = evaluate(generator, tlist[begin])
                evaluate!(op, generator, tlist[end]; vals_dict)
            catch exc
                quiet || @error(
                    "$(px)`evaluate!(op, generator, t; vals_dict)` must be defined.",
                    exception = (exc, catch_abbreviated_backtrace())
                )
                success = false
            end
        end
    end


    if (generator isa Generator) && _check_amplitudes
        try
            for (i, ampl) in enumerate(generator.amplitudes)
                valid_ampl = check_amplitude(
                    ampl;
                    tlist,
                    for_parameterization = false,
                    # We already check the parametrization for the entire
                    # generator, so it would be redundant to do it again for
                    # the amplitudes
                    quiet,
                    _message_prefix = "On ampl $i ($(typeof(ampl))) in `generator`: "
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
