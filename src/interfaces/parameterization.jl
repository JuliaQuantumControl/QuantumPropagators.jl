using ..Controls: get_parameters, ParameterizedFunction
using LinearAlgebra: norm


"""Check a [`ParameterizedFunction`](@ref) instance.

```julia
@test check_parameterized_function(f; tlist; quiet=false)
```

verifies that the given `f`:

* is an instance of [`ParameterizedFunction`](@ref).
* has a field `parameters` that is an `AbstractVector{Float64}`.
* is a [callable](@extref Julia Function-like-objects) as `f(t)` for values of
  `t` in `tlist`, returning a `Float64`.
* [`get_parameters`](@ref) provides access to the `parameters` field.
* passes [`check_parameterized`](@ref)

# See also

* [`check_parameterized`](@ref) for objects that have parameters
  ([`get_parameters`](@ref)), but are not instances of
  [`ParameterizedFunction`](@ref)
"""
function check_parameterized_function(
    f;
    tlist,
    quiet = false,
    _message_prefix = "",  # for recursive calling
)
    px = _message_prefix
    success = true

    if !(f isa ParameterizedFunction)
        quiet || @error "$(px)`f` must be an instance of ParameterizedFunction"
        success = false
    end

    if !hasfield(typeof(f), :parameters)
        quiet || @error "$(px)`f` must have a `parameters` field"
        success = false
    end

    try
        success &= (get_parameters(f) === getfield(f, :parameters))
        success &= check_parameterized(f)
    catch exc
        quiet || @error(
            "$(px)`get_parameters(f)` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    t = tlist[begin]
    try
        v = f(t)
        if !(v isa Float64)
            quiet || @error "$(px)`f(t)` must return a Float64, not $(typeof(v))"
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`f(t)` must be defined for t=$t.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    return success

end


@doc raw"""
Check that that the object supports the parameterization interface.

```julia
@test check_parameterized(object; name="::$typeof(object))", quiet=false)
```

verifies that the given `object`:

* can be passed to [`get_parameters`](@ref), which must return an
  `AbstractVector` of `Float64`
* is mutated by mutating the `parameters` obtained by `get_parameters`

# See also

* [`check_parameterized_function`](@ref) is `object` is a
  [`ParameterizedFunction`](@ref)
"""
function check_parameterized(
    object;
    name = "::$(typeof(object))",
    quiet = false,
    _message_prefix = "",  # for recursive calling
)
    px = _message_prefix
    success = true

    try
        parameters = get_parameters(object)
        try
            # we only check `eltype` to allow for AbstractVector
            if !(eltype(parameters) == Float64)
                quiet ||
                    @error "$(px)`get_parameters($name)` must return a vector of Float64, not $(typeof(parameters))"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`get_parameters($name)` must return a vector of Float64.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`get_parameters($name)` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    if success && !isempty(get_parameters(object))
        parameters = get_parameters(object)
        orig_parameters = copy(Array(parameters))
        try
            parameters .= 1.0
            if !all(isone, get_parameters(object))
                quiet ||
                    @error "$(px)writing to the array returned by `get_parameters($name)` must mutate the the controls inside the object"
                success = false
            end
            copyto!(parameters, orig_parameters)
            Δ = norm(get_parameters(object) - orig_parameters)
            if Δ > 1e-12
                msg = "$(px)`copyto!` the array returned by `get_parameters($name)` must mutate the controls inside the object"
                quiet || @error msg get_parameters(object) orig_parameters Δ
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)Cannot mutate `get_parameters($name)`.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end
    end

    return success

end
