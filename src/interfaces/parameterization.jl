# SPDX-FileCopyrightText: © 2021 Michael Goerz <mail@michaelgoerz.net>
#
# SPDX-License-Identifier: MIT

using ..Controls: get_controls, get_parameters, ParameterizedFunction
using RecursiveArrayTools: ArrayPartition
using LinearAlgebra: norm


"""Check a [`ParameterizedFunction`](@ref) instance.

```julia
@test check_parameterized_function(f; tlist; quiet=false)
```

verifies that the given `f`:

* is an instance of [`ParameterizedFunction`](@ref).
* has a field `parameters` that is an `AbstractVector`.
* is a [callable](@extref Julia Function-like-objects) as `f(t)` for values of
  `t` in `tlist`, returning a `Float64`.
* passes [`check_parameterized`](@ref), which verifies that the `parameters`
  field is a single [parameter block](@ref ParameterBlocks) that meets the
  requirements on a block.
* has the `parameters` field as the only block in
  [`get_parameters(f)`](@ref get_parameters).

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

    if hasfield(typeof(f), :parameters)
        parameters = getfield(f, :parameters)
        if !(parameters isa AbstractVector)
            quiet ||
                @error "$(px)`f.parameters` must be an AbstractVector, not $(typeof(parameters))"
            success = false
        end
    else
        quiet || @error "$(px)`f` must have a `parameters` field"
        success = false
    end

    try
        valid_parameters = check_parameterized(f; quiet, _message_prefix = px)
        success &= valid_parameters
        p = getfield(f, :parameters)
        # If `valid_parameters`, `get_parameters(f)` is an `ArrayPartition`
        if valid_parameters && (get_parameters(f).x !== (isempty(p) ? () : (p,)))
            quiet ||
                @error "$(px)`f.parameters` must be the only parameter block in `get_parameters(f)`"
            success = false
        end
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

verifies that [`get_parameters(object)`](@ref get_parameters) is defined and
returns a [`RecursiveArrayTools.ArrayPartition`](@extref), that a second call
returns an identical (`===`) result, and that the
[parameter blocks](@ref ParameterBlocks) in it (its field `x`) meet the
requirements on a block:

1. It is an `AbstractVector{Float64}` with one-based indexing. Verified for
   every block.
2. It can be mutated in place, and every holder reads its values when it is
   evaluated. Verified partially: writing the values `1:N` into
   `get_parameters(object)` element by element (for `N` parameters) must
   change what a subsequent call to `get_parameters(object)` returns, and so
   must writing the original values back with `copyto!`. Whether the holders
   read the values is not verified.
3. It is not a [`RecursiveArrayTools.ArrayPartition`](@extref). Verified for
   the `parameters` field of `object` if it is a
   [`ParameterizedFunction`](@ref), and of every `ParameterizedFunction` among
   the [controls](@ref get_controls) of `object`. That field must be a single
   block.
4. Any two blocks are either identical (`===`) or share no memory. Verified by
   the element-by-element writes: an overlap between blocks changes one of the
   values written earlier.

The original values of the parameters are restored afterwards.

The function returns `true` if all checks pass and `false` otherwise. Unless
`quiet=true`, it logs an error for every failed check.

# See also

* [`check_parameterized_function`](@ref) if `object` is a
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

    local parameters
    try
        parameters = get_parameters(object)
        if !(parameters isa ArrayPartition)
            quiet ||
                @error "$(px)`get_parameters($name)` must return an ArrayPartition, not $(typeof(parameters))"
            success = false
        elseif get_parameters(object) !== parameters
            quiet ||
                @error "$(px)repeated calls to `get_parameters($name)` must return identical (`===`) results"
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`get_parameters($name)` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end
    success || return false

    # Requirement 1, for every block
    for (k, block) in enumerate(parameters.x)
        if !(block isa AbstractVector) || (eltype(block) != Float64)
            quiet ||
                @error "$(px)parameter block $k of `get_parameters($name)` must be a vector of Float64, not $(typeof(block))"
            success = false
        elseif Base.has_offset_axes(block)
            quiet ||
                @error "$(px)parameter block $k of `get_parameters($name)` must use one-based indexing, not axes $(axes(block))"
            success = false
        end
    end

    # Requirement 3, for the parameterized functions in `object`
    for f in _parameterized_functions(object)
        if hasfield(typeof(f), :parameters)
            p = getfield(f, :parameters)
            if p isa ArrayPartition
                quiet ||
                    @error "$(px)the `parameters` of $(typeof(f)) must be a single parameter block, not an `ArrayPartition` ($(typeof(p)))"
                success = false
            end
        end
    end

    # Requirements 2 and 4
    if success && !isempty(parameters)
        orig_parameters = _values(parameters)
        try
            for (k, i) in enumerate(eachindex(parameters))
                parameters[i] = Float64(k)
            end
            values = _values(get_parameters(object))
            if values != 1:length(parameters)
                msg = "$(px)writing to the array returned by `get_parameters($name)` element by element must mutate the parameter blocks inside the object. Parameter blocks must not share memory."
                changed = [k for (k, v) in enumerate(values) if v != k]
                quiet || @error msg changed
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)Cannot mutate `get_parameters($name)`.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        finally
            try
                copyto!(parameters, orig_parameters)
                Δ = norm(_values(get_parameters(object)) - orig_parameters)
                if Δ > 1e-12
                    msg = "$(px)`copyto!` the array returned by `get_parameters($name)` must mutate the parameter blocks inside the object"
                    quiet || @error msg get_parameters(object) orig_parameters Δ
                    success = false
                end
            catch exc
                quiet || @error(
                    "$(px)Cannot restore the original values of `get_parameters($name)` with `copyto!`.",
                    exception = (exc, catch_abbreviated_backtrace())
                )
                success = false
            end
        end
    end

    return success

end


# The `object` itself if it is a `ParameterizedFunction`, and every
# `ParameterizedFunction` among its controls
function _parameterized_functions(object)
    functions = ParameterizedFunction[]
    (object isa ParameterizedFunction) && push!(functions, object)
    if applicable(get_controls, object)
        for control in get_controls(object)
            if (control isa ParameterizedFunction) && !any(f -> f === control, functions)
                push!(functions, control)
            end
        end
    end
    return functions
end


# The values of a vector, as a `Vector{Float64}` (which, e.g., a typed
# comprehension over a `ComponentVector` does not return)
function _values(v)
    values = Vector{Float64}(undef, length(v))
    for (k, x) in enumerate(v)
        values[k] = x
    end
    return values
end
