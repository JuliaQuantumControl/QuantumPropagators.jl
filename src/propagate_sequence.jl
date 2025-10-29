"""Wrapper around the parameters of a call to [`propagate`](@ref).

```julia
Propagation(
    generator, tlist;
    pre_propagation=nothing, post_propagation=nothing,
    kwargs...
)

Propagation(
    propagator;
    pre_propagation=nothing, post_propagation=nothing,
    kwargs...
)
```

is a wrapper around the arguments for [`propagate`](@ref) /
[`init_prop`](@ref), for use within [`propagate_sequence`](@ref).

The positional and keyword arguments are those accepted by the above
propagation routines, excluding the initial state. A `Propagation` may in
addition include the `pre_propagation` and `post_propagation` keyword arguments
recognized by [`propagate_sequence`](@ref).
"""
struct Propagation
    args::Vector{Any}
    kwargs::Dict{Symbol,Any}
    function Propagation(args...; kwargs...)
        new(Any[args...], kwargs)
    end
end


"""Propagate a state through a sequence of generators.

```
states = propagate_sequence(
    state,
    propagations;
    storage=nothing,
    pre_propagation=nothing,
    post_propagation=nothing,
    kwargs...
)
```

takes an `initial` state and performs a sequence of [`propagate`](@ref) calls
using the parameters in `propagations`. The initial state for each step in the
sequence is the state resulting from the previous step. Optionally, before and
after each step, a `pre_propagation` and `post_propagation` function may modify
the state instantaneously, e.g., to perform a frame transformation. Return the
vector of states at the end of each step (after any `post_propagation`, before
any next `pre_propagation` of the next step).

# Arguments

* `state`: The initial state
* `propagations`: A vector of [`Propagation`](@ref) instances, one per step in
  the sequence, each containing the arguments for the call to
  [`propagate`](@ref) for that step. The [`Propagation`](@ref) contains the
  generator and time grid for each step as positional parameters, or
  alternatively a pre-initialized [`Propagator`](@ref AbstractPropagator), and
  any keyword arguments for [`propagate`](@ref) that are specific to that step.
  Note that [`propagate`](@ref) keyword arguments that are common to all steps
  can be given directly to `propagate_sequence`.
* `storage`: If `storage=true`, return a vector of storage objects as returned
  by `propagate(…, storage=true)` for each propagation step, instead of the
  state after each step. To use a pre-initialized `storage`, each
  [`Propagation`](@ref) in `propagations` should have a `storage` keyword
  argument instead.
* `pre_propagation`: If not `nothing`, must be a function that receives the
  same arguments as [`propagate`](@ref) and returns a state. Called immediately
  before the [`propagate`](@ref) of each step, and the state returned by
  `pre_propagation` will become the initial state for the subsequent call to
  [`propagate`](@ref). Generally, `pre_propagation` would be different in each
  step of the sequence, and should be given as a keyword argument in a
  particular [`Propagation`](@ref).
* `post_propagation`: If not `nothing`, a function that receives the same
  arguments as [`propagate`](@ref) and returns a state, see `pre_propagation`.
  The returned state becomes the initial state for the next step in the
  sequence (and may be further processed by the following `pre_propagation`).
  Like `pre_propagation`, this will generally be set as a keyword argument for
  a particular [`Propagation`](@ref), not as a global keyword argument to
  `propagate_sequence`.

All other keyword arguments are forwarded to `propagate`. Thus, keyword
arguments that are common to all steps in the sequence should be given as
keyword arguments to `propagate_sequence` directly.
"""
function propagate_sequence(
    state,
    propagations::Vector{Propagation};
    pre_propagation = nothing,
    post_propagation = nothing,
    kwargs...,
)
    Ψ = state
    results = []
    for prop in propagations
        prop_kwargs = Dict(prop.kwargs...)
        _pre_propagation = pop!(prop_kwargs, :pre_propagation, pre_propagation)
        _post_propagation = pop!(prop_kwargs, :post_propagation, post_propagation)
        merge!(prop_kwargs, kwargs)
        @assert !haskey(prop_kwargs, :pre_propagate)  # from prototype version
        @assert !haskey(prop_kwargs, :post_propagate)  # from prototype version
        storage = get(prop_kwargs, :storage, nothing)
        if storage === true
            observables = get(prop_kwargs, :observables, _StoreState())
            if prop.args[end] isa AbstractPropagator
                tlist = prop.args[end].tlist
            else
                tlist = prop.args[end]
            end
            prop_kwargs[:storage] =
                QuantumPropagators.Storage.init_storage(Ψ, tlist, observables)
        end
        if !isnothing(_pre_propagation)
            Ψ = _pre_propagation(Ψ, prop.args...; prop_kwargs...)
        end
        Ψ = propagate(Ψ, prop.args...; prop_kwargs...)
        if !isnothing(_post_propagation)
            Ψ = _post_propagation(Ψ, prop.args...; prop_kwargs...)
        end
        if storage === true
            push!(results, prop_kwargs[:storage])
        else
            push!(results, copy(Ψ))
        end
    end
    return results
end
