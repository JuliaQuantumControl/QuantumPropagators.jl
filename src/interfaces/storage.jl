# SPDX-FileCopyrightText: © 2021 Michael Goerz <mail@michaelgoerz.net>
#
# SPDX-License-Identifier: MIT

using LinearAlgebra: norm
using Random: AbstractRNG, default_rng

using ..Storage:
    init_storage, write_to_storage!, get_from_storage, get_from_storage!, map_observables


# Deviation between the stored and the expected data. Falls back to an exact
# comparison for data that does not support `norm` and `-` (e.g., tuples of
# mixed observable values).
function _storage_mismatch(a, b)
    try
        return norm(a - b)
    catch
        return isequal(a, b) ? 0.0 : Inf
    end
end


# Read back every slot of `storage` and compare it to `expected`. Returns
# `true` if all slots match within `atol`.
function _check_storage_roundtrip(storage, expected, atol, quiet, msg)
    for n in eachindex(expected)
        local delta
        try
            delta = _storage_mismatch(get_from_storage(storage, n), expected[n])
        catch exc
            quiet || @error(
                "$msg: `get_from_storage(storage, $n)` must return the data written to slot $n.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            return false
        end
        if !(delta <= atol)
            quiet || @error "$msg: slot $n does not match the written data" delta atol
            return false
        end
    end
    return true
end


"""Check the storage implementation for a `state` or for arbitrary `data`.

```julia
@test check_storage(state, tlist; rng=Random.default_rng(), atol=1e-14, quiet=false)
```

verifies that storing `state` objects in a `storage` object created by 
[`init_storage`](@ref QuantumPropagators.Storage.init_storage)
fulfills [the storage contract](@ref storage-contract) for every point
of `tlist`.

```julia
@test check_storage(state, tlist, observables; rng=Random.default_rng(), atol=1e-14, quiet=false)
```

verifies the same for the data produced by applying `observables` to `state`,
i.e., the full
[`map_observables`](@ref QuantumPropagators.Storage.map_observables) →
[`write_to_storage!`](@ref QuantumPropagators.Storage.write_to_storage!) →
[`get_from_storage`](@ref QuantumPropagators.Storage.get_from_storage) pipeline
that [`propagate`](@ref QuantumPropagators.propagate) relies on.

```julia
@test check_storage(data, nt; transform=(data -> data), atol=1e-14, quiet=false)
```

verifies the write/read round trip for `nt` slots of arbitrary `data`, where
successive values are obtained as `data = transform(data)`. This form checks
ownership only if `transform` mutates its argument in place and returns it;
with the default `transform`, the same object is written to every slot.

The different methods of the check mirror the methods of
[`init_storage`](@ref QuantumPropagators.Storage.init_storage).
The two state-based forms assume that `state` is a valid state (it passes
[`check_state`](@ref QuantumPropagators.Interfaces.check_state)) and that
`tlist` is a valid time grid (it passes
[`check_tlist`](@ref QuantumPropagators.Interfaces.check_tlist)).
They verify the following:

1. Round trip: Each slot is filled with a *fresh* state object of known
   value, and every slot read back with
   [`get_from_storage`](@ref QuantumPropagators.Storage.get_from_storage) must
   reproduce the value written to it.
2. Ownership (only if [`supports_inplace(state)`](@ref QuantumPropagators.Interfaces.supports_inplace)):
   Each slot is filled from a *single* state object that is mutated in place
   between the writes. Reading the slots back, via both
   [`get_from_storage`](@ref QuantumPropagators.Storage.get_from_storage) and
   [`get_from_storage!`](@ref QuantumPropagators.Storage.get_from_storage!),
   must reproduce the values as of the time of each write. This is the part
   that detects a storage that keeps a reference to mutable data: under
   aliasing, every slot reads back as the final value.

The function returns `true` for a valid storage implementation and `false`
otherwise. Unless `quiet=true`, it will log an error to indicate which of the
conditions failed.
"""
function check_storage(
    state,
    tlist;
    rng::AbstractRNG = default_rng(),
    atol = 1e-14,
    quiet = false
)

    success = true
    nt = length(tlist)
    state_initial = state

    inplace = false
    try
        inplace = supports_inplace(state_initial)
    catch exc
        quiet || @error(
            "check_storage: The `QuantumPropagators.Interfaces.supports_inplace` method must be defined for type `$(typeof(state_initial))`.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        return false
    end

    # 1. Round trip: a fresh state object for every slot.
    #
    # Every value is derived from the original state with a factor in
    # [0.5, 1.5) so that all slots stay at order-1 values for any `nt`.
    # Cumulative rescaling would make the norm drift below `atol` after a few
    # dozen slots, leaving late slots effectively unchecked.
    try
        storage = init_storage(state_initial, tlist)
        expected = Any[state_initial]
        write_to_storage!(storage, 1, state_initial)
        for n = 2:nt
            state = (0.5 + rand(rng)) * state_initial  # fresh object; `c * state`
            write_to_storage!(storage, n, state)
            push!(expected, state)
        end
        success &= _check_storage_roundtrip(
            storage,
            expected,
            atol,
            quiet,
            "check_storage (round trip)"
        )
    catch exc
        quiet || @error(
            "check_storage: `init_storage(state, tlist)` and `write_to_storage!` must be defined and compatible.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        return false
    end

    # 2. Ownership: a single state object, mutated in place after each write.
    #    The storage must retain the value as of the time of the write.
    if inplace
        ok = false
        expected = Any[]
        storage = nothing
        try
            storage = init_storage(state_initial, tlist)
            state = copy(state_initial)
            expected = Any[copy(state)]  # the records must be copies!
            write_to_storage!(storage, 1, state)
            for n = 2:nt
                copyto!(state, (0.5 + rand(rng)) * state_initial)  # in-place mutation
                write_to_storage!(storage, n, state)
                push!(expected, copy(state))
            end
            # One more mutation, with no write following it: without this,
            # the final slot still agrees with `state`, and a storage that
            # aliases only that slot would pass.
            copyto!(state, (0.5 + rand(rng)) * state_initial)
            ok = _check_storage_roundtrip(
                storage,
                expected,
                atol,
                quiet,
                "check_storage (ownership): the storage must own its data, not keep a reference to it"
            )
            success &= ok
        catch exc
            quiet || @error(
                "check_storage: `write_to_storage!` must be defined for a mutable state.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            return false
        end
        if ok  # also check in-place extraction
            try
                state = similar(state_initial)
                for n = 1:nt
                    get_from_storage!(state, storage, n)
                    delta = _storage_mismatch(state, expected[n])
                    if !(delta <= atol)
                        quiet ||
                            @error "check_storage: `get_from_storage!` at slot $n does not match the written data" delta atol
                        success = false
                        break
                    end
                end
            catch exc
                quiet || @error(
                    "check_storage: `get_from_storage!(state, storage, i)` must be defined for a mutable state.",
                    exception = (exc, catch_abbreviated_backtrace())
                )
                success = false
            end
        end
    end

    return success

end


function check_storage(
    state,
    tlist,
    observables;
    rng::AbstractRNG = default_rng(),
    atol = 1e-14,
    quiet = false
)

    success = true
    nt = length(tlist)
    state_initial = state

    inplace = false
    try
        inplace = supports_inplace(state_initial)
    catch exc
        quiet || @error(
            "check_storage: The `QuantumPropagators.Interfaces.supports_inplace` method must be defined for type `$(typeof(state_initial))`.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        return false
    end

    try
        storage = init_storage(state_initial, tlist, observables)
        data = map_observables(observables, tlist, 1, state_initial)
        expected = Any[deepcopy(data)]
        write_to_storage!(storage, 1, data)
        for n = 2:nt
            state = (0.5 + rand(rng)) * state_initial
            data = map_observables(observables, tlist, n, state)
            write_to_storage!(storage, n, data)
            push!(expected, deepcopy(data))
        end
        success &= _check_storage_roundtrip(
            storage,
            expected,
            atol,
            quiet,
            "check_storage (observables, round trip)"
        )
    catch exc
        quiet || @error(
            "check_storage: `init_storage(state, tlist, observables)`, `map_observables`, and `write_to_storage!` must be defined and compatible.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        return false
    end

    if inplace
        try
            storage = init_storage(state_initial, tlist, observables)
            state = copy(state_initial)
            data = map_observables(observables, tlist, 1, state)
            expected = Any[deepcopy(data)]
            write_to_storage!(storage, 1, data)
            for n = 2:nt
                copyto!(state, (0.5 + rand(rng)) * state_initial)
                data = map_observables(observables, tlist, n, state)
                write_to_storage!(storage, n, data)
                push!(expected, deepcopy(data))
            end
            # Mutate once more without writing, so that observables holding
            # on to `state` or to a reused buffer expose an alias in the final
            # slot as well.
            copyto!(state, (0.5 + rand(rng)) * state_initial)
            map_observables(observables, tlist, nt, state)
            success &= _check_storage_roundtrip(
                storage,
                expected,
                atol,
                quiet,
                "check_storage (observables, ownership): the storage must own its data, not keep a reference to it"
            )
        catch exc
            quiet || @error(
                "check_storage: `map_observables` and `write_to_storage!` must be defined for a mutable state.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            return false
        end
    end

    return success

end


function check_storage(
    data,
    nt::Integer;
    transform = (data -> data),
    atol = 1e-14,
    quiet = false
)

    try
        storage = init_storage(data, nt)
        expected = Any[deepcopy(data)]
        write_to_storage!(storage, 1, data)
        for n = 2:nt
            data = transform(data)
            write_to_storage!(storage, n, data)
            push!(expected, deepcopy(data))
        end
        # One more transform, with no write following it: for an in-place
        # `transform`, this exposes a storage that aliases only the final slot.
        transform(data)
        return _check_storage_roundtrip(
            storage,
            expected,
            atol,
            quiet,
            "check_storage (data)"
        )
    catch exc
        quiet || @error(
            "check_storage: `init_storage(data, nt)` and `write_to_storage!` must be defined and compatible.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        return false
    end

end
