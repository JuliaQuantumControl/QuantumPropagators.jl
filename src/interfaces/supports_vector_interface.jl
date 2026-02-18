"""Indicate whether a type implements the 1D array interface.

```julia
supports_vector_interface(::Type{T})
```

returns `true` if `T` implements the
[Abstract Array interface](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array)
for one-dimensional arrays. This is `true` for all subtypes of
`AbstractVector`, but may also be `true` for types that implement an array
interface (`size`, `getindex`, etc.) without declaring themselves subtypes
of `AbstractVector`. Calling `supports_vector_interface` on an instance
`x` also works via a convenience fallback that forwards to
`supports_vector_interface(typeof(x))`.

Depending the value of [`supports_inplace`](@ref), `T` may implement a
read-write vector (`setindex!` etc.) or a read-only vector.

The vector interface is encouraged for quantum states, and the
specific conditions of the required interface in that context are checked via
[`QuantumPropagators.Interfaces.check_state`](@ref).

"""
supports_vector_interface(::Type{<:AbstractVector}) = true

supports_vector_interface(::Type{T}) where {T} = false

supports_vector_interface(x) = supports_vector_interface(typeof(x))
