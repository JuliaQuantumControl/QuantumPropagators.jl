"""Indicate whether a type implements the 2D array interface.

```julia
supports_matrix_interface(::Type{T})
```

returns `true` if `T` implements the
[Abstract Array interface](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array)
for two-dimensional arrays. This is `true` for all subtypes of
`AbstractMatrix`, but may also be `true` for types that implement an array
interface (`size`, `getindex`, etc.) without declaring themselves subtypes
of `AbstractMatrix`. Calling `supports_matrix_interface` on an instance
`x` also works via a convenience fallback that forwards to
`supports_matrix_interface(typeof(x))`.

Depending the value of [`supports_inplace`](@ref), `T` may implement a
read-write matrix (`setindex!` etc.) or a read-only matrix.

The matrix interface is encouraged for [operators](@ref Operators), and the
specific conditions of the required interface in that context are checked via
[`QuantumPropagators.Interfaces.check_operator`](@ref).
"""
supports_matrix_interface(::Type{<:AbstractMatrix}) = true

supports_matrix_interface(::Type{T}) where {T} = false

supports_matrix_interface(x) = supports_matrix_interface(typeof(x))
