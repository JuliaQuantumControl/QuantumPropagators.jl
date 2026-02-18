import ..Operator
import ..ScaledOperator
import LinearAlgebra
import SparseArrays: SparseMatrixCSC

"""Indicate whether a type supports in-place operations.

```julia
supports_inplace(::Type{T})
```

where `T` is the type of a state or operator, dispatches on the type `T` to
determine in-place support. Calling `supports_inplace(x)` on an instance also
works via a convenience fallback that forwards to
`supports_inplace(typeof(x))`.

For states, a `true` result indicates that propagators can assume that the
in-place requirements defined in
[`QuantumPropagators.Interfaces.check_state`](@ref) hold. States with in-place
support must also fulfill specific properties when interacting with operators,
see [`QuantumPropagators.Interfaces.check_operator`](@ref).

For operators, a `true` result indicates that the operator can be evaluated
in-place with [`evaluate!`](@ref), see
[`QuantumPropagators.Interfaces.check_generator`](@ref).

Note that `supports_inplace` is not quite the same as
[`Base.ismutabletype`](@extref) and/or [`Base.ismutable`](@extref): When using
[custom structs](@extref Julia :label:`Mutable-Composite-Types`) for states or
operators, even if those structs are not defined as `mutable`, they may still
define the in-place interface (typically because their *components* are
mutable).
"""
supports_inplace(::Type{<:Vector{ComplexF64}}) = true
supports_inplace(::Type{T}) where {T<:AbstractVector} = ismutabletype(T)  # fallback
# The fallback doesn't actually guarantee that the required interface implied
# by `supports_inplace` is fulfilled, but it's a reasonable expectation to
# have, and the `check_state` function will test it.

supports_inplace(::Type{<:Matrix}) = true
supports_inplace(::Type{<:Operator}) = true
supports_inplace(::Type{<:LinearAlgebra.Diagonal}) = true
supports_inplace(::Type{<:SparseMatrixCSC}) = true
supports_inplace(::Type{<:ScaledOperator{<:Any,OT}}) where {OT} = supports_inplace(OT)
supports_inplace(::Type{T}) where {T<:AbstractMatrix} = ismutabletype(T)  # fallback

# Generic catch-all for types without a specific method (prevents StackOverflow
# from the value→type fallback below)
supports_inplace(::Type{T}) where {T} = throw(MethodError(supports_inplace, (T,)))

# Convenience fallback: forward from values to types
supports_inplace(x) = supports_inplace(typeof(x))

# Note: methods for StaticArrays are defined in package extension
