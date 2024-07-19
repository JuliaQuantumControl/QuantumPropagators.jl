import ..Operator
import ..ScaledOperator
import LinearAlgebra

"""Indicate whether a given state or operator supports in-place operations

```julia
supports_inplace(state)
```

Indicates that propagators can assume that the in-place requirements defined
in [`QuantumPropagators.Interfaces.check_state`](@ref) hold. States with
in-place support must also fulfill specific properties when interacting with
operators, see [`QuantumPropagators.Interfaces.check_operator`](@ref).

```julia
supports_inplace(op)
```

Indicates that the operator can be evaluated in-place with [`evaluate!`](@ref),
see [`QuantumPropagators.Interfaces.check_generator`](@ref)

Note that `supports_inplace` is not quite the same as
[`Base.ismutable`](@extref): When using [custom structs](@extref Julia
:label:`Mutable-Composite-Types`) for states or operators, even if those
structs are not defined as `mutable`, they may still define the in-place
interface (typically because their *components* are mutable).
"""
supports_inplace(state::Vector{ComplexF64}) = true
supports_inplace(state::AbstractVector) = ismutable(state)  # fallback
# The fallback doesn't actually guarantee that the required interface implied
# by `supports_inplace` is fulfilled, but it's a reasonable expectation to
# have, and the `check_state` function will test it.

supports_inplace(op::Matrix) = true
supports_inplace(op::Operator) = true
supports_inplace(op::LinearAlgebra.Diagonal) = true
supports_inplace(op::ScaledOperator) = supports_inplace(op.operator)
supports_inplace(op::AbstractMatrix) = ismutable(op)  # fallback

# Note: methods for StaticArrays are defined in package extension
