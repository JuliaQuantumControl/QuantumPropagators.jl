module QuantumPropagatorsStaticArraysExt

import QuantumPropagators.Interfaces: supports_inplace
using StaticArrays: SArray, MArray

supports_inplace(::Type{<:SArray}) = false
supports_inplace(::Type{<:MArray}) = true

end
