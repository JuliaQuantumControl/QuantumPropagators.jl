module QuantumPropagatorsStaticArraysExt

import QuantumPropagators.Interfaces: supports_inplace
using StaticArrays: SArray, MArray

supports_inplace(::SArray) = false
supports_inplace(::MArray) = true

end
