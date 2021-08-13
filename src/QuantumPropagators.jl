module QuantumPropagators

# low-level routines
include("./cheby.jl")
export cheby_coeffs, cheby_coeffs!, ChebyWrk, cheby!

include("./newton.jl")
export NewtonWrk, newton!

include("./expprop.jl")
export ExpPropWrk, expprop!

# high-level interface
include("./propagate.jl")
export initpropwrk, init_storage, propstep!, propagate

end
