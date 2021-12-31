module QuantumPropagators

# low-level routines

include("./arnoldi.jl")

include("./specrad.jl")
export specrange

include("./cheby.jl")
export cheby_coeffs, cheby_coeffs!, ChebyWrk, cheby!, cheby

include("./newton.jl")
export NewtonWrk, newton!

include("./expprop.jl")
export ExpPropWrk, expprop!

include("./storage.jl")
export init_storage, map_observables, map_observable, write_to_storage!
export get_from_storage!

# high-level interface
include("./propagate.jl")
export initpropwrk, init_storage, propstep!, propstep, propagate

end
