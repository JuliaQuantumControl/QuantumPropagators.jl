module QuantumPropagators

include("./arnoldi.jl")  # submodule Arnoldi
include("./specrad.jl")  # submodule SpectralRange
include("./cheby.jl")    # submodule Cheby
include("./newton.jl")   # submodule Newton
include("./expprop.jl")  # submodule ExpProp
include("./storage.jl")  # submodule Storage

# high-level interface
include("./propagate.jl")
export initpropwrk, init_storage, propstep!, propstep, propagate

using .Storage
export init_storage, write_to_storage!, get_from_storage!

end
