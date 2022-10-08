module QuantumPropagators

include("./arnoldi.jl")   # submodule Arnoldi
include("./specrad.jl")   # submodule SpectralRange
include("./cheby.jl")     # submodule Cheby
include("./newton.jl")    # submodule Newton
include("./expprop.jl")   # submodule ExpProp

include("./storage.jl")   # submodule Storage

using .Storage
export init_storage, write_to_storage!, get_from_storage!

include("./generators.jl")  # submodule Generators

using .Generators
export liouvillian, hamiltonian

include("./propagator.jl")
export initprop, reinitprop!, propstep!, set_state!
# not exported: set_t!, choose_propmethod

include("./pwc_utils.jl")
include("./cheby_propagator.jl")
include("./newton_propagator.jl")
include("./exp_propagator.jl")

# high-level interface
include("./propagate.jl")
export propagate

end
