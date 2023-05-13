module QuantumPropagators

include("./arnoldi.jl")   # submodule Arnoldi
include("./specrad.jl")   # submodule SpectralRange
include("./cheby.jl")     # submodule Cheby
include("./newton.jl")    # submodule Newton
include("./expprop.jl")   # submodule ExpProp

include("./storage.jl")   # submodule Storage

include("./shapes.jl")  # submodule Shapes

include("./controls.jl")  # submodule Controls

include("./amplitudes.jl")  # submodule Amplitudes

include("./generators.jl")  # submodule Generators


using .Generators
export liouvillian, hamiltonian

include("./propagator.jl")
export init_prop, reinit_prop!, prop_step!
# not exported: set_t!, set_state!, choose_propmethod

include("./pwc_utils.jl")
include("./cheby_propagator.jl")
include("./newton_propagator.jl")
include("./exp_propagator.jl")

#! format: off
module Interfaces
    export check_operator, check_state, check_amplitude, check_control
    export check_generator
    include(joinpath("interfaces", "state.jl"))
    include(joinpath("interfaces", "operator.jl"))
    include(joinpath("interfaces", "amplitude.jl"))
    include(joinpath("interfaces", "control.jl"))
    include(joinpath("interfaces", "generator.jl"))
end
#! format: on

# high-level interface
include("./propagate.jl")
export propagate



end
