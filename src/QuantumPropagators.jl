module QuantumPropagators


const VERSION = let
    project = joinpath(dirname(dirname(pathof(QuantumPropagators))), "Project.toml")
    Base.include_dependency(project) # Retrigger precompilation when Project.toml changes
    toml = read(project, String)
    m = match(r"(*ANYCRLF)^version\s*=\s\"(.*)\"$"m, toml)
    VersionNumber(m[1])
end


include("arnoldi.jl")   # submodule Arnoldi
include("specrad.jl")   # submodule SpectralRange
include("cheby.jl")     # submodule Cheby
include("newton.jl")    # submodule Newton
include("expprop.jl")   # submodule ExpProp

include("storage.jl")   # submodule Storage

include("shapes.jl")  # submodule Shapes

include("controls.jl")  # submodule Controls

include("parameterized_functions.jl")  # submodule ParameterizedFunctions

include("amplitudes.jl")  # submodule Amplitudes

include("generators.jl")  # submodule Generators


using .Generators
export liouvillian, hamiltonian

include("propagator.jl")
export init_prop, reinit_prop!, prop_step!
# not exported: set_t!, set_state!

#! format: off
module Interfaces
    export supports_inplace, supports_vector_interface, supports_matrix_interface
    export check_operator, check_state, check_tlist, check_amplitude
    export check_control, check_generator, check_propagator
    export check_parameterized_function, check_parameterized
    include("interfaces/supports_inplace.jl")
    include("interfaces/supports_vector_interface.jl")
    include("interfaces/supports_matrix_interface.jl")
    include("interfaces/utils.jl")
    include("interfaces/state.jl")
    include("interfaces/tlist.jl")
    include("interfaces/operator.jl")
    include("interfaces/amplitude.jl")
    include("interfaces/control.jl")
    include("interfaces/generator.jl")
    include("interfaces/propagator.jl")
    include("interfaces/parameterization.jl")
end
#! format: on

include("pwc_utils.jl")
include("exponential_utilities_propagator.jl")
include("cheby_propagator.jl")
include("newton_propagator.jl")
include("exp_propagator.jl")

include("ode_function.jl")


include("timings.jl")

# high-level interface
include("./propagate.jl")
export propagate

include("./propagate_sequence.jl")
export Propagation, propagate_sequence

end
