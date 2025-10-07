module ParameterizedFunctions

using ..Controls: ParameterizedFunction
export ParameterizedFunction

include("parameterized_functions/crab.jl")
export CRABFunction, VariedFrequencyCRABFunction

end
