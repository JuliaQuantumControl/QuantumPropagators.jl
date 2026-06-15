# SPDX-FileCopyrightText: © 2021 Michael Goerz <mail@michaelgoerz.net>
#
# SPDX-License-Identifier: MIT

module ParameterizedFunctions

using ..Controls: ParameterizedFunction
export ParameterizedFunction

include("parameterized_functions/crab.jl")
export CRABFunction, VariedFrequencyCRABFunction

end
