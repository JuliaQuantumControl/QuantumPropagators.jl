# SPDX-FileCopyrightText: © 2021 Michael Goerz <mail@michaelgoerz.net>
#
# SPDX-License-Identifier: MIT

module QuantumPropagatorsRecursiveArrayToolsExt

using RecursiveArrayTools: ArrayPartition
import QuantumPropagators.Controls: _combine_parameter_arrays


function _combine_parameter_arrays(parameter_arrays)
    return ArrayPartition(parameter_arrays...)
end

end
