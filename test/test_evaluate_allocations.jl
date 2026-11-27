# SPDX-FileCopyrightText: © 2026 Michael Goerz <mail@michaelgoerz.net>
#
# SPDX-License-Identifier: MIT

using Test
using QuantumPropagators.Controls: evaluate, StaticParameter
using QuantumPropagators.Amplitudes: ShapedAmplitude

# `C` type parameter: Julia does not specialize on a `::Function` argument that
# is only passed through, which would add a spurious allocation
allocs_tn(c::C, tlist, n; kw...) where {C} = @allocated evaluate(c, tlist, n; kw...)
allocs_t(c::C, t; kw...) where {C} = @allocated evaluate(c, t; kw...)

@testset "Non-allocating evaluate" begin
    tlist = collect(range(0, 1; length = 11))
    f = t -> 2.0 * t
    sp = StaticParameter([1.0])
    v = collect(range(0, 1; length = 10))
    ampl = ShapedAmplitude(sp; shape = t -> 1.0)
    empty_dict = IdDict{Function,Float64}()
    for c in (f, sp, v, ampl)
        allocs_tn(c, tlist, 2)  # compile
        @test allocs_tn(c, tlist, 2) == 0
    end
    for c in (f, sp, ampl)
        allocs_t(c, 0.5)
        @test allocs_t(c, 0.5) == 0
        allocs_t(c, 0.5; vals_dict = empty_dict)
        @test allocs_t(c, 0.5; vals_dict = empty_dict) == 0
    end
end
