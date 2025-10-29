using Test
using QuantumPropagators.Controls
using QuantumPropagators.Shapes: blackman
using OffsetArrays: Origin
using IOCapture: IOCapture


@testset "get_tlist_midpoints" begin

    tlist = Origin(0)([1.0, 3.0, 5.0, 6.0, 7.0])
    @test get_tlist_midpoints(tlist) == [1.0, 4.0, 5.5, 7.0]
    @test get_tlist_midpoints(tlist; preserve_start = false) == [2.0, 4.0, 5.5, 7.0]
    @test get_tlist_midpoints(tlist; preserve_end = false) == [1.0, 4.0, 5.5, 6.5]
    @test get_tlist_midpoints(tlist; preserve_start = false, preserve_end = false) ==
          [2.0, 4.0, 5.5, 6.5]

    c = IOCapture.capture(rethrow = Union{}) do
        get_tlist_midpoints([1.0, 2.0])
    end
    @test c.value isa ArgumentError

    c = IOCapture.capture(rethrow = Union{}) do
        get_tlist_midpoints([0.0, 0.0, 0.0], preserve_start = false, preserve_end = false)
    end
    @test c.value isa AssertionError

    c = IOCapture.capture(rethrow = Union{}) do
        get_tlist_midpoints([0.0, 1.0, 1.0, 0.0])
    end
    @test c.value isa AssertionError

end


@testset "discretize/discretize_on_midpoints" begin
    tlist = collect(range(0, 10, length = 20))

    t₀ = 0.0
    T = 10.0

    function control_func(t)
        return blackman(t, t₀, T)
    end

    control_vals1 = discretize(control_func, tlist; via_midpoints = true)
    pulse_vals1 = discretize_on_midpoints(control_func, tlist)

    control_vals2 = discretize(pulse_vals1, tlist)
    pulse_vals2 = discretize_on_midpoints(control_vals1, tlist)

    control_vals3 = discretize(control_func, tlist; via_midpoints = false)

    control_vals4 = discretize(control_vals1, tlist)
    pulse_vals3 = discretize_on_midpoints(pulse_vals1, tlist)

    # "controls' are on tlist, "pulses" are on intervals
    @test length(control_vals1) == length(control_vals2) == length(tlist)
    @test length(pulse_vals1) == length(pulse_vals2) == length(tlist) - 1

    # discretizing and converting to/from intervals should be interchangeable
    @test maximum(abs.(control_vals1 .- control_vals2)) < 1e-14
    @test maximum(abs.(pulse_vals1 .- pulse_vals2)) < 1e-14

    # discretizing via midpoints or directly should give different results
    @test 1e-3 < maximum(abs.(control_vals1 .- control_vals3)) < 1e-1
    @test control_vals3[13] == control_func(tlist[13])
    @test control_vals1[13] ≠ control_func(tlist[13])
    @test 1e-3 < abs(control_vals1[13] - control_func(tlist[13])) < 1e-1

    # discretizing vectors that are already discretized should return a copy
    @test control_vals4 ≢ control_vals1
    @test maximum(abs.(control_vals4 .- control_vals1)) < 1e-14
    @test pulse_vals3 ≢ pulse_vals1
    @test maximum(abs.(pulse_vals3 .- pulse_vals1)) < 1e-14

end
