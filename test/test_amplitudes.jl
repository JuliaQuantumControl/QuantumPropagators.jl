using Test
using QuantumPropagators.Interfaces: check_amplitude, check_control
using QuantumPropagators.Shapes: flattop
using QuantumPropagators.Amplitudes: LockedAmplitude, ShapedAmplitude
using QuantumPropagators.Controls: evaluate, get_controls, _t
using QuantumPropagators.Controls: discretize_on_midpoints # DEBUG

@testset "LockedAmplitude" begin

    tlist = collect(range(0, 10, length=101))
    S(t) = flattop(t, T=10, t_rise=2, func=:blackman)

    ampl = LockedAmplitude(S)
    @test startswith("$ampl", "LockedAmplitude(")
    @test check_amplitude(ampl; tlist)
    @test length(get_controls(ampl)) == 0
    t = 3.1342
    @test evaluate(ampl, t) ≈ S(t)
    t = _t(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ S(t)

    ampl = LockedAmplitude(S, tlist)
    @test startswith("$ampl", "LockedAmplitude(")
    @test check_amplitude(ampl; tlist)
    @test length(get_controls(ampl)) == 0
    t = _t(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ S(t)

end


@testset "ShapedAmplitude" begin
    tlist = collect(range(0, 10, length=101))
    S(t) = flattop(t, T=10, t_rise=2, func=:blackman)
    ϵ(t) = 1.0

    @test check_control(ϵ; tlist)

    ampl = ShapedAmplitude(ϵ; shape=S)
    @test startswith("$ampl", "ShapedAmplitude(")
    controls = get_controls(ampl)
    @test length(controls) == 1
    @test check_control(controls[1]; tlist)
    @test check_amplitude(ampl; tlist)
    t = 3.1342
    @test evaluate(ampl, t) ≈ ϵ(t) * S(t)
    @test ampl(t) ≈ ϵ(t) * S(t)
    t = _t(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ ϵ(t) * S(t)


    ampl = ShapedAmplitude(ϵ, tlist; shape=S)
    @test startswith("$ampl", "ShapedAmplitude(")
    controls = get_controls(ampl)
    @test length(controls) == 1
    @test check_control(controls[1]; tlist)
    @test check_amplitude(ampl; tlist)
    t = _t(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ ϵ(t) * S(t)

end
