using Test
using IOCapture: IOCapture
using QuantumPropagators.Interfaces: check_amplitude, check_control
using QuantumPropagators.Shapes: flattop
using QuantumPropagators.Amplitudes: LockedAmplitude, ShapedAmplitude
using QuantumPropagators.Controls: evaluate, get_controls, substitute, t_mid
using QuantumPropagators.Controls: discretize_on_midpoints

@testset "LockedAmplitude" begin

    tlist = collect(range(0, 10, length = 101))
    S(t) = flattop(t, T = 10, t_rise = 2, func = :blackman)

    ampl = LockedAmplitude(S)
    @test startswith("$ampl", "LockedAmplitude(")
    @test check_amplitude(ampl; tlist)
    @test length(get_controls(ampl)) == 0
    t = 3.1342
    @test evaluate(ampl, t) ≈ S(t)
    t = t_mid(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ S(t)

    ampl = LockedAmplitude(S, tlist)
    @test startswith("$ampl", "LockedAmplitude(")
    @test check_amplitude(ampl; tlist)
    @test length(get_controls(ampl)) == 0
    t = t_mid(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ S(t)
    captured = IOCapture.capture(rethrow = Union{}) do
        evaluate(ampl, t)
    end
    @test captured.error
    @test contains(captured.value.msg, "can only be evaluated with (tlist, n)")

end


@testset "ShapedAmplitude" begin
    tlist = collect(range(0, 10, length = 101))
    S(t) = flattop(t, T = 10, t_rise = 2, func = :blackman)
    ϵ(t) = 1.0

    @test check_control(ϵ; tlist)

    ampl = ShapedAmplitude(ϵ; shape = S)
    @test startswith("$ampl", "ShapedAmplitude(")
    controls = get_controls(ampl)
    @test length(controls) == 1
    @test check_control(controls[1]; tlist)
    @test check_amplitude(ampl; tlist)
    t = 3.1342
    @test evaluate(ampl, t) ≈ ϵ(t) * S(t)
    @test ampl(t) ≈ ϵ(t) * S(t)
    t = t_mid(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ ϵ(t) * S(t)

    ampl = ShapedAmplitude(ϵ, tlist; shape = S)
    @test startswith("$ampl", "ShapedAmplitude(")
    controls = get_controls(ampl)
    @test length(controls) == 1
    @test check_control(controls[1]; tlist)
    @test check_amplitude(ampl; tlist)
    t = t_mid(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ ϵ(t) * S(t)
    captured = IOCapture.capture(rethrow = Union{}) do
        evaluate(ampl, t)
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "vector control and shape can only be evaluated with (tlist, n)"
    )

end


@testset "Vector{Float64} conversion" begin
    tlist = collect(range(0, 10, length = 101))
    S(t) = flattop(t, T = 10, t_rise = 2, func = :blackman)
    ϵ(t) = 1.0

    # LockedAmplitude: callable returning Int is discretized to Vector{Float64}
    ampl = LockedAmplitude(t -> 1, tlist)
    @test ampl.shape isa Vector{Float64}
    @test check_amplitude(ampl; tlist)

    # LockedAmplitude: Vector{Int} is converted to Vector{Float64}
    S_int = ones(Int, length(tlist) - 1)
    ampl = LockedAmplitude(S_int)
    @test ampl.shape isa Vector{Float64}
    @test check_amplitude(ampl; tlist)

    # ShapedAmplitude: callable shape returning Int
    ampl = ShapedAmplitude(ϵ; shape = t -> 1)
    @test check_amplitude(ampl; tlist)
    t = t_mid(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ ϵ(t) * 1

    # ShapedAmplitude: Vector{Int} control and shape are converted
    ϵ_int = ones(Int, length(tlist) - 1)
    S_int = ones(Int, length(tlist) - 1)
    ampl = ShapedAmplitude(ϵ_int; shape = S_int)
    @test ampl.control isa Vector{Float64}
    @test ampl.shape isa Vector{Float64}
    @test check_amplitude(ampl; tlist)

    # ShapedAmplitude tlist constructor: callable returning Int is discretized to Vector{Float64}
    ampl = ShapedAmplitude(t -> 1, tlist; shape = t -> 1)
    @test ampl.control isa Vector{Float64}
    @test ampl.shape isa Vector{Float64}

end


@testset "ShapedAmplitude mixed" begin
    tlist = collect(range(0, 10, length = 101))
    S(t) = flattop(t, T = 10, t_rise = 2, func = :blackman)
    ϵ(t) = 1.0
    S_vec = discretize_on_midpoints(S, tlist)
    ϵ_vec = discretize_on_midpoints(ϵ, tlist)

    # callable control, vector shape
    ampl = ShapedAmplitude(ϵ; shape = S_vec)
    @test startswith("$ampl", "ShapedAmplitude(")
    controls = get_controls(ampl)
    @test length(controls) == 1
    @test check_amplitude(ampl; tlist)
    t = t_mid(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ ϵ(t) * S(t)
    captured = IOCapture.capture(rethrow = Union{}) do
        evaluate(ampl, t)
    end
    @test captured.error
    @test contains(captured.value.msg, "vector shape can only be evaluated with (tlist, n)")

    # vector control, callable shape
    ampl = ShapedAmplitude(ϵ_vec; shape = S)
    @test startswith("$ampl", "ShapedAmplitude(")
    controls = get_controls(ampl)
    @test length(controls) == 1
    @test check_amplitude(ampl; tlist)
    t = t_mid(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ ϵ(t) * S(t)
    captured = IOCapture.capture(rethrow = Union{}) do
        evaluate(ampl, t)
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "vector control can only be evaluated with (tlist, n)"
    )

end


@testset "Invalid constructions" begin
    tlist = collect(range(0, 10, length = 101))
    S(t) = flattop(t, T = 10, t_rise = 2, func = :blackman)
    ϵ(t) = 1.0
    n = length(tlist) - 1

    # LockedAmplitude: complex vector is not convertible to Vector{Float64}
    captured = IOCapture.capture(rethrow = Union{}) do
        LockedAmplitude(fill(1.0 + 1.0im, n))
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "shape that is a vector must be convertible to Vector{Float64}"
    )

    # LockedAmplitude: non-vector, non-callable shape
    captured = IOCapture.capture(rethrow = Union{}) do
        LockedAmplitude(42)
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "shape must either be a Vector{Float64} or a callable"
    )

    # ShapedAmplitude: complex vector control
    captured = IOCapture.capture(rethrow = Union{}) do
        ShapedAmplitude(fill(1.0 + 1.0im, n); shape = S)
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "control that is a vector must be convertible to Vector{Float64}"
    )

    # ShapedAmplitude: complex vector shape
    captured = IOCapture.capture(rethrow = Union{}) do
        ShapedAmplitude(ϵ; shape = fill(1.0 + 1.0im, n))
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "shape that is a vector must be convertible to Vector{Float64}"
    )

    # ShapedAmplitude: non-vector, non-callable control
    captured = IOCapture.capture(rethrow = Union{}) do
        ShapedAmplitude(42; shape = S)
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "control must either be a Vector{Float64} or a callable"
    )

    # ShapedAmplitude: non-vector, non-callable shape
    captured = IOCapture.capture(rethrow = Union{}) do
        ShapedAmplitude(ϵ; shape = 42)
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "shape must either be a Vector{Float64} or a callable"
    )

    # ShapedAmplitude: control and shape vectors of different lengths
    captured = IOCapture.capture(rethrow = Union{}) do
        ShapedAmplitude(ones(n); shape = ones(n - 1))
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "control and shape vectors must have the same length"
    )

end


@testset "ShapedAmplitude substitute" begin
    tlist = collect(range(0, 10, length = 101))
    S(t) = flattop(t, T = 10, t_rise = 2, func = :blackman)
    ϵ(t) = 1.0
    ϵ_vec = discretize_on_midpoints(ϵ, tlist)

    ampl = ShapedAmplitude(ϵ; shape = S)
    ampl2 = substitute(ampl, IdDict(ϵ => ϵ_vec))
    @test ampl2 isa ShapedAmplitude
    @test get_controls(ampl2)[1] === ϵ_vec
    @test check_amplitude(ampl2; tlist)
    t = t_mid(tlist, 20)
    @test evaluate(ampl2, tlist, 20) ≈ ϵ(t) * S(t)
    captured = IOCapture.capture(rethrow = Union{}) do
        evaluate(ampl2, t)
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "vector control can only be evaluated with (tlist, n)"
    )

end
