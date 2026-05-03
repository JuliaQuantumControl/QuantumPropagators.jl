using Test
using IOCapture: IOCapture
using QuantumPropagators.Interfaces: check_amplitude, check_control
using QuantumPropagators.Shapes: flattop
using QuantumPropagators.Amplitudes: LockedAmplitude, ShapedAmplitude, GuidedAmplitude
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


@testset "GuidedAmplitude" begin
    tlist = collect(range(0, 10, length = 101))
    S(t) = flattop(t, T = 10, t_rise = 2, func = :blackman)
    G(t) = 0.1 * flattop(t, T = 10, t_rise = 2, func = :blackman)
    ϵ(t) = 1.0
    @test typeof(S) <: Function
    @test typeof(G) <: Function
    @test typeof(ϵ) <: Function

    ampl = GuidedAmplitude(ϵ; guide = G, shape = S)
    @test startswith("$ampl", "GuidedAmplitude(")
    controls = get_controls(ampl)
    @test length(controls) == 1
    @test check_control(controls[1]; tlist)
    @test check_amplitude(ampl; tlist)
    t = 3.1342
    @test evaluate(ampl, t) ≈ G(t) + ϵ(t) * S(t)
    @test ampl(t) ≈ G(t) + ϵ(t) * S(t)
    t = t_mid(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ G(t) + ϵ(t) * S(t)

    ampl = GuidedAmplitude(ϵ, tlist; guide = G, shape = S)
    @test startswith("$ampl", "GuidedAmplitude(")
    controls = get_controls(ampl)
    @test length(controls) == 1
    @test check_control(controls[1]; tlist)
    @test check_amplitude(ampl; tlist)
    t = t_mid(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ G(t) + ϵ(t) * S(t)
    captured = IOCapture.capture(rethrow = Union{}) do
        evaluate(ampl, t)
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "vector control can only be evaluated with (tlist, n)"
    )

end


@testset "GuidedAmplitude Number guide" begin
    tlist = collect(range(0, 10, length = 101))
    S(t) = flattop(t, T = 10, t_rise = 2, func = :blackman)
    ϵ(t) = 1.0
    ϵ_vec = discretize_on_midpoints(ϵ, tlist)

    # Number guide with callable control → constant function
    ampl = GuidedAmplitude(ϵ; guide = 0.1, shape = S)
    @test check_amplitude(ampl; tlist)
    t = 3.1342
    @test evaluate(ampl, t) ≈ 0.1 + ϵ(t) * S(t)
    t = t_mid(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ 0.1 + ϵ(t) * S(t)

    # Number guide with vector control → constant vector
    ampl = GuidedAmplitude(ϵ_vec; guide = 0.1, shape = S)
    @test ampl.guide isa Vector{Float64}
    @test all(==(0.1), ampl.guide)
    @test check_amplitude(ampl; tlist)
    t = t_mid(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ 0.1 + ϵ(t) * S(t)
    captured = IOCapture.capture(rethrow = Union{}) do
        evaluate(ampl, t)
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "vector control can only be evaluated with (tlist, n)"
    )

    # Integer Number guide is converted to Float64
    ampl = GuidedAmplitude(ϵ; guide = 0, shape = S)
    @test check_amplitude(ampl; tlist)
    t = 3.1342
    @test evaluate(ampl, t) ≈ ϵ(t) * S(t)

    # guide = 0 equivalence with ShapedAmplitude
    ampl_shaped = ShapedAmplitude(ϵ; shape = S)
    ampl_guided = GuidedAmplitude(ϵ; guide = 0.0, shape = S)
    t = 3.1342
    @test evaluate(ampl_shaped, t) ≈ evaluate(ampl_guided, t)
    t = t_mid(tlist, 20)
    @test evaluate(ampl_shaped, tlist, 20) ≈ evaluate(ampl_guided, tlist, 20)

    # guide = 0 with tlist constructor
    ampl_shaped = ShapedAmplitude(ϵ, tlist; shape = S)
    ampl_guided = GuidedAmplitude(ϵ, tlist; guide = 0.0, shape = S)
    t = t_mid(tlist, 20)
    @test evaluate(ampl_shaped, tlist, 20) ≈ evaluate(ampl_guided, tlist, 20)

end


@testset "GuidedAmplitude mixed" begin
    tlist = collect(range(0, 10, length = 101))
    S(t) = flattop(t, T = 10, t_rise = 2, func = :blackman)
    G(t) = 0.1 * flattop(t, T = 10, t_rise = 2, func = :blackman)
    ϵ(t) = 1.0
    S_vec = discretize_on_midpoints(S, tlist)
    G_vec = discretize_on_midpoints(G, tlist)
    ϵ_vec = discretize_on_midpoints(ϵ, tlist)

    # callable control (Function subtype), vector shape, vector guide
    @test typeof(ϵ) <: Function
    ampl = GuidedAmplitude(ϵ; guide = G_vec, shape = S_vec)
    @test check_amplitude(ampl; tlist)
    t = t_mid(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ G(t) + ϵ(t) * S(t)
    captured = IOCapture.capture(rethrow = Union{}) do
        evaluate(ampl, t)
    end
    @test captured.error
    @test contains(captured.value.msg, "vector shape can only be evaluated with (tlist, n)")

    # callable control, callable shape, vector guide
    ampl = GuidedAmplitude(ϵ; guide = G_vec, shape = S)
    @test check_amplitude(ampl; tlist)
    t = t_mid(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ G(t) + ϵ(t) * S(t)
    captured = IOCapture.capture(rethrow = Union{}) do
        evaluate(ampl, t)
    end
    @test captured.error
    @test contains(captured.value.msg, "vector guide can only be evaluated with (tlist, n)")

    # callable control, vector shape, callable guide
    ampl = GuidedAmplitude(ϵ; guide = G, shape = S_vec)
    @test check_amplitude(ampl; tlist)
    t = t_mid(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ G(t) + ϵ(t) * S(t)
    captured = IOCapture.capture(rethrow = Union{}) do
        evaluate(ampl, t)
    end
    @test captured.error
    @test contains(captured.value.msg, "vector shape can only be evaluated with (tlist, n)")

    # vector control, callable shape, callable guide
    ampl = GuidedAmplitude(ϵ_vec; guide = G, shape = S)
    @test check_amplitude(ampl; tlist)
    t = t_mid(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ G(t) + ϵ(t) * S(t)
    captured = IOCapture.capture(rethrow = Union{}) do
        evaluate(ampl, t)
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "vector control can only be evaluated with (tlist, n)"
    )

    # callable struct (not a Function subtype) as guide — exercises general evaluate fallback
    struct ConstantGuide
        val::Float64
    end
    (g::ConstantGuide)(t) = g.val
    @test !(typeof(ConstantGuide(0.1)) <: Function)
    ampl = GuidedAmplitude(ϵ; guide = ConstantGuide(0.1), shape = S)
    @test check_amplitude(ampl; tlist)
    t = 3.1342
    @test evaluate(ampl, t) ≈ 0.1 + ϵ(t) * S(t)
    t = t_mid(tlist, 20)
    @test evaluate(ampl, tlist, 20) ≈ 0.1 + ϵ(t) * S(t)

end


@testset "GuidedAmplitude invalid constructions" begin
    tlist = collect(range(0, 10, length = 101))
    S(t) = flattop(t, T = 10, t_rise = 2, func = :blackman)
    G(t) = 0.1 * flattop(t, T = 10, t_rise = 2, func = :blackman)
    ϵ(t) = 1.0
    n = length(tlist) - 1

    # complex vector control
    captured = IOCapture.capture(rethrow = Union{}) do
        GuidedAmplitude(fill(1.0 + 1.0im, n); guide = G, shape = S)
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "control that is a vector must be convertible to Vector{Float64}"
    )

    # complex vector shape
    captured = IOCapture.capture(rethrow = Union{}) do
        GuidedAmplitude(ϵ; guide = G, shape = fill(1.0 + 1.0im, n))
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "shape that is a vector must be convertible to Vector{Float64}"
    )

    # complex vector guide
    captured = IOCapture.capture(rethrow = Union{}) do
        GuidedAmplitude(ϵ; guide = fill(1.0 + 1.0im, n), shape = S)
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "guide that is a vector must be convertible to Vector{Float64}"
    )

    # non-vector, non-callable control
    captured = IOCapture.capture(rethrow = Union{}) do
        GuidedAmplitude(42; guide = G, shape = S)
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "control must either be a Vector{Float64} or a callable"
    )

    # non-vector, non-callable shape
    captured = IOCapture.capture(rethrow = Union{}) do
        GuidedAmplitude(ϵ; guide = G, shape = 42)
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "shape must either be a Vector{Float64} or a callable"
    )

    # non-vector, non-callable guide (not a Number)
    captured = IOCapture.capture(rethrow = Union{}) do
        GuidedAmplitude(ϵ; guide = "bad", shape = S)
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "guide must either be a Vector{Float64} or a callable"
    )

    # mismatched vector lengths (control vs shape)
    captured = IOCapture.capture(rethrow = Union{}) do
        GuidedAmplitude(ones(n); guide = G, shape = ones(n - 1))
    end
    @test captured.error
    @test contains(captured.value.msg, "must all have the same length")

    # mismatched vector lengths (control vs guide)
    captured = IOCapture.capture(rethrow = Union{}) do
        GuidedAmplitude(ones(n); guide = ones(n - 1), shape = S)
    end
    @test captured.error
    @test contains(captured.value.msg, "must all have the same length")

    # mismatched vector lengths (shape vs guide, callable control)
    captured = IOCapture.capture(rethrow = Union{}) do
        GuidedAmplitude(ϵ; guide = ones(n - 1), shape = ones(n))
    end
    @test captured.error
    @test contains(captured.value.msg, "must all have the same length")

end


@testset "GuidedAmplitude substitute" begin
    tlist = collect(range(0, 10, length = 101))
    S(t) = flattop(t, T = 10, t_rise = 2, func = :blackman)
    G(t) = 0.1 * flattop(t, T = 10, t_rise = 2, func = :blackman)
    ϵ(t) = 1.0
    ϵ_vec = discretize_on_midpoints(ϵ, tlist)

    ampl = GuidedAmplitude(ϵ; guide = G, shape = S)
    ampl2 = substitute(ampl, IdDict(ϵ => ϵ_vec))
    @test ampl2 isa GuidedAmplitude
    @test get_controls(ampl2)[1] === ϵ_vec
    @test ampl2.guide === G
    @test ampl2.shape === S
    @test check_amplitude(ampl2; tlist)
    t = t_mid(tlist, 20)
    @test evaluate(ampl2, tlist, 20) ≈ G(t) + ϵ(t) * S(t)
    captured = IOCapture.capture(rethrow = Union{}) do
        evaluate(ampl2, t)
    end
    @test captured.error
    @test contains(
        captured.value.msg,
        "vector control can only be evaluated with (tlist, n)"
    )

end
