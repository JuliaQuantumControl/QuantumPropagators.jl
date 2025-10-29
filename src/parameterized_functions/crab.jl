using ..Controls: ParameterizedFunction
import Random


@doc raw"""CRAB control function.

```julia
f = CRABFunction(N; kwargs...)
```

initializes a [`ParameterizedFunction`](@ref) representing a CRAB (Chopped
Random Basis) [DoriaPRL2011, CanevaPRA2011, RachPRA2015, MuellerRPP2022](@cite)
pulse of the form

```math
f(t)
= c_0 g(t) + S(t) \left(
\sum_{i=1}^{N} c^{(+)}_{i} \cos(ω_i t)
+ \sum_{i=1}^{N} c^{(-)}_{i} \sin(ω_i t)
\right)\,,
```

where ``g(t)`` is a guess function, ``S(t)`` is a shape, ``ω_i`` are `N`
randomly chosen frequencies, and ``c_0``, ``c^{(+)}_i``,
``c^{(-)}_i`` are coefficients ("parameters").

# Keyword arguments

* `frequencies`: The frequencies ``ω_i``. Defaults to
  `sort(max_frequency * rand(rng, N))`. Stored in the `frequencies` property.
* `max_frequency`: The maximum possible frequency for the default
  `frequencies`. This is a required keyword argument, unless `frequencies` is
  given manually. When choosing `max_frequency`, consider the time grid on
  which the `CRABFunction` is eventually going to be evaluated, and the general
  properties of the Fourier grid: ``ω_{max} ≈ π/dt`` and ``dω ≈ 2π/T``. The
  `max_frequency` given here should be at least ``N dω`` (the lowest `N`
  principal frequencies), and at most some fraction of ``ω_{max}``.
* `parameters`: A vector of parameters. By default, set by
  [`crab_initial_parameters`](@ref) with forwarded arguments `guess`,
  `scale_guess`, `random_amplitude`, `parity`, `rng` and
  `vary_frequencies=false`. If given by hand, must have the same size and
  layout as specified in [`crab_initial_parameters`](@ref). Stored in the
  `parameters` property. It is recommended to keep all parameters on the order
  of magnitude of 1, and to use a `shape` to scale ``f(t)`` up to a
  physically relevant amplitude.
* `guess`: The function ``g(t)``. Must be callable as `guess(t)`. Stored in
  the `guess` property. If `nothing` (default), behave like ``g(t) = 0`` and
  force `scale_guess=false`.
* `scale_guess`: Whether to include ``c_0=1`` in `parameters`, see
  [`crab_initial_parameters`](@ref). If `guess` is nothing, this is
  automatically reset to `false`. Stored in the `scale_guess` property.
  Setting `scale_guess=false` in conjunction with a `guess` function can be
  used to to enforce non-zero boundary conditions: ``f(t) = g(t)`` where
  ``S(t) = 0``.
* `random_amplitude`: If `true` (default) initialize ``c^{(+)}_i``,
  ``c^{(-)}_i`` randomly, see [`crab_initial_parameters`](@ref). If `false`,
  with the default `parameters`, ``c^{(+)}_i = c^{(-)}_i = 0``, so that
  ``f(t) ≡ g(t)``. This is useful when `g(t)` is (as the name suggests) a guess
  control for an optimization.
* `rng`: A random number generator. The resulting function is guaranteed to be
  reproducible when passed identically seeded RNGs.
* `shape`: The function ``S(t)``. Must be callable as `shape(t)`. Stored in
  the `shape` property. If `nothing` (default), behave like ``S(t) = 1``. Can
  be used to ensure boundary conditions of ``f(t)`` and to scale the overall
  amplitude of ``f(t)``.
* `parity`: One of `:even`, `:odd`, `:evenodd` (default). Which of
  ``c^{(+)}_i``, ``c^{(-)}_i`` to include in `parameters`. Stored in the
  `parity` property.

# Parameters

After instantiation, the coefficients ``c_0``, ``c^{(+)}_i``,
``c^{(-)}_i`` (if applicable, see [`crab_initial_parameters`](@ref), and
concatenated in that order) are accessible
via the `parameters` property. Specifically, `f.parameters[f.i_cos + i]` is the
parameter ``c^{(+)}_i`` and `f.parameters[f.i_sin + i]` is the parameter
``c^{(-)}_i``, with the offset properties `i_cos` and `i_sin`.
"""
struct CRABFunction{G,S} <: ParameterizedFunction
    parameters::Vector{Float64}
    frequencies::Vector{Float64}
    guess::G
    shape::S
    scale_guess::Bool
    parity::Symbol
    i_sin::Int64
    i_cos::Int64
end


function CRABFunction(
    N;
    max_frequency = 0.0,
    rng = Random.GLOBAL_RNG,
    frequencies = sort(max_frequency * rand(rng, N)),
    guess = nothing,
    shape = nothing,
    parity = :evenodd,
    scale_guess = true,
    random_amplitude = true,
    parameters = crab_initial_parameters(
        N;
        guess,
        scale_guess,
        random_amplitude,
        parity,
        rng
    )
)
    return _instantiate_crab_function(
        CRABFunction{typeof(guess),typeof(shape)},
        N;
        rng,
        frequencies,
        guess,
        shape,
        parity,
        scale_guess,
        parameters,
        vary_frequencies = false,
    )
end


"""Generate random parameters for a CRAB function.

```julia
parameters = crab_initial_parameters(
    N;
    guess=nothing,
    scale_guess=true,
    random_amplitude=false,
    vary_frequencies=false,
    parity=:evenodd,
    rng=Random.GLOBAL_RNG
)
```

generates a vector of parameters in the range ``[-1, 1]``. The resulting
coefficients are, in order,

* ``c_0 = 1`` if `guess` is given as a function (not `nothing`) and
  `scale_guess=true`.
* ``c^{(+)}_i ∈ [-1, 1]`` for ``i ∈ [1, N]`` if `parity` is `:even` or
  `:evenodd`
* ``c^{(-)}_i ∈ [-1, 1]`` for ``i ∈ [1, N]`` if `parity` is `:odd` or
  `:evenodd`
* ``r_i = 1`` for ``i ∈ [1, N]`` if `vary_frequencies=true`.

Thus, the length of the returned parameters could be ``N``, ``N + 1``, ``2N``,
``2 N + 1``, ``3N``, or ``3N+1``, depending on the given parameters.

If `random_amplitude` is `true`, choose ``c^{(+)}_i`` and ``c^{(-)}_i``
randomly from the interval ``[-1, 1]``. If `false`, these coefficients are
initialized as zero.

The parameters are used in the equations of [`CRABFunction`](@ref) and
[`VariedFrequencyCRABFunction`](@ref): ``c_0`` as the coefficient of the guess,
``c^{(+)}_i`` as the coefficients for the ``cos`` functions,
``c^{(-)}_i``  as the coefficients for the ``sin`` functions, and ``r_i`` to
scale the frequencies themselves ([`VariedFrequencyCRABFunction`](@ref)
only).

All random numbers are generated via the given `rng`.
"""
function crab_initial_parameters(
    N;
    guess = nothing,
    scale_guess = true,
    random_amplitude = false,
    vary_frequencies = false,
    parity = :evenodd,
    rng = Random.GLOBAL_RNG
)
    isnothing(guess) && (scale_guess = false)
    guess_weight = scale_guess ? Float64[1.0] : Float64[]
    freq_weights = zeros((parity ∈ (:odd, :even)) ? N : 2N)
    if random_amplitude
        freq_weights .= 1.0 .- 2.0 .* rand(rng, length(freq_weights))
    end
    freq_scales = vary_frequencies ? ones(N) : Float64[]
    return vcat(guess_weight, freq_weights, freq_scales)
end


function _instantiate_crab_function(
    F,  # type of CRAB function to instantiate
    N::Int64;
    frequencies::Vector{Float64},
    guess,
    shape,
    parity::Symbol,
    scale_guess::Bool,
    vary_frequencies::Bool,
    parameters::Vector{Float64},
    rng,
)
    (parity == :oddeven) && (parity = :evenodd)
    isnothing(guess) && (scale_guess = false)
    if parity ∉ [:evenodd, :odd, :even]
        error("parity must be one of `:evenodd`, `:odd`, `:even`, not `$(repr(parity))`")
    end
    if all(iszero(frequencies))
        msg = "The `frequencies` in $(nameof(F)) cannot be all zero. Did you forget to pass `max_frequency`?"
        error(msg)
    end
    if length(frequencies) != N
        @error "Invalid frequencies in $(nameof(F)) function" frequencies N
        error("Length of frequencies $(length(frequencies)) must match N=$N")
    end
    if guess isa Vector
        msg = "$F cannot be instantiated with a vector of pulse values as a guess"
        error(msg)
    end
    i_cos = 0 + scale_guess
    i_sin = i_cos
    if parity != :odd
        i_sin += N
    end
    auto_parameters =
        crab_initial_parameters(N; guess, scale_guess, vary_frequencies, parity)
    N_parameters = length(auto_parameters)
    if length(parameters) != N_parameters
        msg = "Unexpected number of $(nameof(F)) parameters"
        @error msg length(parameters) parity scale_guess
        msg = "Number of parameters must be $N_parameters, not $(length(parameters))"
        error(msg)
    end
    return F(parameters, frequencies, guess, shape, scale_guess, parity, i_sin, i_cos)
end


function (crab::CRABFunction)(t)
    f = 0.0
    c₀ = crab.parameters[begin]  # get c₀ early, so we don't have to go back
    if crab.parity == :even || crab.parity == :evenodd
        for (i, ωᵢ) in enumerate(crab.frequencies)
            f += crab.parameters[crab.i_cos+i] * cos(ωᵢ * t)
        end
    end
    if crab.parity == :odd || crab.parity == :evenodd
        for (i, ωᵢ) in enumerate(crab.frequencies)
            f += crab.parameters[crab.i_sin+i] * sin(ωᵢ * t)
        end
    end
    if !isnothing(crab.shape)
        f *= crab.shape(t)
    end
    if !isnothing(crab.guess)
        if crab.scale_guess
            f += c₀ * crab.guess(t)
        else
            f += crab.guess(t)
        end
    end
    return f
end


@doc raw"""CRAB control function with variable frequencies.

```julia
f = VariedFrequencyCRABFunction(N; kwargs...)
```

initializes a [`ParameterizedFunction`](@ref) representing a CRAB (Chopped
Random Basis) pulse of the form

```math
f(t)
= c_0 g(t) + S(t) \left(
\sum_{i=1}^{N} c^{(+)}_{i} \cos(r_i ω_i t)
+ \sum_{i=1}^{N} c^{(-)}_{i} \sin(r_i ω_i t)
\right)\,.
```

The only difference to the equation in [`CRABFunction`](@ref) is the presence
of the parameters `r_i`. The keyword arguments are the same as for
[`CRABFunction`](@ref), but the default `parameters` are initialized with
`vary_frequencies=true`, and thus the last `N` elements of the `parameters`
contain the values ``r_i``.
"""
struct VariedFrequencyCRABFunction{G,S} <: ParameterizedFunction
    parameters::Vector{Float64}
    frequencies::Vector{Float64}
    guess::G
    shape::S
    scale_guess::Bool
    parity::Symbol
    i_sin::Int64
    i_cos::Int64
end


function VariedFrequencyCRABFunction(
    N;
    max_frequency = 0.0,
    rng = Random.GLOBAL_RNG,
    frequencies = sort(max_frequency * rand(rng, N)),
    guess = nothing,
    shape = nothing,
    parity = :evenodd,
    scale_guess = true,
    random_amplitude = true,
    parameters = crab_initial_parameters(
        N;
        guess,
        scale_guess,
        random_amplitude,
        parity,
        rng,
        vary_frequencies = true,
    )
)
    return _instantiate_crab_function(
        VariedFrequencyCRABFunction{typeof(guess),typeof(shape)},
        N;
        rng,
        frequencies,
        guess,
        shape,
        parity,
        scale_guess,
        parameters,
        vary_frequencies = true,
    )
end



function (crab::VariedFrequencyCRABFunction)(t)
    f = 0.0
    c₀ = crab.parameters[begin]  # get c₀ early, so we don't have to go back
    N = length(crab.frequencies)
    for (i, ωᵢ) in enumerate(crab.frequencies)
        rᵢ = crab.parameters[end-N+i]
        if crab.parity == :even || crab.parity == :evenodd
            f += crab.parameters[crab.i_cos+i] * cos(rᵢ * ωᵢ * t)
        end
        if crab.parity == :odd || crab.parity == :evenodd
            f += crab.parameters[crab.i_sin+i] * sin(rᵢ * ωᵢ * t)
        end
    end
    if !isnothing(crab.shape)
        f *= crab.shape(t)
    end
    if !isnothing(crab.guess)
        if crab.scale_guess
            f += c₀ * crab.guess(t)
        else
            f += crab.guess(t)
        end
    end
    return f
end
