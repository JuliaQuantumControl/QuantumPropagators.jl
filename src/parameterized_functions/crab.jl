# SPDX-FileCopyrightText: © 2021 Michael Goerz <mail@michaelgoerz.net>
#
# SPDX-License-Identifier: MIT

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
via the `parameters` property. The properties `cos_indices` and `sin_indices`
are the ranges of indices of the ``c^{(+)}_i`` and ``c^{(-)}_i`` in
`parameters`, that is, `f.parameters[f.cos_indices[i]]` is ``c^{(+)}_i`` and
`f.parameters[f.sin_indices[i]]` is ``c^{(-)}_i``. A range is empty if
`parity` excludes the respective coefficients. If `scale_guess` is `true`,
``c_0`` is `f.parameters[1]`.
"""
struct CRABFunction{G,S} <: ParameterizedFunction
    parameters::Vector{Float64}
    frequencies::Vector{Float64}
    guess::G
    shape::S
    scale_guess::Bool
    parity::Symbol
    cos_indices::UnitRange{Int64}
    sin_indices::UnitRange{Int64}
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
    # The parameters are [c₀; c⁽⁺⁾; c⁽⁻⁾; r], each block only if applicable
    n_cos = (parity == :odd) ? 0 : N
    n_sin = (parity == :even) ? 0 : N
    cos_indices = (1+scale_guess):(scale_guess+n_cos)
    sin_indices = (last(cos_indices)+1):(last(cos_indices)+n_sin)
    r_indices = (last(sin_indices)+1):(last(sin_indices)+N)
    index_ranges = if vary_frequencies
        (cos_indices, sin_indices, r_indices)
    else
        (cos_indices, sin_indices)
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
    return F(parameters, frequencies, guess, shape, scale_guess, parity, index_ranges...)
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
of the parameters ``r_i``. The keyword arguments are the same as for
[`CRABFunction`](@ref), but the default `parameters` are initialized with
`vary_frequencies=true`, and thus the last `N` elements of the `parameters`
contain the values ``r_i``. In addition to `cos_indices` and `sin_indices`,
the property `r_indices` is the range of indices of the ``r_i`` in
`parameters`, that is, `f.parameters[f.r_indices[i]]` is ``r_i``.
"""
struct VariedFrequencyCRABFunction{G,S} <: ParameterizedFunction
    parameters::Vector{Float64}
    frequencies::Vector{Float64}
    guess::G
    shape::S
    scale_guess::Bool
    parity::Symbol
    cos_indices::UnitRange{Int64}
    sin_indices::UnitRange{Int64}
    r_indices::UnitRange{Int64}
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


# The frequency of the i'th term (ωᵢ or rᵢ ωᵢ)
_crab_frequency(crab::CRABFunction, i) = crab.frequencies[i]
_crab_frequency(crab::VariedFrequencyCRABFunction, i) =
    crab.parameters[crab.r_indices[i]] * crab.frequencies[i]


# Evaluate a `CRABFunction` or `VariedFrequencyCRABFunction` at time `t`. For
# `parity=:evenodd`, `sincos` evaluates both basis functions of a frequency at
# about the cost of one. The type parameter `C` forces specialization: since
# both types are subtypes of `Function`, Julia would otherwise not specialize
# on `crab`, which this method passes on without calling it.
function _evaluate_crab(
    crab::C,
    t
) where {C<:Union{CRABFunction,VariedFrequencyCRABFunction}}
    p = crab.parameters
    f = 0.0
    if crab.parity == :evenodd
        for i in eachindex(crab.frequencies)
            (s, c) = sincos(_crab_frequency(crab, i) * t)
            f += p[crab.cos_indices[i]] * c + p[crab.sin_indices[i]] * s
        end
    elseif crab.parity == :even
        for i in eachindex(crab.frequencies)
            f += p[crab.cos_indices[i]] * cos(_crab_frequency(crab, i) * t)
        end
    else  # :odd
        for i in eachindex(crab.frequencies)
            f += p[crab.sin_indices[i]] * sin(_crab_frequency(crab, i) * t)
        end
    end
    isnothing(crab.shape) || (f *= crab.shape(t))
    if !isnothing(crab.guess)
        g = crab.guess(t)
        f += crab.scale_guess ? p[1] * g : g
    end
    return f
end

(crab::CRABFunction)(t) = _evaluate_crab(crab, t)
(crab::VariedFrequencyCRABFunction)(t) = _evaluate_crab(crab, t)
