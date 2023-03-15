module Shapes

export flattop, box, blackman


@doc raw"""Flat shape (amplitude 1.0) with a switch-on/switch-off from zero.

```julia
flattop(t; T, t_rise, t₀=0.0, t_fall=t_rise, func=:blackman)
```

evaluates a shape function that starts at 0 at ``t=t₀``, and ramps to to 1
during the `t_rise` interval. The function then remains at value 1, before
ramping down to 0 again during the interval `t_fall` before `T`. For ``t < t₀``
and ``t > T``, the shape is zero.

The default switch-on/-off shape is half of a Blackman window (see
[`blackman`](@ref)).

For `func=:sinsq`, the switch-on/-off shape is a sine-squared curve.
"""
function flattop(t; T, t_rise, t₀=0.0, t_fall=t_rise, func=:blackman)
    if func == :blackman
        return flattop_blackman(t, t₀, T, t_rise, t_fall)
    elseif func == :sinsq
        return flattop_sinsq(t, t₀, T, t_rise, t_fall)
    else
        throw(
            ArgumentError("Unknown func=$func. Accepted values are :blackman and :sinsq.")
        )
    end
end


function flattop_sinsq(t, t₀, T, t_rise, t_fall=t_rise)
    f::Float64 = 0.0
    if t₀ ≤ t ≤ T
        f = 1.0
        if t < t₀ + t_rise
            f = sin(π * (t - t₀) / (2.0 * t_rise))^2
        elseif t > T - t_fall
            f = sin(π * (t - T) / (2.0 * t_fall))^2
        end
    end
    return f
end


function flattop_blackman(t, t₀, T, t_rise, t_fall=t_rise)
    f::Float64 = 0.0
    if t₀ ≤ t ≤ T
        f = 1.0
        if t < t₀ + t_rise
            f = blackman(t, t₀, t₀ + 2 * t_rise)
        elseif t > T - t_fall
            f = blackman(t, T - 2 * t_fall, T)
        end
    end
    return f
end


@doc raw"""Box shape (Theta-function).

```julia
box(t, t₀, T)
```

evaluates the Heaviside (Theta-) function $\Theta(t) = 1$ for $t_0 \le t \le
T$; and $\Theta(t) = 0$ otherwise.
"""
box(t, t₀, T) = (t₀ ≤ t ≤ T) ? 1.0 : 0.0


@doc raw"""Blackman window shape.

```julia
blackman(t, t₀, T; a=0.16)
```

calculates

```math
B(t; t_0, T) =
    \frac{1}{2}\left(
        1 - a - \cos\left(2π \frac{t - t_0}{T - t_0}\right)
        + a \cos\left(4π \frac{t - t_0}{T - t_0}\right)
    \right)\,,
```

for a scalar `t`, with $a$ = 0.16.

See <http://en.wikipedia.org/wiki/Window_function#Blackman_windows>

A Blackman shape looks nearly identical to a Gaussian with a 6-sigma
interval between `t₀` and `T`.  Unlike the Gaussian, however, it will go
exactly to zero at the edges. Thus, Blackman pulses are often preferable to
Gaussians.
"""
function blackman(t, t₀, T; a=0.16)
    ΔT = T - t₀
    return (
        0.5 *
        box(t, t₀, T) *
        (1.0 - a - cos(2π * (t - t₀) / ΔT) + a * cos(4π * (t - t₀) / ΔT))
    )
end


end
