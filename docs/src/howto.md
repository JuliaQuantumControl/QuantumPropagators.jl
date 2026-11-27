```@meta
# SPDX-FileCopyrightText: © 2021 Michael Goerz <mail@michaelgoerz.net>
# SPDX-License-Identifier: MIT OR CC-BY-4.0
```

# [How-to Guides](@id howto)

```math
\gdef\op#1{\hat{#1}}
\gdef\Liouvillian{\mathcal{L}}
\gdef\Re{\operatorname{Re}}
\gdef\Im{\operatorname{Im}}
```

**Contents**

```@contents
Pages=[
 "howto.md"
]
Depth = 2:2
```

## How to implement a new propagation method

* Define a new sub-type of [`AbstractPropagator`](@ref QuantumPropagators.AbstractPropagator) type that is unique to the propagation method, e.g. `MyNewMethodPropagator`. If appropriate, sub-type [`PiecewisePropagator`](@ref QuantumPropagators.PiecewisePropagator) or [`PWCPropagator`](@ref QuantumPropagators.PWCPropagator).

* The high-level [`propagate`](@ref) and [`init_prop`](@ref) functions have a mandatory `method` keyword argument. That argument should receive a `Module` object for the `module` implementing a propagation method (e.g., `using QuantumPropagators: Cheby; method=Cheby`). This ensures that the module or package implementing the method is fully loaded. Internally, [`init_prop`](@ref) delegates `method=module` to a *positional* argument `method = Val(nameof(module))`

* Thus, if `MyNewMethodPropagator` is implemented in a module `MyNewMethod`, or if it wraps a registered package `MyNewMethod`, implement a new [`init_prop`](@ref) method with the following signature:

  ```
  function init_prop(
      state,
      generator,
      tlist,
      method::Val{:MyNewMethod};
      inplace=true,
      backward=false,
      verbose=false,
      parameters=nothing,
      # ... method-specific keyword arguments
      _...  # ignore other keyword arguments
  )
  ```

  Note the `method::Val{:MyNewMethod}` as the fourth positional (!) parameter. While the *public* interface for [`init_prop`](@ref) takes `method` as a keyword argument, privately [`init_prop`](@ref) dispatches for different methods as above.

* If the propagation method is not associated with a module or package, it is also possible to implement a method [`init_prop`](@ref) with a fourth positional argument with a type of, e.g., `::Val{:my_new_method}`. This would allow to call the high-level [`propagate`](@ref)/[`init_prop`](@ref) with the keyword argument `method=:my_new_method`, i.e., passing a name (`Symbol`) instead of a `Module`.

* Implement the remaining methods in [The Propagator interface](@ref).

* Test the implementation by instantiating a `propagator` and calling [`QuantumPropagators.Interfaces.check_propagator`](@ref) on it.


## How to specify the spectral range for a Chebychev propagation

A propagation with [`method=Cheby`](@ref method_cheby) requires that the dynamic generator ``\op{H}(t)`` be normalized to a spectral range of [-1, 1]. That is, the method needs a (pessimistic) estimate of the "spectral envelope": the minimum and maximum eigenvalue of ``\op{H}(t)`` for any point `t` on the interval of the propagation time grid `tlist`.

By default, the Chebychev propagator uses heuristics to estimate  this spectral envelope. If the spectral envelope is known (either analytically of via a separate numerical exploration of the eigenvalues over the full range of possible controls), the minimum and maximum eigenvalues of ``\op{H}(t)`` can be passed as keyword arguments `E_min` and `E_max` to [`propagate`](@ref) or [`init_prop`](@ref). Since the Chebychev method is only defined for Hermitian generators, `E_min` and `E_max` must be real values. Both values must be given.

Manually specifying `E_min` and `E_max` works with the default `specrange_method=:auto` as well as with the explicit `specrange_method=:manual`. When calculating the Chebychev coefficients, the given values may still be enlarged by the default `specrange_buffer` keyword argument in [`init_prop`](@ref). If `E_min` and `E_max` should be used *exactly*, pass `specrange_buffer=0`.


## [How to define a parameterized control](@id howto_parameterized)

A parameterized control is a [function-like object](@extref Julia Function-like-objects) of `t` that holds a [parameter block](@ref ParameterBlocks), a vector of its parameter values. It is best defined as a subtype of [`ParameterizedFunction`](@ref QuantumPropagators.Controls.ParameterizedFunction), with the block in a field `parameters`. A [`ComponentVector`](https://github.com/SciML/ComponentArrays.jl) makes a convenient block, since it gives names to the parameters. For example,

```@meta
DocTestSetup = quote
    # Ensure there's not pre-compilation from package loading during the actual doctest
    using ComponentArrays
    using UnPack: @unpack
end
```

```jldoctest gaussian
using ComponentArrays
using UnPack: @unpack
using QuantumPropagators.Controls: ParameterizedFunction

struct GaussianControl <: ParameterizedFunction
    parameters::ComponentVector{Float64,Vector{Float64},Tuple{Axis{(A=1, t0=2, sigma=3)}}}
end

function GaussianControl(; A=1.0, t0=0.0, sigma=1.0)
    return GaussianControl(ComponentVector(; A, t0, sigma))
end

function (control::GaussianControl)(t)
    @unpack A, t0, sigma = control.parameters
    return A * exp(- (t - t0)^2 / (2 * sigma^2))
end

# usage

gaussian = GaussianControl(A=2.0, sigma=0.5)
gaussian.parameters.t0 = 5  # shift center from original 0.0

@show gaussian.parameters
println("gaussian(4.5) = $(round(gaussian(4.5); digits=3))")

# output

gaussian.parameters = (A = 2.0, t0 = 5.0, sigma = 0.5)
gaussian(4.5) = 1.213
```

The control must read its parameters from the block every time it is called, so that writing to the block changes the control. Like for any other object, [`get_parameters(gaussian)`](@ref QuantumPropagators.Controls.get_parameters) returns an `ArrayPartition` that wraps the block, see [How to work with control parameters](@ref howto_parameters).

With some extra effort, a control can give direct property access to its parameters, including Unicode aliases. Since this overrides `getproperty`, the control must access its own block with `getfield`:

```jldoctest
using ComponentArrays
using QuantumPropagators.Controls: ParameterizedFunction

struct GaussianControl <: ParameterizedFunction
    parameters::ComponentVector{Float64,Vector{Float64},Tuple{Axis{(A=1, t0=2, sigma=3)}}}
end

function GaussianControl(; A=1.0, t0=0.0, t₀=t0, sigma=1.0, σ=sigma)
    return GaussianControl(ComponentVector(; A, t0=t₀, sigma=σ))
end

function Base.propertynames(g::GaussianControl, private::Bool=false)
    names = (:A, :t0, :t₀, :sigma, :σ)
    return private ? Tuple(union(names, fieldnames(GaussianControl))) : names
end

function Base.getproperty(g::GaussianControl, name::Symbol)
    unicode_aliases = Dict(:σ => :sigma, :t₀ => :t0)
    getproperty(getfield(g, :parameters), get(unicode_aliases, name, name))
end

function Base.setproperty!(g::GaussianControl, name::Symbol, value)
    unicode_aliases = Dict(:σ => :sigma, :t₀ => :t0)
    setproperty!(
        getfield(g, :parameters),
        get(unicode_aliases, name, name),
        value
    )
end

function (control::GaussianControl)(t)
    A, t₀, σ = getfield(control, :parameters)
    return A * exp(- (t - t₀)^2 / (2 * σ^2))
end

# usage

gaussian = GaussianControl(A=2.0, σ=0.5)
gaussian.t₀ = 5  # shift center from original 0.0

@show gaussian.t0
println("gaussian(4.5) = $(round(gaussian(4.5); digits=3))")

# output

gaussian.t0 = 5.0
gaussian(4.5) = 1.213
```

Use [`QuantumPropagators.Interfaces.check_parameterized_function`](@ref) to verify the implementation of a [`ParameterizedFunction`](@ref QuantumPropagators.Controls.ParameterizedFunction).


## [How to work with control parameters](@id howto_parameters)

A method that tunes control parameters works with the parameters of a generator (or of a control problem in QuantumControl), as obtained from [`get_parameters`](@ref QuantumPropagators.Controls.get_parameters). For a generator with two of the Gaussian controls [defined above](@ref howto_parameterized):

```jldoctest gaussian
julia> using QuantumPropagators: hamiltonian

julia> using QuantumPropagators.Controls: get_parameters

julia> ϵ₁ = GaussianControl(A=2.0, t0=1.0, sigma=0.5);

julia> ϵ₂ = GaussianControl(A=1.0, t0=3.0, sigma=0.5);

julia> H = hamiltonian([1.0 0.0; 0.0 -1.0], ([0.0 1.0; 1.0 0.0], ϵ₁), ([0.0 -1im; 1im 0.0], ϵ₂));

julia> u = get_parameters(H)
((A = 2.0, t0 = 1.0, sigma = 0.5), (A = 1.0, t0 = 3.0, sigma = 0.5))

julia> length(u)
6
```

The result `u` is a [`RecursiveArrayTools.ArrayPartition`](@extref) of the parameter blocks of `ϵ₁` and `ϵ₂`, shown as a tuple of the blocks, in order: `u[1:3]` are the parameters of `ϵ₁`, and `u[4:6]` are those of `ϵ₂`. It stores no values of its own. Writing to `u` writes into the blocks:

```jldoctest gaussian
julia> u[4] = 1.5;

julia> ϵ₂.parameters.A
1.5
```

The names of the parameters are available on the blocks, not on `u`. Access a block through the control that holds it, or through `u.x`, the tuple of blocks:

```jldoctest gaussian
julia> u.x[2].t0
3.0

julia> u.x[1] === ϵ₁.parameters
true
```

An optimizer that works with a `Vector{Float64}` reads the parameters with `collect`, and writes new values with `copyto!`:

```jldoctest gaussian
julia> x = collect(u);

julia> x[1] = 2.5;

julia> copyto!(u, x);

julia> ϵ₁.parameters.A
2.5
```

As long as `H` contains the same controls, `get_parameters(H)` returns the identical object:

```jldoctest gaussian
julia> get_parameters(H) === u
true
```

Still, each call goes through all the controls in `H`, so call `get_parameters` once, when setting up the optimization.


## [How to implement a generator with its own parameters](@id howto_custom_parameters)

The built-in [`Generator`](@ref QuantumPropagators.Generators.Generator) holds no [parameter blocks](@ref ParameterBlocks) of its own. It only inherits the blocks of its controls. A custom generator type can also hold blocks in its own fields. For example, a transmon with a tunable frequency ``ν`` and anharmonicity ``α``, driven by the [Gaussian control defined above](@ref howto_parameterized),

```math
\op{H}(t) = ν \op{b}^† \op{b} + \frac{α}{2} \op{b}^† \op{b}^† \op{b} \op{b} + ϵ(t) (\op{b} + \op{b}^†),
```

can hold ``ν`` and ``α`` in a block `w`:

```jldoctest gaussian
using LinearAlgebra: diagm
using QuantumPropagators.Generators: Operator
using QuantumPropagators.Controls: combine_parameters
import QuantumPropagators.Controls:
    get_controls, get_parameters, evaluate, evaluate!, substitute

struct TransmonGenerator{BT,OT,CT}
    w::BT  # parameter block (ν, α)
    ops::Vector{OT}  # b̂⁺b̂, b̂⁺b̂⁺b̂b̂/2, b̂ + b̂⁺
    control::CT
end

function TransmonGenerator(w, control; N=5)
    b̂ = diagm(1 => sqrt.(1.0:(N-1)))
    ops = [b̂' * b̂, b̂' * b̂' * b̂ * b̂ / 2, b̂ + b̂']
    return TransmonGenerator(w, ops, control)
end

get_controls(G::TransmonGenerator) = (G.control,)

get_parameters(G::TransmonGenerator) = combine_parameters(G.w, get_parameters(G.control))

function evaluate(G::TransmonGenerator, args...; vals_dict=IdDict())
    ν, α = G.w
    return Operator(G.ops, [ν, α, evaluate(G.control, args...; vals_dict)])
end

function evaluate!(op::Operator, G::TransmonGenerator, args...; vals_dict=IdDict())
    ν, α = G.w
    op.coeffs .= (ν, α, evaluate(G.control, args...; vals_dict))
    return op
end

function substitute(G::TransmonGenerator, replacements)
    return TransmonGenerator(G.w, G.ops, substitute(G.control, replacements))
end

# usage

w = ComponentVector(ν=5.0, α=-0.3)
ϵ = GaussianControl(A=0.1, t0=5.0, sigma=1.0)
G = TransmonGenerator(w, ϵ)
get_parameters(G)

# output

((ν = 5.0, α = -0.3), (A = 0.1, t0 = 5.0, sigma = 1.0))
```

Like any generator, the type implements [`get_controls`](@ref QuantumPropagators.Controls.get_controls), [`evaluate`](@ref QuantumPropagators.Controls.evaluate), [`evaluate!`](@ref QuantumPropagators.Controls.evaluate!), and [`substitute`](@ref QuantumPropagators.Controls.substitute). In addition, it implements its own [`get_parameters`](@ref QuantumPropagators.Controls.get_parameters) method. The generic method finds only the blocks of the controls in `get_controls(G)`, and would miss `w`.

The function [`combine_parameters`](@ref QuantumPropagators.Controls.combine_parameters) (public, but not exported) combines `w` with the parameters of the control. It returns a flat `ArrayPartition` of the blocks and drops any block that occurs more than once, so that `get_parameters(G)` follows the same rules as for any other object. The order of the blocks is up to the method. This one puts its own block first:

```jldoctest gaussian
julia> u = get_parameters(G);

julia> u.x[1] === w
true

julia> u.x[2] === ϵ.parameters
true
```

With a control that holds no block, such as a plain function, the parameters contain only `w`:

```jldoctest gaussian
julia> get_parameters(TransmonGenerator(w, t -> 0.1))
((ν = 5.0, α = -0.3),)
```

Verify the implementation with [`check_generator`](@ref QuantumPropagators.Interfaces.check_generator) and `for_parameterization=true`:

```jldoctest gaussian
julia> using QuantumPropagators.Interfaces: check_generator

julia> tlist = collect(range(0, 10; length=101));

julia> Ψ = ComplexF64[1, 0, 0, 0, 0];

julia> check_generator(G; state=Ψ, tlist, for_parameterization=true)
true
```

For a gradient-based optimization of the parameters, the type also needs a method for [`get_parameter_deriv`](@extref `QuantumControl.Controls.get_parameter_deriv`) in QuantumControl.
