module QuantumPropagatorsExponentialUtilitiesExt

using LinearAlgebra
using TimerOutputs: @timeit_debug, TimerOutput
using ExponentialUtilities

using QuantumPropagators:
    QuantumPropagators,
    ExponentialUtilitiesPropagator,
    _pwc_get_max_genop,
    _pwc_get_genop,
    _pwc_set_genop!,
    _pwc_advance_time!,
    _pwc_process_parameters
using QuantumPropagators.Controls: get_controls
using QuantumPropagators.Interfaces: supports_inplace
import QuantumPropagators: init_prop, prop_step!


"""
```julia
using ExponentialUtilities

expv_propagator = init_prop(
    state,
    generator,
    tlist;
    method=ExponentialUtilities,
    inplace=QuantumPropagators.Interfaces.supports_inplace(state),
    backward=false,
    verbose=false,
    parameters=nothing,
    expv_kwargs=(;),
    _...
)
```

initializes an [`ExponentialUtilitiesPropagator`](@ref).

# Method-specific keyword arguments

* `expv_kwargs`: NamedTuple of keyword arguments forwarded to the underlying
  [`ExponentialUtilities`](https://docs.sciml.ai/ExponentialUtilities/stable/)
  routines. For in-place propagation, these are passed to both `arnoldi!`
  (Krylov subspace construction) and `expv!` (exponentiation). For
  out-of-place propagation, they are passed to `expv`, which forwards them
  internally to `arnoldi`. The most useful keyword arguments are:

  * `m::Int`: Maximum Krylov subspace dimension. Defaults to
    `min(30, size(A, 1))`. Larger values improve accuracy but increase memory
    and computation cost per step.
  * `ishermitian::Bool`: Whether the operator is Hermitian. Defaults to
    `LinearAlgebra.ishermitian(A)`. When `true`, the cheaper Lanczos iteration
    is used instead of the full Arnoldi process. Note that the propagator
    passes `-𝕚 dt` as the time argument, so a Hermitian generator ``Ĥ`` still
    results in a non-Hermitian product ``-𝕚 Ĥ dt``; however,
    `ExponentialUtilities` handles the complex time scaling internally.
  * `tol::Real`: Tolerance for happy breakdown. Defaults to `1e-7`. The
    Arnoldi iteration stops early when the norm of the next Krylov vector
    drops below `tol * opnorm`.
  * `opnorm`: Operator norm of `A`. By default, this is computed
    automatically. Supplying it avoids redundant norm computations.
  * `iop::Int`: Incomplete orthogonalization procedure length. Defaults to `0`
    (full orthogonalization). A positive value limits the number of previous
    Krylov vectors used for orthogonalization, reducing cost at the expense of
    numerical stability.
  * `mode::Symbol`: Termination strategy for `expv`. Either
    `:happy_breakdown` (default) or `:error_estimate`. In `:happy_breakdown`
    mode, the iteration relies on early termination when the Krylov subspace
    captures the relevant dynamics. In `:error_estimate` mode, an adaptive
    step-size strategy is used with additional tolerance parameters `rtol`
    (relative tolerance, defaults to `√tol`).
"""
function init_prop(
    state,
    generator,
    tlist,
    method::Val{:ExponentialUtilities};
    inplace = supports_inplace(state),
    backward = false,
    verbose = false,
    parameters = nothing,
    expv_kwargs = (;),
    _...
)
    tlist = convert(Vector{Float64}, tlist)
    controls = get_controls(generator)
    G = _pwc_get_max_genop(generator, controls, tlist)

    parameters = _pwc_process_parameters(parameters, controls, tlist)
    timing_data = TimerOutput()

    n = 1
    t = tlist[1]
    if backward
        n = length(tlist) - 1
        t = float(tlist[n+1])
    end

    if inplace
        A₀ = QuantumPropagators.Controls.evaluate(generator, tlist, n)
        is_herm = get(expv_kwargs, :ishermitian, ishermitian(A₀))
        if haskey(expv_kwargs, :mode) && expv_kwargs[:mode] == :error_estimate && !is_herm
            throw(
                ArgumentError(
                    "`mode=:error_estimate` in `expv_kwargs` requires a " *
                    "Hermitian generator. Use the default " *
                    "`mode=:happy_breakdown` for non-Hermitian generators."
                )
            )
        end
        Ks = ExponentialUtilities.arnoldi(A₀, state; expv_kwargs...)
        if haskey(expv_kwargs, :mode) && expv_kwargs[:mode] == :error_estimate
            cache = ExponentialUtilities.get_subspace_cache(Ks)
        else
            T = promote_type(eltype(A₀), eltype(state))
            cache_type = is_herm ? real(T) : T
            cache = ExponentialUtilities.ExpvCache{cache_type}(Ks.m)
        end
    else
        Ks = nothing
        cache = nothing
    end

    GT = typeof(generator)
    OT = typeof(G)
    ST = typeof(state)
    KST = typeof(Ks)
    CT = typeof(cache)

    return ExponentialUtilitiesPropagator{GT,OT,ST,KST,CT}(
        generator,
        inplace ? copy(state) : state,
        t,
        n,
        tlist,
        parameters,
        controls,
        G,
        Ks,
        cache,
        backward,
        inplace,
        expv_kwargs,
        timing_data,
    )
end


function prop_step!(propagator::ExponentialUtilitiesPropagator)
    @timeit_debug propagator.timing_data "prop_step!" begin
        n = propagator.n
        tlist = getfield(propagator, :tlist)
        (0 < n < length(tlist)) || return nothing
        dt = tlist[n+1] - tlist[n]
        if propagator.backward
            dt = -dt
        end

        if propagator.inplace
            if supports_inplace(propagator.genop)
                _pwc_set_genop!(propagator, n)
                H = propagator.genop
            else
                H = _pwc_get_genop(propagator, n)
            end
            @timeit_debug propagator.timing_data "expv" begin
                if haskey(propagator.expv_kwargs, :mode) &&
                   propagator.expv_kwargs[:mode] == :error_estimate
                    ExponentialUtilities.expv!(
                        propagator.state,
                        -1im * dt,
                        H,
                        propagator.state,
                        propagator.Ks,
                        propagator.cache;
                        propagator.expv_kwargs...,
                    )
                else
                    ExponentialUtilities.arnoldi!(
                        propagator.Ks,
                        H,
                        propagator.state;
                        propagator.expv_kwargs...,
                    )
                    ExponentialUtilities.expv!(
                        propagator.state,
                        -1im * dt,
                        propagator.Ks;
                        cache = propagator.cache,
                    )
                end
            end
        else
            H = _pwc_get_genop(propagator, n)
            @timeit_debug propagator.timing_data "expv" begin
                Ψ = ExponentialUtilities.expv(
                    -1im * dt,
                    H,
                    propagator.state;
                    propagator.expv_kwargs...
                )
            end
            setfield!(propagator, :state, convert(typeof(propagator.state), Ψ))
        end

        _pwc_advance_time!(propagator)
        return propagator.state
    end
end

end
