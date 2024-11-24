module SpectralRange

export specrange

using ..Arnoldi: arnoldi!, diagonalize_hessenberg_matrix, extend_arnoldi!
using LinearAlgebra
using Random: GLOBAL_RNG


"""Calculate the spectral range of a Hamiltonian `H` on the real axis.

```julia
E_min, E_max = specrange(H; method=:auto, kwargs...)
```

calculates the approximate lowest and highest eigenvalues of `H`. Any imaginary
part in the eigenvalues is ignored: the routine is intended for (although not
strictly limited to) a Hermitian `H`.

This delegates to

```julia
specrange(H, method; kwargs...)
```

for the different methods.

The default `method=:auto` chooses the best method for the given `H`. This is
`:diag` for small matrices, and `:arnoldi` otherwise. If both `E_min` and
`E_max` are given in the `kwargs`, those will be returned directly
(`method=:manual`).

Keyword arguments not relevant to the underlying implementation will be
ignored.
"""
function specrange(H; method=Val(:auto), kwargs...)
    return specrange(H, method; kwargs...)
end

function specrange(H, method::Symbol; kwargs...)
    return specrange(H, Val(method); kwargs...)
end


function specrange(H, method::Val{:auto}; kwargs...)
    if haskey(kwargs, :E_min) && haskey(kwargs, :E_max)
        return specrange(H, Val(:manual); kwargs...)
    end
    if hasmethod(size, (typeof(H),)) && hasmethod(Array, (typeof(H),))
        if size(H)[1] <= 32
            # TODO: benchmark what a good cross-over point from exact
            # diagonalization to Arnoldi is
            return specrange(H, Val(:diag); kwargs...)
        end
    else
        @warn "`size(H)` and `Array(H)` are not implemented for $(typeof(H))"
    end
    return specrange(H, Val(:arnoldi); kwargs...)
end


"""
```julia
E_min, E_max = specrange(
    H, :arnoldi;
    rng=Random.GLOBAL_RNG,
    state=random_state(H; rng),
    m_min=20,
    m_max=60,
    prec=1e-3,
    norm_min=1e-15,
    enlarge=true
)
```

uses [Arnoldi iteration](https://en.wikipedia.org/wiki/Arnoldi_iteration) with
`state` as the starting vector. It approximates the eigenvalues of `H` with
between `m_min` and `m_max` Ritz values, until the lowest and highest
eigenvalue are stable to a relative precision of `prec`. The `norm_min`
parameter is passed to the underlying [`arnoldi!`](@ref).

If `enlarge=true` (default) the returned `E_min` and `E_max` will be enlarged
via a heuristic to slightly over-estimate the spectral radius instead of
under-estimating it.
"""
function specrange(
    H,
    method::Val{:arnoldi};
    rng=GLOBAL_RNG,
    state=random_state(H; rng),
    kwargs...
)

    m_max = get(kwargs, :m_max, 60)
    m_min = max(5, min(get(kwargs, :m_min, 25), m_max - 1))
    prec = get(kwargs, :prec, 1e-3)
    norm_min = get(kwargs, :norm_min, 1e-15)
    enlarge = get(kwargs, :enlarge, true)

    R = ritzvals(H, state, m_min, m_max; prec=prec, norm_min=norm_min)
    E_min = real(R[1])
    E_max = real(R[end])
    if enlarge && (length(R) > 1)
        # We want to overestimate the spec.rad., not underestimate it, so we
        # use the distance to the next eigenvalues as a buffer
        E_min = 2 * E_min - real(R[2])
        E_max = 2 * E_max - real(R[end-1])
    end
    return E_min::Float64, E_max::Float64
end


"""
```julia
E_min, E_max = specrange(H, :diag)
```

uses exact diagonalization via the standard `eigvals` function to obtain the
smallest and largest eigenvalue. This should only be used for relatively small
matrices.
"""
function specrange(H, method::Val{:diag}; kwargs...)
    evals = real.(eigvals(Array(H)))
    # "By default, the eigenvalues [...] are sorted [...] by (real(λ),imag(λ))"
    return evals[1]::Float64, evals[end]::Float64
end


"""
```julia
E_min, E_max = specrange(H, :manual; E_min, E_max)
```

directly returns the given `E_min` and `E_max` without considering `H`.
"""
function specrange(H, method::Val{:manual}; E_min, E_max, _...)
    return float(E_min), float(E_max)
end


"""Random normalized quantum state.

```julia
    Ψ = random_state(H; rng=Random.GLOBAL_RNG)
```

returns a random normalized state compatible with the Hamiltonian `H`. This is
intended to provide a starting vector for estimating the spectral radius of `H`
via an Arnoldi method.
"""
function random_state(H; rng=GLOBAL_RNG)
    N = size(H)[2]
    Ψ = rand(rng, N) .* exp.((2π * im) .* rand(rng, N))
    Ψ ./= norm(Ψ)
    return Ψ
end


"""Calculate a vector for Ritz values converged to a given precision.

```julia
R = ritzvals(G, state, m_min, m_max=2*m_min; prec=1e-5, norm_min=1e-15)
```

calculates a complex vector `R` of at least `m_min` (assuming a sufficient
Krylov dimension) and at most `m_max` Ritz values.
"""
function ritzvals(G, state, m_min, m_max=2 * m_min; prec=1e-5, norm_min=1e-15)
    if m_max <= m_min
        throw(ArgumentError("m_max=$m_max must be smaller than m_min=$min"))
    end
    m = max(5, min(m_min, m_max - 1))

    Hess = Matrix{ComplexF64}(undef, m_max, m_max)
    q = [similar(state) for _ = 1:(m_max+1)]

    # Get the Ritz eigenvalues for order m-1, so that we can decide whether
    # the values for order m are already precise enough
    m₀ = m - 1
    m₀ = arnoldi!(Hess, q, m₀, state, G; extended=false, norm_min=norm_min)
    eigenvals = diagonalize_hessenberg_matrix(Hess, m₀)
    v̲r₀ = minimum(real(eigenvals))
    v̄r₀ = maximum(real(eigenvals))
    v̄i₀ = maximum(abs.(imag(eigenvals)))
    if m₀ == m - 1
        # Starting from original order m, increase `m` until the Ritz values
        # reach the desired precision
        extend_arnoldi!(Hess, q, m, G; norm_min=norm_min)
        eigenvals = diagonalize_hessenberg_matrix(Hess, m)
        v̲r = minimum(real(eigenvals))
        v̄r = maximum(real(eigenvals))
        v̄i = maximum(abs.(imag(eigenvals)))
        e̲r = (v̲r₀ ≠ 0.0) ? abs(1.0 - v̲r / v̲r₀) : 0.0
        ēr = (v̄r₀ ≠ 0.0) ? abs(1.0 - v̄r / v̄r₀) : 0.0
        ēi = (v̄i₀ ≠ 0.0) ? abs(1.0 - v̄i / v̄i₀) : 0.0
        while (e̲r > prec) || (ēr > prec) || ((v̄i₀ > 1e-14) && ēi > prec)
            v̲r₀ = v̲r
            v̄r₀ = v̄r
            v̄i₀ = v̄i
            m₀ = m
            m = m + 1
            extend_arnoldi!(Hess, q, m, G; norm_min=norm_min)
            (m == m₀) && break  # dimensionality exhausted
            eigenvals = diagonalize_hessenberg_matrix(Hess, m)
            v̲r = minimum(real(eigenvals))
            v̄r = maximum(real(eigenvals))
            v̄i = maximum(abs.(imag(eigenvals)))
            e̲r = abs(1.0 - (v̲r / v̲r₀))
            ēr = abs(1.0 - (v̄r / v̄r₀))
            ēi = abs(1.0 - (v̄i / v̄i₀))
            if (m == m_max)
                #@warn "Ritz values did not converge within m_max=$m_max"
                break
            end
        end
    end
    return eigenvals
end

end
