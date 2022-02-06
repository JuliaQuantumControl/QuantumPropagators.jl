# Implementation of Newton-with-restarted-Arnoldi propagation routine
module Newton

export NewtonWrk, newton!

using ..Arnoldi: arnoldi!, diagonalize_hessenberg_matrix
using LinearAlgebra
using OffsetArrays


"""
```julia
NewtonWrk(v0, m_max=10)
```

Workspace for the Newton-with-restarted-Arnoldi propagation routine.

Initializes the workspace for the propagation of a vector `v0`, using a
maximum Krylov dimension of `m_max` in each restart iteration. Note that
`m_max` should be smaller than the length of `v0`.
"""
mutable struct NewtonWrk{T}
    # arnoldi_vecs and v are pre-allocated for efficiency. Other values and
    # vectors are mainly for debugging, so that we can inspect internal
    # parameters like the Newton coefficients (a) after a call to newton!
    arnoldi_vecs::Array{T}
    v::T
    a::OffsetVector{ComplexF64}
    leja::OffsetVector{ComplexF64}
    radius::Float64
    n_a::Int64
    n_leja::Int64
    restarts::Int64
    function NewtonWrk(v0::T; m_max::Int64=10) where {T}
        if m_max >= length(v0)
            m_max = length(v0) - 1
        end
        new{T}(
            T[similar(v0) for _ = 1:m_max+1], # arnoldi_vecs
            similar(v0),                       # v
            OffsetVector(zeros(ComplexF64, 10 * m_max + 1), 0:10*m_max),  # a
            OffsetVector(zeros(ComplexF64, 10 * m_max + 1), 0:10*m_max),  # leja
            0.0,                               # radius
            0,                                 # n_a
            0,                                 # n_leja
            0,                                 # restarts
        )
    end
end


lbound(array::OffsetArray, dim=1) = first(axes(array)[dim])
ubound(array::OffsetArray, dim=1) = last(axes(array)[dim])


function leja_radius(z)
    r_max = maximum(map(abs, z))
    return 1.2 * r_max
end


"""
```julia
extend_leja!(leja, n, newpoints, n_use)
```

Given an array of `n` (ordered) Leja points, extract `n_use`
points from `newpoints`, and append them to the existing Leja points. The
array `leja` should be sufficiently large to hold the new Leja points,
which are appended after index `n_old`. It will be re-allocated if
necessary and may have a size of up to `2*(n+n_use)`.

# Arguments

- `leja`: Array of leja values. Must contain the "old" leja values to be kept
   in `leja(0:n-1)`. On output, `n_use` new leja points will be in
   `leja(n+:n+n_use-1)`, for the original value of `n`.  The `leja` array must
   use zero-based indexing.
- `n`: On input, number of "old" leja points in `leja`. On output, total number
  of leja points (i.e. `n=n+n_use`)
- `newpoints`: On input, candidate points for new leja points.  The `n_use`
  best values will be chosen and added to `leja`. On output, the values of
  `new_points` are undefined.
- `n_use`: Number of points that should be added to `leja`
"""
function extend_leja!(
    leja::OffsetVector{ComplexF64},
    n,
    newpoints::OffsetVector{ComplexF64},
    n_use
)
    @assert leja.offsets[1] == -1 "leja must be a zero-based vector"
    @assert newpoints.offsets[1] == -1 "newpoints must be a zero-based vector"
    if length(leja) < n + n_use
        # we allocate twice the required size, so we don't have to re-allocate
        # on every call
        resize!(leja, 2 * (n + n_use))
        leja[n:end] .= 0
    end
    u = ubound(newpoints, 1)
    i_add_start = lbound(leja, 1)
    if n == 0
        # Use point of largest absolute value as starting point by moving it to
        # the end of the `newpoints` array
        z_last::ComplexF64 = newpoints[u]
        for i = lbound(newpoints, 1):(u-1)
            if abs(newpoints[i]) > abs(z_last)
                newpoints[u] = newpoints[i]
                newpoints[i] = z_last
                z_last = newpoints[u]
            end
        end
        leja[0] = newpoints[end]
        i_add_start = 1
    end
    exponent = 1.0 / (n + n_use)
    for i_add = i_add_start:n_use-1
        p_max::Float64 = 0
        i_max = 0
        for i = lbound(newpoints, 1):(u-i_add)
            p::Float64 = 1
            for j = 0:(n+i_add-1)  # existing Leja points
                p = p * (abs(newpoints[i] - leja[j])^exponent)
            end
            if p > p_max
                p_max = p
                i_max = i
            end
        end
        # TODO: check that p_max did not overflow or underflow
        leja[n+i_add] = newpoints[i_max]
        # remove the used point by replacing it with the last unused point
        newpoints[i_max] = newpoints[u-i_add]
    end
    return n + n_use
end


"""
```julia
extend_newton_coeffs!(a, n_a, leja, func, n_leja, radius)
```

Extend the array `a` of existing Newton coefficients for the expansion of
the `func` from `n_a` coefficients to `n_leja` coefficients. Return a new value
`n_a=n_a+n_leja` with the total number of Newton coefficients in the updated
`a`.

# Arguments

- `a`: On input, a zero-based array of length `n_a` or greater, containing
  Newton coefficients. On output, array containing a total `n_leja`
  coefficients. The array `a` will be resized if necessary, and may have a
  length greater than `n_leja` on output
- `n_a`:  The number of Newton coefficients in `a`, on input. Elements of `a`
   beyond the first `n_a` elements will be overwritten.
- `leja`: Array of normalized Leja points, containing at least `n_leja`
  elements.
- `func`: Function for which to calcluate Newton coeffiecients
- `n_leja`: The number of elements in `leja` to use for calculating new
  coefficients, and the total number of Newton coefficients on output
- `radius`: Normalization radius for divided differences
"""
function extend_newton_coeffs!(
    a::OffsetVector{ComplexF64},
    n_a::Int64,
    leja::OffsetVector{ComplexF64},
    func,
    n_leja::Int64,
    radius::Float64
)
    @assert a.offsets[1] == -1 "a must be a zero-based vector"
    m = n_leja - n_a
    n0 = n_a
    if length(a) < n_a + m
        # we allocate twice the required size, so we don't have to re-allocate
        # on every call
        resize!(a, 2 * n_leja)
        a[n_a:end] .= 0
    end
    @assert length(leja) >= n_leja
    @assert radius > 0
    if (n_a == 0)
        a[0] = func(leja[0])
        n0 = 1
    end
    for k = n0:(n_a-1+m)
        d::ComplexF64 = 1
        pn::ComplexF64 = 0
        for n = 1:k-1
            zd = leja[k] - leja[n-1]
            d = d * zd / radius
            pn = pn + a[n] * d
        end
        zd = leja[k] - leja[k-1]
        d = d * zd / radius
        @assert abs(d) > 1e-200 "Divided differences too small"
        a[k] = (func(leja[k]) - a[0] - pn) / d
    end
    n_a = n_a + m
    return n_a
end


"""
```julia
newton!(Ψ, H, dt, wrk; func=(z -> exp(-1im*z)), norm_min=1e-14, relerr=1e-12,
        max_restarts=50)
```

Evaluate `Ψ = func(H*dt) Ψ` using a Newton-with-restarted-Arnoldi scheme.

# Arguments

- `Ψ`: The state to propagate, will be overwritten in-place with the propagated
  state
- `H`: Operator acting on `Ψ`. Together with `dt`, this is the argument to
  `func`
- `dt`: Implicit time step. Together with `H`, this is the argument to `func`
- `wkr`: Work array, initialized with [`NewtonWrk`](@ref)
- `func`: The function to apply to `H dt`, taking a single (scalar)
  complex-valued argument `z` in place of `H dt`. The default `func`
  is to evaluate the time evoluation operator for the Schrödinger equation
- `norm_min`: the minium norm at which to consider a state similar to `Ψ` as
  zero
- `relerr`: The relative error defining the convergence condition for the
  restart iteration. Propagation stops when the norm of the accumulated `Ψ`
  is stable up to the given relative error
- `max_restart`: The maximum number of restart iterations. Exceeding
  `max_restart` will throw an `AssertionError`.
"""
function newton!(Ψ, H, dt, wrk; kwargs...)
    func = get(kwargs, :func, z -> exp(-1im * z))
    norm_min::Float64 = get(kwargs, :norm_min, 1e-14)
    relerr::Float64 = get(kwargs, :relerr, 1e-12)
    max_restarts::Int64 = get(kwargs, :max_restarts, 50)

    m_max::Int64 = length(wrk.arnoldi_vecs) - 1
    m::Int64 = m_max
    fill!(wrk.a, 0)
    fill!(wrk.leja, 0)
    R = Array{ComplexF64,1}(undef, m + 1)
    P = Array{ComplexF64,1}(undef, m + 1)
    R_abs = Array{Float64,1}(undef, m + 1)
    Hess = zeros(ComplexF64, m_max + 1, m_max + 1)
    # there seems to be a problem declaring dt as Float64 directly (if you pass
    # it an Int, the kwargs seems to prevent Julia from auto-converting it
    _dt::Float64 = dt

    @assert _dt ≠ 0.0

    n_a::Int64 = 0  # number of Newton coefficients (accumulated)
    n_leja::Int64 = 0  # number of Leja points (accumulated)
    copyto!(wrk.v, Ψ)

    s = 0  # counter for the restart loop
    β::Float64 = norm(wrk.v)
    lmul!(1 / β, wrk.v)  # normalize

    while true  # restart loop (s → s+1)

        m = arnoldi!(
            Hess,
            wrk.arnoldi_vecs,
            m,
            wrk.v,
            H,
            _dt;
            extended=true,
            norm_min=norm_min
        )
        if m == 1 && s == 0
            λ = β * Hess[1, 1]
            # wrk.v is an eigenvec of H with eigenval Hess[1, 1],
            # so Ψ is an eigenvec with eigenval λ, and f(H) Ψ = f(λ) Ψ
            lmul!(func(λ), Ψ)
            break
        end
        ritz = diagonalize_hessenberg_matrix(Hess, m, accumulate=true)

        # In the first iteration, the radius will be determined
        if s == 0
            wrk.radius = leja_radius(ritz)
        end

        # Get the Leja points (i.e. Ritz values in the proper order
        n_s = n_leja  # we need to keep the old n_leja
        n_leja = extend_leja!(wrk.leja, n_leja, OffsetVector(ritz, 0:length(ritz)-1), m)

        # Extend the Newton coefficients
        n_a = extend_newton_coeffs!(wrk.a, n_a, wrk.leja, func, n_leja, wrk.radius)
        @assert n_a == n_leja

        # allocate array R for Newton basis polynomicals in the extended
        # Hessenberg matrix and P for the full Newton series
        if length(R) ≠ m + 1
            # we always treat R, R_abs, and P together
            resize!(R, m + 1)
            resize!(R_abs, m + 1)
            resize!(P, m + 1)
        end

        # Evaluate Newton Polynomial in the extended Hessenberg matrix
        P .= 0.0
        R .= 0.0
        R[1] = β
        P[1] = wrk.a[n_s] * β
        for k = 1:m-1
            R[1:m+1] = (
                (Hess[1:m+1, 1:m+1] * R[1:m+1] - wrk.leja[n_s+k-1] * R[1:m+1]) / wrk.radius
            )
            # P += a_{n_s+k} R
            axpy!(wrk.a[n_s+k], view(R, 1:m+1), view(P, 1:m+1))
        end

        # Calculate the new solution
        if s == 0
            fill!(Ψ, 0)
        end
        # Ψ += P[i] * arnoldi_vecs[i]
        for i = 1:m
            axpy!(P[i], wrk.arnoldi_vecs[i], Ψ)
        end

        # Calculate the starting vector v_{s+1} for the next iteration
        R[1:m+1] =
            (((Hess[1:m+1, 1:m+1] * R[1:m+1]) - wrk.leja[n_s+m-1] * R[1:m+1]) / wrk.radius)
        R_abs[1:m+1] = abs.(R[1:m+1])
        β = norm(R_abs[1:m+1])
        lmul!(1 / β, view(R, 1:m+1))
        copyto!(wrk.arnoldi_vecs[1], wrk.v)
        lmul!(R[1], wrk.v)
        for i = 2:m+1
            axpy!(R[i], wrk.arnoldi_vecs[i], wrk.v)
        end

        # Convergence check (relative error of Newton series)
        Ψ_relerr = β * abs(wrk.a[n_a-1]) / (1 + norm(Ψ))
        if Ψ_relerr < relerr
            break
        else
            s += 1
            @assert s <= max_restarts
        end

    end

    wrk.restarts = max(0, s - 1)
    wrk.n_leja = n_leja
    wrk.n_a = n_a
    return Ψ

end

end
