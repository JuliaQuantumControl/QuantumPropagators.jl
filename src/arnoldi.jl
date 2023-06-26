module Arnoldi

using LinearAlgebra

using TimerOutputs: @timeit_debug, TimerOutput

const _DEFAULT_TIMING_DATA = TimerOutput()
# Fallback if no _timing_data is passed to arnoldi!


"""
```julia
m = arnoldi!(Hess, q, m, Ψ, H, dt=1.0; extended=true, norm_min=1e-15)
```

Calculate the Hessenberg matrix and Arnoldi vectors of `H dt`, from `Ψ`.

For a given order `m`, the `m×m` Hessemberg matrix is calculated and stored in
in the pre-allocated `Hess`. Further  an array of `m` normalized Arnoldi
vectors is stored in in the pre-allocated `q`, plus one additional unnormalized
Arnoldi vector.  The unnormalized `m+1`st vector could be used to easily
extend a given `m×m` Hessenberg matrix to a `(m+1)×(m+1)` matrix.

If the extended Hessenberg matrix is requested (`extended=true`, default), the
`m+1`st Arnoldi vector is also normalized, and it's norm will be stored in
`m+1, m` entry of the (extended) Hessenberg matrix, which is an `(m+1)×(m+1)`
matrix.

Return the size `m` of the calculated Hessenberg matrix. This will usually be
the input `m`, except when the Krylov dimension of `H` starting from `Ψ` is
less then `m`. E.g., if `Ψ` is an eigenstate of `H`, the returned `m` will
be 1.

See <http://en.wikipedia.org/wiki/Arnoldi_iteration> for a description of
the algorithm.

# Arguments

- `Hess::Matrix{ComplexF64}`: Pre-allocated storage for the Hessemberg matrix.
   Can be uninitialized on input. The matrix must be at least of size `m×m`, or
   `(m+1)×(m+1)` if `extended=true`. On output, the `m×m` sub-matrix of `Hess`
   (with the returned output `m`) will contain the Hessenberg matrix, and all
   other elements of `Hess` be be set to zero.
- `q`: Pre-allocated array of states similar to `Ψ`, as storage for the
  calculated Arnoldi vectors. These may be un-initialized on input. Must be at
  least of length `m+1`
- `m`: The requested dimensions of the output Hessenberg matrix.
- `Ψ`: The starting vector for the Arnoldi procedure. This can be of any type,
   as long as `Φ = H * Ψ` results in a vector similar to `Ψ`, there is an inner
   products of `Φ` and `Ψ` (`Ψ⋅Φ` is defined), and `norm(Ψ)` is defined.
- `H`: The operator (up to `dt`) for which to calculate the Arnoldi procedure.
  Can be of any type, as long as `H * Ψ` is defined.
- `dt`: The implicit time step; the total operator for which to calculate the
  Arnoldi procedure is `H * dt`
- `extended`: If `true` (default), calculate the extended Hessenberg matrix,
  and normalized the final Arnoldi vector
- `norm_min`: the minimum value of the norm of `Ψ` at which `Ψ` should be
   considered the zero vector
"""
function arnoldi!(
    Hess::Matrix{ComplexF64},
    q::Array{T},
    m::Int64,
    Ψ::T,
    H,
    dt::Float64=1.0;
    extended=true,
    norm_min=1e-15,
    _timing_data=_DEFAULT_TIMING_DATA  # undocumented (internal use)
) where {T}
    if extended
        dim_hess = m + 1
    else
        dim_hess = m
    end
    @assert size(Hess, 1) >= dim_hess && size(Hess, 2) >= dim_hess
    @assert length(q) >= m + 1
    fill!(Hess, 0)
    copyto!(q[1], Ψ)
    for j = 1:m
        @timeit_debug _timing_data "matrix-vector product" begin
            mul!(q[j+1], H, q[j])
        end
        for i = 1:j  # Orthogonalization with Gram-Schmidt
            Hess[i, j] = dt * (q[i] ⋅ q[j+1]) # = dt ⟨qᵢ|qⱼ₊₁⟩
            axpy!(-Hess[i, j] / dt, q[i], q[j+1])  # qⱼ₊₁ += -(Hessᵢⱼ/dt) qᵢ
        end
        if (j < m) || extended
            h = norm(q[j+1])
            Hess[j+1, j] = dt * h
            if h < norm_min
                # dimensionality exhausted. Returning reduced m
                m = j
                break
            end
            lmul!(1 / h, q[j+1])
        end
    end
    return m
end


""" Extend dimension of Hessenberg matrix by one.

```julia
extend_arnoldi!(Hess, q, m, H, dt; norm_min=1e-15)
```

extends the entries in `Hess` from size (m-1)×(m-1) to size m×m, and the list
`q` of Arnoldi vectors from m to (m+1). It is assumed that the input `Hess` was
created by a call to [`arnoldi!`](@ref) with `extended=false` or a previous
call to `extend_arnoldi!`. Note that `Hess` itself is not resized, so it must
be allocated to size m×m or greater on input.
"""
function extend_arnoldi!(Hess, q, m, H, dt::Float64=1.0; norm_min=1e-15)
    h = norm(q[m])
    (h < norm_min) && return m
    Hess[m, m-1] = dt * h
    lmul!(1 / h, q[m])
    mul!(q[m+1], H, q[m])
    for i = 1:m
        Hess[i, m] = dt * (q[i] ⋅ q[m+1])
        axpy!(-Hess[i, m] / dt, q[i], q[m+1])
    end
    # everything below the first sub-diagonal should be zero. We'll check the
    # last row (previous rows were checked in earlier extend_arnoldi!)
    @assert all(Hess[m, 1:m-2] .== 0.0)
    return Hess
end


"""
```julia
diagonalize_hessenberg_matrix(Hess, m; accumulate=false)
```

Diagonalize the m × m top left submatrix of the given Hessenberg matrix.

If `accumulate` is `true`, return the concatenated eigenvalues for
`Hess[1:1,1:1]` to `Hess[1:m,1:m]`, that is, all sumatrices of size 1 through
`m`.
"""
function diagonalize_hessenberg_matrix(Hess, m; accumulate=false)
    j_min = m
    j_max = m
    if accumulate
        j_min = 1
        eigenvals = zeros(ComplexF64, (m * (m + 1)) ÷ 2)
    else
        eigenvals = zeros(ComplexF64, m)
    end
    offset = 0
    for j = j_min:j_max
        if j == 1
            eigenvals[1] = Hess[1, 1]
        elseif j == 2
            a = Hess[1, 1]
            c = Hess[2, 1]
            b = Hess[1, 2]
            d = Hess[2, 2]
            s = sqrt(a^2 + 4 * b * c - 2 * a * d + d^2)
            eigenvals[offset+1] = 0.5 * (a + d - s)
            eigenvals[offset+2] = 0.5 * (a + d + s)
        else
            # TODO: Wrap in UpperHessenberg?
            # https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#Special-matrices
            # TODO: use a view?
            eigenvals[offset+1:offset+j] .= eigvals(Hess[1:j, 1:j])
        end
        offset += j
    end
    return eigenvals
end

end
