# These functions are meant to be included in individual tests or benchmarks
using Random
using Distributions
using LinearAlgebra
using SparseArrays


function random_complex_matrix(N, ρ)
    σ = 1/√N
    d = Normal(0.0, σ)
    H = ρ * (rand(d, (N, N)) + rand(d, (N, N)) * 1im) / √2
end


function random_real_matrix(N, ρ)
    σ = 1/√N
    d = Normal(0.0, σ)
    H = ρ * rand(d, (N, N))
end


function random_hermitian_matrix(N, ρ)
    σ = 1/√N
    d = Normal(0.0, σ)
    X = rand(d, (N, N))
    H = ρ * (X + X') / (2*√2)
end


function random_complex_sparse_matrix(N, ρ, density)
    σ = 1/√(density * N)
    d = Normal(0.0, σ)
    Hre = sprand(N, N, density, (dims...) -> rand(d, dims...))
    Him = sprand(N, N, density, (dims...) -> rand(d, dims...))
    H = ρ * (Hre + Him * 1im) / √2
end


function random_real_sparse_matrix(N, ρ, density)
    σ = 1/√(density * N)
    d = Normal(0.0, σ)
    H = ρ * sprand(N, N, density, (dims...) -> rand(d, dims...))
end


function random_hermitian_sparse_matrix(N, ρ, density)
    σ = 1/√(density * N)
    d = Normal(0.0, σ)
    H = sprand(N, N, density, (dims...) -> rand(d, dims...))
    return 0.5ρ * (H + H') / √2
end


function random_state_vector(N)
    Ψ = rand(N) .* exp.((2π * im) .* rand(N))
    Ψ ./= norm(Ψ)
    return Ψ
end
