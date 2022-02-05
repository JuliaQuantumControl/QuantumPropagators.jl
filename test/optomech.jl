using SparseArrays
using LinearAlgebra

# Parameters
ω_mech = 10.0
Δ = -ω_mech

# Constants
g = 1.0
η = 2.0
κ = 1.0

N_cav = 4
N_mech = 10

function destroy(N)
    spdiagm(1 => complex.(sqrt.(collect(1:N))))
end

function create(N)
    spdiagm(-1 => complex.(sqrt.(collect(1:N))))
end

function Id(N)
    spdiagm(0 => ones(ComplexF64, N + 1))
end

function ket(n, N)
    Ψ = zeros(ComplexF64, N + 1)
    Ψ[n+1] = 1.0
    return Ψ
end

⊗(A, B) = kron(A, B)

a = destroy(N_cav) ⊗ Id(N_mech)
at = create(N_cav) ⊗ Id(N_mech)

b = Id(N_cav) ⊗ destroy(N_mech)
bt = Id(N_cav) ⊗ create(N_mech)

H_cav = -Δ * at * a + η * (a + at)
H_mech = ω_mech * bt * b
H_int = -g * (bt + b) * at * a
