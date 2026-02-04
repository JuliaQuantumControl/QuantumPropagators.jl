using QuantumControl
using GRAPE
using LinearAlgebra
using QuantumPropagators.Controls: get_controls
using ExponentialUtilities
using Test

# Scale
const GHz = 1.0
const MHz = 0.001 * GHz
const ns = 1.0

function qubit_lvn(eps)
    # Basis vectors
    sigma_x = [0 1; 1 0]
    sigma_z = [1 0; 0 -1]

    H0 = 1 / 2 * sigma_z

    return liouvillian((H0, (sigma_x, eps[1])); convention=:TDSE)
end

@testset "GRAPE ExponentialUtilities Liouvillian" begin
    # Pulse Parameters (Same as unitary)
    tf = 20ns
    Nsteps = 100
    tgrid = collect(range(0, tf; length=Nsteps + 1))

    # Target gate on the qubit subspace (embed 2×2 into 3×3)
    U_tgt = ComplexF64[0 1; 1 0]

    # Basis for process verification (3×3 density operators)
    basis = Matrix{ComplexF64}[]
    push!(basis, [1 0; 0 0])
    push!(basis, [0 0; 0 1])
    push!(basis, 1 / 2 * [1 1; 1 1])
    push!(basis, 1 / 2 * [1 im; -im 1])
    basis_tgt = [U_tgt * b * U_tgt' for b in basis]

    eps = [
        0.2 * (1 + 0.05 * rand()) *
        QuantumControl.Shapes.flattop.(tgrid, T=tf, t_rise=0.3, func=:blackman)
        for _ in 1:2
    ]

    lvn = qubit_lvn(eps)

    trajectories = [
        Trajectory(
            initial_state=reshape(basis[i], :),
            generator=lvn,
            target_state=reshape(basis_tgt[i], :)
        ) for i in eachindex(basis)
    ]

    problem = ControlProblem(
        trajectories=trajectories,
        tlist=tgrid,
        iter_stop=10,
        prop_method=:expv,
        gradient_method=:taylor,
        prop_expv_kwargs=(; ishermitian=false),
        J_T=QuantumControl.Functionals.J_T_re
    )

    println("Starting optimization with QuantumPropagators: ExponentialUtilities...")
    result = QuantumControl.optimize(problem; method=GRAPE, print_iters=true)
    display(result)

    @test isfinite(result.J_T)
end

@testset "GRAPE ExponentialUtilities Gradgen Error" begin
    tf = 1.0
    Nsteps = 10
    tgrid = collect(range(0, tf; length=Nsteps + 1))

    basis = Matrix{ComplexF64}[]
    push!(basis, [1 0; 0 0])
    push!(basis, [0 0; 0 1])
    push!(basis, 1 / 2 * [1 1; 1 1])
    push!(basis, 1 / 2 * [1 im; -im 1])

    U_tgt = ComplexF64[0 1; 1 0]
    basis_tgt = [U_tgt * b * U_tgt' for b in basis]

    eps = [0.2 * QuantumControl.Shapes.flattop.(tgrid, T=tf, t_rise=0.3) for _ in 1:2]
    lvn = qubit_lvn(eps)

    trajectories = [
        Trajectory(
            initial_state=reshape(basis[i], :),
            generator=lvn,
            target_state=reshape(basis_tgt[i], :),
        ) for i in eachindex(basis)
    ]

    problem = ControlProblem(
        trajectories=trajectories,
        tlist=tgrid,
        iter_stop=1,
        prop_method=:expv,
        prop_expv_kwargs=(; ishermitian=false),
        grad_prop_method=:expv,
        gradient_method=:gradgen,
        J_T=QuantumControl.Functionals.J_T_re
    )

    res = QuantumControl.optimize(problem; method=GRAPE, print_iters=false)
    @test occursin("gradient_method=:gradgen", res.message)
end
