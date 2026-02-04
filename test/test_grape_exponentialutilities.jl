using Test
using GRAPE
using QuantumControl
using QuantumControl.Functionals: J_T_sm
using QuantumControlTestUtils.DummyOptimization: dummy_control_problem
using StableRNGs: StableRNG
using ExponentialUtilities

@testset "GRAPE with ExponentialUtilities" begin
    rng = StableRNG(20260204)
    problem = dummy_control_problem(
        N=4,
        n_trajectories=1,
        n_controls=1,
        n_steps=10,
        dt=0.1,
        density=1.0,
        hermitian=true,
        complex_operators=true,
        prop_method=ExponentialUtilities,
        J_T=J_T_sm,
        iter_stop=1,
        rng=rng
    )

    res = QuantumControl.optimize(problem; method=GRAPE, iter_stop=1, verbose=false)

    @test res.iter >= res.iter_start
    @test isfinite(res.J_T)
end
