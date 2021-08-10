using QuantumPropagators
using Test

@testset "QuantumPropagators.jl" begin
    include("test_cheby.jl")
    include("test_newton.jl")
end
