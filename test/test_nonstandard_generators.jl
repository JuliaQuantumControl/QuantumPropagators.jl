using Test

@testset "PseudoGenerator" begin

    # It can be useful to override the propagation for specific custom
    # states/generators by adding a custom `init_prop`, independent of a
    # specific propagation method. For this to work correctly, it is important
    # that # `propagate` calls `init_prop` with `method` as a keyword argument.
    #
    # This test ensures that a custom `init_prop` works as expected with
    # `propagate`.

    import QuantumPropagators
    using QuantumPropagators: propagate

    struct MyCustomPseudoGeneratorXXX end

    # Does not implement the full interface (`check_propagator`), but works for
    # a simple `propagate` call
    struct MyCustomPseudoPropagatorXXX
        state
        tlist
        backward
    end

    function QuantumPropagators.init_prop(
        state,
        generator::MyCustomPseudoGeneratorXXX,
        tlist;
        kwargs...
    )
        backward = false
        return MyCustomPseudoPropagatorXXX(state, tlist, backward)
    end

    function QuantumPropagators.prop_step!(p::MyCustomPseudoPropagatorXXX)
        return p.state
    end

    Ψ = ComplexF64[0, 0]
    H = MyCustomPseudoGeneratorXXX()
    tlist = [0.0, 1.0]
    Ψ_out = QuantumPropagators.propagate(Ψ, H, tlist; method=:Cheby, check=false)
    # The `method=:Cheby` is irrelevant here: The dispatch for the custom type
    # of `H` should take precedence.
    @test Ψ_out ≡ Ψ

end
