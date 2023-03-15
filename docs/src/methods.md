# Propagation Methods

## Chebychev Propagation

For `method=:cheby`, the time evolution operator of the piecewise-constant Schrödinger equation ``|Ψ(t)⟩ = e^{-i Ĥ dt} |Ψ(0)⟩`` is evaluated by an expansion into Chebychev polynomials [Tal-EzerJCP1984](@cite)[KosloffJCP1988](@cite). This requires ``Ĥ`` to be Hermitian (have real eigenvalues) and to have a known spectral range, so that it can be normalized to the domain ``[-1, 1]`` on which the Chebychev polynomials are defined.

… (TODO) …

## Newton Propagation

For `method=:newton`, the time evolution operator of the piecewise-constant Schrödinger equation ``|Ψ(t)⟩ = e^{-i Ĥ dt} |Ψ(0)⟩`` is evaluated by an expansion into Newton polynomials [BermanJPA1992](@cite)[AshkenaziJCP1995](@cite)[Tal-EzerSJSC2007](@cite). Unlike for Chebychev polynomials, this expansion does not require ``Ĥ`` to be Hermitian or to have a known spectral radius. This makes the Newton propagation applicable to open quantum systems, where ``Ĥ`` is replaced by a Liouvillian to calculate the time evolution of the density matrix.

… (TODO) …

## Propagation with Explicit Matrix Exponentiation

TODO
