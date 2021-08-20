"""Test script to compare test_newton.jl with the Python newtonprop reference
implementation.

This relies on test_newton.jl having been run with SERIALIZATION=1.

To compare test_newton.jl and test_newton.py step-by-step, test_newton.jl must
then run with SERIALIZATION=2
"""

import numpy as np
from newtonprop import newton

# Load files generated in test_newton.jl (SERIALIZATION=1)
H = np.load("test/H.npy")
Ψ0 = np.load("test/psi0.npy")
dt = 0.5


def make_apply_H(H):
    """Generate a function that applies H to Ψ."""

    def apply_H(psi):
        return H @ psi

    return apply_H


apply_H = make_apply_H(H)

Ψ_out = newton(apply_H, Ψ0, dt, func='expmi', m_max=5)

assert abs(np.linalg.norm(Ψ_out) - 1) < 1e-10
