"""Regression checks for MPO multiplication orientation and thermal factor."""
import numpy as np
import quimb.operator as qop

from xtrg_quimb import build_spin_mpo
from xtrg_quimb_purified import purified_square
from matrix_free_hamiltonian import SpinBond


def chain_bonds(n):
    coupling = np.diag([0.7, 0.7, 1.1])
    return [SpinBond(i, i + 1, coupling) for i in range(n - 1)]


def main():
    n = 4
    builder, h_mpo = build_spin_mpo(n, chain_bonds(n), longitudinal_field_K=0.2)
    dense_h = np.asarray(builder.build_dense())
    assert np.max(np.abs(np.asarray(h_mpo.to_dense()) - dense_h)) < 1e-12

    # Explicitly validate quimb's index orientation for X^dagger X.
    x = h_mpo + h_mpo.multiply(0.1j)
    dense_x = np.asarray(x.to_dense())
    rho = x.H.apply(x, which_A="upper", which_B="upper", compress=False)
    assert np.max(np.abs(np.asarray(rho.to_dense()) - dense_x.conj().T @ dense_x)) < 1e-11

    eigenvalues = np.linalg.eigvalsh(dense_h)
    rows, diagnostic = purified_square(
        h_mpo,
        n,
        beta_final=1.0,
        squarings=10,
        max_bond=16,
        cutoff=1e-13,
        eigenvalues=eigenvalues,
    )
    # The first-order high-temperature initialization leaves an O(beta_0)
    # error; at ten doublings both observables should still be below 1e-4.
    assert abs(rows[-1]["errors"]["entropy_kB_per_spin"]) < 1e-4
    assert abs(rows[-1]["errors"]["energy_K"]) < 1e-4
    assert diagnostic["final_min_density_eigenvalue"] > -1e-12
    print("quimb MPO orientation and purified thermal regression passed")


if __name__ == "__main__":
    main()
