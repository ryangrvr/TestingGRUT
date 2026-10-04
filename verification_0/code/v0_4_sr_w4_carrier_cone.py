#!/usr/bin/env python3
"""V0-4 second reader, witness W4 (E-9 influence cone; E-22 CA-1 classical-carrier cut).

Fresh code. One ring operator K (Laplacian of an N-ring, plus a small pin), one site-local readout.
Three sectors on the same K (the record's CA-1 list):
  phonon       : q'' = -K q            omega_k = sqrt(lambda_k)   J(w) = pi sum w_k delta(w - omega_k)/(2 omega_k)
  Schroedinger : i psi' = K psi        omega_k = lambda_k         J(w) = pi sum w_k delta(w - lambda_k)
  classical relaxational : x' = -K x   J(w) = sum w_k w lambda_k/(w^2 + lambda_k^2)   (no upper edge)
Noise: quantum vacuum nu = (hbar/2) J (the floor, by construction) for the two quantum carriers;
classical FDT nu = (2T/w) J for the relaxational carrier. hbar = 1.
Cone (E-9): J >= 0 and nu >= J/2.
Checks: (i) cone membership; (ii) the classical carrier fails iff w > 4T (any T, unbounded support);
(iii) the two SURVIVING carriers have different spectral-bottom (w -> 0) edge exponents and different
upper band edges, so the earned elimination (CON) leaves edge data free.
"""
import numpy as np


def main():
    out = []
    P = out.append
    P("W4  E-9 cone / E-22 classical-carrier cut on one ring K (fresh code)")
    P("=" * 100)
    N = 4000
    pin = 0.0
    lam = 2 - 2*np.cos(2*np.pi*np.arange(N)/N) + pin
    wk = np.full(N, 1.0/N)       # site-local weights on a translation-invariant ring
    lam_pos = lam[lam > 1e-12]; w_pos = wk[lam > 1e-12]
    # classical relaxational carrier at T = 0 and T = 0.5
    for T in (0.0, 0.5):
        for w in (1.0, 3.0):
            J = np.sum(w_pos*w*lam_pos/(w**2 + lam_pos**2))
            nu = 2*T*J/w
            P(f"  classical relaxational: T={T} w={w}: J={J:.4f} nu={nu:.4f} nu-J/2={nu-J/2:+.4f}"
              f"  ({'admissible' if nu-J/2 >= 0 else 'VIOLATES cone'}) ; threshold 4T={4*T}")
    # surviving quantum carriers: low-frequency edge exponents from counting measure
    def dos_exponent(omegas, weights, w_lo, w_hi):
        edges = np.logspace(np.log10(w_lo), np.log10(w_hi), 9)
        h, _ = np.histogram(omegas, bins=edges, weights=weights)
        dens = h/np.diff(edges)
        mid = np.sqrt(edges[1:]*edges[:-1])
        m = dens > 0
        return np.polyfit(np.log(mid[m]), np.log(dens[m]), 1)[0]
    om_ph = np.sqrt(lam_pos)
    om_sc = lam_pos
    e_ph = dos_exponent(om_ph, w_pos*np.pi/(2*om_ph), 0.02, 0.2)
    e_sc = dos_exponent(om_sc, w_pos*np.pi, 0.002, 0.04)
    P(f"  phonon       : band (0, {om_ph.max():.3f}],  low-w exponent of J ~ w^{e_ph:.3f} (1D theory: -1); on the floor nu = J/2")
    P(f"  Schroedinger : band (0, {om_sc.max():.3f}],  low-w exponent of J ~ w^{e_sc:.3f} (1D theory: -1/2); on the floor nu = J/2")
    P("  -> both survive the cone (E-9 / E-22) on the same K; their dynamical exponents (z = 1 vs 2),")
    P("     spectral-bottom exponents and upper edges differ. The cut removes one carrier; it fixes no edge datum.")
    # E-9 alone: two cone members with different edge exponents
    w = np.linspace(1e-3, 5, 2000)
    for name, J in (("J1 = w e^-w", w*np.exp(-w)), ("J2 = 7 w^3 e^-w/2", 7*w**3*np.exp(-w/2))):
        nu = 0.5*J/np.tanh(w/(2*0.2))   # thermal at T = 0.2, hbar = 1: nu = (J/2) coth(w/2T) >= J/2
        P(f"  E-9: {name}: min(J) = {J.min():.2e} >= 0, min(nu - J/2) = {np.min(nu - J/2):.2e} >= 0 -> in cone")
    P("  -> same earned cone predicate, edge exponents 1 vs 3: the cone constrains SIGN, not edge values.")
    print("\n".join(out))


if __name__ == "__main__":
    main()
