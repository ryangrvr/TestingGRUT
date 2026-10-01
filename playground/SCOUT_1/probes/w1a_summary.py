"""W1-A summary: g*L at the predicted soft momentum q=2*pi*nu, at the midpoint q=pi*nu (expected hard),
and at the first momentum q=2*pi/L (sound mode), across L. Plus U(1)-breaking hostile (Kitaev/BdG chain)."""
import numpy as np
from w1a_lsm_soft_points import gaps


def run(label, stat, nu, Ls, **kw):
    print(f"\n{label}")
    print("   L   g*L(q=2pi/L)   g*L(q=pi*nu)   g*L(q=2pi*nu)   g*L(q=4pi*nu)")
    for L in Ls:
        N = int(round(nu * L))
        g, Lc, step, E0, m0 = gaps(L, N, stat, **kw)
        def at(qmom):
            idx = int(round(qmom * Lc / (2 * np.pi * step / step) / step)) % Lc if step == 1 else \
                int(round(qmom / (2 * np.pi) * Lc * step)) % Lc
            return g[idx] * L
        # momentum in units of the (cell) Brillouin zone: q_cell = q*step
        def atq(q):
            m = int(round(q * step / (2 * np.pi) * Lc)) % Lc
            return g[m] * L
        print(f"  {L:3d}   {atq(2*np.pi/L*step if step==1 else 2*np.pi/(L)):12.3f}   {atq(np.pi*nu):12.3f}"
              f"   {atq(2*np.pi*nu):13.3f}   {atq(4*np.pi*nu):13.3f}")


if __name__ == "__main__":
    run("free fermions nu=1/2", "F", 0.5, [8, 12, 16, 20])
    run("fermions V=+1 nu=1/2", "F", 0.5, [8, 12, 16, 20], V=1.0)
    run("hard-core bosons nu=1/2", "HCB", 0.5, [8, 12, 16, 20])
    run("Bose-Hubbard U=1 nmax=3 nu=1/2", "BH", 0.5, [8, 10, 12, 14], U=1.0, nmax=3)
    run("Bose-Hubbard U=4 nmax=3 nu=1/2", "BH", 0.5, [8, 10, 12, 14], U=4.0, nmax=3)
    run("free fermions nu=1/4", "F", 0.25, [8, 12, 16, 20, 24])
    run("hard-core bosons nu=1/4", "HCB", 0.25, [8, 12, 16, 20, 24])
    run("Bose-Hubbard U=4 nmax=3 nu=1/4", "BH", 0.25, [8, 12, 16])
    run("HCB V=1 nu=1/3", "HCB", 1 / 3, [9, 12, 15, 18, 21])
    run("fermions V=3 nu=1/3 (strong NN repulsion)", "F", 1 / 3, [9, 12, 15, 18, 21], V=3.0)
    print("\nHOSTILE: translation broken to period 2 (staggered 0.5); momenta now in the cell zone (mod pi)")
    run("F nu=1/2 stag (cell filling 1)", "F", 0.5, [8, 12, 16, 20], stag=0.5)
    run("HCB nu=1/2 stag (cell filling 1)", "HCB", 0.5, [8, 12, 16, 20], stag=0.5)
    run("F nu=1/4 stag (cell filling 1/2)", "F", 0.25, [8, 12, 16, 20, 24], stag=0.5)

    # U(1) breaking: Kitaev/BdG chain, E(k) = sqrt((2t cos k + mu)^2 + (2 Delta sin k)^2), t=1
    print("\nHOSTILE: U(1) broken (pairing Delta); soft momenta = zeros of E(k)")
    ks = np.linspace(-np.pi, np.pi, 200001)
    for Delta in (0.0, 0.3):
        for mu in (-1.0, 0.0, 0.7, 2.0):
            E = np.sqrt((2 * np.cos(ks) + mu) ** 2 + (2 * Delta * np.sin(ks)) ** 2)
            gap = E.min()
            kz = ks[np.argmin(E)]
            nu_free = np.arccos(np.clip(-mu / 2, -1, 1)) / np.pi  # filling the Delta=0 chain would have
            print(f"   Delta={Delta:.1f} mu={mu:+.1f}: min E = {gap:.4f} at |k| = {abs(kz):.4f}"
                  + (f"   (Delta=0: k_F = pi*nu, nu={nu_free:.4f})" if Delta == 0 else ""))
    print("   Delta != 0: gapless only at mu = +-2 with soft k in {0, pi} = fixed points of k -> -k")
