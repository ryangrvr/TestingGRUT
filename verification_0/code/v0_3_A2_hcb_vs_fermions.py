"""A2: hard-core bosons vs free fermions (SF-1), exact sector blocks.
HCB: real-space ED, no signs.  Fermions: real-space ED with canonical anticommutation signs (independent of JW),
plus the free-fermion form-factor formula.  All odd N at L in {6,10,12}; even-N counter-checks."""
import numpy as np
from v0_3_lib import *

hops = {1: -1.0, -1: -1.0}
eps = eps_from_hops(hops)


def sector(L, N, kind, hops=hops):
    B, idx, H = hamiltonian(L, N, kind, ring_bonds(L, hops))
    ev, V = np.linalg.eigh(H)
    return B, ev, V


def supports(B, ev, V, L, N):
    E0 = ev[0]
    deg = int(np.sum(ev - E0 < 1e-9))
    return deg, [ed_support(ev, V, V[:, 0], E0, rho_diag(B, range(L), 2 * np.pi * m / L), N) for m in range(L)]


allok = True
for L in (6, 10, 12):
    for N in range(1, L, 2):
        Bh, evh, Vh = sector(L, N, 'hcb')
        Bf, evf, Vf = sector(L, N, 'fermion')
        dh, Sh = supports(Bh, evh, Vh, L, N)
        df, Sf = supports(Bf, evf, Vf, L, N)
        occ, uniq, gap = fermion_occupation(eps, L, N)
        Sff = [ff_support(eps, L, occ, m) for m in range(L)]
        spec_eq = np.allclose(evh, evf, atol=1e-9)
        s1 = all(same_support(a, b) for a, b in zip(Sh, Sf))
        s2 = all(same_support(a, b) for a, b in zip(Sf, Sff))
        nlines = sum(len(s) for s in Sh)
        ok = dh == 1 and df == 1 and spec_eq and s1 and s2
        allok &= ok
        print(f"L={L:2d} N={N:2d}: GSdeg HCB={dh} F={df}  E0={evh[0]:+.10f} (sum of N lowest eps={np.sort(eps(kgrid(L)))[:N].sum():+.10f})"
              f"  full spectra equal={spec_eq}  supp(HCB)==supp(F-ED)={s1}  supp(F-ED)==form-factor={s2}"
              f"  #(q,w) points={nlines}  ok={ok}")
print("A2 (all odd N, L in {6,10,12}):", "PASS" if allok else "FAIL")

print("\nEven-N counter-checks (HCB = antiperiodic fermions, NOT SF-1 periodic):")
for L, N in ((6, 2), (10, 4), (12, 6), (12, 4)):
    Bh, evh, Vh = sector(L, N, 'hcb')
    Bf, evf, Vf = sector(L, N, 'fermion')
    dh, Sh = supports(Bh, evh, Vh, L, N)
    df = int(np.sum(evf - evf[0] < 1e-9))
    # antiperiodic free fermions: k = 2pi(j+1/2)/L
    ka = 2 * np.pi * (np.arange(L) + 0.5) / L
    E_ap = np.sort(-2 * np.cos(ka))[:N].sum()
    print(f"  L={L} N={N}: E0(HCB)={evh[0]:+.10f}  E0(periodic F)={evf[0]:+.10f} (GS deg {df})  "
          f"E0(antiperiodic F)={E_ap:+.10f}  HCB GS deg={dh}  spectra equal(HCB vs periodic F)={np.allclose(evh, evf)}")

print("\nScope check: HCB with an added next-nearest-neighbour hop t2 (JW string no longer cancels):")
hops2 = {1: -1.0, -1: -1.0, 2: -0.5, -2: -0.5}
eps2 = eps_from_hops(hops2)
for L, N in ((10, 3), (12, 5)):
    Bh, evh, Vh = sector(L, N, 'hcb', hops2)
    Bf, evf, Vf = sector(L, N, 'fermion', hops2)
    dh, Sh = supports(Bh, evh, Vh, L, N)
    df, Sf = supports(Bf, evf, Vf, L, N)
    eqs = all(same_support(a, b) for a, b in zip(Sh, Sf))
    print(f"  L={L} N={N}: E0(HCB)={evh[0]:+.10f}  E0(F)={evf[0]:+.10f}  spectra equal={np.allclose(evh, evf)}  "
          f"supports equal={eqs}")
