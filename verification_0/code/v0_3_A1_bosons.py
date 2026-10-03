"""A1: free bosons, exact sector blocks in real-space Fock basis.  L in {6,10}, N in {1,3,5}."""
import numpy as np
from v0_3_lib import *

hops = {1: -1.0, -1: -1.0}
eps = eps_from_hops(hops)
allok = True
for L in (6, 10):
    for N in (1, 3, 5):
        B, idx, H = hamiltonian(L, N, 'boson', ring_bonds(L, hops))
        ev, V = np.linalg.eigh(H)
        E0 = ev[0]
        degGS = int(np.sum(ev - E0 < 1e-9))
        gap = ev[degGS] - E0
        psi0 = V[:, 0]
        print(f"L={L:2d} N={N}  dim={len(B):5d}  E0={E0:+.12f}  N*eps(0)={-2*N:+d}  GS degeneracy={degGS}  "
              f"gap={gap:.12f}  4sin^2(pi/L)={4*np.sin(np.pi/L)**2:.12f}")
        # ground state = (b_0^dag)^N|0>/sqrt(N!) : in real space all amplitudes prop. to sqrt(N!/prod n_j!)
        from math import factorial
        cond = np.array([np.sqrt(factorial(N) / np.prod([factorial(n) for n in b])) for b in B])
        cond /= np.linalg.norm(cond)
        ov = abs(np.vdot(cond, psi0))
        print(f"      |<condensate|psi0>| = {ov:.15f}")
        ok = degGS == 1 and abs(E0 + 2 * N) < 1e-10 and ov > 1 - 1e-10
        for m in range(L):
            q = 2 * np.pi * m / L
            sup = ed_support(ev, V, psi0, E0, rho_diag(B, range(L), q), N)
            pred = 4 * np.sin(q / 2) ** 2
            line = (len(sup) == 1 and abs(sup[0][0] - pred) < 1e-9)
            wexp = N if m == 0 else 1.0
            wok = abs(sup[0][1] - wexp) < 1e-9 if len(sup) == 1 else False
            ok = ok and line and wok
            print(f"      m={m}: support={[(round(float(e),10), round(float(W),10)) for e,W in sup]}  "
                  f"4sin^2(q/2)={pred:.10f}  single-line={line}  weight={'N' if m==0 else 1} ok={wok}")
        print(f"      => A1 checks pass for (L,N)=({L},{N}): {ok}")
        allok = allok and ok
print("A1 overall:", "PASS" if allok else "FAIL")
print("Note: at m=0, rho_0 = N (number operator) so S(0,0) = N: N-independence of S holds only for q != 0.")
