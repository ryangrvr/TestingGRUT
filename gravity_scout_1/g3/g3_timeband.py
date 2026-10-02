"""GRAVITY-SCOUT-1 / G3-2 -- time-band conditioning control.  ILLUSTRATION ONLY.

Independent code path, not independent reviewer.

Model: the l = 0 sector of a free scalar in global AdS_4 (d = 3), Delta = 3, truncated to N normal modes,
omega_n = Delta + 2n (standard global-AdS spectrum).  The boundary operator restricted to this sector is
    O(t) = sum_{n<N} c_n ( a_n e^{-i omega_n t} + a_n^dag e^{+i omega_n t} ),
with the mode normalisations c_n set to 1 (they are polynomial in n, rescale the answer by c_m^{-1}, and do not change
the epsilon-dependence studied here).  To extract a_m from O(t) on the time band [0, eps] one needs a smearing f with
    int_0^eps f e^{-i omega_n t} dt = delta_nm ,   int_0^eps f e^{+i omega_n t} dt = 0     (all n < N).
The minimum-L2 solution has ||f||^2 = (G^{-1})_{mm}, G_kj = int_0^eps e^{-i(s_k - s_j)t} dt over the 2N exponents
s in {+omega_n} u {-omega_n}.  At eps = pi all exponent differences are multiples of 2 (Delta odd + odd = even), so
G = pi*I: perfectly conditioned.  For every eps > 0 the system is solvable (exact algebraic access); the question is
only how ill-conditioned it becomes.

NOT identified with a physical measurement precision (owner firewall G3-2).  The construction contains no gravity: it
is the kinematic Fourier part of any time-band reconstruction, gravitational or not (G3-9 control reading).
"""
import mpmath as mp

mp.mp.dps = 160
Delta = 3


def gram(N, eps):
    s = [Delta + 2*n for n in range(N)] + [-(Delta + 2*n) for n in range(N)]
    G = mp.matrix(2*N, 2*N)
    for k in range(2*N):
        for j in range(2*N):
            w = s[k] - s[j]
            G[k, j] = eps if w == 0 else (1 - mp.exp(-1j*w*eps))/(1j*w)
    return G


def analyse(N, eps):
    G = gram(N, eps)
    Gi = G**-1
    ev = sorted(mp.re(x) for x in mp.eighe(G, eigvals_only=True))   # G is Hermitian positive definite
    if ev[0] <= mp.mpf(10)**(-mp.mp.dps + 20):       # beyond working precision: treat as singular
        return mp.inf, mp.inf, mp.inf
    cond = ev[-1]/ev[0]
    f0 = mp.sqrt(mp.re(Gi[0, 0]))            # lowest mode a_0
    ftop = mp.sqrt(mp.re(Gi[N-1, N-1]))      # highest kept mode a_{N-1}
    # benchmark: the same norm on the full period is 1/sqrt(pi)
    return cond, f0*mp.sqrt(mp.pi), ftop*mp.sqrt(mp.pi)


if __name__ == "__main__":
    print("columns: N, eps/pi, cond(G), ||f_0||/||f_0||_(eps=pi), ||f_top||/||f_top||_(eps=pi)")
    for N in [2, 3, 4, 6, 8]:
        for k in [1, 2, 4, 8, 16, 32]:
            eps = mp.pi/k
            c, a0, at = analyse(N, eps)
            print(f"N={N}  eps=pi/{k:<3d} cond={mp.nstr(c, 4):>10}  f0={mp.nstr(a0, 4):>10}  ftop={mp.nstr(at, 4):>10}")
        print()
    print("time-band / cutoff trade-off: smallest eps (bisection, 1e-4 relative) with cond(G) <= threshold")
    for N in [2, 3, 4, 6, 8]:
        om = Delta + 2*(N - 1)
        row = []
        for thr in [mp.mpf(10)**6, mp.mpf(10)**12]:
            lo, hi = mp.pi/4096, mp.pi                 # cond decreases monotonically in eps on this range (checked above)
            if analyse(N, lo)[0] <= thr:
                row.append(lo); continue
            while (hi - lo)/hi > mp.mpf("1e-4"):
                mid = (lo + hi)/2
                if analyse(N, mid)[0] <= thr: hi = mid
                else: lo = mid
            row.append(hi)
        print(f"N={N}  omega_max={om:2d}  eps_min/pi: cond<=1e6 -> {mp.nstr(row[0]/mp.pi, 4)},"
              f" cond<=1e12 -> {mp.nstr(row[1]/mp.pi, 4)}   omega_max*eps_min(1e6) = {mp.nstr(om*row[0], 4)}")
