"""SCOUT-1 W1-P: does complete passivity force a single bath temperature (KMS), collapsing Q2's
mode-temperature profile T(omega) to one number?

Bath = independent two-level modes (frequencies w_i), state = product of mode Gibbs states at T_i.
States are diagonal in the energy basis, so ergotropy of rho^{(x)N} is computed exactly by sorting:
  W_N = E - sum_k r_k^(desc) e_k^(asc).
Pusz-Woronowicz / Lenard: rho completely passive (W_N = 0 for all N) <=> Gibbs at a single beta in [0, inf].
"""
import itertools
import numpy as np
from scipy.optimize import brentq

np.set_printoptions(precision=6, suppress=True)


def mode_pops(w, T):
    if T == 0:
        return np.array([1.0, 0.0])
    z = 1 + np.exp(-w / T)
    return np.array([1 / z, np.exp(-w / T) / z])


def joint(ws, Ts):
    """populations and energies of the product state over 2^m levels"""
    p = np.array([1.0]); e = np.array([0.0])
    for w, T in zip(ws, Ts):
        pm = mode_pops(w, T)
        p = np.outer(p, pm).ravel()
        e = np.add.outer(e, np.array([0.0, w])).ravel()
    return p, e


def ergotropy(p, e):
    E = np.sum(p * e)
    return E - np.sum(np.sort(p)[::-1] * np.sort(e))


def copies(p, e, N):
    P, Eng = p.copy(), e.copy()
    for _ in range(N - 1):
        P = np.outer(P, p).ravel()
        Eng = np.add.outer(Eng, e).ravel()
    return P, Eng


def entropy(p):
    q = p[p > 0]
    return -np.sum(q * np.log(q))


def bound_limit(p, e):
    """asymptotic ergotropy per copy: E(rho) - E(gibbs_beta*) with S(gibbs_beta*) = S(rho)"""
    S = entropy(p)
    def gib(beta):
        x = np.exp(-beta * (e - e.min())); return x / x.sum()
    f = lambda beta: entropy(gib(beta)) - S
    if abs(f(0.0)) < 1e-14:
        b = 0.0
    else:
        b = brentq(f, 1e-9, 1e4)
    return np.sum(p * e) - np.sum(gib(b) * e), b


def report(label, ws, Ts, Nmax=9):
    p, e = joint(ws, Ts)
    single = ergotropy(p, e)
    lim, bstar = bound_limit(p, e)
    print(f"\n{label}: w = {ws}, T = {Ts}, beta_i*w_i = {[round(w / T, 4) if T else np.inf for w, T in zip(ws, Ts)]}")
    print(f"   single-copy ergotropy W_1 = {single:.3e}  -> {'PASSIVE' if single < 1e-13 else 'NOT passive'}")
    row = []
    for N in range(1, Nmax + 1):
        if 2 ** (len(ws) * N) > 2 ** 22:
            break
        P, E = copies(p, e, N)
        row.append(ergotropy(P, E) / N)
    print("   W_N / N, N = 1..:", " ".join(f"{x:.3e}" for x in row))
    print(f"   asymptotic per-copy limit E - E(gibbs at equal entropy) = {lim:.4e}  (beta* = {bstar:.4f})")
    return row, lim


if __name__ == "__main__":
    print("=== KMS controls (single temperature) ===")
    report("KMS T=1", [1.0, 2.0], [1.0, 1.0])
    report("KMS T=0.3, three modes", [0.5, 1.0, 1.7], [0.3, 0.3, 0.3], Nmax=6)
    report("ground state (T=0, beta=inf)", [1.0, 2.0], [0, 0])

    print("\n=== non-KMS but single-copy passive (beta_i w_i non-decreasing with w) ===")
    report("T(w): T1=1.0, T2=1.5", [1.0, 2.0], [1.0, 1.5])   # w/T = 1, 1.333: passive
    report("T(w): T1=1.0, T2=0.8", [1.0, 2.0], [1.0, 0.8])   # w/T = 1, 2.5: passive
    report("three modes, T rising slowly", [0.5, 1.0, 1.7], [0.3, 0.35, 0.4], Nmax=6)

    print("\n=== non-KMS and not even single-copy passive (population inversion across modes) ===")
    report("T1=1.0, T2=3.0", [1.0, 2.0], [1.0, 3.0])         # w/T = 1, 0.667: not passive

    print("\n=== single-copy passivity condition (two qubits): passive iff w1/T1 <= w2/T2 for w1<w2 ===")
    ok = True
    rng = np.random.default_rng(3)
    for _ in range(2000):
        w1, w2 = np.sort(rng.uniform(0.1, 3, 2)); T1, T2 = rng.uniform(0.1, 3, 2)
        p, e = joint([w1, w2], [T1, T2])
        pred = w1 / T1 <= w2 / T2 + 1e-12
        ok &= (ergotropy(p, e) < 1e-12) == pred
    print("   2000 random cases agree with the analytic condition:", ok)
    print("   => for TWO modes single-copy passivity leaves a profile freedom (w/T(w) non-decreasing);")
    print("      the three-mode case above already fails: level w1+w2=1.5 < w3=1.7 but carries less population")

    print("\n=== many modes, ONE copy: how much temperature-profile freedom survives single-copy passivity? ===")
    print("   profile T(w) = T0*(1 + eps*(w-1)), T0 = 0.5, modes equally spaced in [0.5, 1.5]")
    def passive_single(ws, Ts):
        p, e = joint(ws, Ts)
        return ergotropy(p, e) < 1e-13 * max(1, np.sum(p * e))
    for m in (2, 3, 4, 6, 8, 12, 16, 20):
        ws = list(np.linspace(0.5, 1.5, m))
        lo, hi = 0.0, 2.0
        if passive_single(ws, [0.5 * (1 + hi * (w - 1)) for w in ws]):
            print(f"   m={m:2d}: passive up to eps >= {hi}"); continue
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if passive_single(ws, [0.5 * (1 + mid * (w - 1)) for w in ws]):
                lo = mid
            else:
                hi = mid
        print(f"   m={m:2d}: single-copy passive iff |eps| <~ {lo:.4f} (eps>0 side)")
    print("   dense-spectrum argument: levels {a,b} vs {c} with c ~ a+b force f(a)+f(b) = f(a+b), f(w) = w/T(w)")
    print("   (Cauchy) => f linear => single beta. A macroscopic bath acts as its own copies.")

    print("\n=== what complete passivity leaves: the FDT ratio shape fixed, its scale T free ===")
    for T in (0.1, 1.0, 10.0):
        w = np.array([0.5, 1.0, 2.0])
        print(f"   T={T:5.1f}: noise/dissipation ratio coth(w/2T) at w=0.5,1,2 ->", np.round(1 / np.tanh(w / (2 * T)), 4),
              " ; in scaled variable w/T it is the same function")
