#!/usr/bin/env python3
"""V0-4 second reader, witness W1 (E-13 / SF-1 against the corrected edge-data quotient).

Fresh code. Nothing imported from the record or from any earlier audit.

Question: does the SF-1 earned predicate ("one fixed free-fermion parent, different exact conserved
particle-number scaling sectors -> inequivalent IR effective-law classes") FIX the edge-data set?

Method: build three one-band free (exclusion / CAR) parents on a periodic ring; for each, compute the
ground-state occupation set O in a dense sector (fixed density nu) and in a dilute sector (fixed N),
then the particle-hole lower edge
    omega_minus(q) = min_{k in O, k+q not in O} [eps(k+q) - eps(k)]
(free-fermion density-response support, same object as SF-1's I-FF support).
From it extract the edge-data set (edge momenta, edge velocities, first non-zero derivative order r_e),
the dynamical exponent z (log-slope of omega_minus at small q) and the soft set (q in [0, pi], quotient
q ~ -q ~ q + 2 pi, where omega_minus -> 0).

Parents:
  A  eps = -2 cos k                     (SF-1's parent; single pocket, quadratic bottom)
  B  eps = -2 cos k + 0.5 cos 2k        (single pocket, QUARTIC bottom: eps = -1.5 + k^4/4 + ...)
  C  eps = -2 cos k + 2 cos 2k          (two degenerate minima at +-arccos(1/4): two pockets)

Each parent satisfies the SF-1 predicate (dense class != dilute class). The edge data differ.
"""
import numpy as np

PARENTS = {
    "A": (lambda k: -2*np.cos(k), lambda k: 2*np.sin(k)),
    "B": (lambda k: -2*np.cos(k) + 0.5*np.cos(2*k), lambda k: 2*np.sin(k) - np.sin(2*k)),
    "C": (lambda k: -2*np.cos(k) + 2*np.cos(2*k), lambda k: 2*np.sin(k) - 4*np.sin(2*k)),
}


def ground_occupation(eps, L, N):
    """N lowest single-particle levels on the L-ring (momenta 2 pi m / L). Returns occupied index set
    and a flag saying whether the N-th / (N+1)-th levels are split (unique ground state)."""
    k = 2*np.pi*np.arange(L)/L
    e = eps(k)
    order = np.argsort(e, kind="stable")
    occ = np.zeros(L, dtype=bool)
    occ[order[:N]] = True
    unique = (e[order[N]] - e[order[N-1]]) > 1e-12 if N < L else True
    return occ, e, unique


def omega_minus(occ, e, L, m):
    """lower edge of the particle-hole support at q = 2 pi m / L."""
    idx = np.nonzero(occ)[0]
    tgt = (idx + m) % L
    ok = ~occ[tgt]
    if not ok.any():
        return np.inf
    return float(np.min(e[tgt[ok]] - e[idx[ok]]))


def edge_momenta(occ, L):
    """momenta (in (-pi, pi]) of occupied levels adjacent to an empty level: the occupation edges."""
    k = 2*np.pi*np.arange(L)/L
    k = np.where(k > np.pi, k - 2*np.pi, k)
    edges = []
    for i in np.nonzero(occ)[0]:
        if (not occ[(i+1) % L]) or (not occ[(i-1) % L]):
            edges.append(k[i])
    return np.array(sorted(edges))


def quotient(q):
    q = np.mod(q, 2*np.pi)
    return np.minimum(q, 2*np.pi - q)


def soft_set(occ, e, L, vmax):
    """soft momenta: local minima of omega_minus(q) on [0, pi] whose value is O(1/L)
    (below 3 * vmax * 2 pi / L); q = 0 is soft by definition (the support always starts at 0)."""
    ms = np.arange(0, L//2 + 1)
    w = np.array([omega_minus(occ, e, L, m) if m > 0 else 0.0 for m in ms])
    q = 2*np.pi*ms/L
    tol = 3.0*max(vmax, 1e-3)*2*np.pi/L
    reps = [0.0]
    for m in range(2, len(ms) - 1):
        if w[m] <= w[m-1] and w[m] <= w[m+1] and w[m] < tol:
            if q[m] - reps[-1] > 0.05:
                reps.append(float(q[m]))
    if w[-1] < tol and w[-1] <= w[-2] and q[-1] - reps[-1] > 0.05:
        reps.append(float(q[-1]))
    mask = np.ones_like(q, dtype=bool)
    for r in reps:
        mask &= np.abs(q - r) > 0.15
    away = float(np.min(w[mask])) if mask.any() else float("nan")
    return reps, away


def soft_set_dilute_limit(eps):
    """dilute (fixed-N) sector, prescription P: L -> infinity first, so every particle sits at a band
    minimum k_a; omega_inf(q) = min_a [eps(k_a + q) - eps_min]. Soft set = its zeros on [0, pi]."""
    kk = np.linspace(-np.pi, np.pi, 400001)
    ee = eps(kk)
    emin = ee.min()
    mins = kk[ee < emin + 1e-12]
    # cluster minima
    ka = []
    for x in mins:
        if not ka or x - ka[-1] > 1e-3:
            ka.append(x)
    q = np.linspace(0, np.pi, 200001)
    w = np.min(np.stack([eps(a + q) - emin for a in ka] + [eps(a - q) - emin for a in ka]), axis=0)
    reps = []
    for i in range(len(q)):
        lo = w[i-1] if i > 0 else np.inf
        hi = w[i+1] if i + 1 < len(q) else np.inf
        if w[i] <= lo and w[i] <= hi and w[i] < 1e-8:
            if not reps or q[i] - reps[-1] > 0.05:
                reps.append(float(q[i]))
    mask = np.ones_like(q, dtype=bool)
    for r in reps:
        mask &= np.abs(q - r) > 0.15
    return reps, float(np.min(w[mask]))


def z_exponent_dense(eps, nu, even=False, L=200001):
    """prescription P at finite density: L large, small q; z = log-slope of omega_minus vs q."""
    N = int(round(nu*L))
    if (N % 2 == 0) != even:
        N += 1
    occ, e, uniq = ground_occupation(eps, L, N)
    ms = np.array([20, 40, 80, 160])
    w = np.array([omega_minus(occ, e, L, m) for m in ms])
    q = 2*np.pi*ms/L
    s = np.polyfit(np.log(q), np.log(w), 1)[0]
    return s, uniq


def z_exponent_dilute(eps, N, L=2**18):
    """prescription P at fixed N: k_F -> band bottom first, then small q."""
    occ, e, uniq = ground_occupation(eps, L, N)
    qs = np.array([0.02, 0.04, 0.08, 0.16])
    ms = np.round(qs*L/(2*np.pi)).astype(int)
    w = np.array([omega_minus(occ, e, L, m) for m in ms])
    q = 2*np.pi*ms/L
    s = np.polyfit(np.log(q), np.log(w), 1)[0]
    return s, uniq


def deriv_order(eps, k0, h=1e-2):
    """first non-zero derivative order of eps at k0 (numerical, orders 1..4)."""
    from math import factorial
    # finite-difference Taylor coefficients via polyfit on a small window
    x = np.linspace(-h, h, 41)
    c = np.polyfit(x, eps(k0 + x) - eps(k0), 6)[::-1]  # c[0] + c[1] x + ...
    for r in range(1, 5):
        if abs(c[r]) > 1e-4:
            return r, c[r]
    return None, None


def main():
    out = []
    P = out.append
    P("W1  E-13 / SF-1 versus the corrected edge-data quotient (fresh code, numpy %s)" % np.__version__)
    P("=" * 100)
    Ls = 2002  # soft-set scan size (even L; N chosen per family below)
    for name, (eps, deps) in PARENTS.items():
        P(f"\nPARENT {name}")
        # dense family: nu = 1/2 (N odd for single pocket, as in SF-1; even for the 2-pocket parent)
        fams = {"A": (1001, 501, 1), "B": (1001, 501, 1), "C": (1000, 500, 2)}[name]
        for fam, N in zip(("D (nu=1/2)", "D-1/4 (nu=1/4)", "E (fixed N)"), fams):
            occ, e, uniq = ground_occupation(eps, Ls, N)
            edges = edge_momenta(occ, Ls)
            vel = np.abs(deps(edges))
            if fam.startswith("E"):
                reps, away = soft_set_dilute_limit(eps)
            else:
                reps, away = soft_set(occ, e, Ls, float(np.max(vel)) if len(vel) else 1.0)
            if fam.startswith("E"):
                z, u2 = z_exponent_dilute(eps, N)
                kmin = edges
            else:
                nu = N/Ls
                z, u2 = z_exponent_dense(eps, nu, even=(name == "C"))
            orders = []
            for k0 in np.unique(np.round(np.abs(edges), 6)):
                if fam.startswith("E"):
                    # dilute: the L -> infinity edge sits at the band minimum; order of eps there
                    from scipy.optimize import minimize_scalar
                    kk = np.linspace(0, np.pi, 20001)
                    k_guess = kk[np.argmin(eps(kk))]
                    res = minimize_scalar(eps, bracket=(max(k_guess-0.01, -0.01), k_guess, k_guess+0.01)) if k_guess > 0 else None
                    kb = float(res.x) if res is not None else 0.0
                    r, c = deriv_order(eps, kb)
                    orders.append((round(float(kb), 4), r))
                    break
                r, c = deriv_order(eps, k0)
                orders.append((round(float(k0), 4), r))
            P(f"  {fam:16s} N={N:5d} L={Ls} unique_GS={bool(uniq)}")
            P(f"     edge momenta (finite L)   : {np.round(edges, 4).tolist()}")
            P(f"     |edge velocities|         : {np.round(vel, 4).tolist()}")
            P(f"     (edge, first nonzero deriv order r_e) in L->inf : {orders}")
            P(f"     soft set reps in [0,pi]   : {np.round(reps, 4).tolist()}   (n_soft = {len(reps)});"
              f"  min omega_minus away from soft set = {away:.4f}")
            P(f"     z (prescription P, log-slope) = {z:.4f}")
    P("\nREADING")
    P("  * Every parent satisfies the SF-1 predicate: class(dense) != class(dilute).")
    P("  * Parent A reproduces SF-1: D -> (z=1, soft {0, pi}); D-1/4 -> (1, {0, pi/2}); E -> (2, {0}).")
    P("    D and D-1/4 share the class but have different edge positions and velocities (2 vs sqrt2):")
    P("    the class is a strict compression of the edge-data set, so it cannot fix those values.")
    P("  * Parent B: same predicate, but the dilute edge has r_e = 4 -> z = 4, class (4, 1) != (2, 1).")
    P("  * Parent C: two pockets. The dense sector has 4 occupation edges with two distinct |velocities|")
    P("    (one v_b cannot describe it) and more soft momenta; the dilute sector has n_soft = 2 (0 and 2k0).")
    P("  => same earned predicate, different edge data, different effective-law class: E-13 does not FIX.")
    text = "\n".join(out)
    print(text)


if __name__ == "__main__":
    main()
