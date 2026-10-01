"""SCOUT-1 W1-A: LSM/Oshikawa — are the soft-point momenta fixed at 2*pi*nu*m (mod 2pi)
independent of statistics and interaction, while velocities stay free?

Momentum-resolved exact diagonalization on periodic rings.
Models (all U(1) x translation symmetric):
  F   spinless fermions, H = -t sum (c+_i c_{i+1} + h.c.) + V sum n_i n_{i+1}
  HCB hard-core bosons (same H, no signs)
  BH  Bose-Hubbard, H = -t sum (b+_i b_{i+1} + h.c.) + U/2 sum n_i(n_i-1), occupation cutoff nmax
Diagnostic: g(q) = E_min(K0+q) - E0 in each momentum sector. A soft point is a q with g(q)*L bounded as L grows;
a hard q has g(q)*L growing ~ L.
Hostile: staggered potential (translation -> period 2), U(1) breaking pairing/field, filling change.
"""
import itertools
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla

np.set_printoptions(precision=4, suppress=True)


def fock_basis(L, N, nmax):
    out = []
    def rec(prefix, left, sites):
        if sites == 0:
            if left == 0:
                out.append(tuple(prefix))
            return
        for n in range(min(nmax, left) + 1):
            rec(prefix + [n], left - n, sites - 1)
    rec([], N, L)
    return out


def translate(s, stat):
    """T: site x -> x+1. returns (state, sign)."""
    t = (s[-1],) + s[:-1]
    sign = 1
    if stat == "F" and s[-1] == 1:
        N = sum(s)
        sign = (-1) ** (N - 1)
    return t, sign


def hop_terms(s, L, stat, t, nmax):
    """yield (new_state, amplitude) for -t sum_<ij> (a+_i a_j + h.c.) on the PBC ring"""
    for i in range(L):
        j = (i + 1) % L
        for (a, b) in ((i, j), (j, i)):  # a+_a a_b
            if s[b] == 0 or s[a] >= nmax:
                continue
            new = list(s)
            amp = -t
            if stat == "BH":
                amp *= np.sqrt(s[b]) * np.sqrt(s[a] + 1)
            new[b] -= 1
            new[a] += 1
            if stat == "F":
                lo, hi = min(a, b), max(a, b)
                amp *= (-1) ** sum(s[lo + 1:hi])
            yield tuple(new), amp


def diag_energy(s, L, stat, V, U, stag):
    e = 0.0
    for i in range(L):
        if stat in ("F", "HCB"):
            e += V * s[i] * s[(i + 1) % L]
        else:
            e += 0.5 * U * s[i] * (s[i] - 1)
        e += stag * (-1) ** i * s[i]
    return e


def sector_spectra(L, N, stat, t=1.0, V=0.0, U=0.0, nmax=None, stag=0.0, nev=2):
    """lowest eigenvalues in each momentum sector (translation by `step` sites; step=1 or 2)"""
    if nmax is None:
        nmax = 1 if stat in ("F", "HCB") else N
    step = 2 if stag != 0 else 1
    Lc = L // step
    basis = fock_basis(L, N, nmax)

    def tr(s):
        sign = 1
        for _ in range(step):
            s, sg = translate(s, stat)
            sign *= sg
        return s, sign

    # orbits
    rep_of = {}
    reps = []
    for s in basis:
        if s in rep_of:
            continue
        orb = [(s, 1)]
        cur, sg = s, 1
        while True:
            cur, g = tr(cur)
            sg *= g
            if cur == s:
                break
            orb.append((cur, sg))
        r = min(o[0] for o in orb)
        # phases relative to representative r: state = sign * T^j r
        # recompute starting from r
        orb_r = {r: (0, 1)}
        cur, sg = r, 1
        j = 0
        while True:
            cur, g = tr(cur)
            sg *= g
            j += 1
            if cur == r:
                period, wrap_sign = j, sg
                break
            orb_r[cur] = (j, sg)
        for st, (jj, ss) in orb_r.items():
            rep_of[st] = (r, jj, ss)
        reps.append((r, period, wrap_sign))
    out = {}
    for m in range(Lc):
        k = 2 * np.pi * m / Lc
        # allowed reps: exp(-i k p) * wrap_sign == 1   (T^p r = wrap_sign r)
        allowed = [(r, p, w) for (r, p, w) in reps if abs(np.exp(1j * k * p) * w - 1) < 1e-9]
        if not allowed:
            continue
        idx = {r: n for n, (r, p, w) in enumerate(allowed)}
        norm = {r: np.sqrt(p) for (r, p, w) in allowed}
        rows, cols, vals = [], [], []
        for (r, p, w) in allowed:
            a = idx[r]
            rows.append(a); cols.append(a); vals.append(diag_energy(r, L, stat, V, U, stag))
            for new, amp in hop_terms(r, L, stat, t, nmax):
                rr, jj, ss = rep_of[new]
                if rr not in idx:
                    continue
                # |r,k> = (1/sqrt p) sum_j e^{-ikj} T^j r ; T^j rr = ss * new  =>  new = ss * T^jj rr
                b = idx[rr]
                rows.append(b); cols.append(a)
                vals.append(amp * ss * np.exp(1j * k * jj) * norm[r] / norm[rr])
        D = len(allowed)
        H = sps.csr_matrix((vals, (rows, cols)), shape=(D, D), dtype=complex)
        H = 0.5 * (H + H.getH())
        if D <= 400:
            ev = np.linalg.eigvalsh(H.toarray())[:nev]
        else:
            ev = np.sort(spla.eigsh(H, k=nev, which="SA", return_eigenvectors=False))
        out[m] = ev
    return out, Lc, step


def full_spectrum_check(L, N, stat, **kw):
    """sanity: union of sector spectra (lowest) vs direct full-Fock diagonalization ground energy"""
    nmax = kw.get("nmax") or (1 if stat in ("F", "HCB") else N)
    basis = fock_basis(L, N, nmax)
    ix = {s: n for n, s in enumerate(basis)}
    rows, cols, vals = [], [], []
    for s in basis:
        a = ix[s]
        rows.append(a); cols.append(a)
        vals.append(diag_energy(s, L, stat, kw.get("V", 0), kw.get("U", 0), kw.get("stag", 0)))
        for new, amp in hop_terms(s, L, stat, 1.0, nmax):
            rows.append(ix[new]); cols.append(a); vals.append(amp)
    H = sps.csr_matrix((vals, (rows, cols)), shape=(len(basis),) * 2)
    return np.sort(spla.eigsh(H, k=1, which="SA", return_eigenvectors=False))[0] if len(basis) > 50 \
        else np.linalg.eigvalsh(H.toarray())[0]


def gaps(L, N, stat, **kw):
    sec, Lc, step = sector_spectra(L, N, stat, **kw)
    E0 = min(v[0] for v in sec.values())
    m0 = min(sec, key=lambda m: sec[m][0])
    g = {}
    for m, v in sec.items():
        q = (m - m0) % Lc
        e = v[0] if m != m0 else v[1]
        g[q] = e - E0
    return g, Lc, step, E0, m0


def table(label, stat, N_of_L, Ls, **kw):
    print(f"\n--- {label} ---")
    rows = {}
    for L in Ls:
        N = N_of_L(L)
        g, Lc, step, E0, m0 = gaps(L, N, stat, **kw)
        nu = N / L
        qs = sorted(g)
        soft_pred = {(round(2 * np.pi * nu * mm * step / (2 * np.pi) * Lc)) % Lc for mm in range(1, 8)}
        # g*L listing
        gl = {q: g[q] * L for q in qs}
        # nontrivial soft candidates: q != 0 with g*L smallest
        nz = [q for q in qs if q != 0]
        best = sorted(nz, key=lambda q: gl[q])[:3]
        print(f"L={L:2d} N={N:2d} nu={nu:.3f} cell={step} K0={m0} E0={E0:.6f}")
        print("   q*cell/(2pi*nu) -> g*L : " + ", ".join(
            f"{(2*np.pi*q/Lc)/(2*np.pi*nu*step):.3f}:{gl[q]:.3f}" for q in qs))
        print(f"   lowest three q!=0 by g*L (in units of 2pi*nu): "
              + ", ".join(f"{(2*np.pi*q/Lc)/(2*np.pi*nu*step):.3f}" for q in best)
              + f"   predicted soft q indices {sorted(soft_pred)}")
        rows[L] = (gl, Lc)
    return rows


def velocity(L, N, stat, **kw):
    g, Lc, step, E0, m0 = gaps(L, N, stat, **kw)
    return g[1] * L / (2 * np.pi)  # E(2pi/L)/(2pi/L): sound velocity estimate


if __name__ == "__main__":
    # sanity: sector decomposition reproduces the full ground state
    for stat, L, N, kw in (("F", 8, 4, {}), ("HCB", 8, 3, {}), ("BH", 6, 3, {"U": 2.0, "nmax": 3}),
                           ("F", 8, 4, {"V": 1.0, "stag": 0.5})):
        g, Lc, step, E0, m0 = gaps(L, N, stat, **kw)
        print(f"check {stat} L={L} N={N} {kw}: sector E0={E0:.10f}  full E0={full_spectrum_check(L, N, stat, **kw):.10f}")

    # ---- core test: nu = 1/2 and 1/4 across statistics and interactions ----
    Ls_half = [8, 12, 16]
    Ls_quarter = [8, 12, 16, 20]
    R = {}
    R["F nu=1/2"] = table("free fermions nu=1/2", "F", lambda L: L // 2, Ls_half)
    R["F V=1 nu=1/2"] = table("fermions V=1 nu=1/2", "F", lambda L: L // 2, Ls_half, V=1.0)
    R["HCB nu=1/2"] = table("hard-core bosons nu=1/2", "HCB", lambda L: L // 2, Ls_half)
    R["BH U=1 nu=1/2"] = table("Bose-Hubbard U=1 (nmax=3) nu=1/2", "BH", lambda L: L // 2, [8, 10, 12], U=1.0, nmax=3)
    R["BH U=4 nu=1/2"] = table("Bose-Hubbard U=4 (nmax=3) nu=1/2", "BH", lambda L: L // 2, [8, 10, 12], U=4.0, nmax=3)
    R["F nu=1/4"] = table("free fermions nu=1/4", "F", lambda L: L // 4, Ls_quarter)
    R["HCB nu=1/4"] = table("hard-core bosons nu=1/4", "HCB", lambda L: L // 4, Ls_quarter)
    R["BH U=4 nu=1/4"] = table("Bose-Hubbard U=4 (nmax=3) nu=1/4", "BH", lambda L: L // 4, [8, 12, 16], U=4.0, nmax=3)
    R["F nu=1/3"] = table("free fermions nu=1/3", "F", lambda L: L // 3, [9, 12, 15, 18])
    R["HCB V=1 nu=1/3"] = table("HCB V=1 nu=1/3", "HCB", lambda L: L // 3, [9, 12, 15, 18], V=1.0)

    # ---- velocities are free ----
    print("\n--- velocity estimate v = E(2pi/L)*L/2pi at L=12, nu=1/2 (free quantity) ---")
    for lab, stat, kw in (("F", "F", {}), ("F V=1", "F", {"V": 1.0}), ("F V=-1", "F", {"V": -1.0}),
                          ("HCB", "HCB", {}), ("BH U=1", "BH", {"U": 1.0, "nmax": 3}),
                          ("BH U=4", "BH", {"U": 4.0, "nmax": 3})):
        print(f"   {lab:8s} v = {velocity(12, 6, stat, **kw):.4f}")

    # ---- hostile 1: staggered potential (translation -> period 2) ----
    print("\n--- hostile: staggered potential delta=0.5 (unit cell 2) ---")
    table("F nu=1/2 staggered (cell filling 1: integer)", "F", lambda L: L // 2, Ls_half, stag=0.5)
    table("HCB nu=1/2 staggered (cell filling 1)", "HCB", lambda L: L // 2, Ls_half, stag=0.5)
    table("F nu=1/4 staggered (cell filling 1/2: LSM still applies with doubled cell)", "F", lambda L: L // 4,
          [8, 12, 16, 20], stag=0.5)
