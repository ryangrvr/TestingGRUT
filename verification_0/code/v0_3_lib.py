"""
v0_3_lib.py -- independent exact-diagonalisation and free-particle utilities for the V0-3 / P-02 blind reproduction.

Conventions (fixed once, used everywhere):
  * sites j = 0..L-1 on a ring, H = sum_j sum_r h_r a_j^dag a_{j+r}   (indices mod L), h_{-r} = conj(h_r).
  * c_k = L^{-1/2} sum_j e^{-ikj} c_j  =>  H = sum_k eps(k) a_k^dag a_k with eps(k) = sum_r h_r e^{ikr}.
  * rho_q = sum_j e^{-iqj} n_j = sum_k a_{k-q}^dag a_k   (q>0 moves a particle from k to k-q).
    The spec writes sum_k a_{k+q}^dag a_k; this is the same physics in the opposite Fourier convention
    (eps(k) -> eps(-k)).  Irrelevant for reflection-symmetric bands, decisive for asymmetric ones.
  * S(q,w) weights |<m|rho_q|0>|^2 / N, summed over each degenerate eigenspace; support = weight > 1e-12.
"""
import itertools
import numpy as np

WTOL = 1e-12      # support threshold (spec)
ETOL = 1e-8       # degeneracy grouping tolerance


# ---------------------------------------------------------------- many-body bases
def basis(nsites, N, kind):
    if kind == 'boson':
        out = []
        for bars in itertools.combinations(range(N + nsites - 1), nsites - 1):
            occ, prev = [], -1
            for b in bars:
                occ.append(b - prev - 1)
                prev = b
            occ.append(N + nsites - 1 - prev - 1)
            out.append(tuple(occ))
        return out
    out = []
    for c in itertools.combinations(range(nsites), N):
        occ = [0] * nsites
        for s in c:
            occ[s] = 1
        out.append(tuple(occ))
    return out


def apply_hop(occ, i, l, kind):
    """a_i^dag a_l |occ>, i != l.  Returns (amplitude, new_occ) or None."""
    if occ[l] == 0:
        return None
    if kind == 'boson':
        n = list(occ)
        a = np.sqrt(n[l]); n[l] -= 1
        a *= np.sqrt(n[i] + 1); n[i] += 1
        return a, tuple(n)
    if occ[i] == 1:
        return None
    n = list(occ)
    if kind == 'hcb':
        n[l] = 0; n[i] = 1
        return 1.0, tuple(n)
    # fermion: canonical ordering by site index
    s1 = (-1) ** sum(n[:l]); n[l] = 0
    s2 = (-1) ** sum(n[:i]); n[i] = 1
    return float(s1 * s2), tuple(n)


def hamiltonian(nsites, N, kind, bonds):
    """bonds: list of (i, l, amp) meaning amp * a_i^dag a_l (i != l)."""
    B = basis(nsites, N, kind)
    idx = {b: t for t, b in enumerate(B)}
    H = np.zeros((len(B), len(B)), dtype=complex)
    for t, b in enumerate(B):
        for (i, l, amp) in bonds:
            r = apply_hop(b, i, l, kind)
            if r is None:
                continue
            a, nb = r
            H[idx[nb], t] += amp * a
    assert np.allclose(H, H.conj().T), "H not hermitian"
    return B, idx, H


def ring_bonds(L, hops):
    """hops: dict r -> h_r (must contain both r and -r)."""
    bonds = []
    for j in range(L):
        for r, h in hops.items():
            if h != 0:
                bonds.append((j, (j + r) % L, h))
    return bonds


def eps_from_hops(hops):
    return lambda k: np.real(sum(h * np.exp(1j * k * r) for r, h in hops.items()))


def rho_diag(B, xs, q):
    """diagonal of rho_q = sum_s e^{-i q x_s} n_s in basis B; xs[s] = position of site s."""
    xs = np.asarray(xs, dtype=float)
    occ = np.array(B, dtype=float)
    return occ @ np.exp(-1j * q * xs)


def group(evals, weights, tol=ETOL):
    """Group eigenvalues into degenerate eigenspaces, sum weights. Returns list of (E, W)."""
    if len(evals) == 0:
        return []
    order = np.argsort(evals)
    out = []
    cur_E, cur_W, first = None, 0.0, None
    for t in order:
        e = evals[t]
        if first is None or e - first > tol:
            if first is not None:
                out.append((cur_E, cur_W))
            first, cur_E, cur_W = e, e, 0.0
        cur_W += weights[t]
    out.append((cur_E, cur_W))
    return out


def ed_support(evals, evecs, psi0, E0, rhod, N):
    amp = evecs.conj().T @ (rhod * psi0)
    w = np.abs(amp) ** 2 / N
    g = group(evals - E0, w)
    return [(e, W) for e, W in g if W > WTOL]


# ---------------------------------------------------------------- free particles
def kgrid(L):
    return 2 * np.pi * np.arange(L) / L


def wrap(k):
    """map to (-pi, pi]"""
    k = np.mod(k + np.pi, 2 * np.pi) - np.pi
    return np.where(np.isclose(k, -np.pi), np.pi, k)


def fermion_occupation(eps, L, N):
    """indices m (k=2pi m/L) of the N lowest levels; returns (occ_bool, unique_flag, gap)"""
    k = kgrid(L)
    e = eps(k)
    order = np.argsort(e, kind='stable')
    occ = np.zeros(L, bool)
    occ[order[:N]] = True
    gap = (e[order[N]] - e[order[N - 1]]) if N < L else np.inf
    return occ, gap > 1e-10, gap


def ff_support(eps, L, occ, m, tol=ETOL):
    """free-fermion grouped support at q=2pi m/L (convention k -> k-q). occ: bool array over m-index."""
    k = kgrid(L)
    e = eps(k)
    N = occ.sum()
    if m % L == 0:
        return [(0.0, float(N))]
    ws = []
    for a in range(L):
        b = (a - m) % L
        if occ[a] and not occ[b]:
            ws.append(e[b] - e[a])
    ws = np.array(ws)
    return group(ws, np.ones(len(ws)) / N, tol)


def ff_lower_edge(eps, L, occ, m):
    k = kgrid(L)
    e = eps(k)
    a = np.nonzero(occ)[0]
    b = (a - m) % L
    ok = ~occ[b]
    if not ok.any():
        return np.nan
    return np.min(e[b[ok]] - e[a[ok]])


def fb_support(eps, L, nocc, m, tol=ETOL):
    """free-boson support from a Fock state with occupation numbers nocc[m] (any state: 'anything admissible')."""
    k = kgrid(L)
    e = eps(k)
    N = nocc.sum()
    if m % L == 0:
        return [(0.0, float(N))]
    ws, wt = [], []
    for a in range(L):
        if nocc[a] > 0:
            b = (a - m) % L
            ws.append(e[b] - e[a]); wt.append(nocc[a] * (nocc[b] + 1) / N)
    return group(np.array(ws), np.array(wt), tol)


def same_support(s1, s2, tol=1e-8, check_weights=True):
    if len(s1) != len(s2):
        return False
    for (e1, w1), (e2, w2) in zip(sorted(s1), sorted(s2)):
        if abs(e1 - e2) > tol:
            return False
        if check_weights and abs(w1 - w2) > 1e-8:
            return False
    return True


# ---------------------------------------------------------------- thermodynamic lower edge (continuum)
def _norm_intervals(ivs):
    """intervals on circle [0,2pi): list of (a,b) with a<b, after splitting."""
    out = []
    for a, b in ivs:
        a0 = np.mod(a, 2 * np.pi)
        w = b - a
        if w >= 2 * np.pi - 1e-15:
            return [(0.0, 2 * np.pi)]
        if a0 + w <= 2 * np.pi:
            out.append((a0, a0 + w))
        else:
            out.append((a0, 2 * np.pi)); out.append((0.0, a0 + w - 2 * np.pi))
    return out


def _complement(ivs):
    ivs = sorted(_norm_intervals(ivs))
    out, cur = [], 0.0
    for a, b in ivs:
        if a > cur:
            out.append((cur, a))
        cur = max(cur, b)
    if cur < 2 * np.pi:
        out.append((cur, 2 * np.pi))
    return out


def _intersect(A, B):
    out = []
    for a1, b1 in _norm_intervals(A):
        for a2, b2 in _norm_intervals(B):
            lo, hi = max(a1, a2), min(b1, b2)
            if hi > lo:
                out.append((lo, hi))
    return out


def thermo_lower_edge(eps, O, q, ns=801):
    """inf over k in O with k-q not in O of eps(k-q)-eps(k).  O: list of (a,b) intervals (occupied set)."""
    from scipy.optimize import minimize_scalar
    U = _complement(O)
    Uq = [(a + q, b + q) for a, b in U]          # k-q in U  <=>  k in U+q
    adm = _intersect(O, Uq)
    best = np.inf
    f = lambda k: eps(k - q) - eps(k)
    for a, b in adm:
        ks = np.linspace(a, b, ns)
        v = f(ks)
        t = int(np.argmin(v))
        best = min(best, v[t])
        lo, hi = ks[max(t - 1, 0)], ks[min(t + 1, ns - 1)]
        if hi > lo:
            r = minimize_scalar(f, bounds=(lo, hi), method='bounded', options={'xatol': 1e-14})
            best = min(best, r.fun)
    return best


def occupied_intervals(eps, mu, ngrid=200001):
    """{k in [-pi,pi): eps(k) < mu} as intervals (roots refined by brentq)."""
    from scipy.optimize import brentq
    ks = np.linspace(-np.pi, np.pi, ngrid)
    v = eps(ks) - mu
    ivs, start = [], None
    if v[0] < 0:
        start = -np.pi
    for t in range(ngrid - 1):
        if v[t] < 0 <= v[t + 1]:
            r = brentq(lambda x: eps(x) - mu, ks[t], ks[t + 1], xtol=1e-15)
            ivs.append((start, r)); start = None
        elif v[t] >= 0 > v[t + 1]:
            start = brentq(lambda x: eps(x) - mu, ks[t], ks[t + 1], xtol=1e-15)
    if start is not None:
        if ivs and np.isclose(ivs[0][0], -np.pi):
            a0, b0 = ivs.pop(0)
            ivs.append((start, b0 + 2 * np.pi))
        else:
            ivs.append((start, np.pi))
    return ivs


def soft_set(eps, O, nq=3001, thr=1e-7):
    """Q_soft representatives in [0,pi] (quotient q~-q~q+2pi) from the thermodynamic lower edge on q in (0,pi]."""
    from scipy.optimize import minimize_scalar
    qs = np.linspace(1e-6, np.pi, nq)
    w = np.array([thermo_lower_edge(eps, O, q) for q in qs])
    soft = []
    if w[0] < 1e-4:
        soft.append(0.0)
    for t in range(1, nq):
        left = w[t - 1]
        right = w[t + 1] if t + 1 < nq else np.inf
        if w[t] <= left and w[t] <= right and t > 1:
            lo, hi = qs[t - 1], qs[min(t + 1, nq - 1)]
            r = minimize_scalar(lambda x: thermo_lower_edge(eps, O, x), bounds=(lo, hi), method='bounded',
                                options={'xatol': 1e-12})
            if r.fun < thr:
                rep = min(np.mod(r.x, 2 * np.pi), 2 * np.pi - np.mod(r.x, 2 * np.pi))
                if all(abs(rep - s) > 1e-5 for s in soft):
                    soft.append(rep)
    return sorted(soft), qs, w


def local_exponent(f, qs):
    """d log f / d log q at each q (central, multiplicative step 1.01) and log f/log q."""
    out = []
    for q in qs:
        a, b = f(q / 1.01), f(q * 1.01)
        out.append((q, f(q), np.log(b / a) / np.log(1.01 ** 2), np.log(f(q)) / np.log(q)))
    return out


def nearest_m(L, q):
    """nearest grid index to q; ties go to the lower m (spec)."""
    x = q * L / (2 * np.pi)
    lo = int(np.floor(x))
    return lo if (x - lo) <= (lo + 1 - x) + 1e-12 else lo + 1
