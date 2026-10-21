"""Tree-level scalar amplitudes of the unchanged action in the +v vacuum (Berends-Giele recursion).

Every tree diagram built from the cubic tensor L3 = d^3U and quartic tensor L4 = d^4U contributes
M = - prod(lambda) prod 1/(p^2 - m^2) (mostly-minus metric), so the recursion
  J_a(P) = [sum_{2-partitions} L3_abc J_b J_c + sum_{3-partitions} L4_abcd J_b J_c J_d] / (p_P^2 - m_a^2)
and M = - e^{(n)}_a [same bracket for P = all but the last leg] reproduce the full tree amplitude.
All legs are incoming; an outgoing particle is an incoming antiparticle with reversed momentum.
"""
import itertools
import math

import numpy as np

from common import MASS, vertex_tensors

L3, L4 = vertex_tensors()
S2 = 1.0 / math.sqrt(2.0)
POL = {  # <0| q_a |particle incoming>
    "h": np.array([1, 0, 0, 0], dtype=complex),
    "X": np.array([0, S2, -1j * S2, 0], dtype=complex),
    "Xb": np.array([0, S2, 1j * S2, 0], dtype=complex),
    "D": np.array([0, 0, 0, 1], dtype=complex),
}
CONJ = {"h": "h", "D": "D", "X": "Xb", "Xb": "X"}
PMASS = {"h": MASS[0], "X": MASS[1], "Xb": MASS[1], "D": MASS[3]}
CHARGE = {"h": 0, "D": 0, "X": 1, "Xb": -1}
Z2D = {"h": 0, "D": 1, "X": 0, "Xb": 0}


def mdot(a, b):
    return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]


def _partitions2(S):
    low = S & -S
    sub = (S - 1) & S
    while sub:
        if sub & low:
            yield sub, S ^ sub
        sub = (sub - 1) & S


def _partitions3(S):
    """Unordered partitions of bitmask S into three nonempty blocks (first block holds the lowest bit)."""
    low = S & -S
    blocks = []
    rest = S ^ low
    extra = rest
    while True:
        A = low | extra
        R = S ^ A
        if R:
            for B, C in _partitions2(R):
                blocks.append((A, B, C))
        if extra == 0:
            break
        extra = (extra - 1) & rest
    return blocks


def amplitude(types, moms, L3=L3, L4=L4, mass=MASS):
    """All-incoming tree amplitude. types: list of particle names; moms: list of 4-momenta (incoming)."""
    n = len(types)
    pols = [POL[t] for t in types]
    moms = [np.asarray(p, dtype=float) for p in moms]
    m2 = mass ** 2
    J = {}
    psum = {}
    for i in range(n - 1):
        J[1 << i] = pols[i]
        psum[1 << i] = moms[i]
    full = (1 << (n - 1)) - 1
    # build currents by increasing subset size
    subsets = sorted(range(1, full + 1), key=lambda s: bin(s).count("1"))

    def bracket(S):
        acc = np.zeros(4, dtype=complex)
        for A, B in _partitions2(S):
            acc += np.einsum("abc,b,c->a", L3, J[A], J[B])
        for A, B, C in _partitions3(S):
            acc += np.einsum("abcd,b,c,d->a", L4, J[A], J[B], J[C])
        return acc

    for S in subsets:
        if S in J:
            continue
        A = S & -S
        psum[S] = psum[A] + psum[S ^ A]
        if S == full:
            continue
        p2 = mdot(psum[S], psum[S])
        J[S] = bracket(S) / (p2 - m2)
    return -np.dot(pols[n - 1], bracket(full))


# ----------------------------------------------------------------------------- kinematics
def two_body_final(sqrt_s, m1, m2, nhat):
    """CM-frame momenta of a two-body final state along unit vector nhat."""
    s = sqrt_s ** 2
    lam = (s - (m1 + m2) ** 2) * (s - (m1 - m2) ** 2)
    p = math.sqrt(max(lam, 0.0)) / (2 * sqrt_s)
    E1 = math.sqrt(p * p + m1 * m1)
    E2 = math.sqrt(p * p + m2 * m2)
    k1 = np.array([E1, *(p * nhat)])
    k2 = np.array([E2, *(-p * nhat)])
    return k1, k2, p


def amp_process(initial, final, pin, pout):
    """Physical amplitude for initial -> final, momenta given as physical (positive-energy) vectors."""
    types = list(initial) + [CONJ[t] for t in final]
    moms = list(pin) + [-k for k in pout]
    return amplitude(types, moms)


def sym_factor(lst):
    f = 1
    for t in set(lst):
        f *= math.factorial(list(lst).count(t))
    return f


def allowed(initial, final):
    if sum(CHARGE[t] for t in initial) != sum(CHARGE[t] for t in final):
        return False
    if sum(Z2D[t] for t in initial) % 2 != sum(Z2D[t] for t in final) % 2:
        return False
    return True


def sigma_v2_threshold(initial, final):
    """3->2 at threshold (all initial at rest). Returns |M|^2, p_f, sigma v^2 [GeV^-5] incl. 1/S_f."""
    ms = [PMASS[t] for t in initial]
    sqrt_s = sum(ms)
    mf = [PMASS[t] for t in final]
    if sqrt_s <= sum(mf):
        return None
    pin = [np.array([m, 0.0, 0.0, 0.0]) for m in ms]
    k1, k2, p = two_body_final(sqrt_s, mf[0], mf[1], np.array([0.0, 0.0, 1.0]))
    M = amp_process(initial, final, pin, [k1, k2])
    M2 = abs(M) ** 2
    sv2 = M2 / (8.0 * ms[0] * ms[1] * ms[2]) * p / (4 * math.pi * sqrt_s) / sym_factor(final)
    return dict(M2=M2, p_f=p, sqrt_s=sqrt_s, sigma_v2=sv2)


def sigma_22(initial, final, sqrt_s, nq=48):
    """2->2 total cross section [GeV^-2] in the CM frame, with 1/S_f for identical final particles."""
    m1, m2 = PMASS[initial[0]], PMASS[initial[1]]
    m3, m4 = PMASS[final[0]], PMASS[final[1]]
    if sqrt_s <= max(m1 + m2, m3 + m4):
        return 0.0
    s = sqrt_s ** 2
    pi_ = math.sqrt((s - (m1 + m2) ** 2) * (s - (m1 - m2) ** 2)) / (2 * sqrt_s)
    E1 = math.sqrt(pi_ ** 2 + m1 ** 2)
    E2 = math.sqrt(pi_ ** 2 + m2 ** 2)
    p1 = np.array([E1, 0.0, 0.0, pi_])
    p2 = np.array([E2, 0.0, 0.0, -pi_])
    x, w = np.polynomial.legendre.leggauss(nq)
    acc = 0.0
    pf = None
    for c, wt in zip(x, w):
        sn = math.sqrt(max(0.0, 1 - c * c))
        k1, k2, pf = two_body_final(sqrt_s, m3, m4, np.array([sn, 0.0, c]))
        M = amp_process(initial, final, [p1, p2], [k1, k2])
        acc += wt * abs(M) ** 2
    # dsigma/dOmega = |M|^2 p_f / (64 pi^2 s p_i); integrate over phi -> 2 pi
    return acc * 2 * math.pi * pf / (64 * math.pi ** 2 * s * pi_) / sym_factor(final)


def boost_to_cm(moms):
    P = np.sum(moms, axis=0)
    M = math.sqrt(mdot(P, P))
    beta = P[1:] / P[0]
    b2 = float(beta @ beta)
    if b2 < 1e-30:
        return [m.copy() for m in moms], M
    gam = 1.0 / math.sqrt(1 - b2)
    out = []
    for p in moms:
        bp = float(beta @ p[1:])
        E = gam * (p[0] - bp)
        vec = p[1:] + ((gam - 1) * bp / b2 - gam * p[0]) * beta
        out.append(np.array([E, *vec]))
    return out, M


def _mb_sampler(m, T, rng, size):
    """Exact (gridded inverse-CDF) sampling of |p| from the relativistic MB weight p^2 exp(-(E-m)/T)."""
    pmax = math.sqrt((m + 60 * T) ** 2 - m * m)
    grid = np.linspace(0.0, pmax, 20001)
    E = np.sqrt(grid ** 2 + m * m)
    pdf = grid ** 2 * np.exp(-(E - m) / T)
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(grid))])
    cdf /= cdf[-1]
    p = np.interp(rng.random(size), cdf, grid)
    u = rng.normal(size=(size, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    vec = p[:, None] * u
    return np.column_stack([np.sqrt(p * p + m * m), vec])


def mc_sigma_v2(initial, final, T, nsample=4000, seed=1):
    """Thermal average <(1/(8E1E2E3)) int dPhi2 |M|^2 / S_f> over relativistic MB momenta (mu = 0).

    Returns (mean, standard error) in GeV^-5. The flux factor is evaluated in the plasma frame;
    the two-body phase space and |M|^2 in the CM frame of the initial triple (invariant)."""
    rng = np.random.default_rng(seed)
    ms = [PMASS[t] for t in initial]
    mf = [PMASS[t] for t in final]
    samples = [_mb_sampler(m, T, rng, nsample) for m in ms]
    vals = np.zeros(nsample)
    for i in range(nsample):
        pin = [samples[j][i] for j in range(3)]
        cm, sqrt_s = boost_to_cm(pin)
        if sqrt_s <= sum(mf):
            continue
        u = rng.normal(size=3)
        u /= np.linalg.norm(u)
        k1, k2, pf = two_body_final(sqrt_s, mf[0], mf[1], u)
        M = amp_process(initial, final, cm, [k1, k2])
        flux = 1.0 / (8 * pin[0][0] * pin[1][0] * pin[2][0])
        vals[i] = flux * abs(M) ** 2 * pf / (4 * math.pi * sqrt_s) / sym_factor(final)
    return float(vals.mean()), float(vals.std() / math.sqrt(nsample))


def enumerate_processes():
    """All charge- and Z2_D-conserving 2->2 conversions and 3->2 processes open at threshold."""
    types = ["X", "Xb", "h", "D"]
    conv, three = [], []
    for ini in itertools.combinations_with_replacement(types, 2):
        for fin in itertools.combinations_with_replacement(types, 2):
            if sorted(ini) == sorted(fin) or not allowed(ini, fin):
                continue
            # number-conserving conversion between different species pairs; keep heavier-initial direction
            mi = sum(PMASS[t] for t in ini)
            mf = sum(PMASS[t] for t in fin)
            if mi > mf or (mi == mf and ini < fin):
                # elastic-equivalent (same species content) excluded above; X Xb <-> X Xb trivially excluded
                key = (tuple(ini), tuple(fin))
                if set(ini) <= {"X", "Xb"} and set(fin) <= {"X", "Xb"}:
                    continue
                conv.append(key)
    for ini in itertools.combinations_with_replacement(types, 3):
        for fin in itertools.combinations_with_replacement(types, 2):
            if not allowed(ini, fin):
                continue
            if sum(PMASS[t] for t in ini) <= sum(PMASS[t] for t in fin):
                continue
            three.append((tuple(ini), tuple(fin)))
    return conv, three
