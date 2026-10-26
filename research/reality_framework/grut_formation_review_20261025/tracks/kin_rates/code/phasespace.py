#!/usr/bin/env python3
"""kin_rates (Claude, own code): Lorentz-invariant n-body phase space by sequential two-body splitting.

Measure: dPhi_n(P) = prod_i d^3p_i/((2pi)^3 2E_i) (2pi)^4 delta^4(P - sum p_i), so that Phi_2 = |p*|/(4 pi sqrt s).
Recursion: dPhi_n(P; m_1..m_n) = dPhi_2(P; q, m_n) (ds_q/2pi) dPhi_{n-1}(q; m_1..m_{n-1}), applied as a binary tree
    n=2: P -> (1,2)
    n=3: P -> (q12, 3), q12 -> (1,2)
    n=4: P -> (q12, q34), q12 -> (1,2), q34 -> (3,4)
with the invariant masses sampled either uniformly in their kinematic ranges or with a tan (Breit-Wigner)
map concentrated at a chosen mass^2 (importance sampling of a resonance in the (1,2) pair); isotropic decay angles.
Returns CM-frame momenta (N, n, 4) and weights w (N,) such that mean(w f) estimates int dPhi_n f.
Checks (main): massless Phi_3, Phi_4 against (2pi)^{4-3n} (pi/2)^{n-1} s^{n-2}/((n-1)!(n-2)!); massive Phi_3 against
a direct Dalitz-plot quadrature; invariance of the weighted integral under the choice of map.
"""
import math
import numpy as np
from amp import kallen


def two_body(P, m1, m2, rng, N):
    """Decay P (N,4) -> p1, p2 isotropic in the rest frame of P, boosted back; returns p1, p2, |p*|/(4 pi M)."""
    M2 = P[:, 0]**2 - P[:, 1]**2 - P[:, 2]**2 - P[:, 3]**2
    Mm = np.sqrt(M2)
    pst = np.sqrt(np.clip(kallen(M2, m1 * m1, m2 * m2), 0, None)) / (2 * Mm)
    c = rng.uniform(-1, 1, N); ph = rng.uniform(0, 2 * np.pi, N)
    sn = np.sqrt(1 - c * c)
    E1 = (M2 + m1 * m1 - m2 * m2) / (2 * Mm)
    p1 = np.stack([E1, pst * sn * np.cos(ph), pst * sn * np.sin(ph), pst * c], 1)
    p1 = boost(p1, P)
    p2 = P - p1
    return p1, p2, pst / (4 * np.pi * Mm)


def boost(p, P):
    """Boost p (N,4) from the rest frame of P (N,4) to the frame where P is given."""
    M = np.sqrt(np.clip(P[:, 0]**2 - P[:, 1]**2 - P[:, 2]**2 - P[:, 3]**2, 1e-300, None))
    b = P[:, 1:] / P[:, [0]]
    g = P[:, 0] / M
    bp = np.sum(b * p[:, 1:], 1)
    b2 = np.sum(b * b, 1)
    gam2 = np.where(b2 > 0, (g - 1) / np.where(b2 > 0, b2, 1), 0.0)
    E = g * (p[:, 0] + bp)
    pv = p[:, 1:] + b * (gam2 * bp + g * p[:, 0])[:, None]
    return np.column_stack([E, pv])


def sample_mass2(lo, hi, rng, N, peak=None, width=None, mix=0.5):
    """Sample s in [lo, hi]: uniform (weight = hi-lo), or a mixture (fraction mix) of a tan (Breit-Wigner) map
    around `peak` with scale `width` and (1-mix) uniform; weight = 1/density (multichannel importance sampling)."""
    lo = lo * np.ones(N); hi = hi * np.ones(N)
    if peak is None:
        s = lo + (hi - lo) * rng.uniform(0, 1, N)
        return s, (hi - lo)
    ua = np.arctan((lo - peak) / width); ub = np.arctan((hi - peak) / width)
    use_tan = rng.uniform(0, 1, N) < mix
    u = ua + (ub - ua) * rng.uniform(0, 1, N)
    s_tan = peak + width * np.tan(u)
    s_flat = lo + (hi - lo) * rng.uniform(0, 1, N)
    s = np.where(use_tan, s_tan, s_flat)
    p_tan = width / ((ub - ua) * ((s - peak)**2 + width * width))
    p_flat = 1.0 / (hi - lo)
    return s, 1.0 / (mix * p_tan + (1 - mix) * p_flat)


def generate(sqrt_s, masses, N, rng, peak12=None, width12=None, peak34=None, width34=None):
    """n-body phase space at total momentum (sqrt_s, 0) in the CM.  Returns (mom (N,n,4), w (N,))."""
    n = len(masses)
    s = sqrt_s**2
    P = np.zeros((N, 4)); P[:, 0] = sqrt_s
    if n == 2:
        p1, p2, w = two_body(P, masses[0], masses[1], rng, N)
        return np.stack([p1, p2], 1), w * np.ones(N)
    if n == 3:
        lo = (masses[0] + masses[1])**2; hi = (sqrt_s - masses[2])**2
        s12, w12 = sample_mass2(lo, hi, rng, N, peak12, width12)
        q, p3, wA = two_body(P, np.sqrt(s12), masses[2], rng, N)
        p1, p2, wB = two_body(q, masses[0], masses[1], rng, N)
        return np.stack([p1, p2, p3], 1), w12 / (2 * np.pi) * wA * wB
    if n == 4:
        lo12 = (masses[0] + masses[1])**2; hi12 = (sqrt_s - masses[2] - masses[3])**2
        s12, w12 = sample_mass2(lo12, hi12, rng, N, peak12, width12)
        lo34 = (masses[2] + masses[3])**2 * np.ones(N); hi34 = (sqrt_s - np.sqrt(s12))**2
        s34, w34 = sample_mass2(lo34, hi34, rng, N, peak34, width34)
        q12, q34, wA = two_body(P, np.sqrt(s12), np.sqrt(s34), rng, N)
        p1, p2, wB = two_body(q12, masses[0], masses[1], rng, N)
        p3, p4, wC = two_body(q34, masses[2], masses[3], rng, N)
        return np.stack([p1, p2, p3, p4], 1), w12 * w34 / (2 * np.pi)**2 * wA * wB * wC
    raise ValueError('n must be 2, 3 or 4')


def volume_massless(n, s):
    return (2 * np.pi)**(4 - 3 * n) * (np.pi / 2)**(n - 1) * s**(n - 2) / (math.factorial(n - 1) * math.factorial(n - 2))


def phi3_dalitz(sqrt_s, m1, m2, m3, n=400):
    """Direct Phi_3 = (1/(2pi)^3) (1/(16 s))... use dPhi_3 = (1/(2 pi)^3) (1/8) ds12 ds23 / (4 s)  -> 1/(128 pi^3 s) int ds12 ds23."""
    s = sqrt_s**2
    lo, hi = (m1 + m2)**2, (sqrt_s - m3)**2
    tot = 0.0
    s12s = np.linspace(lo, hi, n + 1); s12s = 0.5 * (s12s[1:] + s12s[:-1]); d = (hi - lo) / n
    for s12 in s12s:
        E2 = (s12 - m1 * m1 + m2 * m2) / (2 * np.sqrt(s12)); E3 = (s - s12 - m3 * m3) / (2 * np.sqrt(s12))
        p2 = np.sqrt(max(E2 * E2 - m2 * m2, 0)); p3 = np.sqrt(max(E3 * E3 - m3 * m3, 0))
        s23_max = (E2 + E3)**2 - (p2 - p3)**2; s23_min = (E2 + E3)**2 - (p2 + p3)**2
        tot += (s23_max - s23_min) * d
    return tot / (128 * np.pi**3 * s)


if __name__ == '__main__':
    rng = np.random.default_rng(11)
    sq = 1.0
    for n in (3, 4):
        mom, w = generate(sq, [0.0] * n, 400000, rng)
        print(f'massless Phi_{n}: MC {w.mean():.6e} +- {w.std()/np.sqrt(len(w)):.1e}  exact {volume_massless(n, sq*sq):.6e}')
        # momentum conservation and on-shell
        P = mom.sum(1); print('   sum p - P max:', np.abs(P - np.array([sq, 0, 0, 0])).max())
    ms = [0.1, 0.25, 0.3]
    mom, w = generate(sq, ms, 400000, rng)
    print(f'massive Phi_3: MC {w.mean():.6e} +- {w.std()/np.sqrt(len(w)):.1e}  Dalitz {phi3_dalitz(sq, *ms):.6e}')
    mom, w = generate(sq, ms, 400000, rng, peak12=0.2, width12=0.01)
    print(f'massive Phi_3 with tan map: MC {w.mean():.6e} +- {w.std()/np.sqrt(len(w)):.1e}')
    m2 = (mom[..., 0]**2 - (mom[..., 1:]**2).sum(-1))
    print('   on-shell check max |p^2 - m^2|:', np.abs(m2 - np.array(ms)**2).max())
    ms4 = [0.1, 0.1, 0.2, 0.15]
    mom, w = generate(sq, ms4, 400000, rng)
    mom2, w2 = generate(sq, ms4, 400000, rng, peak12=0.3, width12=0.02, peak34=0.4, width34=0.05)
    print(f'massive Phi_4: flat {w.mean():.6e} +- {w.std()/np.sqrt(len(w)):.1e}; mapped {w2.mean():.6e} +- {w2.std()/np.sqrt(len(w2)):.1e}')
