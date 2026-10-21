"""Verifier v_residual_bias, part 1: independent zero-momentum width of the soft phi mode.

Independent of the module code: different kinematic variables (E_b, E_c instead of E_c, E_d), statistical
factor written out as f_b(1+f_c)(1+f_d) - (1+f_b) f_c f_d (no KMS simplification), own vertex derivation, own
thermal masses.  numpy/scipy only.

Width of a phi quantum at rest, a + b -> c + d, p_a = 0, p_d = p_b - p_c, E_d = m_a + E_b - E_c:
  Gamma = sum_ch mult |M|^2 S /(2 m_a) Int d^3p_b/((2pi)^3 2E_b) d^3p_c/((2pi)^3 2E_c) (2pi) delta(E) /(2E_d) F
        = sum_ch mult |M|^2 S /(64 pi^3 m_a) Int dE_b dE_c F Theta,
Theta: |cos(theta_bc)| <= 1, cos = (p_b^2 + p_c^2 + m_d^2 - E_d^2)/(2 p_b p_c).

Vertices from U (real basis X = (q1 + i q2)/sqrt2):
  (lambda/4) phi^4                     -> d^4/dphi^4 = 6 lambda
  g^2 phi^2 |X|^2 = (g^2/2) phi^2 q_i^2 -> d^4/dphi^2 dq_i^2 = 2 g^2
  (h^2/2) phi^2 D^2                     -> 2 h^2
Thermal masses at phi = 0 (leading order, mu = 0):
  m_phi^2 = -lambda v^2 + Pi_phi T^2, Pi_phi = lambda/4 + g^2/6 + h^2/12 (written as Pi_phi (T^2 - Tc^2))
  m_q^2 = (g^2/12 + lambda_X/3) T^2;  m_D^2 = m_D0^2 + (h^2/12 + lambda_D/4) T^2
"""
import json
import math
import os
import time

import numpy as np
from scipy.integrate import quad

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(OUT, exist_ok=True)

v, lam, g2, lamX, h2, lamD, mD0 = 1000.0, 0.01, 0.01, 0.001, 1e-4, 0.001, 100.0
PI_P = lam / 4 + g2 / 6 + h2 / 12
PI_Q = g2 / 12 + lamX / 3
PI_D = h2 / 12 + lamD / 4
TC = 1557.5            # module's adopted T_c (resummed range 1557.1-1564.2 GeV, quoted)
TC_LO = v * math.sqrt(lam / PI_P)

# (name, b, c, d, |M|^2, S, multiplicity)
CH = [("pp->pp", "p", "p", "p", (6 * lam) ** 2, 0.5, 1),
      ("pq->pq", "q", "p", "q", (2 * g2) ** 2, 1.0, 2),
      ("pp->qq", "p", "q", "q", (2 * g2) ** 2, 0.5, 2),
      ("pD->pD", "D", "p", "D", (2 * h2) ** 2, 1.0, 1),
      ("pp->DD", "p", "D", "D", (2 * h2) ** 2, 0.5, 1)]


def masses_over_T(T, Tc=TC):
    mp2 = PI_P * (T * T - Tc * Tc)
    return dict(p=math.sqrt(max(mp2, 0.0)) / T, q=math.sqrt(PI_Q), D=math.sqrt(mD0 ** 2 + PI_D * T * T) / T)


def nB(E):
    return 1.0 / math.expm1(E)


def Fstat(Eb, Ec, Ed, mode):
    if mode == "full":
        fb, fc, fd = nB(Eb), nB(Ec), nB(Ed)
        return fb * (1 + fc) * (1 + fd) - (1 + fb) * fc * fd
    if mode == "classical":
        fb, fc, fd = 1 / Eb, 1 / Ec, 1 / Ed
        return fb * fc + fb * fd - fc * fd
    if mode == "nsv":
        return nB(Eb)
    raise ValueError


def Ec_window(Eb, ma, mb, mc, md):
    """Allowed E_c for fixed E_b: quadratic (alpha^2 - pb^2) Ec^2 + 2 alpha beta Ec + beta^2 + pb^2 mc^2 <= 0."""
    pb2 = Eb * Eb - mb * mb
    c1 = 0.5 * (md * md - ma * ma - mb * mb - mc * mc)
    al = Eb + ma
    be = c1 - ma * Eb
    A = al * al - pb2
    B = 2 * al * be
    C = be * be + pb2 * mc * mc
    disc = B * B - 4 * A * C
    if disc <= 0:
        return None
    s = math.sqrt(disc)
    lo, hi = (-B - s) / (2 * A), (-B + s) / (2 * A)
    lo = max(lo, mc)
    hi = min(hi, ma + Eb - md)
    return (lo, hi) if hi > lo else None


def chan_int(ma, mb, mc, md, mode, rtol=1e-8):
    def inner(Eb):
        w = Ec_window(Eb, ma, mb, mc, md)
        if w is None:
            return 0.0
        f = lambda Ec: Fstat(Eb, Ec, ma + Eb - Ec, mode)
        return quad(f, w[0], w[1], epsabs=0, epsrel=rtol, limit=200)[0]
    # split the outer range: soft region resolved on the scale of the masses
    pts = sorted(set([mb, mb + 0.1 * ma + 1e-12, mb + ma, mb + 10 * ma, 1.0, 5.0, 60.0]))
    pts = [p for p in pts if p >= mb]
    tot = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        tot += quad(inner, a, b, epsabs=0, epsrel=rtol, limit=400)[0]
    if mode == "classical":      # f = 1/E has a power-law tail: integrate the outer variable to infinity
        tot += quad(inner, pts[-1], np.inf, epsabs=0, epsrel=rtol, limit=400)[0]
    return tot


def width_over_T(T, mode="full", Tc=TC, per=False):
    x = masses_over_T(T, Tc)
    rows = {}
    for nm, b, c, d, M2, S, mult in CH:
        rows[nm] = mult * M2 * S / (64 * math.pi ** 3 * x["p"]) * chan_int(x["p"], x[b], x[c], x[d], mode)
    tot = sum(rows.values())
    return (tot, rows, x) if per else tot


def main():
    t0 = time.time()
    out = {}
    # (1) Parwani check: one field, all masses m; classical integral Int F_cl = pi^2/6 (i.e. Int 1/(EbEcEd) =
    # pi^2/(6m)).  Then Gamma = |M|^2 S T^2/(384 pi m) and for M = g^2 (g^2 phi^4/4!) S=1/2:
    # Gamma = g^4 T^2/(768 pi m) = 2 * g^4 T^2/(1536 pi m)  (Parwani gamma = Gamma/2).
    par = {}
    for m in (0.003, 0.03, 0.3, 1.0):
        I = chan_int(m, m, m, m, "classical")
        par[str(m)] = dict(I=I, pi2_over_6=math.pi ** 2 / 6, rel=I / (math.pi ** 2 / 6) - 1)
    gg = 6 * lam
    mT = math.sqrt(gg / 24)
    If = chan_int(mT, mT, mT, mT, "full")
    Gf = gg ** 2 * 0.5 / (64 * math.pi ** 3 * mT) * If
    par["pure_phi4_full_over_LO"] = Gf / (2 * gg ** 2 / (1536 * math.pi * mT))
    out["parwani_check"] = par
    # (2) high-T and table
    tab = []
    for eps in (6.917e-4, 2e-3, 5e-3, 7e-3, 1e-2, 3e-2, 0.08, 0.3, 0.93, 3.0, 10.0, 31.1, 100.0, 1e3):
        T = TC * (1 + eps)
        gf, rows, x = width_over_T(T, "full", per=True)
        gc = width_over_T(T, "classical")
        gn = width_over_T(T, "nsv")
        tab.append(dict(eps=eps, mphi_over_T=x["p"], full=gf, classical=gc, nsv=gn, Gamma_over_mphi=gf / x["p"],
                        channels=rows))
        print(eps, x["p"], gf, gc, gn, flush=True)
    out["table"] = tab
    # one-to-two kinematics: m_phi/(2 m_q) at T -> inf
    out["mphi_over_2mq_highT"] = math.sqrt(PI_P) / (2 * math.sqrt(PI_Q))
    # sensitivity of the high-T width to T_c choice (LO zero vs adopted): irrelevant at T >> Tc
    out["TC_LO"] = TC_LO
    out["runtime_s"] = time.time() - t0
    with open(os.path.join(OUT, "width.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out["parwani_check"], indent=1))


if __name__ == "__main__":
    main()
