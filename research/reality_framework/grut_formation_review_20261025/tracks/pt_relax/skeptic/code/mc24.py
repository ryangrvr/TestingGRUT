"""pt_relax skeptic (Claude, own code): Monte Carlo of the thermal 2<->4 chemical-equilibration rate.

Purpose: calibrate pt_relax's parametric bracket Gamma_{2<->4} = c V^4 T with c = (10 - 1e3)/(4 pi)^5
(3.2e-5 - 3.2e-3), V = Feynman vertex.  Independent method: explicit tree-level collision integral.

Model: one real scalar, L_int = -(lambda/4!) phi^4 (vertex V = lambda), thermal (Hartree) mass m on every line,
Bose statistics at mu = 0, temperature T = 1.  Tree 2 -> 4 amplitude: two vertices, one propagator, all 10
partitions of the six legs into two triples:  M = lambda^2 sum_{triples} 1/(q^2 - m^2),  q = sum of signed momenta.
With all external masses m, |q^2 - m^2| >= m^2 (the 1 -> 3 splitting is closed), so there is no on-shell pole.

  R_{2->4} = (1/2!) (1/4!) Int dPi_1 dPi_2 f_1 f_2  Int dPhi_4 |M|^2 prod_{i=3..6} (1 + f_i)        (events/vol/time)
  dn/dt    = 2 (R_{2->4} - R_{4->2}) = -4 R_eq mu/T  (linear in a number chemical potential mu)
  chi      = dn/dmu = Int d^3p/(2pi)^3 f (1 + f)/T
  Gamma_chem = 4 R_eq / (T chi);     c_eff = Gamma_chem / (lambda^4 T).

Phase space: RAMBO (Kleiss-Stirling-Ellis) massless 4-body in the pair rest frame, rescaled to mass m with the
standard massive weight; momenta boosted to the plasma frame for the Bose factors (|M|^2 is invariant).
Checks: (a) massive 2-body weight reproduces beta = sqrt(1 - 4m^2/s) exactly; (b) the same machinery for 2 -> 2
reproduces the analytic massless Boltzmann rate.
"""
import math, json, os, sys, time
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'mc24.json')
rng = np.random.default_rng(20261101)


def bose(E):
    return 1.0 / np.expm1(E)


def make_sampler(m, stat='bose', pmax=40.0, ngrid=200001):
    p = np.linspace(0, pmax, ngrid); E = np.sqrt(p * p + m * m)
    f = bose(E) if stat == 'bose' else np.exp(-E)
    w = p * p * f
    cdf = np.concatenate([[0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(p))])
    n = cdf[-1] / (2 * math.pi ** 2)
    fp = f * (1 + f) if stat == 'bose' else f
    chi = np.trapezoid(p * p * fp, p) / (2 * math.pi ** 2)       # T = 1
    cdf /= cdf[-1]
    def draw(k):
        u = rng.random(k); pm = np.interp(u, cdf, p)
        c = 2 * rng.random(k) - 1; ph = 2 * math.pi * rng.random(k); s = np.sqrt(1 - c * c)
        v = np.stack([pm * s * np.cos(ph), pm * s * np.sin(ph), pm * c], -1)
        return np.concatenate([np.sqrt(pm * pm + m * m)[:, None], v], -1)
    return draw, n, chi


def mdot(a, b):
    return a[..., 0] * b[..., 0] - np.sum(a[..., 1:] * b[..., 1:], -1)


def rambo(rs, nb, k):
    """massless RAMBO in the rest frame of total energy rs (array, shape k).  Returns (k, nb, 4) momenta."""
    r = rng.random((k, nb, 4))
    c = 2 * r[..., 0] - 1; ph = 2 * math.pi * r[..., 1]; q0 = -np.log(r[..., 2] * r[..., 3]); s = np.sqrt(1 - c * c)
    q = np.stack([q0, q0 * s * np.cos(ph), q0 * s * np.sin(ph), q0 * c], -1)
    Q = q.sum(1); M = np.sqrt(np.maximum(mdot(Q, Q), 1e-300))
    b = -Q[:, 1:] / M[:, None]; x = rs / M; g = Q[:, 0] / M; a = 1 / (1 + g)
    bq = np.sum(b[:, None, :] * q[..., 1:], -1)
    p0 = x[:, None] * (g[:, None] * q[..., 0] + bq)
    pv = x[:, None, None] * (q[..., 1:] + b[:, None, :] * q[..., 0:1] + a[:, None, None] * bq[..., None] * b[:, None, :])
    return np.concatenate([p0[..., None], pv], -1)


def massive(p, rs, m):
    """rescale massless momenta to mass m; return momenta and the massive/massless weight ratio."""
    pa = p[..., 0]                                         # |p| of massless momenta
    xi = np.sqrt(np.maximum(1 - (p.shape[1] * m / rs) ** 2, 1e-12))
    for _ in range(60):
        E = np.sqrt(m * m + (xi[:, None] * pa) ** 2)
        F = E.sum(1) - rs; dF = np.sum(xi[:, None] * pa * pa / E, 1)
        xi = xi - F / dF
    kv = xi[:, None, None] * p[..., 1:]; ka = xi[:, None] * pa; E = np.sqrt(m * m + ka * ka)
    n = p.shape[1]
    w = (ka.sum(1) / rs) ** (2 * n - 3) * np.prod(ka / E, 1) * rs / np.sum(ka * ka / E, 1)
    return np.concatenate([E[..., None], kv], -1), w


def phi_massless(nb, s):
    return (2 * math.pi) ** (4 - 3 * nb) * (math.pi / 2) ** (nb - 1) * s ** (nb - 2) / (math.factorial(nb - 1) * math.factorial(nb - 2))


def boost(k, P):
    """boost momenta k (rest frame of P) to the frame where the pair has 4-momentum P."""
    M = np.sqrt(mdot(P, P)); bv = P[:, 1:] / P[:, 0:1]; g = P[:, 0] / M
    bk = np.sum(bv[:, None, :] * k[..., 1:], -1); b2 = np.sum(bv * bv, -1)
    E = g[:, None] * (k[..., 0] + bk)
    coef = ((g - 1) / np.maximum(b2, 1e-300))[:, None] * bk + g[:, None] * k[..., 0]
    v = k[..., 1:] + coef[..., None] * bv[:, None, :]
    return np.concatenate([E[..., None], v], -1)


TRIPLES = [(0, a, b) for a in range(1, 6) for b in range(a + 1, 6)]   # 10 partitions of 6 legs
SIGN = np.array([1, 1, -1, -1, -1, -1], float)


def rate(m, nb_final=4, stat='bose', nev=2_000_000, batch=200_000, lam=1.0):
    draw, n, chi = make_sampler(m, stat)
    acc, acc2, ntot = 0.0, 0.0, 0
    thr = (nb_final * m) ** 2
    while ntot < nev:
        k = batch
        p1, p2 = draw(k), draw(k); P = p1 + p2; s = mdot(P, P)
        val = np.zeros(k)
        ok = s > thr * (1 + 1e-9)
        if ok.any():
            Po, so = P[ok], s[ok]; rs = np.sqrt(so)
            q = rambo(rs, nb_final, ok.sum())
            km, wm = massive(q, rs, m) if m > 0 else (q, np.ones(ok.sum()))
            kl = boost(km, Po)
            El = kl[..., 0]
            fin = np.prod(1 + bose(El), 1) if stat == 'bose' else np.ones(ok.sum())
            if nb_final == 4:
                legs = np.concatenate([p1[ok][:, None], p2[ok][:, None], kl], 1) * SIGN[None, :, None]
                amp = np.zeros(ok.sum())
                for t in TRIPLES:
                    qq = legs[:, t[0]] + legs[:, t[1]] + legs[:, t[2]]
                    amp += 1.0 / (mdot(qq, qq) - m * m)
                M2 = lam ** 4 * amp * amp
                sym = 1.0 / math.factorial(4)
            else:
                M2 = lam ** 2 * np.ones(ok.sum()); sym = 0.5
            ph = phi_massless(nb_final, so) * wm
            val[ok] = sym * ph * M2 * fin / (4 * p1[ok][:, 0] * p2[ok][:, 0])
        acc += val.sum(); acc2 += (val * val).sum(); ntot += k
    mean = acc / ntot; err = math.sqrt(max(acc2 / ntot - mean * mean, 0) / ntot)
    R = 0.5 * n * n * mean; dR = 0.5 * n * n * err
    nb = nb_final
    G = (2 * (nb - 2)) * R / chi if nb == 4 else None    # 4 R / (T chi) for 2<->4  (T = 1)
    return dict(m_over_T=m, stat=stat, n=n, chi=chi, R=R, R_err=dR, Gamma_chem=G,
                Gamma_chem_err=(4 * dR / chi if nb == 4 else None), c_eff=(G / lam ** 4 if G is not None else None),
                events=ntot)


def checks():
    out = {}
    # (a) massive 2-body weight equals beta
    rs = np.array([3.0, 5.0, 10.0]); m = 1.0
    q = rambo(rs, 2, 3); _, w = massive(q, rs, m)
    out['two_body_weight_vs_beta'] = dict(weight=w.tolist(), beta=np.sqrt(1 - 4 * m * m / rs ** 2).tolist())
    # (b) massless Boltzmann 2 -> 2 rate: R = (1/2) n^2 lambda^2/(64 pi) <1/E>^2, <1/E> = 1/2, n = 1/pi^2
    r = rate(0.0, nb_final=2, stat='boltzmann', nev=400_000)
    n = 1 / math.pi ** 2
    out['two_to_two_boltzmann'] = dict(mc=r['R'], mc_err=r['R_err'], analytic=0.5 * n * n / (64 * math.pi) * 0.25)
    # (c) massless 4-body phase-space volume vs RAMBO normalisation (trivially exact; report value at s = 1)
    out['phi4_massless_s1'] = phi_massless(4, 1.0)
    return out


def main():
    t0 = time.time()
    res = dict(checks=checks())
    print('CHECKS', json.dumps(res['checks']), flush=True)
    rows = []
    nev = int(sys.argv[1]) if len(sys.argv) > 1 else 2_000_000
    for m in [1.0, 0.5, 0.36, 0.16, 0.08]:
        for stat in ['bose', 'boltzmann']:
            r = rate(m, 4, stat, nev=nev); rows.append(r)
            print('RATE m/T=%.3f %-9s Gamma_chem/(lam^4 T) = %.3e +- %.1e   (pt_relax c bracket 3.2e-5 - 3.2e-3)  t=%.0fs'
                  % (m, stat, r['c_eff'], r['Gamma_chem_err'], time.time() - t0), flush=True)
    res['rows'] = rows
    res['pt_relax_bracket'] = [10 / (4 * math.pi) ** 5, 1e3 / (4 * math.pi) ** 5]
    json.dump(res, open(OUT, 'w'), indent=1)


if __name__ == '__main__':
    main()
