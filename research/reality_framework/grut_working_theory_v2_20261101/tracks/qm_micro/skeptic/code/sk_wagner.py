"""qm_micro SKEPTIC (Claude, 2 Nov 2026): is qm_micro's self-similar ripening coefficient R_c^2 = 0.372 K_W t a
property of the kinetics or of the initial distribution's edge?  COPY of qm_micro coarsen.wagner_test() with the
initial distribution parametrised; the kinetics (exchange_rate, refine, Q_min) are imported from qm_micro's own
coarsen.py unchanged.

LSW/Wagner theory [R: Lifshitz-Slyozov 1961, Wagner 1961; Niethammer-Pego 1999]: for dR/dt = K(1/R_c - 1/R) the
self-similar solutions form a one-parameter family R_c^2 = alpha K t with alpha <= 1/2; which one is reached depends
on the behaviour of the initial number density near its largest radius ("weak selection").  If alpha moves with the
initial edge, qm_micro's <Q> (set by ripening: N_obj ~ (alpha K t)^{-3/2}) carries that freedom (<Q> ~ alpha^{3/2}).
Shapes: 'trunc4' (qm_micro: lognormal sigma 0.7 cut at +-4 sigma, a jump edge), 'trunc2' (cut at +-2 sigma),
'smooth2' and 'smooth6' (N(R) ~ (1 - R/R_max)^p with p = 2, 6: vanishing edges), 'exp' (lognormal to +-8 sigma).
"""
import json, math, os, sys, time
import numpy as np
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__))
QMC = os.path.join(HERE, '..', '..', 'code')
sys.path.insert(0, QMC)
import coarsen as cq            # qm_micro's own module, imported read-only
import qm_vapour as qv
import bessel_vapour_copy as bv


def run(shape, nit=20000):
    T = 0.5
    mc = 0.775905
    FX, dFX = qv.flux(mc, T, 'X', True, True); FB, dFB = qv.flux(-mc, T, 'Xbar', True, True)
    v = bv.vapour(mc, T)
    KW = 2 * cq.SIGMA * (dFX + dFB) / cq.NB ** 2
    nQ = 1.0
    s = dict(T=T, Phi_X=FX, Phi_Xbar=FB, n_X=v['nX'], n_Xbar=v['nXb'], nQ=nQ)
    Q0 = 1e12
    if shape in ('trunc4', 'trunc2', 'exp'):
        cut = dict(trunc4=2.8, trunc2=1.4, exp=5.6)[shape]
        lq = np.linspace(math.log(Q0) - cut, math.log(Q0) + cut, 201)
        N = np.exp(-0.5 * ((lq - math.log(Q0)) / 0.7) ** 2)
    else:
        p = dict(smooth2=2.0, smooth6=6.0)[shape]
        R0 = cq.radius(Q0)
        Rmax = 2.0 * R0
        R = np.linspace(0.2 * R0, Rmax, 202)[:-1]
        N = (1 - R / Rmax) ** p * (R / R0) ** 2      # number per unit R (arbitrary interior shape), smooth zero at R_max
        Q = 4 * math.pi / 3 * cq.NB * R ** 3
        lq = np.log(Q)
        N = N * np.gradient(R)                    # cohort numbers
    eta0 = 1 - v['n'] / nQ
    N *= eta0 / (N * np.exp(lq)).sum()
    t = 0.0
    dt = 0.002 * cq.radius(Q0) ** 2 / KW
    rows = []
    for it in range(nit):
        Q = np.exp(lq)
        def g(d):
            Qn = Q + cq.exchange_rate(Q, d, s) * dt
            Qn = np.where(Qn < cq.QMIN, 0.0, Qn)
            return (N * Qn).sum() + (s['n_X'] * math.exp(d / T) - s['n_Xbar'] * math.exp(-d / T)) / nQ - 1
        d = brentq(g, -1e-3, 1e-3, xtol=1e-20, rtol=1e-15, maxiter=400)
        Qn = Q + cq.exchange_rate(Q, d, s) * dt
        keep = Qn > cq.QMIN
        lq, N = np.log(Qn[keep]), N[keep]
        wq = N * np.exp(lq); o = np.argsort(-wq); cum = np.cumsum(wq[o]) / wq.sum()
        if np.searchsorted(cum, 0.999) + 1 < 80:
            lq, N = cq.refine(lq, N)
        t += dt
        dt *= 1.003
        if it % 500 == 0 and d > 0:
            Rc = 2 * cq.SIGMA / (cq.NB * d)
            Qq = np.exp(lq); Rq = cq.radius(Qq)
            rows.append(dict(it=it, t=t, Rc2_over_KWt=Rc ** 2 / (KW * t), Rmax_over_Rc=float(Rq.max() / Rc),
                             nobj=float(N.sum()), Q_charge_weighted=float((N * Qq * Qq).sum() / (N * Qq).sum())))
    late = [r['Rc2_over_KWt'] for r in rows[-10:]]
    return dict(shape=shape, K_W=KW, rows=rows, alpha_late=float(np.median(late)), alpha_spread=float(np.ptp(late)),
                Rmax_over_Rc_late=float(np.median([r['Rmax_over_Rc'] for r in rows[-10:]])))


def main(out, shapes):
    res = dict(notes=__doc__, runs=[])
    if os.path.exists(out):
        res = json.load(open(out))
    for sh in shapes.split(','):
        t0 = time.time()
        r = run(sh)
        res['runs'] = [x for x in res['runs'] if x['shape'] != sh] + [r]
        print('%s: alpha_late = %.4f (spread over last 10 records %.4f), Rmax/Rc = %.3f, %.0fs' % (
            sh, r['alpha_late'], r['alpha_spread'], r['Rmax_over_Rc_late'], time.time() - t0), flush=True)
        json.dump(res, open(out, 'w'), indent=None, default=float, separators=(',', ':'))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
