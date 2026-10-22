"""v_dissipation, script 6: the K'(q) != 0 version of the KV16a drain (track NOTES sec 4.5, 'QV2-diss' with QV1's
K(q)), recomputed independently, plus a test of the quasi-static closed form derived in v1_symbolic.py (check G).

Units of 0806.2805/KSV16 (sigma = 1, chi0 q0^2 = 1): K = f/3, eps = lambda + (1/2)(-f^2 + f^4/3), radiation w = 1/3.
Drain (the track's transcription of KV16a sec 7, second remark): u-dot = s with f s = gamma |h| (hdot + 2 h^2)^2.
State (f, fdot, h, hdot, r_M, ln a), same independent formulation as v4_ksv16.py (fddot from (g)_ii, hddot from u-dot = s).
Quasi-static prediction (v1 check G, first order in the drift, all orders in H/m_q):
    d rV/dN = -4 gamma rV^2 eps''(f) / (f^2 eps''(f) - 2 rV),   rV = f eps'(f) - eps(f) on the static branch,
    = -(4 gamma/f0^2) rV^2 [1 - eps'''(f0) df^2/(f0 eps''(f0)) + O(df^3)]  (the O(df) = O(H^2/m_q^2) term cancels).
Track values [NOTES sec 4.5, results/qv1diss.json].  Output: ../results/v6_qv1diss.json.  Run: python3 -I v6_qv1diss.py
"""
import json, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

W = 1 / 3


def eps(f, lam): return lam + 0.5 * (-f * f + f**4 / 3)
def e1(f): return -f + 2 * f**3 / 3
def e2(f): return -1 + 2 * f * f
def e3(f): return 4 * f


def rhs(tau, y, lam, gam):
    f, fd, h, hd, rM, lna = y
    u = e1(f) - hd - 2 * h * h
    k = hd + 2 * h * h
    s = gam * abs(h) * k * k / f
    fdd = -3 * W * rM + 3 * eps(f, lam) - 3 * f * h * h - 3 * f * u - 2 * f * hd - 2 * h * fd
    hdd = -4 * h * hd + fd * e2(f) - s
    return [fd, fdd, hd, hdd, -4 * h * rM + f * s, h]


def rVt(y, lam):
    f, fd, h, hd, rM, lna = y
    return eps(f, lam) - (e1(f) - hd - 2 * h * h) * f


def constraint(y, lam):
    f, fd, h, hd, rM, lna = y
    return f * h * h + h * fd - rVt(y, lam) - rM


def f0_of(lam):          # rV_flat(f) = lam + f^2/2 - f^4/2 = 0
    return np.sqrt((1 + np.sqrt(1 + 8 * lam)) / 2)


def static_branch(u, lam):
    """f > f0 with u = 2 eps/f - eps' (static de Sitter branch); h^2 = eps' - eps/f."""
    f0 = f0_of(lam)
    U = lambda f: 2 * eps(f, lam) / f - e1(f) - u
    f = brentq(U, f0 * (1 + 1e-12), 5.0)
    return f, np.sqrt(e1(f) - eps(f, lam) / f)


def qs_rate(r, lam, gam):
    """quasi-static d rV/dN at given rV (branch f(rV) with rV = f eps' - eps)."""
    f0 = f0_of(lam)
    f = brentq(lambda f: f * e1(f) - eps(f, lam) - r, f0, 5.0)
    return -4 * gam * r * r * e2(f) / (f * f * e2(f) - 2 * r)


out = {}
gam = 0.01
track_dev = {0.0: [1.11e-2, 7.29e-4, 2.54e-5, 1.02e-6, 2.85e-8], 1.0: [6.19e-4, 2.62e-5, 7.43e-7, 4.67e-8, 2.28e-9]}
rows = []
for lam in (0.0, 1.0):
    f0 = f0_of(lam)
    u0 = 2 * eps(f0, lam) / f0 - e1(f0)
    mq = np.sqrt(2 * f0 * e2(f0))
    for j, du in enumerate((-0.3, -0.1, -0.03, -0.01, -0.003)):
        fi, hi = static_branch(u0 + du, lam)
        y0 = [fi, 0.0, hi, 0.0, 0.0, 0.0]
        rV1 = rVt(y0, lam)

        def eN(tau, y, *a):
            return y[5] - 6.0
        eN.terminal = True
        s = solve_ivp(rhs, (0.0, 1e7), y0, args=(lam, gam), method='DOP853', rtol=1e-12, atol=1e-15, events=[eN],
                      dense_output=True)
        Ns = np.linspace(0, 6, 61)
        taus = []
        # invert ln a(tau) on the dense output
        tg = np.linspace(0, s.t[-1], 4001)
        lna = s.sol(tg)[5]
        devs, qsdevs = [], []
        # quasi-static law integrated in N
        qs = solve_ivp(lambda N, r: [qs_rate(r[0], lam, gam)], (0, 6), [rV1], rtol=1e-12, atol=1e-18, dense_output=True)
        for N in Ns[1:]:
            i_ = int(np.searchsorted(lna, N))
            a_, b_ = tg[max(i_ - 1, 0)], tg[min(i_ + 1, len(tg) - 1)]
            tN = brentq(lambda t_: s.sol(t_)[5] - N, a_, b_) if N < 6.0 else s.t[-1]
            yN = s.sol(tN)
            r = rVt(yN, lam)
            law = rV1 / (1 + 4 * gam * rV1 * N / f0**2)
            devs.append(abs(r / law - 1))
            qsdevs.append(abs(r / qs.sol(N)[0] - 1))
        dfi = fi - f0
        pred_rate = -e3(f0) * dfi**2 / (f0 * e2(f0))
        x6 = 4 * gam * rV1 * 6 / f0**2
        rows.append({'lambda': lam, 'du': du, 'f_i': fi, 'H_i/m_q': hi / mq, 'rV_i': rV1,
                     'max_rel_dev_from_constG_law': max(devs), 'track': track_dev[lam][j],
                     'max_rel_dev_from_quasistatic_law': max(qsdevs),
                     'leading_rate_correction -e3 df^2/(f0 e2)': pred_rate,
                     'x(6)=4 gamma rV N/f0^2': x6, 'predicted dev ~ x |rate corr|/(1+x)': x6 * abs(pred_rate) / (1 + x6),
                     'max|constraint|': float(max(abs(constraint(s.sol(t_), lam)) for t_ in np.linspace(0, s.t[-1], 50)))})
        print(json.dumps(rows[-1]))
out['reduction'] = rows
# slopes of log dev vs log(H/m_q)
for lam in (0.0, 1.0):
    rr = [r_ for r_ in rows if r_['lambda'] == lam]
    sl = [np.log(rr[i]['max_rel_dev_from_constG_law'] / rr[i + 1]['max_rel_dev_from_constG_law']) /
          np.log(rr[i]['H_i/m_q'] / rr[i + 1]['H_i/m_q']) for i in range(len(rr) - 1)]
    out['local_slopes_lambda=%g' % lam] = sl

# KV nonconstant-G remark: {h, hdot, f, r_M}(1) = {1/2, 0, 1/2, 1/20}
rem = {}
for lam in (0.0, 1.0):
    f0 = f0_of(lam)
    h1, hd1, f1, rM1 = 0.5, 0.0, 0.5, 0.05
    u1 = e1(f1) - hd1 - 2 * h1 * h1
    fd1 = (eps(f1, lam) - u1 * f1 + rM1 - f1 * h1 * h1) / h1
    s = solve_ivp(rhs, (1.0, 1e5), [f1, fd1, h1, hd1, rM1, 0.0], args=(lam, gam), method='DOP853', rtol=1e-11,
                  atol=1e-15, t_eval=[1e3, 1e4, 1e5])
    y = s.y[:, -1]
    rem['lambda=%g' % lam] = {'u_i': u1, 'tau^(2/3) rV(1e5)': float(1e5**(2 / 3) * rVt(y, lam)),
                              'tau^(1/3) h(1e5)': float(1e5**(1 / 3) * y[2]), 'f(1e5)': float(y[0]),
                              'pred f0^(5/3)(6g)^(-2/3)': f0**(5 / 3) * (6 * gam)**(-2 / 3),
                              'pred f0^(1/3)(6g)^(-1/3)': f0**(1 / 3) * (6 * gam)**(-1 / 3), 'f0': f0,
                              'track': {'tau^(2/3) rV': [6.54, 11.64][int(lam)], 'tau^(1/3) h': [2.554, 2.867][int(lam)],
                                        'f': [1.003, 1.4155][int(lam)]},
                              'constraint(1e5)': float(constraint(y, lam)), 'status': int(s.status)}
out['KV_remark'] = rem
print(json.dumps(rem, indent=1))
print(out['local_slopes_lambda=0'], out['local_slopes_lambda=1'])
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'v6_qv1diss.json'), 'w') as fh:
    json.dump(out, fh, indent=1)
