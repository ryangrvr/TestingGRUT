"""Track `dissipation`, script 5: the dissipative law with the QV1 structure K'(q) != 0 (a propagating q-mode).

Our transcription (KV 1601.00601 sec 7 states the result for nonconstant G but prints no equations):
units of 0806.2805 / 1601.04676 (sigma = 1, chi0 q0^2 = 1), K = f/3, eps = lambda + (1/2)(-f^2 + f^4/3), radiation
w = 1/3, and the drain  d u_eff/dtau = s,  f s = gamma |h| (hdot + 2 h^2)^2  (energy into radiation: r_M' += f s).
(diss_models.ksv_rhs with source='pp').

Checks:
  1. Reduction to the constant-G law at leading order in H/m_q.  Constant-G prediction with K0 = f0/3:
       k = hdot + 2h^2 = 2 r_V/f0  =>  d r_V/dN = -(4 gamma/f0^2) r_V^2,  r_V(N) = r_V1/(1 + 4 gamma r_V1 dN/f0^2).
     m_q^2 = 2 f0 eps''(f0) in these units (equations track section 4.1 with K = f/3).  Runs start on the static de
     Sitter branch (fdot = hdot = r_M = 0) at several u below u0 and record the deviation from the law vs H/m_q.
  2. KV's statement (1601.00601 sec 7, second remark) that the nonconstant-G model with an added lambda shows the same
     tau^(-1/3), tau^(-2/3) asymptote: run lambda = 0, +1 from boundary data {h, hdot, f, r_M} = {1/2, 0, 1/2, 1/20}
     (0806.2805 Figs 1-2 data with u_eff free).  For this eps, r_V(f) = lambda + f^2/2 - f^4/2 has no zero for
     lambda < -1/8, so no lambda below -1/8 can be cancelled with this eps [D].
Output: ../results/qv1diss.json.  Run: python3 -I s5_qv1diss.py  (about 1-3 min)
"""
import json, os, sys, time
import numpy as np
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from diss_models import *  # noqa

GAM = 0.01


def static_dS(lam, du):
    u0, f0 = ksv_u0(lam)
    # static branch: u = eps'(f) - 2h^2, h^2 = eps'(f) - eps/f  =>  u(f) = eps/f... solve u(f) = u0 + du near f0
    from scipy.optimize import brentq
    def g(f):
        h2 = ksv_epsp(f) - ksv_eps(f, lam) / f
        return ksv_epsp(f) - 2 * h2 - (u0 + du)
    f = brentq(g, f0 * 0.5, f0 * 0.99999999) if g(f0 * 0.5) * g(f0 * 0.99999999) < 0 else brentq(g, f0 * 1.00000001, f0 * 1.5)
    h2 = ksv_epsp(f) - ksv_eps(f, lam) / f
    return f, np.sqrt(h2), u0, f0


def run_reduction(args):
    lam, du, Nspan = args
    f, h, u0, f0 = static_dS(lam, du)
    y0 = np.array([f, 0.0, h, u0 + du, 0.0, 0.0])
    rV1 = ksv_rVtilde(y0, lam)
    mq = np.sqrt(2 * f0 * ksv_epspp(f0))

    def ev_N(t, y, *a):
        return y[5] - Nspan
    ev_N.terminal = True
    sol = integrate(ksv_rhs, 0.0, 1e9, y0, (0.0, 0.0, 1 / 3, lam, 'pp', GAM), rtol=1e-10, atol=1e-14,
                    events=[ev_N], dense=True)
    # sample in N
    taus = np.linspace(0, sol.t[-1], 4000)
    Y = sol.sol(taus)
    N = Y[5]
    rV = np.array([ksv_rVtilde(Y[:, i], lam) for i in range(Y.shape[1])])
    law = rV1 / (1 + 4 * GAM * rV1 * N / f0**2)
    dev = np.max(np.abs(rV / law - 1))
    return {'lambda': lam, 'du': du, 'f_i': f, 'h_i': h, 'r_V,i': rV1, 'm_q': mq, 'H_i/m_q': h / mq,
            'N_span': float(N[-1]), 'r_V_end': float(rV[-1]), 'law_end': float(law[-1]),
            'max_rel_dev_from_constG_law': float(dev), 'constraint_end': float(ksv_constraint(sol.y[:, -1], lam))}


def run_kv_bc(args):
    lam, tend = args
    y0 = ksv_initial(0.5, 0.0, 0.5, 1 / 20, lam)
    te = [1e1, 1e2, 1e3, 1e4, 1e5][: [1e1, 1e2, 1e3, 1e4, 1e5].index(tend) + 1]
    sol = integrate(ksv_rhs, 1.0, tend, y0, (0.0, 0.0, 1 / 3, lam, 'pp', GAM), rtol=1e-10, atol=1e-14, t_eval=te)
    u0, f0 = ksv_u0(lam)
    rows = []
    for tt, yy in zip(sol.t, sol.y.T):
        rows.append({'tau': float(tt), 'h': float(yy[2]), 'f': float(yy[0]), 'u_eff': float(yy[3]),
                     'r_V': float(ksv_rVtilde(yy, lam)), 'r_M': float(yy[4]),
                     'tau^(2/3) r_V': float(tt**(2 / 3) * ksv_rVtilde(yy, lam)), 'tau^(1/3) h': float(tt**(1 / 3) * yy[2]),
                     'constraint': float(ksv_constraint(yy, lam))})
    # constant-G-law prediction of the late asymptote: r_V = f0^2/(4 gamma N) and h^2 ~ r_V/f0
    return {'lambda': lam, 'u_eff,i': float(y0[3]), 'u0': u0, 'f0': f0, 'r_V,i': float(ksv_rVtilde(y0, lam)),
            'rows': rows, 'status': sol.message}


if __name__ == '__main__':
    T0 = time.time()
    pool = Pool(4)
    out = {'labels': '[C] computed here; transcription [D] (ours)'}
    jobs = [(lam, du, 6.0) for lam in (0.0, 1.0) for du in (-0.3, -0.1, -0.03, -0.01, -0.003)]
    out['reduction_to_constant_G_law'] = pool.map(run_reduction, jobs)
    out['KV_remark_nonconstant_G'] = pool.map(run_kv_bc, [(0.0, 1e5), (1.0, 1e5)])
    out['lambda_range_note_[D]'] = ('eps = lambda + (1/2)(-f^2 + f^4/3): r_V(f) = eps - f eps\' = lambda + f^2/2 - f^4/2, '
                                    'max lambda + 1/8 at f^2 = 1/2: no Minkowski state for lambda < -1/8.')
    pool.close()
    out['runtime_s'] = time.time() - T0
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../results/qv1diss.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    for r in out['reduction_to_constant_G_law']:
        print('lam %.0f du %.3f H/m %.4f Nspan %.1f dev %.3e  rVend %.4e law %.4e C %.1e' % (
            r['lambda'], r['du'], r['H_i/m_q'], r['N_span'], r['max_rel_dev_from_constG_law'], r['r_V_end'],
            r['law_end'], r['constraint_end']))
    for r in out['KV_remark_nonconstant_G']:
        print(r['lambda'], r['u_eff,i'], r['u0'], [(x['tau'], round(x['tau^(2/3) r_V'], 5), round(x['tau^(1/3) h'], 5),
                                                    round(x['f'], 5)) for x in r['rows']])
    print(time.time() - T0)
