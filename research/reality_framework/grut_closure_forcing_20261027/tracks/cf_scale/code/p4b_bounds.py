#!/usr/bin/env python3
"""cf_scale P4b: 'no Landau pole / perturbative and bounded below up to M_P' applied to GRUT couplings (own code).

Criteria (no gravity terms; one loop; betas from beta_derive.py):
  PT-4pi : no |lambda_i| reaches 4 pi below M_P = 2.435e18 GeV;
  PT-1   : (stricter) every |lambda_i(M_P)| <= 1;
  STAB   : V4 bounded below (minimum of V4 on the unit sphere > 0) at every scale from the matching scale to M_P.
Candidate 1 (X-phi sector, matched at m_X = 92.4 MeV; the hidden U(1) is broken at h = 5.33e12 GeV, so e is switched
on only above h, at e = 1e-7). Heavy spot point (matched at M_Y = 4.0e12 GeV, hw_gates START without the D field).
Also: the kappa boundedness ratio kappa/kappa_max (hw_gates gate 4 criterion) along the flow, and the bound on e.
Usage: python3 -B p4b_bounds.py <out.json>
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import rg_common as R

NAMES, NC = R.NAMES, R.NC
MP = R.MP_RED
M_Y, H_SIG, M_X = 4.0e12, 5.33e12, 0.0924
SPOT = dict(lam_X=1e-3, lam_Y=1.0, lam_S=1.0, kap=0.0034586847211042523, g2=0.01, lam_phi=0.03, e=1e-7)
C1 = dict(lam_X=1e-3, g2=0.01, lam_phi=0.03, e=1e-7)


def flow(y0, mu0, mu1=MP, gauge_on_above=None, dense=False):
    """Run from mu0 up to mu1; gauge switched on only above gauge_on_above (GeV) if given. Stops at 4 pi."""
    ev = lambda t, y: np.max(np.abs(y[:NC])) - 4 * np.pi
    ev.terminal = True

    def rhs(t, y):
        mu = mu0 * np.exp(t)
        on = gauge_on_above is None or mu > gauge_on_above
        return R.beta(y, gauge_active=on)
    sol = solve_ivp(rhs, (0.0, np.log(mu1 / mu0)), y0, method='LSODA', rtol=1e-10, atol=1e-30, events=ev,
                    dense_output=dense)
    pole = float(mu0 * np.exp(sol.t_events[0][0])) if len(sol.t_events[0]) else None
    return sol, pole


def ok(y0, mu0, crit, gauge_on_above=None):
    sol, pole = flow(y0, mu0, gauge_on_above=gauge_on_above)
    if pole is not None:
        return False
    if crit == 'PT-1':
        return bool(np.max(np.abs(sol.y[:NC, -1])) <= 1.0)
    return True


def max_param(base, name, mu0, crit, lo, hi, gauge_on_above=None, also=None):
    def f(x):
        d = dict(base); d[name] = x
        if also:
            for n in also:
                d[n] = x
        return 1.0 if ok(R.vec(d), mu0, crit, gauge_on_above) else -1.0
    if f(hi) > 0:
        return None
    if f(lo) < 0:
        return 0.0
    a, b = lo, hi
    for _ in range(60):
        m = np.sqrt(a * b)
        if f(m) > 0:
            a = m
        else:
            b = m
        if b / a < 1 + 1e-6:
            break
    return a


def kap_ratio(c):
    lx, ly, k = c['lam_X'], c['lam_Y'], c['kap']
    if lx <= 0 or ly <= 0:
        return float('inf')
    return k / ((lx / (1.5 * 2**(-1 / 3)))**0.75 * ly**0.25)


def main(out):
    res = dict(method=__doc__)
    # ---------------- Candidate 1
    c1 = {}
    sol, pole = flow(R.vec(C1), M_X, gauge_on_above=H_SIG)
    c1['benchmark_at_MP'] = R.asdict(sol.y[:, -1]); c1['benchmark_pole'] = pole
    c1['benchmark_V4min_at_MP'] = R.V4_unit_min(sol.y[:NC, -1], n=20)[0]
    for crit in ('PT-4pi', 'PT-1'):
        c1[crit] = dict(lam_X_max=max_param(C1, 'lam_X', M_X, crit, 1e-4, 10.0, H_SIG),
                        lam_phi_max=max_param(C1, 'lam_phi', M_X, crit, 1e-4, 10.0, H_SIG),
                        g2_max=max_param(C1, 'g2', M_X, crit, 1e-4, 10.0, H_SIG))
        c1[crit]['g_max'] = np.sqrt(c1[crit]['g2_max'])
        c1[crit]['ms_over_mX_max_at_g0.1'] = np.sqrt(2 * c1[crit]['lam_phi_max']) / 0.1
    # derived Candidate-1 quantities at m_X (tree level): m_s/m_X = sqrt(2 lam_phi)/g
    c1['ms_over_mX_benchmark'] = np.sqrt(2 * 0.03) / 0.1
    res['candidate1'] = c1
    print('C1:', json.dumps(c1, default=float), flush=True)
    # ---------------- heavy spot point
    hs = {}
    sol, pole = flow(R.vec(SPOT), M_Y, dense=True)
    hs['spot_pole_4pi_GeV'] = pole
    for crit in ('PT-4pi', 'PT-1'):
        hs[crit] = dict(common_lamY_lamS_max=max_param(SPOT, 'lam_Y', M_Y, crit, 1e-3, 5.0, also=['lam_S']),
                        lamY_max_lamS1=max_param(SPOT, 'lam_Y', M_Y, crit, 1e-3, 5.0),
                        lamS_max_lamY1=max_param(SPOT, 'lam_S', M_Y, crit, 1e-3, 5.0),
                        kap_max_PT=max_param(SPOT, 'kap', M_Y, crit, 1e-4, 5.0))
    # at the PT-4pi-compatible common value (and at 0.5, 0.3): kappa ratio and V4 minimum along M_Y -> M_P
    rows = []
    lam_pt = hs['PT-4pi']['common_lamY_lamS_max']
    for lam0 in (round(0.999 * lam_pt, 4), 0.5, 0.3, 0.2):
        d = dict(SPOT); d['lam_Y'] = d['lam_S'] = lam0
        sol, pole = flow(R.vec(d), M_Y, dense=True)
        prof = []
        for mu in (M_Y, 1.07e13, 1e14, 1e15, 1e16, 1e17, MP):
            if pole is not None and mu > pole:
                continue
            y = sol.sol(np.log(mu / M_Y)); c = R.asdict(y)
            prof.append(dict(mu=mu, kappa_over_max=kap_ratio(c), V4min=R.V4_unit_min(y[:NC], n=16)[0],
                             lam_X=c['lam_X'], lam_Y=c['lam_Y'], lam_S=c['lam_S'], kap=c['kap'], l_XY=c['l_XY']))
        rows.append(dict(lamY_eq_lamS_at_MY=lam0, pole=pole, profile=prof,
                         kappa_over_max_at_MY=kap_ratio(dict(lam_X=1e-3, lam_Y=lam0, kap=SPOT['kap']))))
        print('HS lam0=%.4f pole=%s k/kmax(MY)=%.3f max over flow=%.3f minV4=%.2e' % (
            lam0, pole, rows[-1]['kappa_over_max_at_MY'], max(p['kappa_over_max'] for p in prof), min(p['V4min'] for p in prof)), flush=True)
    hs['reduced_profiles'] = rows
    # kappa_max needed: largest kappa with kappa/kappa_max < 1 at all scales M_Y..M_P given lam_Y(M_Y) = lam0
    res['heavy_spot'] = hs
    print('HS:', json.dumps({k: v for k, v in hs.items() if k != 'reduced_profiles'}, default=float), flush=True)
    # ---------------- the hidden gauge coupling
    gb = {}
    # pure gauge Landau pole: 1/e^2(MP) = 1/e^2(h) - 2 b ln(MP/h)/16pi^2 > 1/(4pi) (alpha < 1) ; or e(MP) < 1
    L = np.log(MP / H_SIG)
    for tag, emp in (('e(MP)<sqrt(4pi)', np.sqrt(4 * np.pi)), ('e(MP)<1', 1.0)):
        gb[tag] = 1 / np.sqrt(1 / emp**2 + 2 * R.B_E * L / R.K16)
    # with the scalar sector at lam_Y = lam_S = 0.3: e max for PT-4pi (charge-6 Sigma quartic is driven by 7776 e^4)
    d = dict(SPOT); d['lam_Y'] = d['lam_S'] = 0.3
    gb['e_max_PT4pi_with_lam0.3'] = max_param(d, 'e', M_Y, 'PT-4pi', 1e-7, 2.0)
    gb['e_max_PT1_with_lam0.3'] = max_param(d, 'e', M_Y, 'PT-1', 1e-7, 2.0)
    res['hidden_gauge'] = gb
    print('E:', gb, flush=True)
    Path(out).write_text(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main(sys.argv[1])
