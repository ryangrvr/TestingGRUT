"""v_dissipation, script 5: K22 (Klinkhamer, arXiv:2207.03453, 4D-brane q) recomputed independently.

Formulation (not the track's): state (xi, r_M, ln a); h is NOT evolved but taken from the first Friedmann equation
h = +sqrt((8 pi/3)(r_M + r_V)) (expanding branch; integration stops if h^2 -> 0, i.e. turnaround), and
hdot = -(3/2) h^2 - 4 pi (w r_M - r_V) from the second (only needed for the pp source).  E_P units (G = 1) as K22 (7).
    r_M-dot = -4 h r_M + Gamma,    xi-dot = -Gamma / r_V'(xi),   r_V = lambda + xi^2 - mu xi,   w = 1/3
    Gamma_adhoc = g~ |h| h^4 (K22 eq 17);  Gamma_pp = gamma |h| (hdot + 2 h^2)^2 (K22 eq 8)
Paper values [Q]: 2207.03453v4 Figs 1-3 captions, Table III, eq (20); track values from ../../dissipation/NOTES.md
secs 3 and 4.4.  Output: ../results/v5_k22.json.  Run: python3 -I v5_k22.py
"""
import json, os
import numpy as np
from scipy.integrate import solve_ivp
from multiprocessing import Pool

W = 1 / 3
MU = 4.0


def rV(xi, lam):
    return lam + xi * xi - MU * xi


def hubble(xi, rM, lam):
    return np.sqrt(max(8 * np.pi / 3 * (rM + rV(xi, lam)), 0.0))


def rhs(tau, y, lam, g, src):
    xi, rM, lna = y
    h = hubble(xi, rM, lam)
    if src == 'adhoc':
        G = g * h**5
    else:
        hd = -1.5 * h * h - 4 * np.pi * (W * rM - rV(xi, lam))
        G = g * h * (hd + 2 * h * h)**2
    return [-G / (2 * xi - MU), -4 * h * rM + G, h]


def make_events(lam):
    def e_neg(tau, y, *a):          # r_V < 0: no return (Gamma >= 0)
        return rV(y[0], lam)
    e_neg.terminal = True; e_neg.direction = -1

    def e_sing(tau, y, *a):         # xi reaches mu/2: xi-dot diverges
        return abs(2 * y[0] - MU) - 1e-7
    e_sing.terminal = True; e_sing.direction = -1

    def e_turn(tau, y, *a):
        return y[1] + rV(y[0], lam)
    e_turn.terminal = True; e_turn.direction = -1
    return [e_neg, e_sing, e_turn]


def run(lam, h1, xi1, g, src, tmax, teval=None, rtol=1e-11, events=True):
    rM1 = 3 * h1 * h1 / (8 * np.pi) - rV(xi1, lam)
    return solve_ivp(rhs, (1.0, tmax), [xi1, rM1, 0.0], args=(lam, g, src), method='DOP853', rtol=rtol, atol=1e-16,
                     t_eval=teval, events=make_events(lam) if events else None)


def classify(args):
    lam, xi1, src, g = args
    h1 = 9.0
    rM1 = 3 * h1 * h1 / (8 * np.pi) - rV(xi1, lam)
    if rM1 < 0:
        return 'invalid'
    if rV(xi1, lam) < 0:
        return 'crunch'
    s = run(lam, h1, xi1, g, src, 1e5)
    if s.status == 1:
        if len(s.t_events[0]):
            return 'crunch'
        if len(s.t_events[1]):
            return 'singular'
        return 'crunch'
    xi, rM, _ = s.y[:, -1]
    r = rV(xi, lam)
    if r > 0 and rM < r:
        return 'mink'
    return 'undecided'


if __name__ == '__main__':
    out = {}
    gt = 1e-3
    # Figs 1-2 and 3 boundary data
    figs = {}
    for lam, rMpaper in ((1.0, 6.16616), (-1.0, 8.16616), (4.0, 3.16616)):
        s = run(lam, 9.0, 4.55, gt, 'adhoc', 1e6, teval=[1e2, 1e3, 1e4, 1e5, 1e6], events=False)
        xib = (MU + np.sqrt(MU * MU - 4 * lam)) / 2
        rows = []
        for i, t in enumerate(s.t):
            xi, rM, _ = s.y[:, i]
            rows.append({'tau': float(t), 'h': float(hubble(xi, rM, lam)), 'xi': float(xi), 'r_M': float(rM), 'r_V': float(rV(xi, lam))})
        figs['lambda=%g' % lam] = {'r_M(1)': 3 * 81 / (8 * np.pi) - rV(4.55, lam), 'paper r_M(1) [Q]': rMpaper,
                                   'xibar': xib, 'rows': rows}
    out['figs'] = figs
    out['TableIII_paper_[Q]'] = {'h': [0.926, 0.430, 0.200, 0.0927, 0.0430], 'xi': [2.32, 2.15, 2.07, 2.03, 2.01],
                                 'r_M': [1.85e-4, 8.55e-6, 3.97e-7, 1.84e-8, 8.56e-10], 'r_V': [0.102, 0.0221, 0.00476, 0.00102, 0.000221]}
    cR = (1.5)**(-2 / 3) * (8 * np.pi / 3)**(-5 / 3)
    out['asymptote'] = {'r_V coeff': cR, 'h coeff': np.sqrt(8 * np.pi / 3 * cR),
                        'r_V(1e5) asym': cR * (gt * 1e5)**(-2 / 3), 'r_V(1e6) asym': cR * (gt * 1e6)**(-2 / 3)}
    # rectangle (20) corners
    cor = []
    for lam in (1.0, -1.0):
        for h1 in (8.9, 9.1):
            for xi1 in (4.54, 4.56):
                s = run(lam, h1, xi1, gt, 'adhoc', 1e5, teval=[1e5], events=False)
                cor.append({'lambda': lam, 'h1': h1, 'xi1': xi1, 'r_V(1e5)': float(rV(s.y[0, -1], lam)), 'xi(1e5)': float(s.y[0, -1])})
    out['rectangle_corners'] = cor
    # ad hoc overshoot example: lambda = 1, xi = 3.8, h = 9
    s = run(1.0, 9.0, 3.8, gt, 'adhoc', 10.0)
    out['adhoc_example_lambda1_xi3.8'] = {'r_V,i': rV(3.8, 1.0), 'stopped_at_tau': float(s.t[-1]),
                                          'reason': ['r_V<0', 'xi=mu/2', 'h^2=0'][[len(e) > 0 for e in s.t_events].index(True)] if s.status == 1 else 'none',
                                          'track': 'crunches at tau = 1.005'}
    s2 = run(1.0, 9.0, 3.8, gt, 'adhoc', 10.0, events=False)
    out['adhoc_example_lambda1_xi3.8']['continued_status'] = int(s2.status)
    out['adhoc_example_lambda1_xi3.8']['continued_last_tau'] = float(s2.t[-1])
    print(json.dumps(out, indent=1)[:4000])
    # basin grids
    lams = np.round(np.arange(-1.0, 5.0 + 1e-9, 0.25), 10)
    xis = np.round(np.arange(0.1, 6.0 + 1e-9, 0.1), 10)
    grids = {}
    for src, g in (('adhoc', 1e-3), ('pp', 1e-3)):
        jobs = [(float(l), float(x), src, g) for l in lams for x in xis]
        with Pool(4) as P:
            res = P.map(classify, jobs)
        res = np.array(res).reshape(len(lams), len(xis))
        valid = res != 'invalid'
        low = (lams <= MU * MU / 4)[:, None] & valid
        high = (lams > MU * MU / 4)[:, None] & valid
        grids[src] = {'basin_frac_valid': float((res == 'mink').sum() / valid.sum()),
                      'basin_frac_lambda<=4': float((res == 'mink')[low].sum() / low.sum()),
                      'basin_frac_lambda>4': float((res == 'mink')[high].sum() / high.sum()),
                      'outcomes_lambda>4': sorted(set(res[high].tolist())),
                      'counts': {k: int((res == k).sum()) for k in ('mink', 'crunch', 'singular', 'invalid', 'undecided')}}
        # analytic expectation for pp: basin = valid & lambda <= 4 & r_V,i >= 0
        LL, XX = np.meshgrid(lams, xis, indexing='ij')
        ana = valid & (LL <= 4) & (rV(XX, LL) >= 0)
        grids[src]['analytic_rV>=0_frac_valid'] = float(ana.sum() / valid.sum())
        print(src, json.dumps(grids[src]))
    grids['track'] = {'adhoc': [0.309, 0.363, 0.0], 'pp': [0.392, 0.460, 0.0]}
    out['basin'] = grids
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'v5_k22.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
