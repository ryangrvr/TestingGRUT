"""v_dissipation, script 4: KSV16 (Klinkhamer-Savelainen-Volovik, arXiv:1601.04676) recomputed independently.

Formulation (not the track's): state (f, fdot, h, hdot, r_M, ln a); hddot from (17a), fddot from the space-space
component of the BRIEF (g) with K = f/3 (derived in v1_symbolic.py, check E: fddot = -3 w r_M + 3 eps - 3 f h^2 - 3 f u
- 2 f hdot - 2 h fdot, u = eps' - hdot - 2 h^2), r_M from (17c).  u_eff is a derived quantity, not a state variable
(the track evolves u and solves hdot from it).  fdot(1) from the constraint (17b).  The constraint
C = f h^2 + h fdot - (eps - u f) - r_M is monitored.  eps = lambda + (1/2)(-f^2 + f^4/3); gamma_f = gamma_h = 1; w = 1/3.
Paper values [Q] (1601.04676v4 sec IV, Figs 1-3, Tables I-II); track values [NOTES secs 3, 4.3].
Output: ../results/v4_ksv16.json.  Run: python3 -I v4_ksv16.py  (grid part uses 4 processes)
"""
import json, os, sys
import numpy as np
from scipy.integrate import solve_ivp
from multiprocessing import Pool

W = 1.0 / 3.0


def epsf(f, lam):
    return lam + 0.5 * (-f * f + f**4 / 3)


def e1f(f):
    return -f + 2 * f**3 / 3


def e2f(f):
    return -1 + 2 * f * f


def rhs(tau, y, lam, gf, gh):
    f, fd, h, hd, rM, lna = y
    u = e1f(f) - hd - 2 * h * h
    s = gf * fd * fd + gh * hd * hd
    fdd = -3 * W * rM + 3 * epsf(f, lam) - 3 * f * h * h - 3 * f * u - 2 * f * hd - 2 * h * fd
    hdd = -4 * h * hd + fd * e2f(f) - s
    rMd = -3 * h * (1 + W) * rM + f * s
    return [fd, fdd, hd, hdd, rMd, h]


def init(h, hd, f, rM, lam):
    u = e1f(f) - hd - 2 * h * h
    fd = (epsf(f, lam) - u * f + rM - f * h * h) / h
    return [f, fd, h, hd, rM, 0.0]


def derived(y, lam):
    f, fd, h, hd, rM, lna = y
    u = e1f(f) - hd - 2 * h * h
    rV = epsf(f, lam) - u * f
    C = f * h * h + h * fd - rV - rM
    return u, rV, C


def ev_h0(tau, y, *a):
    return y[2]
ev_h0.terminal = True; ev_h0.direction = -1


def ev_f0(tau, y, *a):
    return y[0] - 1e-3
ev_f0.terminal = True; ev_f0.direction = -1


def run(h, hd, f, rM, lam, tmax, gf=1.0, gh=1.0, rtol=1e-11, events=(ev_h0, ev_f0), teval=None):
    return solve_ivp(rhs, (1.0, tmax), init(h, hd, f, rM, lam), args=(lam, gf, gh), method='DOP853', rtol=rtol,
                     atol=1e-14, events=list(events) if events else None, t_eval=teval, dense_output=False)


def outcome(args):
    h, hd, f, rM, lam, tmax = args
    s = run(h, hd, f, rM, lam, tmax)
    y = s.y[:, -1]
    u, rV, C = derived(y, lam)
    if s.status == 1:
        return ('collapse' if len(s.t_events[0]) else 'f->0'), float(s.t[-1]), float(u), float(rV)
    hdot = y[3]
    if y[2] > 0 and abs(hdot) < 1e-3 * y[2]**2:
        return 'dS', float(s.t[-1]), float(u), float(rV)
    return 'undecided', float(s.t[-1]), float(u), float(rV)


def singularity(h, hd, f, rM, lam):
    """integrate through turnaround until the curvature singularity (r_M or |h| blows up)."""
    def ev_big(tau, y, *a):
        return y[4] - 1e6
    ev_big.terminal = True
    s = solve_ivp(rhs, (1.0, 400.0), init(h, hd, f, rM, lam), args=(lam, 1.0, 1.0), method='LSODA', rtol=1e-11,
                  atol=1e-14, events=[ev_h0, ev_big])
    t_turn = None
    if len(s.t_events[0]):
        t_turn = float(s.t_events[0][0])
        y0 = s.y_events[0][0]
        s2 = solve_ivp(rhs, (t_turn, 400.0), y0, args=(lam, 1.0, 1.0), method='LSODA', rtol=1e-11, atol=1e-14,
                       events=[ev_big])
        s = s2
    return t_turn, float(s.t[-1]), float(np.exp(s.y[5, -1])), s.status


if __name__ == '__main__':
    out = {}
    # ---------------- Figs 1-3
    s1 = run(0.55, 0.0, 0.2953, 0.0, 0.0, 400.0, gf=0.0, gh=0.0, teval=[1.0, 100.0, 400.0])
    u1 = [derived(s1.y[:, i], 0.0)[0] for i in range(s1.y.shape[1])]
    fin = s1.y[:, -1]
    out['Fig1_no_dissipation'] = {'u_eff(1,100,400)': u1, 'paper_[Q]': -0.883133, 'track': -0.8831328,
                                  'h(400)': float(fin[2]), 'track_h_inf': 0.795511,
                                  'h^2 - (eps\' - eps/f)': float(fin[2]**2 - (e1f(fin[0]) - epsf(fin[0], 0) / fin[0]))}
    s2 = run(0.55, 0.0, 0.2953, 0.0, 0.0, 3000.0, teval=[100.0, 3000.0])
    out['Fig2'] = {'u_eff(100)': derived(s2.y[:, 0], 0.0)[0], 'paper_[Q]': -0.333499, 'track': -0.3334995,
                   'h(3000)': float(s2.y[2, -1]), 'track_h_inf': 0.01243, 'status': int(s2.status),
                   'max|C|': float(max(abs(derived(s2.y[:, i], 0.0)[2]) for i in range(s2.y.shape[1])))}
    s3 = run(0.55, 0.0, 0.2952, 0.0, 0.0, 98.0, teval=[98.0], events=None)
    tt, tsing, a_sing, st = singularity(0.55, 0.0, 0.2952, 0.0, 0.0)
    out['Fig3'] = {'u_eff(98)': derived(s3.y[:, 0], 0.0)[0], 'paper_[Q]': -0.333276, 'track': -0.333276,
                   'turnaround_tau': tt, 'singularity_tau(r_M=1e6)': tsing, 'a_at_singularity': a_sing,
                   'paper_[Q] a': 'O(10) at tau ~ 198', 'track': 'singularity tau = 197.58, a = 10.7'}
    print(json.dumps(out, indent=1))
    sys.stdout.flush()

    # ---------------- Table I (three slices) and Table II
    hs = [0.525, 0.5375, 0.55, 0.60, 0.65, 0.70]
    fs = [0.24, 0.29, 0.34, 0.39, 0.44, 0.49]
    paper = ['NNNYYY', 'NNNYYY', 'NYYYYY', 'YYYYYY', 'YYYYYY', 'YYYYYY']
    slices = [(0.0, 0.0), (-0.01, 0.0), (0.0, 0.02)]
    jobs = [(h, hd, f, rM, 0.0, 2000.0) for (hd, rM) in slices for f in fs for h in hs]
    with Pool(4) as P:
        res = P.map(outcome, jobs)
    tabs = []
    k = 0
    for (hd, rM) in slices:
        rows = []
        for f in fs:
            row = ''
            for h in hs:
                row += 'Y' if res[k][0] == 'dS' else ('N' if res[k][0] in ('collapse', 'f->0') else '?')
                k += 1
            rows.append(row)
        tabs.append({'slice (hdot, r_M)': [hd, rM], 'rows f=0.24..0.49': rows, 'matches_paper': rows == paper})
    out['TableI'] = tabs
    print(json.dumps(tabs))
    sys.stdout.flush()
    t2 = []
    for f in [0.352, 0.353, 0.354, 0.355, 0.36, 0.37, 0.38, 0.39]:
        s = run(0.525, 0.0, f, 0.0, 0.0, 200.0, teval=[200.0])
        if s.status == 1:
            te = [e[0] for e in s.t_events if len(e)]
            t2.append({'f1': f, 'outcome': 'turnaround (h=0) at tau=%.2f' % te[0]})
        else:
            y = s.y[:, -1]
            t2.append({'f1': f, 'h(200)': float(y[2]), 'f(200)': float(y[0]), 'rV(200)': float(derived(y, 0.0)[1]), 'rM(200)': float(y[4])})
    out['TableII'] = t2
    out['TableII_paper_[Q]'] = {'0.352': '-', '0.353': [0.0362, 1.00131, 0.00132, 5e-12], '0.354': [0.0510, 1.00259, 0.00260, 3e-14],
                                '0.355': [0.0615, 1.00376, 0.00379, 5e-15], '0.36': [0.0952, 1.00895, 0.00915, 0],
                                '0.37': [0.136, 1.0180, 0.0188, 0], '0.38': [0.164, 1.0260, 0.0277, 0], '0.39': [0.187, 1.0334, 0.0363, 0]}
    print(json.dumps(t2))
    sys.stdout.flush()

    # ---------------- separatrices and fine-tuning (lambda = 0, 1; h = 0.55 slice)
    def bisect(lam, a, b, tmax=1500.0):
        oa = outcome((0.55, 0.0, a, 0.0, lam, tmax))[0]
        ob = outcome((0.55, 0.0, b, 0.0, lam, tmax))[0]
        assert oa != ob, (lam, a, b, oa, ob)
        for _ in range(60):
            m = 0.5 * (a + b)
            om = outcome((0.55, 0.0, m, 0.0, lam, tmax))[0]
            if om == oa:
                a = m
            else:
                b = m
            if b - a < 1e-11:
                break
        return 0.5 * (a + b), oa, ob
    seps = {}
    for lam, br, trk in ((0.0, (0.2952, 0.2953), 0.295245), (0.0, (1.25, 1.30), 1.286909), (1.0, (1.35, 1.45), 1.400694)):
        fsep, oa, ob = bisect(lam, *br)
        side = +1 if ob == 'dS' else -1
        tune = []
        for df in (1e-3, 1e-4, 1e-5):
            r1500 = outcome((0.55, 0.0, fsep + side * df, 0.0, lam, 1500.0))
            r6000 = outcome((0.55, 0.0, fsep + side * df, 0.0, lam, 6000.0))
            tune.append({'df': df, 'outcome': r6000[0], 'rV(1500)': r1500[3], 'rV(6000)': r6000[3], 'u(6000)': r6000[2]})
        p = [np.log(tune[i]['rV(6000)'] / tune[i + 1]['rV(6000)']) / np.log(10) for i in range(2)]
        seps['lambda=%g near %g' % (lam, trk)] = {'f_sep': fsep, 'track': trk, 'dS side': '+' if side > 0 else '-',
                                                  'fine_tuning': tune, 'local exponents p': p}
        print(json.dumps(seps['lambda=%g near %g' % (lam, trk)]))
        sys.stdout.flush()
    out['separatrices'] = seps
    out['track_fine_tuning_[NOTES 4.3]'] = {'lambda=0 f_sep=0.295245': [1.69e-3, 2.51e-4, 4.29e-5],
                                            'lambda=1': [2.16e-2, 2.98e-3, 4.46e-4], 'p': '0.77-0.85'}

    # ---------------- outcome grid (track: dS 33.0 %, collapse 67.0 %; min dS rV 6.4e-3, median 0.38)
    lams = [-0.12, -0.1, -0.05, 0.0] + [round(0.1 * i, 1) for i in range(1, 11)]
    fgrid = np.round(np.arange(0.10, 1.5001, 0.05), 10)
    jobs = [(0.55, 0.0, float(f), 0.0, lam, 1500.0) for lam in lams for f in fgrid]
    with Pool(4) as P:
        res = P.map(outcome, jobs)
    kinds = [r_[0] for r_ in res]
    rVdS = [r_[3] for r_ in res if r_[0] == 'dS']
    out['grid'] = {'n': len(res), 'dS_frac': kinds.count('dS') / len(res), 'collapse_frac': kinds.count('collapse') / len(res),
                   'f->0_frac': kinds.count('f->0') / len(res), 'undecided': kinds.count('undecided'),
                   'min_rV_dS': float(min(rVdS)), 'median_rV_dS': float(np.median(rVdS)),
                   'track': {'dS': 0.330, 'collapse': 0.670, 'min_rV_dS': 6.4e-3, 'median': 0.38}}
    print(json.dumps(out['grid']))
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'v4_ksv16.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
