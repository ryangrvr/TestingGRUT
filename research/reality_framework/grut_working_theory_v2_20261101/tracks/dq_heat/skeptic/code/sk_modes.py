#!/usr/bin/env python3
"""dq_heat SKEPTIC (Claude, own code, written fresh): longitudinal-vector modes in cosmic time, Stueckelberg variable.

Independent of dq_heat/code/vector_modes.py in formulation, time variable, initial state, background handling and read-out:
  * variable: the Stueckelberg (Goldstone-gauge) angle theta_k (A_L = 0 gauge, A_0 integrated out),
        S = 1/2 int dt [ p(t) theta_dot^2 - q(t) theta^2 ],  p = a^3 F^2 k^2/(k^2 + a^2 e^2 F^2),  q = a F^2 k^2;
  * cosmic time t = tau (H_i units), pair (theta, pi = p theta_dot);
  * initial state: the exact massless de Sitter Goldstone mode u = aF theta = e^{-ik eta}(1 - i/(k eta))/sqrt(2k),
    eta = -1/a (H = 1), started at k/aH = 50 (dq_heat: adiabatic vacuum at k/aH = 1000 in conformal time);
  * background: splines in t of ln a(t) and F^2(t) from the recorded rptri_dq run (dense tau <= 200, coarse to 2000; no
    conformal-time reconstruction), then exact radiation a(t) = (a_end^2 + 2K(t - t_end))^{1/2}, F^2 = 72|Sigma|^2 + 2chi^2/a^2;
  * read-out: (i) the cosmic adiabatic number (|pi|^2/p + q|theta|^2)/(2 Omega) - 1/2, Omega^2 = q/p, and (ii) the conformal
    canonical number of u = sqrt(P_eta) theta (P_eta = p/a), each at two late times.
Usage: python3 -B sk_modes.py e variant [kmin kmax nk | list]  -> ../results/sk_modes_<variant>_<e>.json
"""
import json, math, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
DQ = HERE.parents[1] / 'results'
SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2')
NPZ = SP / 'dq_heat' / 'runs' / 'm1r_background.npz'
BGJ = json.load(open(DQ / 'background.json'))
F_INF = BGJ['prep']['F_inf']
XFRAC = BGJ['prep']['X_part_of_F2']
FV2 = 72.0 * BGJ['late']['absS_late'] ** 2
CHI2 = BGJ['late']['chi2_mean_late']
_B = {}


def bg():
    if _B:
        return _B
    d = np.load(NPZ)
    td, tc = d['d_t'], d['c_t']
    sel = tc > td[-1] + 1e-9
    t = np.concatenate([td, tc[sel]])
    la = np.log(np.concatenate([d['d_a'], d['c_a'][sel]]))
    z = np.concatenate([d['d_z'], d['c_z'][sel]])
    F2 = np.concatenate([d['d_F2'], d['c_F2'][sel]])
    X2 = 2 * np.abs(z[:, 0]) ** 2
    H = np.concatenate([d['d_H'], d['c_H'][sel]])
    _B.update(t=t, la=CubicSpline(t, la), F2full=CubicSpline(t, F2), F2noX=CubicSpline(t, F2 - X2),
              t_end=float(t[-1]), a_end=float(np.exp(la[-1])), H_end=float(H[-1]))
    _B['K'] = _B['a_end'] ** 2 * _B['H_end']
    return _B


def background(t, variant):
    """a, F^2 at cosmic time t (t < 0: de Sitter at H = 1 with the prepared-state F_inf)."""
    B = bg()
    if t < 0:
        a = math.exp(t)
        F2 = F_INF ** 2 if variant == 'full' else F_INF ** 2 * (1 - XFRAC)
        return a, F2
    if t <= B['t_end']:
        a = math.exp(float(B['la'](t)))
        F2 = float(B['F2full'](t)) if variant == 'full' else float(B['F2noX'](t))
        return a, F2
    a = math.sqrt(B['a_end'] ** 2 + 2 * B['K'] * (t - B['t_end']))
    F2 = FV2 + (2 * CHI2 / a ** 2 if variant == 'full' else 0.0)
    return a, F2


def t_of_a_late(a):
    B = bg()
    return B['t_end'] + (a * a - B['a_end'] ** 2) / (2 * B['K'])


def pq(t, k, e, variant):
    a, F2 = background(t, variant)
    M2 = a * a * e * e * F2
    p = a ** 3 * F2 * k * k / (k * k + M2)
    q = a * F2 * k * k
    return a, F2, p, q


def numbers(t, y, k, e, variant):
    a, F2, p, q = pq(t, k, e, variant)
    th = y[0] + 1j * y[1]; pi = y[2] + 1j * y[3]
    Om = math.sqrt(q / p)
    n_cos = (abs(pi) ** 2 / p + q * abs(th) ** 2) / (2 * Om) - 0.5
    # conformal canonical u = sqrt(P) theta, P = p/a ; u' = d/deta = a d/dt
    P = p / a
    h = 1e-5 * max(1.0, abs(t))
    Pp = (pq(t + h, k, e, variant)[2] / background(t + h, variant)[0] - pq(t - h, k, e, variant)[2] / background(t - h, variant)[0]) / (2 * h)
    dsP_eta = a * Pp / (2 * math.sqrt(P))
    thd = pi / p
    u = math.sqrt(P) * th
    up = math.sqrt(P) * a * thd + dsP_eta * th
    w = math.sqrt(k * k + a * a * e * e * F2)
    n_conf = (abs(up) ** 2 + w * w * abs(u) ** 2) / (2 * w) - 0.5
    wr = (th * np.conj(pi) - np.conj(th) * pi).imag          # = 1 for the canonical normalisation
    return n_cos, n_conf, wr, a


def run_mode(args):
    k, e, variant = args
    B = bg()
    a0 = k / 50.0
    t0 = math.log(a0)
    eta0 = -1.0 / a0
    a, F2 = background(t0, variant)
    F = math.sqrt(F2)
    u = np.exp(-1j * k * eta0) * (1 - 1j / (k * eta0)) / math.sqrt(2 * k)
    up = -1j * k * u + np.exp(-1j * k * eta0) * (1j / (k * eta0 ** 2)) / math.sqrt(2 * k)
    th = u / (a * F)
    dth_deta = (up - a * u) / (a * F)          # a'/a = a in eta = -1/a
    thd = dth_deta / a
    _, _, p, q = pq(t0, k, e, variant)
    pi = p * thd
    y = np.array([th.real, th.imag, pi.real, pi.imag])
    m = e * math.sqrt(FV2)
    a_star = math.sqrt(B['K'] / m)
    if k <= 0.3:
        a_reads = [30 * a_star, 60 * a_star]
    else:
        a_reads = [1.3 * B['a_end'], 1.6 * B['a_end']]
    t_reads = [t_of_a_late(x) for x in a_reads]

    def rhs(t, y):
        _, _, p, q = pq(t, k, e, variant)
        return [y[2] / p, y[3] / p, -q * y[0], -q * y[1]]

    segs = [t0, 0.0, 20.0, 200.0, B['t_end']] + t_reads
    segs = sorted(set(s for s in segs if s >= t0))
    out = {}
    wr0 = numbers(t0, y, k, e, variant)[2]
    for s1, s2 in zip(segs[:-1], segs[1:]):
        ms = 0.005 if 0 <= s1 < 20 else (0.05 if 20 <= s1 < 200 else (0.5 if s1 < B['t_end'] else np.inf))
        sc = float(np.max(np.abs(y)))
        sol = solve_ivp(rhs, (s1, s2), y, method='DOP853', rtol=1e-10, atol=1e-14 * sc, max_step=ms)
        if not sol.success:
            raise RuntimeError(sol.message)
        y = sol.y[:, -1]
        if s2 in t_reads:
            out[s2] = numbers(s2, y, k, e, variant)
    r1, r2 = out[t_reads[0]], out[t_reads[1]]
    return dict(k=k, e=e, variant=variant, n_cos=[r1[0], r2[0]], n_conf=[r1[1], r2[1]], wronskian=[wr0, r1[2], r2[2]],
                a_reads=a_reads)


def main():
    e = float(sys.argv[1]); variant = sys.argv[2]
    if len(sys.argv) > 3 and sys.argv[3] == 'dq':
        ks = json.load(open(DQ / 'vector_modes.json'))['kgrid']
        ks = [k for k in ks if k >= float(sys.argv[4])] if len(sys.argv) > 4 else ks
    elif len(sys.argv) > 3 and sys.argv[3] == 'list':
        ks = [float(x) for x in sys.argv[4].split(',')]
    else:
        ks = list(np.geomspace(float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])))
    t0 = time.time()
    bg()
    with Pool(2) as pool:
        rows = pool.map(run_mode, [(k, e, variant) for k in ks], chunksize=1)
    tag = sys.argv[6] if len(sys.argv) > 6 else (sys.argv[3] if len(sys.argv) > 3 else 'grid')
    fn = RES / ('sk_modes_%s_%.1e_%s.json' % (variant, e, tag))
    json.dump(dict(method=__doc__, e=e, variant=variant, K=bg()['K'], rows=rows, seconds=time.time() - t0), open(fn, 'w'), indent=1)
    for r in rows:
        print('k=%.4g n_cos=%.6e/%.6e n_conf=%.6e/%.6e wr=%s' % (r['k'], r['n_cos'][0], r['n_cos'][1], r['n_conf'][0], r['n_conf'][1],
              ['%.2e' % (w - 1) for w in r['wronskian']]))
    print('seconds %.1f -> %s' % (time.time() - t0, fn))


if __name__ == '__main__':
    main()
