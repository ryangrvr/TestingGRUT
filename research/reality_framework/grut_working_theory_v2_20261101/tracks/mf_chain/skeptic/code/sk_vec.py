#!/usr/bin/env python3
"""mf_chain SKEPTIC: independent check of the longitudinal-vector occupation n_k on the hybrid histories (the "vector heat
with the new histories").  Different code path from vmodes.c:
  * F(t) from the skeptic's own release (sk_release.py), F^2 = 2(|X|^2 + 9|Y|^2 + 36|Sigma|^2), continued as
    72|Sigma_late|^2 + 2 <|X|^2 a^2>/a^2 after t = 100; background a, H from the bench via two-piece splines, then radiation;
  * e-folds N = ln a as the independent variable, scipy DOP853 (vmodes: cosmic time, hand-coded DP5, linear interpolation);
  * the same single-Higgs-direction action S = 1/2 int dt [p thetadot^2 - q theta^2], p = a^3F^2k^2/(k^2 + a^2e^2F^2),
    q = aF^2k^2 (dq_heat / its skeptic), started as the exact massless de Sitter Goldstone at k/aH = 50 with F = F_s;
  * read-out n = (|pi|^2/p + q|theta|^2)/(2 Omega) - 1/2 at Omega/H = 100 and 200.
Usage: python3 -B sk_vec.py BENCH VECJSON e k1 [k2 ...]  -> prints mine vs the track's spectrum (nearest stored k)
Labels: [C] given the model choices.
"""
import json, math, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import sk_release as sr


def background_F(bench_path, T=100.0):
    o, (tt, Y, D) = sr.run(bench_path, T=T, tag='vec')
    x = Y[0] + 1j * Y[3]; y = Y[1] + 1j * Y[4]; s = Y[2] + 1j * Y[5]
    F2 = 2 * (np.abs(x) ** 2 + 9 * np.abs(y) ** 2 + 36 * np.abs(s) ** 2)
    a = D['a']
    late = tt >= 0.75 * T
    chi2 = float(np.mean(np.abs(x[late]) ** 2 * a[late] ** 2)); Fv2 = float(np.mean(72 * np.abs(s[late]) ** 2))
    bench = json.load(open(bench_path))
    bg = sr.Background(bench)
    return dict(t=tt, lna=np.log(a), F2=F2, chi2=chi2, Fv2=Fv2, bg=bg, Fs=math.sqrt(F2[0]))


def mode(k, e, B, KS=50.0, AD=100.0):
    Fs2 = B['F2'][0]
    lnaF = CubicSpline(B['lna'], B['F2'])
    Nd = float(B['lna'][-1])
    bg = B['bg']
    # N -> H map for N > 0 (from the bench background); tabulate H on a fine N grid via t
    tg = np.r_[np.linspace(0, 1, 2001)[:-1], np.geomspace(1.0, 1e12, 6000)]
    lag = np.array([bg.at(t)[0] for t in tg]); Hg = np.array([bg.at(t)[1] for t in tg])
    keep = np.r_[True, np.diff(lag) > 1e-12]
    HN = CubicSpline(lag[keep], np.log(Hg[keep]))

    def pqH(N):
        if N < 0:
            H = 1.0; f2 = Fs2
        else:
            H = math.exp(float(HN(N)))
            f2 = float(lnaF(N)) if N <= Nd else B['Fv2'] + 2 * B['chi2'] * math.exp(-2 * N)
        a = math.exp(N)
        p = a ** 3 * f2 * k * k / (k * k + a * a * e * e * f2)
        q = a * f2 * k * k
        return p, q, H

    def rhs(N, yv):
        p, q, H = pqH(N)
        th = yv[0] + 1j * yv[1]; pi = yv[2] + 1j * yv[3]
        dth = pi / (p * H); dpi = -q * th / H
        return [dth.real, dth.imag, dpi.real, dpi.imag]

    N0 = math.log(k / KS)
    a0 = math.exp(N0); eta = -1.0 / a0
    u = np.exp(-1j * k * eta) * (1 - 1j / (k * eta)) / math.sqrt(2 * k)
    du = np.exp(-1j * k * eta) * (-1j * k * (1 - 1j / (k * eta)) + 1j / (k * eta * eta)) / math.sqrt(2 * k)   # d/deta
    F = math.sqrt(Fs2)
    th = u / (a0 * F)
    dth_deta = (du - u * a0) / (a0 * F)          # d(u/a)/deta = u'/a - u a'/a^2, a' = a^2 (H = 1)
    thdot = dth_deta / a0
    p0, q0, _ = pqH(N0)
    pi = p0 * thdot
    W = (th * np.conj(pi) - np.conj(th) * pi)
    y0 = [th.real, th.imag, pi.real, pi.imag]

    def number(N, yv):
        p, q, H = pqH(N)
        th = yv[0] + 1j * yv[1]; pi = yv[2] + 1j * yv[3]
        Om = math.sqrt(q / p)
        return (abs(pi) ** 2 / p + q * abs(th) ** 2) / (2 * Om) - 0.5, Om / H

    def ev(N, yv, thr):
        p, q, H = pqH(N)
        return math.sqrt(q / p) / H - thr
    e1 = lambda N, yv: ev(N, yv, AD); e1.terminal = False; e1.direction = 1
    e2 = lambda N, yv: ev(N, yv, 2 * AD); e2.terminal = True; e2.direction = 1
    Nend = 40.0
    sol = solve_ivp(rhs, (N0, Nend), y0, method='DOP853', rtol=1e-10, atol=1e-30, events=[e1, e2],
                    max_step=0.02, first_step=1e-4)
    out = dict(k=k, e=e, wronskian=complex(W).imag)
    if len(sol.y_events[0]):
        out['n1'], _ = number(sol.t_events[0][0], sol.y_events[0][0]); out['N1'] = float(sol.t_events[0][0])
    if len(sol.y_events[1]):
        out['n2'], _ = number(sol.t_events[1][0], sol.y_events[1][0]); out['N2'] = float(sol.t_events[1][0])
    out['nfev'] = int(sol.nfev)
    return out


if __name__ == '__main__':
    bench, vj, e = sys.argv[1], sys.argv[2], float(sys.argv[3])
    ks = [float(x) for x in sys.argv[4:]]
    B = background_F(bench)
    v = json.load(open(vj))
    sp = [x for x in v['by_e'] if abs(x['e'] - e) < 1e-12 * e][0]['spectrum']
    res = []
    for k in ks:
        m = mode(k, e, B)
        tr = min(sp, key=lambda r: abs(math.log(r['k'] / k)))
        m['track_k'] = tr['k']; m['track_n'] = tr['n_k']
        res.append(m)
        print('k %.5g: mine n = %.5g (Omega/H=200: %.5g) | track (k=%.5g) %.5g | ratio %.4f | W %.3g | N_read %.2f | nfev %d' % (
            k, m.get('n1', float('nan')), m.get('n2', float('nan')), tr['k'], tr['n_k'], m.get('n2', float('nan')) / tr['n_k'],
            m['wronskian'], m.get('N2', float('nan')), m['nfev']), flush=True)
    print(json.dumps(dict(Fs=B['Fs'], Fv2=B['Fv2'], chi2=B['chi2'])))
