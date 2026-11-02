#!/usr/bin/env python3
"""inf_chain SKEPTIC (Claude, own code): independent longitudinal-vector count for the quiet T-model TB1 (alpha = 0.3)
in the matter era, to test the track's vmodes.c + fluid continuation outside its M1r (radiation) validation.

Method (different code path from the track: Python, own background, own F(t) model):
  * background: my own T-model field from S1 (sk_cmb_fluid.background/oscillate machinery, Gamma = 0) for t <= t1 = 150,
    then an exact matter era a^{3/2} = a1^{3/2} + (3/2)(x0/3)^{1/2}(t - t1) (Gamma_cal = 7.5e-9 is irrelevant before the
    readout, t ~ 1e5-1e6 << 1/Gamma);
  * Higgs norm: F^2 = F_s^2 before the slice, F^2 = F_v^2 + 2 chi^2/a^2 after it (F_v = 240, chi^2 = 6538.8 from the
    track's background summary; the early-release transient of F is ignored: superhorizon the angle is frozen);
  * Stueckelberg-angle action S = 1/2 int dt [p thetadot^2 - q theta^2], p = a^3 F^2 k^2/(k^2 + a^2 e^2 F^2), q = a F^2 k^2;
    exact massless de Sitter Goldstone (H = 1, F_s) at t = 0 ('static' exit, f = 1), adiabatic number read when
    omega/H > 200 (and 100 as a check);
  * N_A(static) = int dlnk k^3 n_k / 2 pi^2, compared with the track's N_A_k_le_1_static (114.07 at e = 1e-6).
Output ../results/sk_vmodes.json.  Label [C].
"""
import json, math, sys, pathlib
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import sk_cmb_fluid as sk

FV, CHI2, FS = 240.0, 6538.8038808057745, 271.2857435012218


def make_background(alpha=0.3, t1=150.0):
    m = sk.TModel(alpha)
    sol, Nend = sk.background(m, m.phi_start(80.0))
    f, fN = sol.sol(Nend - 1.0); Hs, eps = sk.H_of(m, f, fN)
    sc = 1.0 / Hs ** 2

    def rhs(t, y):
        f, p, lna = y
        rho = 0.5 * p * p + sc * m.V(f)
        return [p, -3 * math.sqrt(rho / 3) * p - sc * m.dV(f), math.sqrt(rho / 3)]
    s = solve_ivp(rhs, (0, t1), [f, fN, 0.0], method='DOP853', rtol=1e-11, atol=1e-14, dense_output=True)
    # x0 from the last cycles
    mI = math.sqrt(m.m2() * sc); T = 2 * math.pi / mI
    tt = np.linspace(t1 - 6 * T, t1, 3001); Y = s.sol(tt)
    x0 = float(np.mean((0.5 * Y[1] ** 2 + sc * np.array([m.V(v) for v in Y[0]])) * np.exp(3 * Y[2])))
    a1 = math.exp(s.y[2, -1])

    def aH(t):
        if t <= t1:
            f, p, lna = s.sol(t)
            return math.exp(lna), math.sqrt((0.5 * p * p + sc * m.V(f)) / 3)
        a = (a1 ** 1.5 + 1.5 * math.sqrt(x0 / 3) * (t - t1)) ** (2 / 3)
        return a, math.sqrt(x0 / 3) * a ** -1.5
    return aH, x0, t1


def mode(k, e, aH, t1, wh=200.0):
    def pq(t):
        a, H = aH(t)
        F2 = FV ** 2 + 2 * CHI2 / a ** 2
        m2 = e * e * F2
        p = a ** 3 * F2 * k * k / (k * k + a * a * m2)
        q = a * F2 * k * k
        return p, q, a, H, m2
    # initial state at t = 0: dS Goldstone with F_s
    x = k
    pi0 = np.exp(1j * x) * (x + 1j) / (k * math.sqrt(2 * k))
    pid0 = -1j * x * x * np.exp(1j * x) / (k * math.sqrt(2 * k))
    p0 = FS ** 2 * k * k / (k * k + e * e * FS ** 2)
    th0 = pi0 / FS; P0 = p0 * pid0 / FS

    def rhs(t, y):
        p, q, a, H, m2 = pq(t)
        return [y[2] / p, y[3] / p, -q * y[0], -q * y[1]]

    def num(t, y):
        p, q, a, H, m2 = pq(t)
        w = math.sqrt(q / p)
        return ((y[2] ** 2 + y[3] ** 2) / p + q * (y[0] ** 2 + y[1] ** 2)) / (2 * w) - 0.5, w / H

    y = [th0.real, th0.imag, P0.real, P0.imag]
    t = 0.0; out = {}
    # march in segments, check omega/H
    seg = [t1] + list(np.geomspace(t1 * 1.5, 1e9, 120))
    for tn in seg:
        s = solve_ivp(rhs, (t, tn), y, method='DOP853', rtol=1e-9, atol=1e-30)
        y = s.y[:, -1]; t = tn
        n, wH = num(t, y)
        for lev in (100.0, wh):
            if wH > lev and lev not in out:
                out[lev] = n
        if wh in out:
            break
    return out.get(wh), out.get(100.0)


def NA(e, aH, t1, ks):
    rows = []
    for k in ks:
        n2, n1 = mode(k, e, aH, t1)
        rows.append(dict(k=k, n200=n2, n100=n1, dN_dlnk=k ** 3 * n2 / (2 * math.pi ** 2)))
    lk = np.log(ks); d = np.array([r['dN_dlnk'] for r in rows])
    sel = ks <= 1.0 + 1e-12
    N_le1 = float(np.trapezoid(d[sel], lk[sel])); N_all = float(np.trapezoid(d, lk))
    d1 = np.array([r['k'] ** 3 * r['n100'] / (2 * math.pi ** 2) for r in rows])
    return dict(e=e, N_A_k_le_1_static=N_le1, N_A_all_static=N_all, N_A_k_le_1_static_readout100=float(np.trapezoid(d1[sel], lk[sel])),
                spectrum=rows)


def main():
    aH, x0, t1 = make_background()
    print('x0 =', x0, flush=True)
    track = json.load(open(pathlib.Path(__file__).resolve().parent.parent.parent / 'results' / 'vectors_TB1_g7.512e-09.json'))
    out = dict(method=__doc__, x0=x0, rows=[])
    for e in (1e-6, 1e-5):
        ks = np.geomspace(1e-4, 10.0, 51)
        r = NA(e, aH, t1, ks)
        tr = [x for x in track['by_e'] if abs(x['e'] - e) / e < 1e-9]
        r['track_N_A_k_le_1_static'] = tr[0]['N_A_k_le_1_static'] if tr else None
        r['ratio_mine_over_track'] = (r['N_A_k_le_1_static'] / r['track_N_A_k_le_1_static']) if tr else None
        out['rows'].append(r)
        print('e=%.0e: N_A(k<=1, static) mine %.4g (readout100 %.4g) track %s ratio %s ; all k %.4g' % (
            e, r['N_A_k_le_1_static'], r['N_A_k_le_1_static_readout100'], r['track_N_A_k_le_1_static'], r['ratio_mine_over_track'],
            r['N_A_all_static']), flush=True)
    if len(out['rows']) == 2:
        a, b = out['rows']
        out['slope_dlnN_dlne'] = math.log(b['N_A_k_le_1_static'] / a['N_A_k_le_1_static']) / math.log(10.0)
        print('d ln N_A / d ln e =', out['slope_dlnN_dlne'])
    p = pathlib.Path(__file__).resolve().parent.parent / 'results' / 'sk_vmodes.json'
    json.dump(out, open(p, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
