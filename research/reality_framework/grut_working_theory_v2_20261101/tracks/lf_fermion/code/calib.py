"""Calibrate a local blocking closure to the line fermion runs (non-back-reacting, same recorded background).
Family (model choice): the lf3d local closure P0 with the Pauli capacity multiplied by kappa,
   bl(x) = max(0, 1 - B(x) / (kappa F_sat(E(x)) E(x))),   dB/dt = P - 6 H B,   P = max(G bl, G_pert) phidot^2,
evaluated here in Python along the recorded lines (sites independent; validated against the C tracer at kappa = 1).
kappa = 1 is P0; kappa -> infinity is unblocked tracking (TR).  Fit kappa per species to the fermion energy at
t = 0.3 (homogeneous stage, also on the homogeneous control) and over 0.5-1.1 (fragmented stage).
Claude, lf_fermion (issue #13 context). Usage: python3 -B calib.py RUN -> ../results/calib_<RUN>.json"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
sys.dont_write_bytecode = True
import lines_io, lfcommon as LC
import closure_cmp as CC


def local_tracer(rd, D, lines, kappa, homog=False, t_end=1.1):
    P = CC.par(rd)
    u0, du, ch = CC.tables(rd)
    ug = u0 + du * np.arange(len(ch[0]['G']))
    t, H = D['t'], D['H']
    nt = int(np.searchsorted(t, t_end)) + 1
    F = D['F'][:nt][:, lines]                         # (nt, nl, N, 4)
    if homog:
        F = np.repeat(F.mean(axis=2, keepdims=True), F.shape[2], axis=2)
    s, sd, x, xd = F[..., 0], F[..., 1], F[..., 2], F[..., 3]
    d = s - P['sct']
    V = 0.25 * P['Lv2'] * (x * x - 1) ** 2 + 0.5 * P['Gv2'] * s * s * x * x - P['cbv'] * x + P['c1v'] * d - 0.5 * P['m2c'] * d * d + P['Vc']
    M2 = P['Gv2'] * s * s - P['Lv2']
    with np.errstate(divide='ignore', invalid='ignore'):
        loop = np.where((s > P['sct']) & (M2 != 0), M2 * M2 * (np.log(np.abs(M2) / P['mu2']) - 1.5) / (64 * math.pi ** 2 * P['vH'] ** 2), 0.0)
    V = V + loop
    Et = 0.5 * (sd * sd + xd * xd) + V
    Es = np.maximum(0.5 * sd * sd + 0.5 * P['Gv2'] * x * x * (s - P['s0t']) ** 2, 0.0)
    Ex = np.maximum(Et - Es, 0.0)
    out = np.zeros((nt, 4))
    for j, c in enumerate(ch):
        E = Ex if c['w'] else Es
        vel2 = xd * xd if c['w'] else sd * sd
        with np.errstate(divide='ignore'):
            uE = np.log(np.maximum(E, 1e-300))
        G = np.interp(uE, ug, c['G'])
        cap = kappa * np.interp(uE, ug, c['F']) * E
        B = np.zeros(E.shape[1:])
        TE = np.zeros(E.shape[1:])
        for n in range(nt - 1):
            dt = t[n + 1] - t[n]
            bl = np.where(cap[n] > 0, np.maximum(0.0, 1 - B / np.where(cap[n] > 0, cap[n], 1)), 0.0)
            dE = np.where(E[n] > 1e-300, np.maximum(G[n] * bl, c['Gp']), 0.0) * vel2[n] * dt
            B = B * math.exp(-6 * H[n] * dt) + dE
            TE = TE * math.exp(-4 * H[n] * dt) + dE
            out[n + 1, j] = TE.mean()
    return t[:nt], out


def main():
    run = sys.argv[1] if len(sys.argv) > 1 else 'REC_P0'
    rd = LC.RUNS / run
    D = lines_io.read_lines(rd)
    cmpd = json.loads((LC.TRACK / 'results' / ('cmp_%s.json' % run)).read_text())
    tg = np.array(cmpd['t'])
    SP = {'N1': 0, 'N2': 1, 'N3': 2, 'N4': 3}
    res = dict(run=run, label='C', note=__doc__, species={})
    kgrid = [1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0, 32.0, 64.0, 1e6]
    for nm, d in cmpd['species'].items():
        j = SP[nm]
        lines = d['lines']
        q = np.array(d['Q_mean'])
        rows = {}
        for kap in kgrid:
            tt, te = local_tracer(rd, D, lines, kap)
            rows[kap] = np.interp(tg, tt, te[:, j])
        hrow = {}
        if nm in cmpd['homog']:
            h = cmpd['homog'][nm]
            for kap in kgrid:
                tt, te = local_tracer(rd, D, [h['line']], kap, homog=True)
                hrow[kap] = np.interp(tg, tt, te[:, j])
        # check kappa = 1 against the C tracer (lambda = 0)
        chk = np.array(rows[1.0]) / np.array(d['tracer_mean']['0.0'])

        def kfit(target, table, i):
            ks = np.array(kgrid)
            vals = np.array([table[k][i] for k in kgrid])
            if target <= vals[0]:
                return float('nan') if target < 0.98 * vals[0] else 1.0
            if target >= vals[-1]:
                return float('inf')
            return float(np.exp(np.interp(target, vals, np.log(ks))))   # vals increase with kappa
        kap_t = [kfit(q[i], rows, i) for i in range(len(tg))]
        kap_h = [kfit(cmpd['homog'][nm]['Q_hom'][i], hrow, i) for i in range(len(tg))] if hrow else None
        res['species'][nm] = dict(t=list(tg), kappa_fit_inhom=kap_t, kappa_fit_homog=kap_h,
                                  check_kappa1_over_C_tracer=list(chk),
                                  tracer_by_kappa={str(k): list(v) for k, v in rows.items()},
                                  tracer_by_kappa_homog={str(k): list(v) for k, v in hrow.items()}, Q=list(q))
        print(nm, 'check k=1/C', np.round(chk[[3, 5, 9, 15, 21]], 3))
        print('   kappa(t) inhom:', ' '.join('%.2f:%.2f' % (tg[i], kap_t[i]) for i in range(3, len(tg), 2)))
        if kap_h:
            print('   kappa(t) homog:', ' '.join('%.2f:%.2f' % (tg[i], kap_h[i]) for i in range(3, len(tg), 2)))
        print('   TR (kappa 1e6) / Q:', np.round(np.array(rows[1e6])[[5, 9, 13, 17, 21]] / q[[5, 9, 13, 17, 21]], 2))
    Path(__file__).resolve().parents[1].joinpath('results', 'calib_%s.json' % run).write_text(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()
