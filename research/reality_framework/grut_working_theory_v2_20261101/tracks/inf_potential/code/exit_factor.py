#!/usr/bin/env python3
"""inf_potential (Claude, own code): the hidden vectors' horizon-exit factor f_eff on the benchmark inflaton histories,
by the dq_heat SKEPTIC's method (dq_heat/skeptic/code/sk_exit_path.py), applied to OUR computed paths.

Each mode k < 1 (H_i units, a = 1 at the slice) left the horizon N(k) e-folds before the slice, k = e^{-N} h_r(N).  Its
frozen Goldstone power is (H/2 pi F)^2 at exit, so dq_heat's static-slice per-mode spectrum is multiplied by
    f(N) = [ h_r(N) F(1, 1) / F(s(N), h_r(N)) ]^2 ,
with F re-solved at the exit state: u = mu + s (10i - mu), R = 6 (2 - eps_H) h_r^2 (H_i units), physical M_Y, h, mu, kappa,
lambdas unchanged.  s(N) = phi_I(N)/phi_I(slice) and h_r(N) = H(N)/H(slice) come from our background (inf_common), with
phi_I = phi - phi_min the holder field.  Modes k >= 1 keep f = 1.  S/Q from dq_entropy's closed form
S/Q = 0.190255 (W - 0.97494)^0.75013 - 0.00428 with W = W0 (1 + R_total) [Q dq_entropy, as the skeptic].

Read-only imports: dq_heat/skeptic/code/dq_model_copy.py (stationary-point solver), dq_heat/results/{vector_modes,
heat_table}.json.  Nothing there is modified.
CAVEAT [A]: the per-mode spectrum and the release (N_X, E_X, K, the mediator heat) are the M1r static-slice ones; a real
benchmark changes them (the chain track's job).  This is the exit factor only, on the M1r release.
Labels: [C] given the path and the M1r release; the paths are [C] from our potentials; the transfer is [E].
Usage: python3 exit_factor.py T1 E1 HT6t ...   Output: ../results/exit_factor.json (merged).
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
HERE = Path(__file__).resolve().parent
TRK = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(TRK / 'dq_heat' / 'skeptic' / 'code'))
import inf_common as ic
import dq_model_copy as fc            # read-only import (stationary solver)

DQ = TRK / 'dq_heat' / 'results'
RES = HERE.parent / 'results'
Q = np.array([1.0, 3.0, 6.0])
WTHR = 49.05678346694957               # W at the bulk threshold S/Q = 3.4716 [Q skeptic]


def SQ(W):
    return 0.190255 * (W - 0.97494) ** 0.75013 - 0.00428


def model_from_bench(name, Hi_index=0):
    d = json.load(open(RES / ('bench_%s.json' % name)))
    rec = d['by_Hi'][Hi_index]
    kind = {'T1': 'T', 'E1': 'E', 'HT6t': 'HT', 'HT6': 'HT', 'HT4': 'HT', 'HT8': 'HT', 'PL2': 'PL'}[name]
    return ic.Model(kind, **rec['params']), rec


CONVS = {'S1': (1.0, 1.0), 'Send': (0.0, 0.0), 'Cexit': (3.6, 1.0)}   # (N_cal, N_geo) e-folds before eps_H = 1
K_M1R = 1.98581                          # a^2 H of the M1r release (a = 1 at its slice) [Q dq_heat]


def paths(m, conv, Ngrid):
    """s = phi_I/phi_I(N_cal) (the holder calibration point: mu_eff = 10i there) and h_r = H/H(N_geo) (H_i = H at the
    geometric slice, where a = 1), at N e-folds before the geometric slice."""
    bg = ic.background(m, 75.0)
    Nc, Ngeo = CONVS[conv]
    pmin = m.phi_min()
    fc_, _ = ic.bg_at(bg, Nc)
    fs, fNs = ic.bg_at(bg, Ngeo); hs, _ = ic.h_of(m, fs, fNs)
    out = []
    for N in Ngrid:
        f, fN = ic.bg_at(bg, N + Ngeo)
        h, eps = ic.h_of(m, f, fN)
        out.append(((f - pmin) / (fc_ - pmin), h / hs, eps))
    return out


def main(names):
    p = fc.params(0.5)
    ref = fc.stationary(fc.U_PREP, p, R=12.0)
    F11 = float(np.sqrt(2 * np.sum(Q ** 2 * np.abs(ref['z']) ** 2)))
    vm = json.load(open(DQ / 'vector_modes.json'))
    hb = json.load(open(DQ / 'heat_table.json'))['budgets']['portals_vanish_at_MY_benchmark']
    W0 = hb['W0']; rows_e = {r['e']: r for r in hb['rows']}
    Ng = np.linspace(0.0, 8.0, 65)
    try:
        allres = json.load(open(RES / 'exit_factor.json'))
    except Exception:
        allres = dict(method=__doc__, F11=F11, results={})
    allres['method'] = __doc__
    for name in names:
        m, rec = model_from_bench(name)
        for conv in ('S1', 'Send', 'Cexit'):
            P = paths(m, conv, Ng)
            # K shift [A]: the relic modes sit at k_* ~ (K m_A)^{1/2} with K = a^2 H (a = 1 at the geometric slice), so the
            # M1r per-mode spectrum is moved by k -> k (K_bench/K_M1r)^{1/2}; K_bench = K_end e^{2 N_geo} (Gamma_I = 0.5 H_i)
            Kb = rec['K_end_gamma05'] * math.exp(2 * CONVS[conv][1])     # converged (radiation) value from calibrate.py
            kfac = math.sqrt(Kb / K_M1R)
            fN = []; info = []; guess = None
            for (s, hr, eps), N in zip(P, Ng):
                u = p['mu'] + (fc.U_PREP - p['mu']) * s
                R = 6.0 * (2.0 - eps) * hr * hr
                st = None
                for g in (guess, None):
                    try:
                        st = fc.stationary(u, p, R=R, guess=g)
                        if st['grad'] < 1e-8:
                            break
                    except Exception:
                        st = None
                guess = st['z']
                F = float(np.sqrt(2 * np.sum(Q ** 2 * np.abs(st['z']) ** 2)))
                f = (hr * F11 / F) ** 2
                fN.append(f)
                info.append(dict(N=float(N), s=float(s), h_r=float(hr), eps=float(eps), u_re=float(np.real(u)), u_im=float(np.imag(u)),
                                 amps=[float(x) for x in np.abs(st['z'])], F=F, f=f, grad=float(st['grad']),
                                 m2min_over_H2=float(sorted(st['eig'])[1] / hr ** 2)))
            fN = np.array(fN); hN = np.array([x[1] for x in P])
            lnk_of_N = -Ng + np.log(hN)

            def f_of_k(k):
                if k >= hN[0]:
                    return 1.0
                N = float(np.interp(math.log(k), lnk_of_N[::-1], Ng[::-1]))
                return float(np.interp(N, Ng, fN))
            by_e = []
            for c in vm['cases']:
                if c['variant'] != 'full':
                    continue
                e = c['e']; r = rows_e[e]
                ks = np.array([mm['k'] for mm in c['spectrum']]); dn = np.array([mm['dN_dlnk'] for mm in c['spectrum']])
                lo = ks <= 1.0
                fk = np.array([f_of_k(k * kfac) for k in ks])
                f_eff = float(np.trapezoid((dn * fk)[lo], np.log(ks[lo])) / np.trapezoid(dn[lo], np.log(ks[lo])))
                Rt = r['R_total'] - r['R_vec_k_le_1'] + f_eff * r['R_vec_k_le_1']
                by_e.append(dict(e=e, f_eff=f_eff, R_total=Rt, S_per_Q=SQ(W0 * (1 + Rt))))
            es = sorted(x['e'] for x in by_e); le = np.log(es)
            lW = [math.log(W0 * (1 + [x for x in by_e if x['e'] == e][0]['R_total'])) for e in es]
            fn = lambda x: np.interp(x, le, lW) - math.log(WTHR)
            ec = math.exp(brentq(fn, le[0], le[-1])) if fn(le[0]) > 0 > fn(le[-1]) else None
            b6 = [x for x in by_e if x['e'] == 1e-6][0]
            key = '%s_%s' % (name, conv)
            allres['results'][key] = dict(benchmark=name, H_i_GeV=rec['H_i_GeV'], slice=conv, N_cal=CONVS[conv][0], N_geo=CONVS[conv][1],
                                          K_bench=Kb, k_shift=kfac, path=info, by_e=by_e, e_crit_bulk=ec,
                                          f_eff_e1e6=b6['f_eff'], SQ_e1e6=b6['S_per_Q'],
                                          branch_e1e6='free X' if b6['S_per_Q'] > 3.4716191256166837 else 'Q-matter')
            print('%-10s s(N=1,2,3.8,6)=%s h_r=%s f=%s m2min/H2(N=3.8)=%.2f | e=1e-6: f_eff=%.3f S/Q=%.3f -> %s ; e_crit=%s' % (
                key, ['%.3f' % np.interp(x, Ng, [i['s'] for i in info]) for x in (1, 2, 3.8, 6)],
                ['%.4f' % np.interp(x, Ng, hN) for x in (1, 2, 3.8)], ['%.3f' % np.interp(x, Ng, fN) for x in (1, 2, 3.8, 6)],
                np.interp(3.8, Ng, [i['m2min_over_H2'] for i in info]), b6['f_eff'], b6['S_per_Q'],
                allres['results'][key]['branch_e1e6'], '%.2e' % ec if ec else None), flush=True)
            json.dump(allres, open(RES / 'exit_factor.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main(sys.argv[1:])
