#!/usr/bin/env python3
"""mf_chain SKEPTIC: how robust is "free X at every e, no e_crit" on the hybrid histories?  With the independent release
(sk_release.py) I vary model choices the track held fixed, one at a time, on the track's own histories:
  * the spectator's late mass m_s (HFa: m_s = g_s |chi|, so g_s scales; the brief allows m_s in {2, 3, 4} H_i);
  * the holder held (s frozen, the holder never switches off): separates the background's Y_L drive from the turn-on;
  * the holder amplitude c0 (mu_eff(slice) = mu + (10i - mu) c0), as mf_ideal's pocket scan.
For each variant the late (75-100) means give eps, E/E_c, N_YL, N_YH, N_sigma; the heat ledger uses the base history's own
vector run (N_A x deposit, which depends on K and only weakly on the release: |dN_A| <= 1% across c0 = 0.99-1.01 in the
track's runs), so R_vec = (N_A dep)/E_X; R_Y = (N_YL + N_YH) m_L <a_L>/E_X; R_casc = (N_YH + 2 N_sigma) dep/E_X.
Approximate (late-mean instead of a^-2-frozen eps: -2% in eps for HFa_raw) [C on the inherited ledger]; the branch-changing
variants are re-run through the full chain (sk_chainrun.py).
Usage: python3 -B sk_scan.py BENCH VECTAG ENGDIR OUTTAG [variants...]  variant = ms:<m> | hold | c0:<c>
"""
import json, math, sys
from pathlib import Path
import numpy as np
import sk_release as sr

SQB, SQ2 = 3.4716191256166837, 2.9283778774213727
M_L, G_L = 22.24859546128699, 1.005903273746766e-07


def SQ(W):
    return 0.190255 * (W - 0.97494) ** 0.75013 - 0.00428


def heat(lm, vec):
    K = vec['background']['K']
    EX = lm['EX']; eps = abs(lm['eps'])
    E_Ec = eps ** (-4 / 3)
    aL = math.sqrt(math.pi * K / (2 * G_L))
    RY = (lm['NYL'] + lm['NYH']) * M_L * aL / EX
    es = np.array([x['e'] for x in vec['by_e']])
    ND = np.array([x['N_A'] * x['deposit_per_vector'] for x in vec['by_e']]); DA = np.array([x['deposit_per_vector'] for x in vec['by_e']])
    rows = []
    eg = np.geomspace(1e-7, 1e-2, 101)
    for e in eg:
        nd = math.exp(float(np.interp(math.log(e), np.log(es), np.log(ND), right=None)))
        if e > es[-1]:
            sl = math.log(ND[-1] / ND[-2]) / math.log(es[-1] / es[-2]); nd = ND[-1] * (e / es[-1]) ** sl
        da = math.exp(float(np.interp(math.log(e), np.log(es), np.log(DA))))
        if e > es[-1]:
            da = DA[-1] * (e / es[-1]) ** -0.5
        Rv = nd / EX; Rc = (lm['NYH'] + 2 * lm['Nsig']) * da / EX
        rows.append((e, Rv, Rc, SQ(E_Ec * (1 + Rv + RY + Rc))))
    sq = {('%g' % e): s for e, _, _, s in rows if abs(math.log10(e) - round(math.log10(e))) < 1e-9}
    f = lambda e: float(np.interp(math.log(e), np.log(eg), [r[3] for r in rows]))
    ec = None
    for r1, r2 in zip(rows, rows[1:]):
        if (r1[3] - SQB) * (r2[3] - SQB) < 0:
            ec = math.exp(math.log(r1[0]) + (math.log(r2[0]) - math.log(r1[0])) * (r1[3] - SQB) / (r1[3] - r2[3])); break
    return dict(K=K, E_Ec=E_Ec, R_Y=RY, SQ_1e6=f(1e-6), SQ_3e6=f(3e-6), SQ_1e5=f(1e-5), SQ_inf=SQ(E_Ec * (1 + RY)), e_crit=ec,
                R_vec_1e6=rows[[abs(r[0] - 1e-6) for r in rows].index(min(abs(r[0] - 1e-6) for r in rows))][1])


def main():
    bench, vectag, engdir, outtag = sys.argv[1:5]
    vec = json.load(open(Path(engdir) / 'vectors' / ('vectors_%s.json' % vectag)))
    out = dict(bench=bench, vectors=vectag, variants={})
    outp = Path(__file__).resolve().parent.parent / 'results' / ('sk_scan_%s.json' % outtag)
    try:
        out = json.load(open(outp))
    except Exception:
        pass
    for v in sys.argv[5:]:
        kw = {}
        if v.startswith('ms:'):
            kw['ms_scale'] = float(v[3:]) / 3.0
        elif v == 'hold':
            kw['hold_s'] = True
        elif v.startswith('c0:'):
            kw['c0'] = float(v[3:])
        o, _ = sr.run(bench, T=100.0, tag=v, **kw)
        lm = o['late_mean']
        h = heat(lm, vec)
        out['variants'][v] = dict(late_mean=lm, ledger=h, gauss=o['gauss_max'], seconds=o['seconds'])
        print('%-10s eps %.4f E/Ec %6.2f NYL %8.4g NYH %7.4g Nsig %7.3g | R_Y %7.3f R_vec(1e-6) %6.3f | S/Q 1e-6 %7.3f 3e-6 %7.3f 1e-5 %7.3f inf %7.3f | e_crit %s' % (
            v, abs(lm['eps']), h['E_Ec'], lm['NYL'], lm['NYH'], lm['Nsig'], h['R_Y'], h['R_vec_1e6'], h['SQ_1e6'], h['SQ_3e6'], h['SQ_1e5'],
            h['SQ_inf'], ('%.3g' % h['e_crit']) if h['e_crit'] else '-'), flush=True)
        json.dump(out, open(outp, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
