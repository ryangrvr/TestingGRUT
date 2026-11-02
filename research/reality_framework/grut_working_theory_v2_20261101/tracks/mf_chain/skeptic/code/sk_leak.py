#!/usr/bin/env python3
"""mf_chain SKEPTIC: independent recomputation of the s -> A_L A_L hidden-heat leak (track E1) and its effect on the
skeptic's strong-holder points.

[D] (re-derived here): the holder u = mu + G_eff c (c = s/s_*) induces A (s/s_*)|Sigma|^2 with A = 2 mu Re(G_eff) c_L
(G_eff = c0 (10i - mu), so |A| = 2 mu^2 c0 c_L).  In the linear basis Sigma = (h + sigma + i pi)/sqrt2 the s -> pi pi amplitude
is the direct A_s plus sigma exchange, A_s [1 + m_sig^2/(m_s^2 - m_sig^2)] = A_s m_s^2/(m_s^2 - m_sig^2) (A_s = A/s_*), so
Gamma = A_s^2 m_s^3/[32 pi (m_sig^2 - m_s^2)^2] (identical pions), h-independent; pi = A_L by equivalence (m_s >> m_A).
s quanta: n_s a^3 = rho_s a^3/m_s = 3 M_P^2 K^2 (f_s/a)/m_s in radiation; BR = Gamma/Gamma_s, Gamma_s = 1.3e-6 H_i (H_i/1.3e11 GeV)
(track's rescaling of mf_hybrid's light-s Yukawa); each A_L decays non-relativistically deep in radiation (the leak A_L
are produced at a ~ 2e3 with p = m_s/2 and become non-relativistic at a ~ 1e7, long before their decay at a ~ 1e9), so it
deposits the same m_A <a_A> as an inflationary vector.  R_leak = 2 n_s a^3 BR dep_A/E_X.
Usage: python3 -B sk_leak.py -> ../results/sk_leak.json
"""
import json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TRK = HERE.parents[2]
MFC = TRK / 'mf_chain' / 'results'
SPD = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/mf_chain_skeptic')
MP_GEV = 2.435e18
MU = 961402368973.181 / 134284273491.30815
MSIG, MS = 40.0, 3.0


def SQ(W):
    return 0.190255 * (W - 0.97494) ** 0.75013 - 0.00428


def fs_over_a(bench):
    h = bench['history']
    a = np.array(h['a'], float); fs = np.array([np.nan if v is None else v for v in h['f_s']], float); t = np.array(h['t'], float)
    sel = (t > 200) & np.isfinite(fs)
    return float(np.median(fs[sel] / a[sel]))


def leak(bench_path, vec_path, chain_path, c0=1.0, e=1e-6):
    b = json.load(open(bench_path)); v = json.load(open(vec_path)); ch = json.load(open(chain_path))
    HiG = float(b['params'].get('H_i_GeV', b.get('H_i_GeV')))
    mP = MP_GEV / HiG
    sst = float(b['params']['s_star_GeV']) / HiG
    K = v['background']['K']; EX = v['background']['E_X']
    dep = [x for x in v['by_e'] if abs(x['e'] - e) < 1e-12][0]['deposit_per_vector']
    fsa = fs_over_a(b)
    ns = 3 * mP ** 2 * K ** 2 * fsa / MS
    Gs = 1.3e-6 * HiG / 1.3e11
    lab = 'e1e-6' if abs(e - 1e-6) < 1e-12 else 'e3e-6'
    W0 = ch[lab]['W']; E_Ec = ch['E_over_Ec']; R0 = W0 / E_Ec - 1
    out = dict(bench=str(bench_path), H_i_GeV=HiG, K=K, E_X=EX, f_s_over_a=fsa, n_s_a3=ns, Gamma_s=Gs, dep_A=dep, c0=c0, e=e,
               S_per_Q_noleak=ch[lab]['S_per_Q'], rows={})
    for name, cL in (('1/16pi^2', 1 / (16 * math.pi ** 2)), ('1/4pi^2', 1 / (4 * math.pi ** 2)),
                     ('ln(MP/MY)/4pi^2', math.log(MP_GEV / (30 * HiG)) / (4 * math.pi ** 2))):
        A = 2 * MU * MU * c0 * cL
        As = A / sst
        G = As ** 2 * MS ** 3 / (32 * math.pi * (MSIG ** 2 - MS ** 2) ** 2)
        BR = G / Gs
        R = 2 * ns * BR * dep / EX
        Wl = E_Ec * (1 + R0 + R)
        out['rows'][name] = dict(c_L=cL, A=A, Gamma_leak=G, BR=BR, R_leak=R, S_per_Q_with_leak=SQ(Wl), R_total_noleak=R0)
    return out


if __name__ == '__main__':
    E = MFC / 'engine'; S = SPD / 'engine'
    cases = [
        ('HFa exported, c0 = 1', TRK / 'mf_hybrid' / 'results' / 'bench_HFa.json', E / 'vectors' / 'vectors_HFa_raw.json', E / 'chain' / 'chain_HFa_raw.json', 1.0),
        ('HFa Pauli sc 2.40e11, c0 = 1', MFC / 'histories' / 'bench_itHFa_P0_it2.json', E / 'vectors' / 'vectors_itHFa_P0_it2.json', E / 'chain' / 'chain_itHFa_P0_it2.json', 1.0),
        ('HFa Pauli sc 8.40e10, c0 = 1.40 (skeptic)', SPD / 'hist' / 'bench_skP0_H084.json', S / 'vectors' / 'vectors_skP0H084_C1.400.json', S / 'chain' / 'chain_skP0H084_C1.400.json', 1.40),
        ('HFa Pauli 8.40e10, c0 = 1.36 (skeptic)', SPD / 'hist' / 'bench_skP0_H084.json', S / 'vectors' / 'vectors_skP0H084_C1.360.json', S / 'chain' / 'chain_skP0H084_C1.360.json', 1.36),
    ]
    res = dict(method=__doc__, cases={})
    for lab, bp, vp, cp, c0 in cases:
        o = leak(bp, vp, cp, c0)
        res['cases'][lab] = o
        print('%s: f_s/a %.3g n_s a^3 %.3g Gamma_s %.3g dep/E_X %.3g S/Q(no leak) %.3f' % (lab, o['f_s_over_a'], o['n_s_a3'], o['Gamma_s'], o['dep_A'] / o['E_X'], o['S_per_Q_noleak']))
        for k, r in o['rows'].items():
            print('    c_L %-16s A %.3g BR %.3g R_leak %.4g -> S/Q %.3f' % (k, r['A'], r['BR'], r['R_leak'], r['S_per_Q_with_leak']))
    json.dump(res, open(HERE.parent / 'results' / 'sk_leak.json', 'w'), indent=1, default=float)
