#!/usr/bin/env python3
"""Critic track: recompute the load-bearing numbers on which the review tracks disagree.

Own analysis code (critic). Inputs are read-only:
  - raw histories of the independent_solver track (its own spectral solver; no Codex code),
  - result JSONs of the spectrum and conversion tracks,
  - Codex packet JSONs read via `git archive` into the scratchpad (ff49e5d, bb211d1).

usage: critic_checks.py OUT.json
Paths are set below (scratchpad layout of this review).
"""
import json, sys
from math import erf, sqrt, pi
import numpy as np
from scipy.interpolate import PchipInterpolator

SCR = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review'
REPO = '/home/user/TestingGRUT/research/reality_framework/grut_formation_review_20261025/tracks'
IND = SCR + '/independent_solver/runs/'
FF = SCR + '/critic/ff49e5d/research/reality_framework/grut_causal_receiver_20261024/'
BB = SCR + '/critic/bb211d1/research/reality_framework/grut_core_conversion_20261023/'

RAD = [30, 60, 120, 240, 480]
COLS = ['tau', 'E', 'Q', 'Erec', 'Qrec', 'E_local_sum', 'S0']
for R in RAD:
    COLS += [f'Q{R}', f'Q{R}_trap', f'E{R}', f'Ech_out{R}', f'Ene_out{R}', f'pad{R}']
C = {k: i for i, k in enumerate(COLS)}


def hist(name):
    a = np.fromfile(IND + name + '_history.bin', '<f8')
    return a.reshape(-1, len(COLS))


def at(h, T, key):
    i = int(np.argmin(abs(h[:, 0] - T)))
    assert abs(h[i, 0] - T) < 1e-6
    return h[i]


def escape(h, T):
    r = at(h, T, None)
    return dict(Q=r[C['Q']] - r[C['Q240']] + r[C['Qrec']], E=r[C['E']] - r[C['E240']] + r[C['Erec']])


out = {}

# ---- A. Mechanism of the old-damper bias (independent_solver raw runs, same h=0.5, dt=1/64) ----
box, old = hist('A_h05_dt2'), hist('D_olddamper_dt2')
Q0 = box[0, C['Q']]
n = min(len(box), len(old)); t = box[:n, 0]
dq = 100 * (old[:n, C['Q60']] - box[:n, C['Q60']]) / Q0           # percentage points of Q0
pert, base = hist('A_h05_dt_pert'), hist('A_h05_dt')
m2 = min(len(pert), len(base)); t2 = base[:m2, 0]
dp = 100 * (pert[:m2, C['Q60']] - base[:m2, C['Q60']]) / Q0
A = {'pointwise_pp': {}, 'window_mean_pp': {}, 'perturbation_1e-6_window_mean_pp': {}, 'perturbation_1e-6_window_max_abs_pp': {}}
for T in [1200, 1400, 1600, 1800, 2000, 2500, 3000, 4000, 6000, 8000, 9600, 10800, 12000]:
    A['pointwise_pp'][T] = float(dq[int(np.argmin(abs(t - T)))])
for lo, hi in [(1000, 2000), (2000, 3000), (3000, 4000), (4000, 6000), (6000, 8000), (8000, 9600), (9600, 12000)]:
    w = (t > lo) & (t <= hi); A['window_mean_pp'][f'{lo}-{hi}'] = float(dq[w].mean())
    w2 = (t2 > lo) & (t2 <= hi)
    A['perturbation_1e-6_window_mean_pp'][f'{lo}-{hi}'] = float(dp[w2].mean())
    A['perturbation_1e-6_window_max_abs_pp'][f'{lo}-{hi}'] = float(abs(dp[w2]).max())
# return times for a wave emitted at r=60, reflected near r~300 in the old cubic damper (270-360)
A['return_time_after_emission'] = {f'p={p}': 2 * (300 - 60) / (p / sqrt(1 + p * p)) for p in [0.02, 0.05, 0.1, 0.5, 1.0]}
A['label'] = '[C] critic analysis of independent_solver raw histories (its own solver; old damper profile [Q-code])'
out['A_old_damper_mechanism'] = A

# ---- B. Escape totals through r=240 by tau=12000, and their trajectory sensitivity ----
eb, ep = escape(base, 12000), escape(pert, 12000)
B = {'indep_h05_dt1_64': escape(box, 12000), 'indep_h05_dt1_32': eb, 'indep_h05_dt1_32_amp_plus_1e-6': ep,
     'sensitivity_Q': ep['Q'] - eb['Q'], 'sensitivity_E': ep['E'] - eb['E']}
spec = json.load(open(REPO + '/spectrum/results/spectrum_summary.json'))['convergence_R240_W12000']
qvals = {'indep_solver(h.5,dt1/64)': B['indep_h05_dt1_64']['Q'],
         'spectrum_F2dt': spec['F2dt']['Q0_minus_Qin240_T12000'],
         'spectrum_F3h': spec['F3h']['Q0_minus_Qin240_T12000'],
         'Codex_P0125_trace(spectrum analysis of reproduction data)': spec['P0125']['Qpass_240_T12000'],
         'Codex_P00625_trace(spectrum analysis of reproduction data)': spec['P00625']['Qpass_240_T12000']}
qv = np.array(list(qvals.values()))
B['converged_Q_values'] = qvals
B['Q_consolidated'] = dict(mean=float(qv.mean()), half_range=float((qv.max() - qv.min()) / 2))
ev = np.array([B['indep_h05_dt1_64']['E'], 25478.8])            # spectrum F3h energy flux, quoted from its NOTES
B['E_consolidated'] = dict(mean=float(ev.mean()), half_range=float((ev.max() - ev.min()) / 2), note='spectrum F3h value 25478.8 quoted [Q]')
i025 = json.load(open(FF + 'I025_RESULTS.json'))
B['Codex_I025'] = dict(Q=23680.857132, E=25618.620689, note='quoted from REPORT.md [Q]')
B['Codex_I025_rel_error'] = dict(Q=23680.857132 / qv.mean() - 1, E=25618.620689 / ev.mean() - 1)
B['ratio_mQ_over_E_consolidated'] = float(qv.mean() / ev.mean())
out['B_escape_totals'] = B

# ---- C. Retention 9600-12000: agreement of independent discretizations; late window offset ----
def jl(p):
    return json.load(open(p))
ret = {'Codex_P0125': jl(FF + 'P0125_RESULTS.json')['core60_late_mean_fraction'] * 100,
       'Codex_P00625': jl(FF + 'P00625_RESULTS.json')['core60_late_mean_fraction'] * 100,
       'spectrum_F2dt': spec['F2dt']['Q60_9600_12000'] * 100,
       'spectrum_F3h': spec['F3h']['Q60_9600_12000'] * 100,
       'indep_solver_h05_dt1_64(nodal)': float(np.mean(box[(box[:, 0] > 9600) & (box[:, 0] <= 12000), C['Q60']]) / Q0 * 100)}
rv = np.array(list(ret.values()))
Cc = dict(values_percent=ret, mean=float(rv.mean()), half_range_pp=float((rv.max() - rv.min()) / 2),
          old_executable_remote_damper_percent={'R720_h.125': 35.375, 'R720_h.0625': 35.383, 'R1440_h.125': 35.335},
          note='old-executable values quoted from reproduce NOTES [Q]; they carry the old layer fast-wave reflection')
P025 = jl(FF + 'P025_RESULTS.json')['core60_late_mean_fraction'] * 100
P025r = jl(FF + 'P025_remote_RESULTS.json')['core60_late_mean_fraction'] * 100
off = rv.mean() - P025
Cc['late_window'] = dict(P025_9600_12000=P025, P025_remote_19200_24000=P025r, offset_pp=off,
                         offset_corrected=P025r + off,
                         indep_dt1_32_plus_dt_offset=34.729 + (35.357 - 35.346),
                         spectrum_F1_plus_offset=34.664 + (spec['F3h']['Q60_9600_12000'] - spec['F1']['Q60_9600_12000']) * 100,
                         label='[E] constant-offset assumption; inputs quoted from the tracks/Codex [Q]')
out['C_retention'] = Cc

# ---- D. Antiparticle number across converged runs ----
nx = {'spectrum_F2dt': spec['F2dt']['NXbar'], 'spectrum_F3h': spec['F3h']['NXbar'],
      'Codex_P0125_arrays': spec['P0125']['NXbar'], 'Codex_P00625_arrays': spec['P00625']['NXbar'],
      'indep_solver_projection': 255.0}
nv = np.array(list(nx.values()))
out['D_antiparticles'] = dict(values=nx, mean=float(nv.mean()), half_range=float((nv.max() - nv.min()) / 2),
                              fraction_of_Q0=float(nv.mean() / Q0),
                              breathing_period_from_spectrum_nu=2 * pi / 0.181)

# ---- E. Rate-table resonance: PCHIP (log) of Codex's 32 nodes vs the conversion track's dense scan ----
rows = jl(BB + 'RATE_RESULTS.json')['rows']
p = np.array([r['p_over_m'] for r in rows]); s = np.array([r['sigma_v_conversion_cm3_s'] for r in rows])
o = np.argsort(p); p, s = p[o], s[o]
tab = np.array(jl(REPO + '/conversion/results/rate_scan.json')['dense_scan']['table']); tab = tab[np.argsort(tab[:, 0])]
f = np.linspace(0.001, 0.708, 200001)
lg = np.exp(PchipInterpolator(p, np.log(s))(f)); ln = PchipInterpolator(p, s)(f); dn = np.interp(f, tab[:, 0], tab[:, 1])
E = {}
for lo, hi in [(0.45, 0.56), (0.001, 0.708)]:
    w = (f >= lo) & (f <= hi)
    E[f'{lo}-{hi}'] = dict(ratio_pchip_log=float(np.trapezoid(lg[w], f[w]) / np.trapezoid(dn[w], f[w])),
                           ratio_pchip_lin=float(np.trapezoid(ln[w], f[w]) / np.trapezoid(dn[w], f[w])))
E['dense_peak'] = dict(p=float(tab[np.argmax(tab[:, 1]), 0]), sigma_v=float(tab[:, 1].max()))
E['codex_optical_depth_ratio_reported'] = 54.61 / 49.74
out['E_rate_resonance'] = E

# ---- F. free_soft packet spectral coverage (Codex packet: |u|^2 ~ exp(-(r-c)^2/w^2), k0=0.02, w=100) ----
w, k0 = 100.0, 0.02
sk = 1 / (w * sqrt(2))
Phi = lambda x: 0.5 * (1 + erf(x / sqrt(2)))
out['F_packet_coverage'] = dict(sigma_k=sk, frac_below_0p0163=Phi((0.0163 - k0) / sk), frac_below_0p01=Phi((0.01 - k0) / sk),
                                k_min_return_estimate=(700 + 2 * 90 * 1.5 - 8) / 60000.0,
                                note='path 300->1000 layer, x1.5 mean CFS stretch over 90+90, back to r<992 [E]')

json.dump(out, open(sys.argv[1], 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
