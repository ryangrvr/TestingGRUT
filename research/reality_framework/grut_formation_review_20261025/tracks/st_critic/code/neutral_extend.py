"""st_critic: REPRODUCTION-TYPE extension of Codex's own 2cbc7f1 neutral (delta-phi) control from tau 80 to tau 800.

Uses Codex's modes.c compiled with -DNEUTRAL_SECTOR (via Codex's neutral.run, mass and build, imported read-only from
the b73 tree), Codex's two Gaussian seed conventions and Codex's energy-ratio definition, on the same quiet histories
b3_g0.5 = (3,.5) and b5_g1 = (5,1), with the b73 timing grid to tau 800.  Our own additions: the longer time range,
the threshold bookkeeping and a late growth fit.  Purpose: settle st_physics P8 versus its skeptic (omitted
delta-phi resonance, mu_phi ~ 0.2 at k ~ 6-7) without new physics choices.
Usage: GRUT_MODE_WORKDIR=<scratch> python3 -B neutral_extend.py <reality_framework_dir> <out_json>
"""
import sys, json, math
from pathlib import Path
import numpy as np
RF = Path(sys.argv[1]); OUT = Path(sys.argv[2])
sys.path.insert(0, str(RF / 'grut_reheating_modes_20261026' / 'code'))
import neutral as N  # Codex, read-only
base = json.loads((RF / 'grut_reheating_modes_20261026' / 'MODE_RESULTS.json').read_text())
tgrid = np.load(RF / 'grut_source_timing_20261026' / 'quiet3.npz')['times']
ks = np.unique(np.r_[np.array(base['ks']), np.geomspace(.8, 20, 97)])
exe = N.build(); neutral = exe.parent / 'neutral'
import subprocess
subprocess.run(['gcc', '-O3', '-std=c11', '-DNEUTRAL_SECTOR', '-I', str(exe.parent), str(RF / 'grut_reheating_modes_20261026' / 'code' / 'modes.c'), '-lm', '-o', str(neutral)], check=True)
pubj = json.loads((RF / 'grut_reheating_modes_20261026' / 'NEUTRAL_RESULTS.json').read_text())
pub = {row['tag']: row['intervals'][-1]['maximum_conditional_phi_energy_fraction'] for row in pubj['rows']}
out = {}
for r in base['rows']:
    if r['tag'] not in ('b3_g0.5', 'b5_g1'): continue
    bg, TT, check = N.run(neutral, r, ks, tgrid, suffix='_st_critic_800')
    m0 = N.mass(bg[0], r['beta']); w0 = np.sqrt(ks * ks + m0); d0 = m0 / (2 * w0 * w0)
    covs = [np.array([[[1 / (2 * w), 0], [0, w / 2]] for w in w0]),
            np.array([[[1 / (2 * w), -d / (2 * w)], [-d / (2 * w), (w * w + d * d) / (2 * w)]] for w, d in zip(w0, d0)])]
    ratios = []; nmax = []; kpeak = []; a = np.exp(bg[:, 15])
    for n, v in enumerate(bg):
        w2 = ks * ks + N.mass(v, r['beta']); w = np.sqrt(w2)
        EX = np.exp(4 * v[15]) * ((v[8]**2 + v[9]**2) + .001 * (v[2]**2 + v[3]**2)**2)
        vals = []
        for C0 in covs:
            C = TT[n] @ C0 @ TT[n].transpose(0, 2, 1); vals.append((C[:, 1, 1] + w2 * C[:, 0, 0]) / (2 * w) - .5)
        ratios.append([float(np.trapezoid(ks**3 * w * val / (2 * np.pi**2), np.log(ks)) / EX) for val in vals])
        nmax.append(float(max(vals[0].max(), vals[1].max()))); kpeak.append(float(ks[int(np.argmax(vals[0]))]))
    ratios = np.array(ratios); nmax = np.array(nmax)
    det = float(np.abs(np.linalg.det(TT[-1]) - 1).max())
    eta = np.r_[0, np.cumsum(np.diff(tgrid) * 0.5 * (1 / a[1:] + 1 / a[:-1]))]
    def first(level, j):
        w = np.flatnonzero(ratios[:, j] >= level); return float(tgrid[w[0]]) if len(w) else None
    late = tgrid >= 300
    sl = float(np.polyfit(eta[late], np.log(nmax[late]), 1)[0]) / 2  # occupation ~ e^{2 mu eta}
    row = dict(tag=r['tag'], check=check, symplectic_det_err_at_800=det,
               energy_fraction_at=dict({str(tt): ratios[int(np.argmin(abs(tgrid - tt)))].tolist() for tt in (80, 134, 147.5, 200, 300, 400, 500, 600, 772, 782, 800)}),
               first_sampled_fraction={str(L): [first(L, 0), first(L, 1)] for L in (0.01, 0.1, 1.0)},
               max_occupation_at={str(tt): nmax[int(np.argmin(abs(tgrid - tt)))] for tt in (80, 134, 300, 400, 600, 800)},
               peak_k_at={str(tt): kpeak[int(np.argmin(abs(tgrid - tt)))] for tt in (80, 300, 800)},
               late_mu_phi_from_max_occupation_tau300_800=sl,
               max_fraction_through_80_own=ratios[tgrid <= 80].max(axis=0).tolist(),
               codex_2cbc7f1_max_fraction_through_80=pub[r['tag']])
    out[r['tag']] = row
    print(json.dumps(row), flush=True)
OUT.write_text(json.dumps(dict(label='reproduction-type extension of Codex neutral control (Codex engine/seeds/definitions), our time range and bookkeeping', ks_n=len(ks), k_range=[float(ks[0]), float(ks[-1])], rows=out), indent=1) + '\n')
