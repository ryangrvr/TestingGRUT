#!/usr/bin/env python3
"""st_clocks_skeptic (Claude, own analysis of CODEX DATA): how much do Codex's 'unit modulation' crossing samples depend on
normalising by the INSTANTANEOUS |X| of a very eccentric homogeneous orbit?

Codex's amplitude criterion is (chi_X . dchi_X)/|chi_X|^2 per radian x seed = 1 at the instantaneous radius (code
timing_common.responses).  The homogeneous X orbit has |X|max/|X|min ~ 8-9 within one radial period, so the criterion fires
near pericentre.  For the field equations (cubic force), second-order terms relative to the leading ones scale as
|dchi|/A with A the orbit amplitude [A], so we also report
   full / orbit max:  |dchi_X| (both components, Codex's saved tangent vectors at k = 0.02, 2.13, 5.8) / max|chi_X| over
                      a +-W a window in tau (one radial period ~ 1.8 a; W = 0.9, 1.5, 3 give identical results).
The local-charge criterion dn/n_X is likewise amplified by ~1/eps (orbit charge fraction ~0.15) relative to |dchi|/A.
Nothing here is a new computation of the dynamics; it is a re-reading of Codex's saved linear solutions [C on Codex data].
Output: ../results/criterion_skeptic.json
"""
import json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / 'results' / 'criterion_skeptic.json'
SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad')
CODEX = SP / 'b73/research/reality_framework/grut_source_timing_20261026'
BG = json.loads((HERE.parent / 'results' / 'background_skeptic.json').read_text())['histories']


def first(t, m):
    w = np.flatnonzero(m >= 1)
    return float(t[w[0]]) if len(w) else None


def main():
    out = {}
    for tag in ('quiet3', 'quiet5'):
        z = np.load(CODEX / (tag + '.npz'))
        t = z['times']; bg = z['bg']; a = np.exp(bg[:, 15])
        X = bg[:, 2] + 1j * bg[:, 3]
        q = math.sqrt(2) * a[:, None] * np.c_[X.real, X.imag]
        r = np.linalg.norm(q, axis=1)
        V = z['selected_vectors']; ks = list(z['ks'])
        rmax = np.array([r[(t >= ti - 0.9 * ai) & (t <= ti + 0.9 * ai)].max() for ti, ai in zip(t, a)])
        rmin = np.array([r[(t >= ti - 0.9 * ai) & (t <= ti + 0.9 * ai)].min() for ti, ai in zip(t, a)])
        h = dict(orbit_rmax_over_rmin={str(tt): float(rmax[np.argmin(abs(t - tt))] / rmin[np.argmin(abs(t - tt))])
                                       for tt in (80, 134, 300, 772)})
        j = list(z['selected_vector_ks']).index(2.13)
        radial = np.einsum('ti,ti->t', V[:, j, :2], q) / r
        full = np.linalg.norm(V[:, j, :2], axis=1)
        h['k'] = 2.13
        h['check_radial_equals_codex_amplitude'] = float(np.max(np.abs(radial / r - z['amplitude_per_radian'][:, ks.index(2.13)])))
        rows = {}
        C = BG[tag]['C_own']; ts = BG[tag]['tau_star_own']
        for seed in (0.03, 1e-6, 1e-9):
            inst = np.abs(radial) / r * seed
            om = full / rmax * seed
            row = dict(codex_criterion_first_sample=first(t, inst), orbit_amplitude_criterion_first_sample=first(t, om),
                       orbit_amplitude_running_max_at_800=float(np.max(om)))
            # value of the orbit-amplitude measure at Codex's own crossing sample
            if row['codex_criterion_first_sample'] is not None:
                i = int(np.argmin(abs(t - row['codex_criterion_first_sample'])))
                row['orbit_amplitude_measure_at_codex_sample'] = float(om[i])
                row['running_max_orbit_measure_to_codex_sample'] = float(np.max(om[:i + 1]))
            if row['orbit_amplitude_criterion_first_sample'] is None:
                # [E] extrapolate the running max of the orbit measure vs a over tau 500-800
                run = np.maximum.accumulate(om); s = t >= 500
                p = np.polyfit(a[s], np.log(run[s]), 1)
                ac = -p[1] / p[0]
                row['E_extrapolated'] = dict(slope_per_a=float(p[0]), a_cross=float(ac), tau_cross=float(ac * ac / (2 * C) + ts))
            rows['seed_%g' % seed] = row
        h['seeds'] = rows
        out[tag] = h
        print(tag, json.dumps(h, indent=1))
    OUT.write_text(json.dumps(out, indent=2) + '\n')


if __name__ == '__main__':
    main()
