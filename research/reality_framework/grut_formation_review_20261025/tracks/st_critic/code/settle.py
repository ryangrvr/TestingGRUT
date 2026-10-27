"""st_critic: quick checks to settle two contradictions between review tracks, using Codex's PUBLISHED arrays as data.

(1) Late growth rate per conformal time: st_rerun N2 (running max in 20-tau windows -> 0.45-0.58) versus st_physics
    P1/P2 (0.6335 / 0.5599; per-k fits within 0.5% of Floquet).  Here: windows of one radial period (pericentre to
    pericentre of |chi_X|), versus st_rerun's 20-tau windows, on the same arrays.
(2) Orbit-scale crossing: st_rerun skeptic (radial d|X| / max|X| within +-8 tau -> 196.3/234.0, 1e-6 max 0.41/0.34)
    versus st_clocks skeptic / st_physics (full |d chi_X| / r_max -> 197.5/215.5, 1e-6 0.34/0.19).
Own code; nothing from Codex imported.  Usage: python3 -I settle.py <packet_dir> <out_json>
"""
import sys, json, math
import numpy as np

P, OUT = sys.argv[1], sys.argv[2]
res = {}
for tag in ['quiet3', 'quiet5']:
    d = np.load(f'{P}/{tag}.npz')
    t = d['times']; bg = d['bg']; ks = d['ks']; amp = d['amplitude_per_radian']; ch = d['charge_fraction_per_radian']
    vec = d['selected_vectors']; vks = d['selected_vector_ks']
    a = np.exp(bg[:, 15]); z = bg[:, 2] + 1j * bg[:, 3]; chiX = np.sqrt(2) * a * z; R = np.abs(chiX)
    eta = np.r_[0, np.cumsum(np.diff(t) * 0.5 * (1 / a[1:] + 1 / a[:-1]))]
    late = t >= 80
    # pericentres: local minima of R on the 0.5-tau grid after tau 80
    idx = np.flatnonzero(late)
    peri = [i for i in idx[1:-1] if R[i] < R[i - 1] and R[i] <= R[i + 1]]
    apo = [i for i in idx[1:-1] if R[i] > R[i - 1] and R[i] >= R[i + 1]]
    row = dict(n_pericentres_after_80=len(peri))
    # (1) growth: one-radial-period windows (between consecutive apocentres, each containing one pericentre)
    def slope_period(series, lo, hi):
        px, py = [], []
        for i0, i1 in zip(apo[:-1], apo[1:]):
            if t[i0] < lo or t[i1] > hi: continue
            seg = slice(i0, i1 + 1); j = int(np.argmax(series[seg])); px.append(eta[seg][j]); py.append(math.log(series[seg][j]))
        return float(np.polyfit(px, py, 1)[0]), len(px)
    def slope_20tau(series, lo, hi):
        px, py = [], []
        for e0 in np.arange(lo, hi, 20.):
            m = (t >= e0) & (t < e0 + 20)
            if m.any():
                j = int(np.argmax(series[m])); px.append(eta[m][j]); py.append(math.log(series[m][j]))
        return float(np.polyfit(px, py, 1)[0]), len(px)
    g = {}
    env_all = np.abs(amp).max(axis=1)
    for name, series in [('maxk_amp', env_all)] + [(f'amp_k{ks[j]:g}', np.abs(amp[:, j])) for j in range(len(ks)) if ks[j] in (2.0, 2.13)] \
            + [(f'charge_k{ks[j]:g}', np.abs(ch[:, j])) for j in range(len(ks)) if ks[j] in (2.0, 2.13)]:
        g[name] = {}
        for lo, hi in [(150, 400), (400, 800), (200, 790), (150, 800)]:
            sp, n1 = slope_period(series, lo, hi); s20, n2 = slope_20tau(series, lo, hi)
            g[name][f'{lo}-{hi}'] = dict(one_radial_period=sp, n_periods=n1, twenty_tau=s20, n_windows=n2)
    # radial period in tau near tau 200 and 780
    pt = np.diff(t[apo]); row['radial_period_tau_near_200'] = float(pt[np.argmin(abs(t[apo][:-1] - 200))])
    row['radial_period_tau_near_780'] = float(pt[np.argmin(abs(t[apo][:-1] - 780))])
    # effective rate between the two published unit samples (Codex brackets) for reference
    C = float((a**2 * np.exp(bg[:, 16]) if False else 0) or 0)
    row['growth_ln_per_eta'] = g
    # (2) orbit-scale measures
    def rmax_period(i):
        # max of R over the radial period containing sample i (apocentre-to-apocentre bracket)
        k0 = max([x for x in apo if x <= i], default=None); k1 = min([x for x in apo if x >= i], default=None)
        if k0 is None or k1 is None: return np.nan
        return float(max(R[k0], R[k1]))
    rmaxP = np.array([rmax_period(i) if late[i] else np.nan for i in range(len(t))])
    rmax8 = np.array([R[(t >= ti - 8) & (t <= ti + 8)].max() for ti in t])
    o = {}
    for seed in [0.03, 1e-6]:
        for j in range(len(ks)):
            if ks[j] not in (2.0, 2.13): continue
            radial = seed * np.abs(amp[:, j]) * R
            for nm, den in [('radial_over_rmax_period', rmaxP), ('radial_over_max_pm8tau', rmax8)]:
                m = radial / den; m[~late] = np.nan
                w = np.flatnonzero(m >= 1)
                o[f'seed{seed:g}_k{ks[j]:g}_{nm}'] = dict(first_sample_ge1=float(t[w[0]]) if len(w) else None, max_by_800=float(np.nanmax(m)))
        for jj, kv in enumerate(vks):
            if kv != 2.13: continue
            full = seed * np.hypot(vec[:, jj, 0], vec[:, jj, 1])
            for nm, den in [('full_over_rmax_period', rmaxP), ('full_over_max_pm8tau', rmax8)]:
                m = full / den; m[~late] = np.nan
                w = np.flatnonzero(m >= 1)
                o[f'seed{seed:g}_k2.13_{nm}'] = dict(first_sample_ge1=float(t[w[0]]) if len(w) else None, max_by_800=float(np.nanmax(m)))
    row['orbit_scale'] = o
    row['ratio_rmax8_over_rmaxPeriod_near_800'] = float(np.nanmin((rmax8 / rmaxP)[t > 700]))
    res[tag] = row
open(OUT, 'w').write(json.dumps(res, indent=1) + '\n')
for tag, row in res.items():
    print(tag, 'periods', row['n_pericentres_after_80'], 'P200', row['radial_period_tau_near_200'], 'P780', row['radial_period_tau_near_780'],
          'rmax8/rmaxP min>700', row['ratio_rmax8_over_rmaxPeriod_near_800'])
    for name, v in row['growth_ln_per_eta'].items():
        print('  ', name, {k: (round(x['one_radial_period'], 4), round(x['twenty_tau'], 4)) for k, x in v.items()})
    for name, v in row['orbit_scale'].items():
        print('  ', name, v)
