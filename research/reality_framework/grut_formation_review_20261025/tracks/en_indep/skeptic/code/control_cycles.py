"""Skeptic check (en_indep_skeptic): orbit-phase analysis of ALL seven published Codex curves.

Reads only the published FIELD_CURVES.npz (git archive of 7d9c93b). analyze.py stores every control
at stride 10 (dtau = 0.5) and halfstep48 at stride 2 (dtau = 0.05), N48 at stride 1 (dtau = 0.05).
For each history: tau=800 value, range over tau 600-800, and full-cycle means between successive
maxima of the energy-outside-mean curve (cycle mean weighted uniformly in conformal time s = eta).
"""
import sys, json
import numpy as np

src, outp = sys.argv[1], sys.argv[2]
z = np.load(src)
TAGS = ['N32', 'N48', 'halfstep48', 'N64', 'seed2_48', 'cut10_48', 'box16_64']
KEYS = dict(outside='energy_outside_spatial_mean_fraction', meanQ='spatial_mean_charge_fraction',
            phimix='phi_plus_interaction_energy_fraction', Ephi='Ephi_fraction',
            Xvar='X_field_fluctuation_variance_fraction', contrast='charge_contrast_RMS')
out = dict(note='Published FIELD_CURVES.npz; sampling per tag given; cycle = successive maxima of outside-mean, mean uniform in eta', tags={})
for tag in TAGS:
    tau = z[tag + '_tau']; eta = z[tag + '_eta']
    keep = np.r_[np.diff(tau) > 1e-9, True]          # drop the duplicated tiny final step
    tau, eta = tau[keep], eta[keep]
    d = {k: z[tag + '_' + v][keep] for k, v in KEYS.items()}
    row = dict(dtau=float(np.median(np.diff(tau))), at800={k: float(v[-1]) for k, v in d.items()})
    w = tau >= 600
    row['range_600_800'] = {k: [float(v[w].min()), float(v[w].max())] for k, v in d.items()}
    # cycle detection on outside-mean (smooth enough at dtau 0.5)
    o = d['outside']; sel = np.flatnonzero(tau >= 250)
    pk = [i for i in sel[1:-1] if o[i] > o[i - 1] and o[i] >= o[i + 1]]
    # keep only well separated maxima (> 15 tau apart)
    pk2 = []
    for i in pk:
        if not pk2 or tau[i] - tau[pk2[-1]] > 15: pk2.append(i)
        elif o[i] > o[pk2[-1]]: pk2[-1] = i
    row['peak_taus'] = [float(tau[i]) for i in pk2]
    per = np.diff(eta[pk2]); row['period_eta_mean'] = float(per.mean()); row['period_eta_last3'] = [float(x) for x in per[-3:]]
    cyc = []
    for a, b in zip(pk2[:-1], pk2[1:]):
        ss = eta[a:b + 1]
        cyc.append(dict(tau=[float(tau[a]), float(tau[b])], **{k: float(np.trapezoid(v[a:b + 1], ss) / (ss[-1] - ss[0])) for k, v in d.items()}))
    row['last_cycle'] = cyc[-1]; row['prev_cycle'] = cyc[-2]
    # mean over tau 600-800 uniform in eta (not whole cycles; secondary)
    ss = eta[w]
    row['eta_mean_600_800'] = {k: float(np.trapezoid(v[w], ss) / (ss[-1] - ss[0])) for k, v in d.items()}
    # phase at 800 relative to last maximum
    row['phase_at_800'] = float((eta[-1] - eta[pk2[-1]]) / per[-1])
    out['tags'][tag] = row

def spread(f):
    v = [f(out['tags'][t]) for t in TAGS]; return [min(v), max(v)]
summ = {}
for k in KEYS:
    summ[k] = dict(at800=spread(lambda r: r['at800'][k]), last_cycle_mean=spread(lambda r: r['last_cycle'][k]),
                   N48_swing_600_800=out['tags']['N48']['range_600_800'][k],
                   per_tag_at800={t: out['tags'][t]['at800'][k] for t in TAGS},
                   per_tag_last_cycle={t: out['tags'][t]['last_cycle'][k] for t in TAGS})
out['summary'] = summ
json.dump(out, open(outp, 'w'), indent=1)
for k in ['outside', 'meanQ', 'phimix', 'Ephi']:
    print(k)
    for t in TAGS:
        r = out['tags'][t]
        print(f"  {t:11s} 800={r['at800'][k]:.5g} last={r['last_cycle'][k]:.5g} prev={r['prev_cycle'][k]:.5g} rng600={r['range_600_800'][k][0]:.4g}-{r['range_600_800'][k][1]:.4g} phase800={r['phase_at_800']:.2f} lastcyc={r['last_cycle']['tau']}")
