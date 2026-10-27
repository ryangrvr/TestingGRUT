#!/usr/bin/env python3
"""st_response skeptic (Claude, own code): Y_L rest-clock fractions at the response crossings, recomputed with OUR fine
a(tau) (sk_tangent background, 0.002 grid, trapezoid) and OUR tail (closed form via scipy gammaincc), and at the
corrected first continuous 10% crossing.  gamma = Gamma_L/H_i is quoted [Q] from Codex's SOURCE_CLOCKS.json.
Writes ../results/sk_clock.json."""
import json
from pathlib import Path
import numpy as np
from scipy.special import gammaincc, gamma as G
import sk_analyze as K

D = json.loads((K.CODEX / 'SOURCE_CLOCKS.json').read_text())
out = {}
for tag, job, taus in (('quiet3', 'q3', [56.458, 83.0, 120.0, 134.0, 133.602]), ('quiet5', 'q5', [55.168, 55.292790025819336, 134.5, 147.5, 147.322])):
    M = K.load(job); t, a, H = M['t'], M['a'], M['H']
    h = next(x for x in D['histories'] if x['tag'] == tag)
    gam = h['sources']['YL']['gamma_over_Hi']
    C = a[-1] ** 2 * H[-1]; off = t[-1] - a[-1] ** 2 / (2 * C)
    f = gam * a * np.exp(-gam * t)
    cum = np.r_[0, np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(t))]
    u = gam * (t[-1] - off)
    tail = np.sqrt(2 * C / gam) * np.exp(-gam * off) * G(1.5) * gammaincc(1.5, u)
    total = cum[-1] + tail
    rows = []
    codex = {c['tau']: c for c in h['sources']['YL']['at_response_crossings'] if c['phase_seed'] == 0.03}
    for tau in taus:
        Ef = float(np.interp(tau, t, cum) / total); Nf = float(-np.expm1(-gam * tau))
        c = codex.get(tau)
        rows.append(dict(tau=tau, number_fraction=Nf, energy_fraction=Ef,
                         codex_number=c and c['decayed_number_fraction_rest_clock'],
                         codex_energy=c and c['cumulative_conformal_energy_fraction_prescribed_source']))
    out[tag] = dict(gamma=gam, late_C=float(C), total_mean_a=float(total), codex_total=h['sources']['YL']['mean_a_per_decay'], rows=rows)
(Path(__file__).resolve().parents[1] / 'results' / 'sk_clock.json').write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
