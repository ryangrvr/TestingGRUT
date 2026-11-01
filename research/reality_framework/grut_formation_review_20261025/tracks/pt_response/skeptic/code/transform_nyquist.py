"""SKEPTIC (pt_response_skeptic), INDEPENDENT code, no Codex code.
What is the 0.1-0.4% canonical-sum residual of the published response arrays?  Codex's finite transform (delays 0..6
step .04, cos^2 taper, trapezoid) is applied with
  (a) Codex's omega grid, +-20 step .02 (control: must reproduce pt_response's max |Z-1| = 0.003966),
  (b) the full Nyquist grid of the Delta = .04 sampling, +-pi/Delta = +-78.54, step .02,
to the published normal arrays and to our low-noise shell traces (meas_traces.npz).
On the Nyquist grid the sum over omega of omega*sin(omega t_j) is the exact discrete kernel -(-1)^j/(j Delta^2), so
Z = -(2/Delta) sum_j (-1)^j g_j taper_j trap_j / j is an identity of the SAMPLED data (= 1 for any sampled sine with
|omega| < pi/Delta); a residual on the +-20 grid is therefore spectral weight beyond |omega| = 20, i.e. line width.
Test lines on both grids: undamped sin(wt)/w, Lorentzian-damped exp(-G t/2) sin(wt)/w, Gaussian-damped
exp(-(G t)^2/2) sin(wt)/w.  A Lorentzian tail (rho ~ G/omega^2) gives a residual ~ (2/pi) G/omega_max on the +-20 grid,
so the published residual sizes the effective width:  G_eff ~ (pi/2) * 20 * |Z-1|.
Usage: python3 -I transform_nyquist.py <packet dir> <meas_traces.npz> <out json>
"""
import sys, json
import numpy as np
P = sys.argv[1].rstrip('/') + '/'; MEAS = sys.argv[2]; OUT = sys.argv[3]
D = .04; t = np.arange(151) * D
taper = np.cos(np.pi * t / 12) ** 2; trap = np.full(151, D); trap[[0, -1]] *= .5
grids = dict(codex=np.linspace(-20, 20, 2001), nyquist=np.linspace(-np.pi / D, np.pi / D, 7855))

def Zsum(g, w):
    dw = w[1] - w[0]
    TR = np.exp(-1j * w[:, None] * t[None, :])
    rho = -2 * (TR @ (g * trap * taper)).imag
    return float(np.sum(w * rho) * dw / (2 * np.pi))

def Zexact_discrete(g):
    # closed form of the Nyquist-grid sum (continuous omega integral over one Nyquist period)
    f = (g * trap * taper).real
    j = np.arange(1, 151)
    return float(-(2 / D) * np.sum((-1.) ** j * f[1:] / j) / D)   # trap already carries one factor D

res = dict(lines={}, published={}, lownoise={})
for gname, w in grids.items():
    L = {}
    for om in (1.25, 4.0, 8.0, 15.0):
        L[f'free_{om:g}'] = Zsum(np.sin(om * t) / om, w) - 1
    for om in (4.0, 8.0):
        for G in (.3, 1.):
            L[f'lorentz_{om:g}_G{G:g}'] = Zsum(np.exp(-G * t / 2) * np.sin(om * t) / om, w) - 1
            L[f'gauss_{om:g}_G{G:g}'] = Zsum(np.exp(-(G * t) ** 2 / 2) * np.sin(om * t) / om, w) - 1
    res['lines'][gname] = L
res['closed_form_nyquist_minus_1'] = {f'free_{om:g}': Zexact_discrete(np.sin(om * t) / om) - 1 for om in (1.25, 4.0, 8.0, 15.0)}
for tag in ('s20_seed7', 's20.44_seed7', 's20_seed11'):
    nr = np.load(P + f'RESPONSE_{tag}.npz')['normal']
    row = {}
    for gname, w in grids.items():
        z = np.array([Zsum(nr[:, k], w) for k in range(nr.shape[1])]) - 1
        row[gname] = dict(per_shell=z.tolist(), max_abs=float(np.max(abs(z))))
    row['closed_form_nyquist'] = [Zexact_discrete(nr[:, k]) - 1 for k in range(nr.shape[1])]
    row['implied_lorentzian_width_from_codex_grid'] = (np.pi / 2 * 20 * np.abs(np.array(row['codex']['per_shell']))).tolist()
    res['published'][tag] = row
m = np.load(MEAS); TN = m['tr_normal']; ph = m['phases']
for i, s0 in enumerate(ph):
    row = {}
    for gname, w in grids.items():
        z = np.array([Zsum(TN[i, n], w) for n in range(TN.shape[1])]) - 1
        row[gname] = dict(per_shell=z.tolist(), max_abs=float(np.max(abs(z))))
    res['lownoise'][f'{s0:g}'] = row
res['summary'] = dict(
    published_max_abs_codex_grid=max(r['codex']['max_abs'] for r in res['published'].values()),
    published_max_abs_nyquist_grid=max(r['nyquist']['max_abs'] for r in res['published'].values()),
    lownoise_max_abs_codex_grid=max(r['codex']['max_abs'] for r in res['lownoise'].values()),
    lownoise_max_abs_nyquist_grid=max(r['nyquist']['max_abs'] for r in res['lownoise'].values()),
    lines_max_abs_codex_grid_undamped_and_gauss=max(abs(v) for k, v in res['lines']['codex'].items() if not k.startswith('lorentz')),
    lines_max_abs_codex_grid_lorentz=max(abs(v) for k, v in res['lines']['codex'].items() if k.startswith('lorentz')),
    lines_max_abs_nyquist_grid=max(abs(v) for v in res['lines']['nyquist'].values()))
print(json.dumps(res['lines'], indent=1)); print(json.dumps(res['summary'], indent=1))
for tag, r in res['published'].items():
    print(tag, 'codex grid max %.3e  nyquist max %.3e  implied widths %s' % (r['codex']['max_abs'], r['nyquist']['max_abs'],
          np.round(r['implied_lorentzian_width_from_codex_grid'], 3).tolist()))
for s0, r in res['lownoise'].items():
    print('lownoise', s0, 'codex grid max %.3e  nyquist max %.3e' % (r['codex']['max_abs'], r['nyquist']['max_abs']))
json.dump(res, open(OUT, 'w'), indent=1)
