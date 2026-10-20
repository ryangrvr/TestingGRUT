"""Independent re-analysis of Codex 3D Candidate-1 raw diagnostics (commit 2ee7f82).

Reads only the raw CSV histories + protocol JSONs; does not import the author's analyze.py.
Usage: python -I independent_dynamics.py <raw_dir> <protocol.json> <TRANSVERSE_PROTOCOL.json>
"""
import sys, json
import numpy as np

RAW, PROT, TPROT = sys.argv[1:4]
prot = json.load(open(PROT)); tprot = json.load(open(TPROT))
cases = [c for c in prot['cases'] if c['case'] in ('single', 'pair')] + tprot['cases']
allruns = prot['cases'] + tprot['cases']


def load(name):
    return np.genfromtxt(f'{RAW}/{name}.csv', delimiter=',', names=True, ndmin=1)


def centers(a, j):
    return np.stack([a[f'c{j}_x'], a[f'c{j}_y'], a[f'c{j}_z']], axis=1)


print('=== A. conservation over all 29 science runs ===')
qd = []; ed = []
for c in allruns:
    a = load(c['id'])
    qd.append((np.max(np.abs(a['Q'] / a['Q'][0] - 1)), c['id']))
    ed.append((np.max(np.abs(a['E'] / a['E'][0] - 1)), c['id']))
print('runs', len(allruns), 'max |dQ/Q0| = %.3e (%s)' % max(qd), ' max |dE/E0| = %.3e (%s)' % max(ed))

print('\n=== B. retention recomputed (late window t>=120) ===')
hdr = '%-24s %5s %4s %4s %8s %9s %9s %9s %8s %8s %10s %10s'
print(hdr % ('case', 'h', 'L', 'Lz', 'nlate', 'retain%', 'core0%', 'core1%', 'minsep', 'wallgap', 'slope%/40t', 'Q12late/Qbox'))
summary = {}
for c in cases:
    a = load(c['id']); Q0 = a['Q'][0]; late = a[a['t'] >= 120]
    L = c['L']; Lz = c.get('Lz', L)
    # wall gap: min distance from R12 ball surface to the fixed-vacuum face
    gaps = []
    for j in (0, 1):
        x = centers(late, j)
        g = np.min(np.stack([L - np.abs(x[:, 0]), L - np.abs(x[:, 1]), Lz - np.abs(x[:, 2])], 1), 1) - 12
        gaps.append(g.min())
    if c['case'] == 'single':
        i = int(np.argmax([a[-1]['c0_Q12'], a[-1]['c1_Q12']]))
        r = late[f'c{i}_Q12'] / Q0
        sep = np.nan; c0 = c1 = r.min()
        slope = np.polyfit(late['t'], late[f'c{i}_Q12'] / Q0, 1)[0] * 40
        eq = np.max(late[f'c{i}_E12'] / late[f'c{i}_Q12'])
    else:
        r = (late['c0_Q12'] + late['c1_Q12']) / Q0
        sep = np.min(np.linalg.norm(centers(late, 0) - centers(late, 1), axis=1))
        c0 = np.min(late['c0_Q12']) / Q0; c1 = np.min(late['c1_Q12']) / Q0
        slope = np.polyfit(late['t'], r, 1)[0] * 40
        eq = np.nan
    summary[c['id']] = dict(retain=r.min(), c0=c0, c1=c1, sep=sep, gap=min(gaps), slope40=slope, n=len(late),
                            tmin=late['t'].min(), tmax=late['t'].max(), eq=eq)
    print(hdr.replace('%8s %9s', '%8d %9.5f').replace('%9s %9s %8s %8s %10s %10s', '%9.4f %9.4f %8.2f %8.2f %10.5f %10.5f') %
          (c['id'], c['h'], L, Lz, len(late), 100 * r.min(), 100 * c0, 100 * c1, sep, min(gaps), 100 * slope, r.mean()))
npass = sum(1 for k, v in summary.items() if v['retain'] >= .9 and (np.isnan(v['sep']) or (v['sep'] > 24 and v['c0'] >= .15 and v['c1'] >= .15)))
print('pass count (R12 criterion, excluding single E/Q and component sub-criteria):', npass, '/', len(summary))
print('late windows: t from', {v['tmin'] for v in summary.values()}, 'to', {v['tmax'] for v in summary.values()}, 'samples', {v['n'] for v in summary.values()})
print('single max late E12/Q12:', {k: round(v['eq'], 6) for k, v in summary.items() if k.startswith('single')})


def maxdiff(a, b, keys):
    x = load(a); y = load(b); Q0 = x['Q'][0]
    tt = np.intersect1d(x['t'], y['t']); tt = tt[tt >= 120]
    ix = np.searchsorted(x['t'], tt); iy = np.searchsorted(y['t'], tt)
    return [np.max(np.abs(x[k][ix] - y[k][iy])) / Q0 for k in keys]


def order(d1, d2, h=(0.75, 0.5, 0.375)):
    # solve (h0^p-h1^p)/(h1^p-h2^p) = d1/d2 for p by bisection
    f = lambda p: (h[0] ** p - h[1] ** p) / (h[1] ** p - h[2] ** p) - d1 / d2
    lo, hi = 0.05, 8.0
    if f(lo) * f(hi) > 0:
        return np.nan
    for _ in range(200):
        m = .5 * (lo + hi)
        if f(lo) * f(m) <= 0: hi = m
        else: lo = m
    return .5 * (lo + hi)


def richardson(vals, h=(0.75, 0.5, 0.375)):
    """vals at h coarse/medium/fine -> (p, extrapolated, fine-extrap)"""
    d1 = vals[0] - vals[1]; d2 = vals[1] - vals[2]
    if d1 == 0 or d2 == 0 or np.sign(d1) != np.sign(d2):
        return np.nan, np.nan, np.nan
    p = order(abs(d1), abs(d2), h)
    C = d2 / (h[1] ** p - h[2] ** p)
    f0 = vals[2] - C * h[2] ** p
    return p, f0, vals[2] - f0


print('\n=== C. grid convergence of late core charges (max |dQ12|/Q0 over common late samples) ===')
for fam in ['single_reference', 'single_generic', 'pair_reference', 'pair_generic']:
    cm = max(maxdiff(fam + '_coarse', fam + '_medium', ['c0_Q12', 'c1_Q12']))
    mf = max(maxdiff(fam + '_medium', fam + '_fine', ['c0_Q12', 'c1_Q12']))
    rets = [summary[f'{fam}_{g}']['retain'] for g in ['coarse', 'medium', 'fine']]
    p, f0, err = richardson(rets)
    print('%-18s coarse-medium=%.5f%% medium-fine=%.5f%% ratio=%.2f  diff-order p=%.2f | min retention c/m/f=%s  Richardson p=%.2f extrap=%.5f%% fine-extrap=%.5f%%' %
          (fam, 100 * cm, 100 * mf, cm / mf, order(cm, mf), [round(100 * r, 5) for r in rets], p, 100 * f0 if f0 == f0 else np.nan, 100 * err if err == err else np.nan))
    # per-core late-mean charges and their Richardson extrapolation (pairs)
    if fam.startswith('pair'):
        for j in (0, 1):
            means = [np.mean(load(f'{fam}_{g}')[load(f'{fam}_{g}')['t'] >= 120][f'c{j}_Q12']) for g in ['coarse', 'medium', 'fine']]
            p, f0, err = richardson(means)
            print('    core%d late-mean Q12 c/m/f = %s  p=%.2f extrap=%.4f fine-extrap=%.4f (%.3f%% of Q0)' %
                  (j, [round(m, 4) for m in means], p, f0, err, 100 * err / load(f'{fam}_fine')['Q'][0]))

print('\n=== D. box / dt / rotation controls ===')
for name, a, b in [('single_box', 'single_generic_medium', 'single_generic_box'), ('single_dt', 'single_generic_medium', 'single_generic_dt'),
                   ('pair_box', 'pair_reference_medium', 'pair_reference_box'), ('pair_dt', 'pair_reference_medium', 'pair_reference_dt'),
                   ('rot_medium', 'pair_reference_medium', 'pair_rotated_medium'), ('rot_fine', 'pair_reference_fine', 'pair_rotated_fine'),
                   ('impact_vs_ref_fine', 'pair_reference_fine', 'impact_fine')]:
    d = maxdiff(a, b, ['c0_Q12', 'c1_Q12'])
    print('%-20s max late |dQ12|/Q0 = %.5f%%, %.5f%%' % (name, 100 * d[0], 100 * d[1]))

print('\n=== E. charge outside the two R12 balls vs time (pairs): reflection check, box 24x48 vs 30x60 (medium) ===')
for nm in ['pair_reference_medium', 'pair_reference_box', 'pair_reference_fine']:
    a = load(nm); Q0 = a['Q'][0]
    sep = np.linalg.norm(centers(a, 0) - centers(a, 1), axis=1)
    out = np.where(sep > 24, (a['Q'] - a['c0_Q12'] - a['c1_Q12']) / Q0, np.nan)
    print(nm, ' '.join('%d:%.4f%%' % (t, 100 * o) for t, o in zip(a['t'], out) if t >= 80))

print('\n=== F. direct impact responses (x/delta), reference-subtracted, CM check, growth form ===')
resp = {}
for c in tprot['cases']:
    a = load(c['id']); d = c['extra']['--impact']
    refname = {'impact_coarse': 'pair_reference_coarse', 'impact_medium': 'pair_reference_medium', 'impact_fine': 'pair_reference_fine',
               'impact_half_fine': 'pair_reference_fine', 'impact_box': 'pair_reference_box', 'impact_dt': 'pair_reference_dt'}[c['id']]
    b = load(refname)
    assert np.allclose(a['t'], b['t'])
    xr = np.array([a[f'c{j}_x'] - b[f'c{j}_x'] for j in (0, 1)]) / d
    xraw = np.array([a[f'c{j}_x'] for j in (0, 1)]) / d
    resp[c['id']] = (a['t'], xr)
    late = a['t'] >= 120
    # Q-weighted transverse "centre of charge" of the two cores (should stay ~0 if momentum conserved & rad small)
    qw = (a['c0_Q12'] * a['c0_x'] + a['c1_Q12'] * a['c1_x']) / (a['c0_Q12'] + a['c1_Q12']) / d
    ew = (a['c0_E12'] * a['c0_x'] + a['c1_E12'] * a['c1_x']) / (a['c0_E12'] + a['c1_E12']) / d
    # linear vs exponential fit on t in [80,160]
    m = a['t'] >= 80; t = a['t'][m]
    lin = []; ex = []
    for j in (0, 1):
        y = xr[j][m]
        pl = np.polyfit(t, y, 1); lin.append(np.max(np.abs(np.polyval(pl, t) - y)) / np.max(np.abs(y)))
        # exponential fit y = A exp(g t) + B via scan over g
        best = (1e9, 0)
        for g in np.linspace(-0.05, 0.1, 3001):
            M = np.stack([np.exp(g * (t - 160)), np.ones_like(t)], 1)
            co, *_ = np.linalg.lstsq(M, y, rcond=None)
            r = np.max(np.abs(M @ co - y))
            if r < best[0]: best = (r, g)
        ex.append((best[0] / np.max(np.abs(y)), best[1]))
    print('%-16s endpoint x/d raw=(%.6f, %.6f) ref-subtracted=(%.6f, %.6f)  ref x at end=(%.2e,%.2e)  Qweighted xcm/d end=%.4f Eweighted=%.4f' %
          (c['id'], xraw[0, -1], xraw[1, -1], xr[0, -1], xr[1, -1], b['c0_x'][-1], b['c1_x'][-1], qw[-1], ew[-1]))
    print('                 late slope d(x/d)/dt = (%.5f, %.5f); linear-fit max rel err t>=80: (%.2e, %.2e); best exp fit growth g=(%.4f,%.4f) err=(%.2e,%.2e)' %
          (np.polyfit(a['t'][late], xr[0][late], 1)[0], np.polyfit(a['t'][late], xr[1][late], 1)[0], lin[0], lin[1], ex[0][1], ex[1][1], ex[0][0], ex[1][0]))


def rdiff(a, b):
    x = resp[a][1]; y = resp[b][1]; ta = resp[a][0]; tb = resp[b][0]
    m = ta >= 120
    return np.max(np.abs(x[:, m] - y[:, m]), axis=1) / np.max(np.abs(y[:, m]), axis=1)


print('relative response differences (late, normalised by 2nd arg):')
for nm, a, b in [('coarse_medium', 'impact_coarse', 'impact_medium'), ('medium_fine', 'impact_medium', 'impact_fine'), ('half_seed', 'impact_half_fine', 'impact_fine'),
                 ('box', 'impact_box', 'impact_medium'), ('dt', 'impact_dt', 'impact_medium')]:
    print('   %-14s %s' % (nm, np.round(100 * rdiff(a, b), 5)))
for j in (0, 1):
    vals = [resp[k][1][j, -1] for k in ['impact_coarse', 'impact_medium', 'impact_fine']]
    p, f0, err = richardson(vals)
    print('core%d endpoint x/d c/m/f=%s  Richardson p=%.2f continuum=%.4f fine-continuum=%.4f (%.2f%%)' % (j, np.round(vals, 6), p, f0, err, 100 * err / abs(f0)))
json.dump({k: {kk: float(vv) for kk, vv in v.items()} for k, v in summary.items()}, open(sys.argv[4] if len(sys.argv) > 4 else '/dev/null', 'w'), indent=1)
