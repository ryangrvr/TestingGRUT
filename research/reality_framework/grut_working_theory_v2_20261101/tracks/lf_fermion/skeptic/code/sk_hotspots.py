"""Skeptic S4: hot-spot statistics on the recorded REC_P0 lines (read-only). Claude, lf_fermion skeptic (issue #13).
(i) <mu5^4> over one chi period on the 64 lattice samples vs lf_fermion's band-limited interpolation (finterp) and a
linear interpolation to 8192 points; (ii) fermion energy by |k| bin (from lfline's n_kbin, k^3 dk weight [approx]);
(iii) share of the P0 tracer energy in the top 5% / 1% of line sites.
Usage: python3 -B sk_hotspots.py -> ../results/sk_interp_mu4.json, sk_kbin_fraction.json"""
import sys, json, math, glob
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'code'))
import lines_io, lfcommon as LC
from lfline import finterp
OUT = Path(__file__).resolve().parents[1] / 'results'
rd = LC.RUNS / 'REC_P0'; D = lines_io.read_lines(rd); t = D['t']; a = D['a']
vf = 367046.6168416639 / 2000


def lin(c, Nf):
    n = c.shape[-1]; x = np.arange(Nf) * n / Nf; i = np.floor(x).astype(int); f = x - i
    return c[..., i % n] * (1 - f) + c[..., (i + 1) % n] * f


res = {}
for tt in [0.3, 0.4, 0.5, 0.6, 0.7, 1.0]:
    i0 = np.searchsorted(t, tt); idx = list(range(i0, i0 + 76, 3)); r = []
    for l in [0, 1, 2, 4, 5, 8, 9, 10]:
        Xd = vf * D['F'][idx, l, :, 3]; aa = a[idx, None] / 2
        mu, mf, ml = np.abs(Xd) * aa, np.abs(finterp(Xd, 8192)) * aa, np.abs(lin(Xd, 8192)) * aa
        r.append([np.mean(mu ** 4), np.mean(mf ** 4), np.mean(ml ** 4), mu.max(), mf.max()])
    r = np.array(r)
    res[tt] = dict(bl_over_samples=float(r[:, 1].sum() / r[:, 0].sum()), lin_over_samples=float(r[:, 2].sum() / r[:, 0].sum()),
                   max_samples=float(r[:, 3].max()), max_bl=float(r[:, 4].max()))
    print(tt, res[tt])
(OUT / 'sk_interp_mu4.json').write_text(json.dumps(res, indent=1))
SP = str(LC.SP / 'lfq' / 'REC_P0') + '/'
edges = np.array([0, 0.25, 0.5, 0.75, 1.0]); kb = {}
for sp, names in [('chi', ['N2', 'N3', 'N4']), ('sig', ['N1'])]:
    for nm in names:
        for tt in [0.7, 1.0]:
            fr = []
            for f in sorted(glob.glob(SP + 'base_l*_%s.json' % sp)):
                d = json.load(open(f)); row = min(d['meas'], key=lambda q: abs(q['t'] - tt))
                kny = math.pi * row['N'] / d['L']; nb = np.array(row[nm]['n_kbin'][:4])
                e = nb * (edges[1:] ** 4 - edges[:-1] ** 4) * kny ** 4; fr.append(e / e.sum())
            kb['%s_%.1f' % (nm, tt)] = list(np.mean(fr, axis=0)); print(nm, tt, kb['%s_%.1f' % (nm, tt)])
T = lines_io.read_trlines(rd, D['N'], D['nl'])
for nm, j, lines in [('N2', 1, [0, 1, 2, 4, 5, 8, 9, 10]), ('N4', 3, [0, 1, 2, 4, 5, 8, 9, 10]), ('N1', 0, [3, 6, 7, 11])]:
    for tt in [0.3, 0.7, 1.0]:
        s = np.sort(T['TE'][np.argmin(abs(T['t'] - tt)), 0, j][lines].ravel())[::-1]
        kb['P0_top5pct_%s_%.1f' % (nm, tt)] = float(s[:len(s) // 20].sum() / s.sum())
(OUT / 'sk_kbin_fraction.json').write_text(json.dumps(kb, indent=1))
