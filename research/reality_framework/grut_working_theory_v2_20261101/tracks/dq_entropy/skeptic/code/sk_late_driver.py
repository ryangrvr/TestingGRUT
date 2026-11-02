#!/usr/bin/env python3
"""dq_entropy SKEPTIC (Claude, own code): drive sk_late.c from the rpdq state at H_i t = 1000 (sk_state1000.py), compare its
conformal windows with dq_entropy's rpdq windows (results/freeze.json, base_1e5) and re-extrapolate N_X, N_Sigma, E_X.
Usage: sk_late_driver.py STATE.json ETA_END C_H OUT.json [TAG]"""
import json, math, os, subprocess, sys
from pathlib import Path
import numpy as np
SP = Path(os.environ.get('SK_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/dq_entropy_skeptic'))
HERE = Path(__file__).resolve().parent
BIN = SP / 'build' / 'sk_late'
subprocess.run(['gcc', '-O3', '-std=c11', '-o', str(BIN), str(HERE / 'sk_late.c'), '-lm'], check=True)
st = json.load(open(sys.argv[1])); eta_end = float(sys.argv[2]); ch = float(sys.argv[3])
Y = st['Y']; p = st['params']; last = st['last']
a0 = math.exp(Y[0]); K = a0 * a0 * last['H']
line = [a0, st['eta_end'], K, p['lx'], p['ly'], p['ls'], p['kap'], p['M'], p['h'], p['mu']] + Y[3:9] + Y[9:15]
out = subprocess.run([str(BIN), '%.17g' % eta_end, '%.17g' % ch, '20', '0.25'], input=' '.join('%.17g' % v for v in line),
                     capture_output=True, text=True, check=True).stdout
W, Fb, info, end = [], [], None, None
for l in out.split('\n'):
    if l.startswith('W '):
        W.append([float(v) for v in l.split()[1:]])
    elif l.startswith('F '):
        Fb.append([float(v) for v in l.split()[1:]])
    elif l.startswith('I '):
        info = l
    elif l.startswith('END'):
        end = l
W = np.array(W); Fb = np.array(Fb)
cols = ['eta_lo', 'eta_hi', 'NX', 'NY', 'NS', 'Econf', 'Eeng', 'a', 'NX_finestd']
rows = [dict(zip(cols, r)) for r in W.tolist()]
# rpdq windows for comparison
FR = json.load(open(HERE.parents[1] / 'results' / 'freeze.json'))
ref = {(w['eta_lo'], w['eta_hi']): w for w in FR['base_1e5']['windows']}
cmp_ = []
for r in rows:
    k = (r['eta_lo'], r['eta_hi'])
    if k in ref:
        cmp_.append(dict(eta_lo=r['eta_lo'], a=r['a'], NX_sk=r['NX'], NX_rpdq=ref[k]['NX'], dNX=r['NX'] - ref[k]['NX'],
                         NS_sk=r['NS'], NS_rpdq=ref[k]['NS'], Eeng_sk=r['Eeng'], EX_rpdq=ref[k]['EX'], dE=r['Eeng'] - ref[k]['EX'],
                         Econf_sk=r['Econf']))
def fit(name, amin, model):
    a = np.array([r['a'] for r in rows]); y = np.array([r[name] for r in rows]); s = a >= amin; a, y = a[s], y[s]
    X = {'a-2': np.c_[np.ones_like(a), a ** -2.0], 'a-1,a-2': np.c_[np.ones_like(a), a ** -1.0, a ** -2.0], 'const': np.c_[np.ones_like(a)]}[model]
    c, *_ = np.linalg.lstsq(X, y, rcond=None); res = y - X @ c; dof = max(len(y) - X.shape[1], 1)
    cov = np.linalg.inv(X.T @ X) * float(res @ res) / dof
    return dict(model=model, a_min=amin, n=int(len(y)), N_inf=float(c[0]), err=float(math.sqrt(cov[0, 0])))
amax = rows[-1]['a']
fits = {nm: [fit(nm, f * amax, m) for f in (1 / 6, 1 / 3, 1 / 2) for m in ('a-2', 'a-1,a-2', 'const')] for nm in ('NX', 'NS', 'Econf', 'Eeng')}
# fluctuation law: std of fine bins inside each wide window vs a
aa = np.array([r['a'] for r in rows]); sd = np.array([r['NX_finestd'] for r in rows]); s = aa > 60
slope = float(np.polyfit(np.log(aa[s]), np.log(sd[s]), 1)[0])
# rpdq-matched windows: freeze_analysis.windows() keeps only 0.25-bins with eta1 <= hi + 1e-9, and rpdq's bins close on the
# first accepted step past each edge, so each 'Delta eta = 20' window is really [lo, hi - 0.25].  Rebuild that from the fine bins.
matched = []
for r in rows:
    k = (r['eta_lo'], r['eta_hi'])
    if k not in ref:
        continue
    s_ = (Fb[:, 0] >= r['eta_lo'] - 1e-9) & (Fb[:, 1] <= r['eta_hi'] - 0.25 + 1e-9)
    nx = float(np.mean(Fb[s_, 2]))
    matched.append(dict(eta_lo=r['eta_lo'], NX_sk_19p75=nx, NX_rpdq=ref[k]['NX'], d=nx - ref[k]['NX'], EX_sk_19p75=float(np.mean(Fb[s_, 6])),
                        EX_rpdq=ref[k]['EX'], dE=float(np.mean(Fb[s_, 6])) - ref[k]['EX']))
o = dict(info=info, end=end, eta_end=eta_end, c_h=ch, windows=rows, compare_rpdq=cmp_, compare_rpdq_matched_windows=matched, fits=fits,
         std_slope_vs_a=slope)
json.dump(o, open(sys.argv[4], 'w'), indent=1)
print(info); print(end)
print('max |dNX| vs rpdq over common windows: %.3e (n=%d); max |dE_engine|: %.3e' % (max(abs(c['dNX']) for c in cmp_) if cmp_ else -1, len(cmp_),
      max(abs(c['dE']) for c in cmp_) if cmp_ else -1))
for nm in ('NX', 'NS', 'Econf', 'Eeng'):
    print(nm, [(f['model'], round(f['a_min'], 1), round(f['N_inf'], 4), round(f['err'], 4)) for f in fits[nm]])
print('std slope', slope)
if matched:
    print('matched (19.75) windows: max |dNX| = %.3e, max |dE_engine| = %.3e' % (max(abs(m['d']) for m in matched), max(abs(m['dE']) for m in matched)))
for r in rows[::max(1, len(rows) // 12)]:
    print('%7.1f a=%7.1f NX=%.4f NS=%.4f Econf=%.3f Eeng=%.3f std=%.3g' % (r['eta_lo'], r['a'], r['NX'], r['NS'], r['Econf'], r['Eeng'], r['NX_finestd']))
