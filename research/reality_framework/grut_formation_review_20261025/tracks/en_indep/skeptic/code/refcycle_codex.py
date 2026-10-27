"""Skeptic (en_indep_skeptic): reference projections over the last orbit cycle, computed with Codex's OWN
analyze.spectral() (source text exec'd from the pristine 7d9c93b archive, with symbol/lap/energy copied from en_common.py)
on Codex-ENGINE fields (early3d.cpp restarted from the en_rerun tau=740 snapshot, snapshots every 4 tau).
This is a reproduction-type check, independent of en_indep's code (refproj.py) and of its BM6 fields.
usage: refcycle_codex.py PUB_PACKET_DIR RUN_DIR OUT.json [EN_RERUN_N48_DIR]
"""
import sys, json, struct, glob, ast
import numpy as np
from scipy.fft import rfftn, irfftn
LX, G2, LP = .001, .01, .03
pub, run, outp = sys.argv[1], sys.argv[2], sys.argv[3]
# --- verbatim definitions from en_common.py (7d9c93b) ---
def symbol(N, L):
    a = 2*np.pi*np.fft.fftfreq(N); b = 2*np.pi*np.fft.rfftfreq(N); dx = L/N
    def f(x): return 4*np.sin(x/2)**2*(1+np.sin(x/2)**2/3)/dx**2
    return f(a)[:, None, None]+f(a)[None, :, None]+f(b)[None, None, :]
def lap(q, L): return irfftn(-symbol(q.shape[-1], L)*rfftn(q, axes=(-3, -2, -1)), s=q.shape[-3:], axes=(-3, -2, -1))
def energy(q, p, L):
    r = q[0]**2+q[1]**2
    return float(np.mean(.5*np.sum(p*p-q*lap(q, L), axis=0)+LX/4*r*r+G2/2*q[2]**2*r+LP/4*q[2]**4))
# --- Codex analyze.spectral, exec'd from the archive text ---
src = open(pub + '/code/analyze.py').read()
tree = ast.parse(src); fn = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'spectral'][0]
g = dict(np=np, rfftn=rfftn, symbol=symbol, energy=energy, LX=LX, G2=G2, LP=LP)
exec(compile(ast.Module(body=[fn], type_ignores=[]), 'analyze.spectral', 'exec'), g)
spectral = g['spectral']
AM, C = 8.869218738806, 1.9857962750634526
def load(path):
    with open(path, 'rb') as f:
        N, L, s, bits = struct.unpack('<iddi', f.read(24))
        a = np.fromfile(f, dtype='<f8').reshape(6, N, N, N)
    return N, L, s, a[:3], a[3:]
clock = json.load(open(run + '/clock.json'))
rows = {}
for path in sorted(glob.glob(run + '/state_s*.bin')):
    N, L, s, q, p = load(path)
    S = clock['s'] + s; tau = 20 + AM*S + .5*C*S*S
    sp = spectral(q, p, L)
    r = {x['name']: x['raw_nminus_over_net_current'] for x in sp['reference_projections']}
    m = q.mean(axis=(1, 2, 3))
    E = energy(q, p, L); mp = p.mean(axis=(1, 2, 3)); mr = m[0]**2+m[1]**2
    E0 = .5*np.sum(mp**2)+LX/4*mr**2+G2/2*m[2]**2*mr+LP/4*m[2]**4
    rows[f'{tau:.3f}'] = dict(massless=r['massless'], hartree=r['isotropic_Hartree_reference'],
                              ratio=r['massless']/r['isotropic_Hartree_reference'], hartree_mass2=sp['light_Hartree_mass2'],
                              mean_radius=float(np.hypot(m[0], m[1])), outside_mean=float(1-E0/E))
    print(f'{tau:8.3f}', json.dumps(rows[f'{tau:.3f}']), flush=True)
v = np.array([[x['massless'], x['hartree'], x['ratio']] for k, x in rows.items() if float(k) >= 751.9])
summ = dict(window='tau 752-800 every 4 (13 samples)', massless=[float(v[:, 0].min()), float(v[:, 0].max()), float(v[:, 0].mean())],
            hartree=[float(v[:, 1].min()), float(v[:, 1].max()), float(v[:, 1].mean())], ratio=[float(v[:, 2].min()), float(v[:, 2].max())])
# full cycle 740-788 (one period ~ 48 tau)
w = np.array([[x['massless'], x['hartree'], x['ratio']] for k, x in rows.items() if 739.9 <= float(k) <= 788.1])
summ['cycle_740_788'] = dict(massless=[float(w[:, 0].min()), float(w[:, 0].max()), float(w[:, 0].mean())],
                             hartree=[float(w[:, 1].min()), float(w[:, 1].max()), float(w[:, 1].mean())], ratio=[float(w[:, 2].min()), float(w[:, 2].max())])
out = dict(method=__doc__, rows=rows, summary=summ)
if len(sys.argv) > 4:   # cross-check: published-engine snapshots (en_rerun) at 740/760/780/800 straight from the main run
    cc = {}
    for k in (36, 37, 38, 39):
        N, L, s, q, p = load(f'{sys.argv[4]}/state_s{k:07d}.bin'); sp = spectral(q, p, L); tau = 20 + AM*s + .5*C*s*s
        cc[f'{tau:.3f}'] = [x['raw_nminus_over_net_current'] for x in sp['reference_projections']]
    out['en_rerun_main_run'] = cc; print('main run', cc)
json.dump(out, open(outp, 'w'), indent=1)
print(json.dumps(summ))
