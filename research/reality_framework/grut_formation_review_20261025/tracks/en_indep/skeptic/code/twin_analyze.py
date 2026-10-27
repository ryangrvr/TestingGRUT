"""Skeptic (en_indep_skeptic): growth of twin differences (Codex engine restarts from the same tau=100 state).
For each snapshot pair: ||d q_X||, ||d q_phi||, ||d p_X|| split into Fourier bands (zero mode, 0<|k|<=2.5, 2.5<|k|<=6, |k|>6),
normalized by the base run's fluctuation norm of the same field group; plus diagnostic differences from state.tsv.
usage: twin_analyze.py BASE_DIR TWIN_DIR [TWIN_DIR ...] OUT.json
"""
import sys, json, struct, glob
import numpy as np
AM, C = 8.869218738806, 1.9857962750634526
base, twins, outp = sys.argv[1], sys.argv[2:-1], sys.argv[-1]
def load(path):
    with open(path, 'rb') as f:
        N, L, s, bits = struct.unpack('<iddi', f.read(24))
        a = np.fromfile(f, dtype='<f8').reshape(6, N, N, N)
    return N, L, s, a
def bands(N, L):
    k1 = 2*np.pi*np.fft.fftfreq(N, d=L/N)
    K = np.sqrt(k1[:, None, None]**2 + k1[None, :, None]**2 + k1[None, None, :]**2)
    return dict(zero=K == 0, ir=(K > 0) & (K <= 2.5), mid=(K > 2.5) & (K <= 6), uv=K > 6)
def grp(a, idx): return a[idx]
clock = json.load(open(base + '/clock.json'))
out = dict(base=base, note=__doc__, twins={})
for tw in twins:
    rows = []
    for pb in sorted(glob.glob(base + '/state_s*.bin')):
        pt = tw + '/' + pb.split('/')[-1]
        N, L, s, A = load(pb); _, _, s2, B = load(pt); assert abs(s - s2) < 1e-12
        S = clock['s'] + s; tau = 20 + AM*S + .5*C*S*S
        bd = bands(N, L); row = dict(tau=round(tau, 4), eta=S)
        for name, idx in (('qX', [0, 1]), ('qphi', [2]), ('pX', [3, 4]), ('pphi', [5])):
            a = A[idx]; d = B[idx] - a
            fl = a - a.mean(axis=(1, 2, 3), keepdims=True)
            nrm = np.sqrt(np.sum(fl**2))
            dh = np.fft.fftn(d, axes=(1, 2, 3)); tot = np.sqrt(np.sum(np.abs(dh)**2) / N**3)
            row[name] = dict(total=float(tot/nrm), **{b: float(np.sqrt(np.sum(np.abs(dh[:, m])**2) / N**3)/nrm) for b, m in bd.items()})
        rows.append(row)
    # diagnostics from tsv
    def tsv(d):
        t = np.genfromtxt(d + '/state.tsv', names=True, delimiter='\t'); return t
    ta, tb = tsv(base), tsv(tw); n = min(len(ta), len(tb))
    ta, tb = ta[:n], tb[:n]
    def derived(t):
        mr = t['q0']**2 + t['q1']**2
        E0 = .5*(t['p0']**2+t['p1']**2+t['pphi']**2) + 1e-3/4*mr**2 + 1e-2/2*t['qphi']**2*mr + 3e-2/4*t['qphi']**4
        return dict(outside=1-E0/t['rho'], meanQ=(t['q0']*t['p1']-t['q1']*t['p0'])/t['nQ'],
                    phimix=(t['Kphi']+t['Gphi']+t['Uphi']+t['Umix'])/t['rho'],
                    contrast=np.sqrt(np.maximum(t['charge_RMS']**2-t['nQ']**2, 0))/abs(t['nQ']))
    da, db = derived(ta), derived(tb)
    tau = ta['tau']; seg = [(100, 150), (150, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 700), (700, 800.01)]
    diag = {k: {f'{a}-{b}': float(np.max(np.abs(db[k]-da[k])[(tau >= a) & (tau < b)]/np.abs(da[k][(tau >= a) & (tau < b)]))) for a, b in seg} for k in da}
    out['twins'][tw] = dict(snapshots=rows, diagnostic_max_rel_by_segment=diag)
    print(tw)
    for r in rows:
        print(f"  tau {r['tau']:7.2f}  qX tot {r['qX']['total']:.2e} [0 {r['qX']['zero']:.1e} ir {r['qX']['ir']:.1e} mid {r['qX']['mid']:.1e} uv {r['qX']['uv']:.1e}]"
              f"  qphi {r['qphi']['total']:.2e} [ir {r['qphi']['ir']:.1e} mid {r['qphi']['mid']:.1e} uv {r['qphi']['uv']:.1e}]  pX {r['pX']['total']:.2e}")
    for k, v in diag.items(): print('  diag', k, ' '.join(f'{a}:{b:.1e}' for a, b in v.items()))
json.dump(out, open(outp, 'w'), indent=1)
