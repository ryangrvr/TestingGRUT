"""en_neutral_skeptic item 2b: ensemble mean / std of the same-seed fixed-orbit control at tau 800 for the N32, N48,
N64 FD4 symbols and continuum k^2 (same continuum seed covariance from published MATCH_L12.npz), and first sampled
unit crossing of the ensemble mean (Delta tau = 0.5). OWN code. usage: python3 -I ensemble_grid.py <rundir> <MATCH> <out>"""
import sys, json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sk_lib as K
rd, mf, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
tr = K.read_trace(rd / 'state.tsv'); rho0 = float(tr['rho'][0])
a0 = float(tr['a'][0]); C = float((tr['a'][-1] - tr['a'][0]) / tr['eta'][-1])
N, L, q, p = K.read_initial(rd / 'initial.bin'); V = L ** 3
xm, vm = q[:2].mean(axis=(1, 2, 3)), p[:2].mean(axis=(1, 2, 3))
Mt = np.load(mf); lookup = {int(v): i for i, v in enumerate(Mt['keys'])}
nmax = 16; a = np.arange(-nmax, nmax + 1)
I, J, Kk = np.meshgrid(a, a, a, indexing='ij'); n = np.stack([I.ravel(), J.ravel(), Kk.ravel()], 1)
kc = 2 * np.pi / L * np.sqrt((n ** 2).sum(1))
first = np.where(n[:, 0] != 0, n[:, 0], np.where(n[:, 1] != 0, n[:, 1], n[:, 2]))
n = n[(kc >= 1.8) & (kc <= 8.0) & (first > 0)]          # one representative per +-n pair
ep = np.array([Mt['evolved_phi'][lookup[int(x)]] for x in (n ** 2).sum(1)])
Cens = ep @ ep.transpose(0, 2, 1) / V                     # per complex coefficient: <c c^dag> (real 2x2)
taus = np.arange(600.0, 800.0 + 1e-9, 0.5)
s_out = 2 * (taus - 20) / (a0 + np.sqrt(a0 * a0 + 2 * C * (taus - 20)))
res = dict(label='OWN ensemble of the fixed-orbit control by grid', pairs=int(len(n)))
for name, k2 in (('N32', None), ('N48', None), ('N64', None), ('continuum', (2 * np.pi / L) ** 2 * (n ** 2).sum(1))):
    if k2 is None:
        Ng = int(name[1:]); k2 = sum(K.fd4_k2_axis(n[:, c], Ng, L) for c in range(3))
    vals, ids = np.unique(np.round(k2, 10), return_inverse=True)
    def rhs(s, y):
        x = y[:2]; v = y[2:4]; r = x @ x; T = y[4:].reshape(len(vals), 2, 2)
        o = np.empty_like(y); o[:2] = v; o[2:4] = -K.LX * r * x
        o[4:] = np.stack([T[:, 1, :], -(vals + K.G2 * r)[:, None] * T[:, 0, :]], axis=1).ravel(); return o
    so = solve_ivp(rhs, (0, s_out[-1]), np.r_[xm, vm, np.tile(np.eye(2).ravel(), len(vals))], t_eval=s_out, method='DOP853', rtol=1e-11, atol=1e-12)
    ens = []
    for j in range(len(s_out)):
        y = so.y[:, j]; r = y[:2] @ y[:2]; T = y[4:].reshape(len(vals), 2, 2)[ids]
        Ce = T @ Cens @ T.transpose(0, 2, 1)
        ens.append(float(np.sum(Ce[:, 1, 1] + (k2 + K.G2 * r) * Ce[:, 0, 0]) / rho0))   # pair: 2 x 1/2
    ens = np.array(ens); f1 = np.flatnonzero(ens >= 1)
    res[name] = dict(ensemble_tau800=float(ens[-1]), first_unit=float(taus[f1[0]]) if len(f1) else None)
    print(name, res[name], flush=True)
out.write_text(json.dumps(res, indent=1) + '\n')
