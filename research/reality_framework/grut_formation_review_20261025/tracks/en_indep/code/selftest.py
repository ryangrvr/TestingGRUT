"""Self-tests of indep3d.py: operator = FD4 stencil, force = -grad H, integrator order, charge conservation."""
import sys, json, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import indep3d as I

rng = np.random.default_rng(7)
out = {}
# 1. spectral symbol vs explicit FD4 stencil (own np.roll implementation)
N, L = 16, 4.0; lat = I.Lattice(N, L); dx = L / N
f = rng.normal(size=(3, N, N, N))
st = np.zeros_like(f)
for ax in (1, 2, 3):
    st += (-np.roll(f, 2, ax) + 16 * np.roll(f, 1, ax) - 30 * f + 16 * np.roll(f, -1, ax) - np.roll(f, -2, ax)) / (12 * dx * dx)
out['lap_vs_stencil_maxabs'] = float(np.max(abs(-lat.minus_lap(f) - st)) / np.max(abs(st)))
# gradient energy by Parseval vs real-space
G = lat.grad_energy_density(f)
out['grad_parseval_vs_realspace'] = float(np.max(abs(G - 0.5 * np.mean(f * lat.minus_lap(f), axis=(1, 2, 3))) / G))
# 2. force = -dH/dq (finite difference of total energy)
q = rng.normal(size=(3, N, N, N)) * 5; p = rng.normal(size=(3, N, N, N))
def Etot(q, p):
    r = q[0]**2 + q[1]**2
    return N**3 * (0.5*np.mean(np.sum(p*p, 0)) + np.sum(lat.grad_energy_density(q)) + np.mean(I.LX/4*r*r + I.G2/2*q[2]**2*r + I.LP/4*q[2]**4))
F = lat.force(q); errs = []
for c, x, y, z in [(0, 1, 2, 3), (1, 5, 5, 5), (2, 0, 7, 9)]:
    e = 1e-5; qa = q.copy(); qa[c, x, y, z] += e; qb = q.copy(); qb[c, x, y, z] -= e
    fd = -(Etot(qa, p) - Etot(qb, p)) / (2 * e); errs.append(abs(fd - F[c, x, y, z]) / abs(F[c, x, y, z]))
out['force_vs_fd_energy_rel'] = float(max(errs))
# 3. order of BM6, VV, RK4 on a small nonlinear lattice
q0 = rng.normal(size=(3, 8, 8, 8)) * 30 + np.array([30., -20., 0.])[:, None, None, None]
p0 = rng.normal(size=(3, 8, 8, 8)) * 50
lat8 = I.Lattice(8, 6.0)
def evolve(meth, h, T=1.0):
    q, p = q0.copy(), p0.copy(); n = int(round(T / h))
    q, p, _ = I.STEPPERS[meth](lat8, q, p, h, n)
    return q, p
ref = evolve('rk4', 1/4096)
for meth in ('bm6', 'rk4', 'vv'):
    e = []
    for h in (1/64, 1/128, 1/256):
        qq, pp = evolve(meth, h); e.append(float(np.max(abs(qq - ref[0]))))
    out[f'order_{meth}'] = dict(errors=e, observed_orders=[float(np.log2(e[0]/e[1])), float(np.log2(e[1]/e[2]))])
# energy/charge conservation BM6
def EQ(q, p):
    r = q[0]**2 + q[1]**2
    E = 0.5*np.mean(np.sum(p*p, 0)) + np.sum(lat8.grad_energy_density(q)) + np.mean(I.LX/4*r*r + I.G2/2*q[2]**2*r + I.LP/4*q[2]**4)
    return E, np.mean(q[0]*p[1] - q[1]*p[0])
E0, Q0 = EQ(q0, p0)
for h in (1/64, 1/128):
    qq, pp = evolve('bm6', h, 4.0); E1, Q1 = EQ(qq, pp)
    out[f'bm6_h{h:.5f}_T4_energy_rel'] = float(E1/E0 - 1); out[f'bm6_h{h:.5f}_T4_charge_rel'] = float(Q1/Q0 - 1)
# 4. timing of one force on 48^3
lat48 = I.Lattice(48, 12.0); qq = rng.normal(size=(3, 48, 48, 48))
t = time.perf_counter()
for i in range(50): lat48.force(qq)
out['force_ms_48'] = (time.perf_counter() - t) / 50 * 1e3
print(json.dumps(out, indent=1))
json.dump(out, open(sys.argv[1], 'w'), indent=1)
