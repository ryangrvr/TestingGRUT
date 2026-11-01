"""Verification of the skeptic's integrator: (1) force and dV/ds vs central differences of V; (2) convergence order
of the KL8 / Y6 compositions on the coupled system at fixed M (kappa = 0) against scipy DOP853 rtol 1e-13;
(3) extended-Hamiltonian (energy-work) conservation with the drift on; (4) speed."""
import time, numpy as np
from scipy.integrate import solve_ivp
from sk_common import *

out = {}
y4, info = orbit('tau20')
par = make_par(info, 1e-3)
# (1) finite differences
rng = np.random.default_rng(1)
q = np.array([60.0, -40.0, 3.0]); s = 1234.5
G = np.zeros(3); dv = ctypes.c_double()
LIB.sk_grad_py(ctypes.byref(par), arr(q).ctypes.data_as(DP), s, G.ctypes.data_as(DP), ctypes.byref(dv))
V = lambda qq, ss: LIB.sk_V(ctypes.byref(par), arr(qq).ctypes.data_as(DP), ss)
fd = []
for i in range(3):
    h = 1e-5 * max(1, abs(q[i])); e = np.zeros(3); e[i] = h
    fd.append((V(q + e, s) - V(q - e, s)) / (2 * h))
hs = 1e-4
fds = (V(q, s + hs) - V(q, s - hs)) / (2 * hs)
out['fd_check'] = dict(grad=G.tolist(), grad_fd=fd, rel=[abs(a - b) / abs(a) for a, b in zip(G, fd)], dVds=dv.value, dVds_fd=fds,
                       rel_s=abs(dv.value - fds) / abs(fds))
print('FD check grad rel', out['fd_check']['rel'], 'dVds rel', out['fd_check']['rel_s'])

# (2) order test at fixed M = 3 Mperi, chi started near M with a kick; 6 radial periods
M = 3 * info['Mperi']
parf = SkPar(LX, G2, LP, EPST, M, 0.0, par.om_fix, par.fm, 0, 8)
q0 = np.r_[y4[:2], 0.6 * M]; p0 = np.r_[y4[2:], 0.3 * M]
def rhs(t, y):
    G = np.zeros(3); dv = ctypes.c_double()
    LIB.sk_grad_py(ctypes.byref(parf), arr(y[:3]).ctypes.data_as(DP), t, G.ctypes.data_as(DP), ctypes.byref(dv))
    return [y[3], y[4], y[5], -G[0], -G[1], -G[2]]
T = 6 * info['T_r']
ref = solve_ivp(rhs, (0, T), np.r_[q0, p0], method='DOP853', rtol=1e-13, atol=1e-13).y[:, -1]
out['order'] = {}
for meth in (8, 6):
    parf.method = meth
    rows = []
    for nper in (4, 6, 8, 12, 16, 24, 32, 48):
        S = Sim(parf, q0, p0, nper=nper); S.to(T)
        err = np.max(np.abs(S.y[:6] - ref) / np.maximum(1, np.abs(ref)))
        rows.append((nper, err, S.steps, S.hext_err()[0]))
    out['order'][f'method{meth}'] = rows
    print(f'method {meth}:')
    for i, r in enumerate(rows):
        slope = '' if i == 0 else f'  slope {np.log(rows[i-1][1]/r[1])/np.log(r[0]/rows[i-1][0]):.2f}'
        print(f'  nper {r[0]:3d} err {r[1]:.3e} steps {r[2]} Hext {r[3]:.1e}{slope}')

# (3)+(4) drift on: eps = 1e-3 tau20 to 5 Mapo, energy-work identity and speed; also compare KL8 vs Y6 endpoints
for meth, nper in ((8, 32), (8, 48), (6, 64)):
    par = make_par(info, 1e-3, method=meth)
    r2 = y4[0] ** 2 + y4[1] ** 2
    S = Sim(par, np.r_[y4[:2], chi_start(r2, par.M0)], np.r_[y4[2:], 0.0], nper=nper)
    t0 = time.time(); S.to(s_of_M(par, 5 * info['Mapo'])); dt = time.time() - t0
    he, E = S.hext_err()
    print(f'drift run meth {meth} nper {nper}: steps {S.steps} time {dt:.1f}s ({S.steps/dt:.3g} steps/s)  Hext rel {he:.2e}  final y {S.y[:6]}')
    out[f'drift_m{meth}_n{nper}'] = dict(steps=S.steps, time=dt, hext_rel=he, y=S.y.tolist())
save_json(os.path.join(RESULTS, 'sk_verify.json'), out)
