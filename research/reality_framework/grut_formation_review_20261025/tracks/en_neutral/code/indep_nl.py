"""en_neutral: INDEPENDENT nonlinear X/phi lattice integrator (numpy), used to check Codex's nonlinear neutral
fraction and to run counterfactual initial states.  Own code: real-space FD4 Laplacian with np.roll, Forest-Ruth/
Yoshida fourth-order composition of velocity-Verlet steps with a FIXED step h, own diagnostics.  Same Hamiltonian as
Codex MATHEMATICS.md (conformal, autonomous in s).  Initial state from en_lib (own generator).

usage: python3 -I indep_nl.py <published_dir> <out_json> <N> <h> <end_tau> [variant] [threads]
variant: 'full' (default) | 'noXnoise' (X Fourier fluctuations removed; zero mode and phi seed kept)
"""
import sys, json, time, os
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import en_lib as E

pub, out = Path(sys.argv[1]), Path(sys.argv[2])
N, h, end_tau = int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
variant = sys.argv[6] if len(sys.argv) > 6 else 'full'
L = 12.0; dx = L / N
M = np.load(pub / 'MATCH_L12.npz')
a0, C = float(M['a_match']), float(M['C_match'])
n, coef = E.realization_coefficients(M, L, seed=1, cutoff=8.0)
if variant == 'noXnoise':
    coef = coef.copy(); coef[:, [0, 1, 3, 4]] = 0.0
q, p = E.fields_from_coefficients(n, coef, M['light_mean_q'], M['light_mean_p'], N)
rho0 = E.energy_density(q, p, L)


def lap(f):
    out = np.zeros_like(f)
    for ax in (1, 2, 3):
        out += (4 / 3) * (np.roll(f, 1, ax) + np.roll(f, -1, ax)) - (1 / 12) * (np.roll(f, 2, ax) + np.roll(f, -2, ax)) - 2.5 * f
    return out / dx ** 2


def force(q):
    r = q[0] ** 2 + q[1] ** 2; s2 = q[2] ** 2
    F = lap(q)
    F[0] -= (E.LX * r + E.G2 * s2) * q[0]
    F[1] -= (E.LX * r + E.G2 * s2) * q[1]
    F[2] -= (E.LP * s2 + E.G2 * r) * q[2]
    return F


def diag(q, p):
    r = q[0] ** 2 + q[1] ** 2
    K = 0.5 * np.mean(p * p, axis=(1, 2, 3))
    G = [E.fd4_gradient_energy_density(q[c], dx) for c in range(3)]
    UX = np.mean(E.LX / 4 * r * r); UP = np.mean(E.LP / 4 * q[2] ** 4); UM = np.mean(E.G2 / 2 * q[2] ** 2 * r)
    rho = K.sum() + sum(G) + UX + UP + UM
    nQ = np.mean(q[0] * p[1] - q[1] * p[0])
    return rho, nQ, (K[2] + G[2] + UP) / rho, (K[2] + G[2] + UP + UM) / rho, float(np.mean(r)), float(q[0].mean() ** 2 + q[1].mean() ** 2)


w1 = 1 / (2 - 2 ** (1 / 3)); w0 = 1 - 2 * w1
s_end = 2 * (end_tau - 20) / (a0 + np.sqrt(a0 * a0 + 2 * C * (end_tau - 20)))
nsteps = int(np.ceil(s_end / h)); h = s_end / nsteps
tau_of = lambda s: 20 + a0 * s + 0.5 * C * s * s
F = force(q)
rows = [(0.0, 20.0) + diag(q, p)]
t0 = time.time(); next_out = 0.5
for i in range(nsteps):
    for w in (w1, w0, w1):  # velocity-Verlet substeps (kick-drift-kick), force reused between substeps
        p += 0.5 * w * h * F
        q += w * h * p
        F = force(q)
        p += 0.5 * w * h * F
    s = (i + 1) * h
    if tau_of(s) >= 20 + next_out - 1e-9 or i == nsteps - 1:
        rows.append((s, tau_of(s)) + diag(q, p)); next_out += 0.5
rows = np.array(rows)
res = dict(label='INDEPENDENT numpy FD4 + 4th-order composition, fixed step', variant=variant, N=N, h=h, steps=nsteps, rho0=rho0,
           energy_drift=float(np.max(abs(rows[:, 2] / rows[0, 2] - 1))), charge_drift=float(np.max(abs(rows[:, 3] / rows[0, 3] - 1))),
           at={str(x): dict(Ephi=float(rows[np.argmin(abs(rows[:, 1] - x)), 4]), phi_plus_int=float(rows[np.argmin(abs(rows[:, 1] - x)), 5]),
                            S_over_X2=float(1 - rows[np.argmin(abs(rows[:, 1] - x)), 7] / rows[np.argmin(abs(rows[:, 1] - x)), 6]))
               for x in (20, 60, 70, 80, 100, 150, 200, 300, 400, 500, 600, 700, 800) if x <= end_tau + 1e-9},
           final_phi_plus_int=float(rows[-1, 5]), final_Ephi=float(rows[-1, 4]),
           mean_last_100tau_phi_plus_int=float(rows[rows[:, 1] >= end_tau - 100, 5].mean()),
           max_phi_plus_int=float(rows[:, 5].max()), tau_at_max=float(rows[np.argmax(rows[:, 5]), 1]), seconds=time.time() - t0)
out.write_text(json.dumps(res, indent=1) + '\n')
np.save(out.with_suffix('.npy'), rows)
print(json.dumps(res, indent=1))
