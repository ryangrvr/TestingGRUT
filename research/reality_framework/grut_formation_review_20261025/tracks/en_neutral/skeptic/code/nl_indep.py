"""en_neutral_skeptic item 4: INDEPENDENT nonlinear X/phi lattice run (own numpy code, different from both Codex's
early3d.cpp and en_neutral's indep_nl.py): FD4 Laplacian applied as the Fourier symbol (30-32cos+2cos2)/(12dx^2) via
rfftn, PEFRL fourth-order symplectic integrator (Omelyan-Mryglod-Folk 2002; 4 force evaluations per step, not the
Yoshida/Forest-Ruth triple-jump used by Codex and en_neutral), fixed step h.
Initial state: the tau-20 N48 field of the en_neutral reproduction run, re-sampled EXACTLY on an N-grid by moving
its Fourier coefficients (all band modes have |n_i| <= 15 so they fit N >= 32). Variant 'noXnoise' zeroes all
nonzero X Fourier modes (zero mode and phi seed kept), as in en_neutral's counterfactual.
Records every 0.5 tau: rho, charge, E_phi, E_phi+E_mix, <|chi_X|^2>, |chibar|^2.
usage: python3 -I nl_indep.py <rundir> <out_json> <N> <h> <variant> [end_tau]
"""
import sys, json, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sk_lib as K

rd, out = Path(sys.argv[1]), Path(sys.argv[2])
N, h0, variant = int(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
end_tau = float(sys.argv[6]) if len(sys.argv) > 6 else 800.0
tr = K.read_trace(rd / 'state.tsv')
a0 = float(tr['a'][0]); C = float((tr['a'][-1] - tr['a'][0]) / tr['eta'][-1])
N0, L, q0, p0 = K.read_initial(rd / 'initial.bin')

def resample(f):
    F = np.fft.fftn(f) / N0 ** 3
    m0 = np.rint(np.fft.fftfreq(N0) * N0).astype(int)
    nz = np.argwhere(abs(F) > 1e-13)
    G = np.zeros((N, N, N), complex)
    for i, j, k in nz:
        a, b, c = m0[i], m0[j], m0[k]
        assert max(abs(a), abs(b), abs(c)) < N // 2
        G[a % N, b % N, c % N] = F[i, j, k]
    g = np.fft.ifftn(G) * N ** 3
    assert np.max(abs(g.imag)) < 1e-9 * max(1, np.max(abs(g.real)))
    return g.real

q = np.array([resample(q0[c]) for c in range(3)]); p = np.array([resample(p0[c]) for c in range(3)])
if variant == 'noXnoise':
    for c in (0, 1):
        q[c] = q[c].mean(); p[c] = p[c].mean()
dx = L / N
mm = np.fft.fftfreq(N) * N; mr = np.fft.rfftfreq(N) * N
sym = lambda m: (30 - 32 * np.cos(2 * np.pi * m / N) + 2 * np.cos(4 * np.pi * m / N)) / (12 * dx * dx)
K2r = sym(mm)[:, None, None] + sym(mm)[None, :, None] + sym(mr)[None, None, :]
wr = np.full(N // 2 + 1, 2.0); wr[0] = 1; wr[-1] = 1  # rfft Parseval weights (N even)

def lap(f):
    return np.fft.irfftn(-K2r * np.fft.rfftn(f, axes=(-3, -2, -1)), s=(N, N, N), axes=(-3, -2, -1))

def acc(q):
    r = q[0] ** 2 + q[1] ** 2; s2 = q[2] ** 2
    A = lap(q)
    A[0] -= (K.LX * r + K.G2 * s2) * q[0]
    A[1] -= (K.LX * r + K.G2 * s2) * q[1]
    A[2] -= (K.LP * s2 + K.G2 * r) * q[2]
    return A

def diag(q, p):
    r = q[0] ** 2 + q[1] ** 2
    Kin = 0.5 * np.mean(p * p, axis=(1, 2, 3))
    qh = np.fft.rfftn(q, axes=(-3, -2, -1)) / N ** 3
    Gr = 0.5 * np.sum(wr[None, None, None, :] * K2r[None] * abs(qh) ** 2, axis=(1, 2, 3))
    UX = np.mean(K.LX / 4 * r * r); UP = np.mean(K.LP / 4 * q[2] ** 4); UM = np.mean(K.G2 / 2 * q[2] ** 2 * r)
    rho = Kin.sum() + Gr.sum() + UX + UP + UM
    nQ = np.mean(q[0] * p[1] - q[1] * p[0])
    return [rho, nQ, (Kin[2] + Gr[2] + UP) / rho, (Kin[2] + Gr[2] + UP + UM) / rho, float(np.mean(r)), float(q[0].mean() ** 2 + q[1].mean() ** 2)]

XI, LAM, CHI = 0.1786178958448091, -0.2123418310626054, -0.06626458266981849
s_end = 2 * (end_tau - 20) / (a0 + np.sqrt(a0 * a0 + 2 * C * (end_tau - 20)))
nsteps = int(np.ceil(s_end / h0)); h = s_end / nsteps
tau_of = lambda s: 20 + a0 * s + 0.5 * C * s * s
rows = [[0.0, 20.0] + diag(q, p)]
t0 = time.time(); next_out = 20.5
for i in range(nsteps):
    q += XI * h * p
    p += (1 - 2 * LAM) * 0.5 * h * acc(q)
    q += CHI * h * p
    p += LAM * h * acc(q)
    q += (1 - 2 * (CHI + XI)) * h * p
    p += LAM * h * acc(q)
    q += CHI * h * p
    p += (1 - 2 * LAM) * 0.5 * h * acc(q)
    q += XI * h * p
    s = (i + 1) * h
    if tau_of(s) >= next_out - 1e-9 or i == nsteps - 1:
        rows.append([s, tau_of(s)] + diag(q, p)); next_out += 0.5
rows = np.array(rows)
tau = rows[:, 1]; S = rows[:, 6] - rows[:, 7]
# period-averaged S/<|chi|^2>: pericentres = local minima of |chibar|^2 on the 0.5-tau record
rb = rows[:, 7]
per = np.flatnonzero((rb[1:-1] < rb[:-2]) & (rb[1:-1] <= rb[2:])) + 1
pav = []
for a_, b_ in zip(per[:-1], per[1:]):
    sl = slice(a_, b_ + 1)
    pav.append((float(0.5 * (tau[a_] + tau[b_])), float(np.trapezoid(S[sl], rows[sl, 0]) / np.trapezoid(rows[sl, 6], rows[sl, 0]))))
def at(x, col): return float(rows[int(np.argmin(abs(tau - x))), col])
res = dict(label='INDEPENDENT numpy (spectral FD4 symbol + PEFRL), own resampled initial state', variant=variant, N=N, h=h, steps=nsteps,
           rho0=float(rows[0, 2]), energy_drift=float(np.max(abs(rows[:, 2] / rows[0, 2] - 1))),
           charge_drift=float(np.max(abs(rows[:, 3] / rows[0, 3] - 1))),
           phi_plus_int_at={str(x): at(x, 5) for x in (20, 60, 100, 200, 300, 400, 500, 600, 700, 800) if x <= end_tau},
           Ephi_at={str(x): at(x, 4) for x in (20, 100, 200, 400, 600, 800) if x <= end_tau},
           S_over_X2_instantaneous_at={str(x): float(1 - at(x, 7) / at(x, 6)) for x in (100, 200, 300, 400, 500, 600, 700, 800) if x <= end_tau},
           S_over_X2_period_avg=[dict(tau_mid=a_, value=b_) for a_, b_ in pav],
           final_phi_plus_int=float(rows[-1, 5]), final_Ephi=float(rows[-1, 4]),
           max_phi_plus_int=float(rows[:, 5].max()), tau_at_max=float(tau[int(np.argmax(rows[:, 5]))]),
           mean_phi_plus_int_last_100tau=float(rows[tau >= end_tau - 100, 5].mean()),
           seconds=time.time() - t0)
out.write_text(json.dumps(res, indent=1) + '\n')
np.save(out.with_suffix('.npy'), rows)
print(json.dumps({k: v for k, v in res.items() if k != 'S_over_X2_period_avg'}, indent=1))
print('period-avg S/X2:', [(round(a_, 1), round(b_, 3)) for a_, b_ in pav][::2])
