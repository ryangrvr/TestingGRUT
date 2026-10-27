"""en_neutral_skeptic item 2: independent checks of en_neutral's fixed-orbit side claims.
  (i)   Floquet exponent of delta-phi on the matched reduced orbit: DOP853 monodromy (rtol 1e-12) over one exact
        radial period found by an event, fine k grid (en_neutral: own RK4 on a 20000-point grid).
  (ii)  Realization scatter of the fixed-orbit control: exact mean and variance of the tau-800 energy over the
        Gaussian ensemble (sum over independent +-n pairs of z^dag H z, E = tr H, Var = tr H^2), to test the
        statement "the single realization sits 4% below the ensemble = realization scatter".
  (iii) Grid dependence of the same-seed fixed-orbit control: identical continuum seed coefficients (all band modes
        fit every grid) evolved with the N32, N48, N64 FD4 symbols and with continuum k^2.
  (iv)  Machinery check of the depleted-orbit construction in factorial.py: concatenated free periods with CONSTANT
        (E, Q) must reproduce control A.
OWN code (sk_lib.py). Inputs: initial.bin of the en_neutral reproduction run (tau-20 field), published MATCH_L12.npz
(git 7d9c93b) for the ensemble covariance, state.tsv for rho0 and the clock.
usage: python3 -I controls_extra.py <rundir> <MATCH_L12.npz> <out_json>
"""
import sys, json, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sk_lib as K

rundir, matchf, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
t0 = time.time()
tr = K.read_trace(rundir / 'state.tsv'); rho0 = float(tr['rho'][0])
N, L, q, p = K.read_initial(rundir / 'initial.bin'); V = L ** 3
a0 = float(tr['a'][0]); C = float((tr['a'][-1] - tr['a'][0]) / tr['eta'][-1])
s800 = float(tr['eta'][-1])
res = dict(label='OWN independent checks (sk_lib)', rho0=rho0)

# ---------- (i) Floquet on the matched orbit ----------
xm, vm = q[:2].mean(axis=(1, 2, 3)), p[:2].mean(axis=(1, 2, 3))
E0 = 0.5 * vm @ vm + K.LX / 4 * (xm @ xm) ** 2
Q0 = xm[0] * vm[1] - xm[1] * vm[0]
sol, T, rmin, rmax = K.one_free_period(E0, Q0)
def mono(k):
    def rhs(s, y):
        z = sol.sol(s); r = z[0] ** 2 + z[1] ** 2; m = k * k + K.G2 * r
        return [y[1], -m * y[0], y[3], -m * y[2]]
    o = solve_ivp(rhs, (0, T), [1, 0, 0, 1], method='DOP853', rtol=1e-12, atol=1e-12)
    M = np.array([[o.y[0, -1], o.y[2, -1]], [o.y[1, -1], o.y[3, -1]]])
    tr_ = np.trace(M); lam = abs(tr_) / 2 + np.sqrt(max(tr_ * tr_ / 4 - 1, 0))
    return np.log(max(lam, 1.0)) / T, np.linalg.det(M)
kg = np.arange(6.50, 7.501, 0.005)
mus = np.array([mono(k)[0] for k in kg])
i = int(np.argmax(mus))
kf = np.arange(kg[max(i - 2, 0)], kg[min(i + 2, len(kg) - 1)], 0.0005)
muf = np.array([mono(k)[0] for k in kf]); j = int(np.argmax(muf))
half = kg[mus > 0.5 * muf[j]]
res['floquet_matched_orbit'] = dict(E0=E0, Q0=Q0, rmin=rmin, rmax=rmax, rmax_over_rmin=rmax / rmin, radial_period=T,
                                    mu_max=float(muf[j]), k_at_max=float(kf[j]), half_max_band=[float(half.min()), float(half.max())],
                                    det_err=float(abs(mono(kf[j])[1] - 1)),
                                    method='DOP853 rtol 1e-12, exact period by event; k grid 0.005 then 0.0005')
print('Floquet', res['floquet_matched_orbit'], flush=True)

# ---------- shared: seed coefficients per mode (band modes, from the N48 field) ----------
ph = np.fft.fftn(q[2]) / N ** 3; pp = np.fft.fftn(p[2]) / N ** 3
use = (abs(ph) + abs(pp)) > 1e-12
m = np.fft.fftfreq(N) * N
I, J, Kk = np.meshgrid(m, m, m, indexing='ij')
nvec = np.rint(np.stack([I[use], J[use], Kk[use]], -1)).astype(int)  # rint: fftfreq*N is not exactly integral
cq, cp = ph[use], pp[use]
def k2_grid(Ng):
    return K.fd4_k2_axis(nvec[:, 0], Ng, L) + K.fd4_k2_axis(nvec[:, 1], Ng, L) + K.fd4_k2_axis(nvec[:, 2], Ng, L)
kcont2 = (2 * np.pi / L) ** 2 * (nvec ** 2).sum(1)

def control(k2modes, s_list):
    vals, ids = np.unique(np.round(k2modes, 10), return_inverse=True)
    def rhs(s, y):
        x = y[:2]; v = y[2:4]; r = x @ x; Tm = y[4:].reshape(len(vals), 2, 2)
        o = np.empty_like(y); o[:2] = v; o[2:4] = -K.LX * r * x
        o[4:] = np.stack([Tm[:, 1, :], -(vals + K.G2 * r)[:, None] * Tm[:, 0, :]], axis=1).ravel()
        return o
    y0 = np.r_[xm, vm, np.tile(np.eye(2).ravel(), len(vals))]
    so = solve_ivp(rhs, (0, s_list[-1]), y0, t_eval=s_list, method='DOP853', rtol=1e-11, atol=1e-12)
    assert so.success
    return vals, ids, so

def tau_of(s): return 20 + a0 * s + 0.5 * C * s * s
s_list = np.array([2 * (t - 20) / (a0 + np.sqrt(a0 * a0 + 2 * C * (t - 20))) for t in (200, 400, 600, 700, 800)])
s_list[-1] = s800

# ---------- (iii) grid dependence ----------
grid = {}
for name, k2m in (('N32', k2_grid(32)), ('N48', k2_grid(48)), ('N64', k2_grid(64)), ('continuum', kcont2)):
    vals, ids, so = control(k2m, s_list)
    row = {}
    for jj, s in enumerate(s_list):
        y = so.y[:, jj]; r = y[:2] @ y[:2]; Tm = y[4:].reshape(len(vals), 2, 2)[ids]
        Qn = Tm[:, 0, 0] ** 2 * abs(cq) ** 2 + 2 * Tm[:, 0, 0] * Tm[:, 0, 1] * np.real(np.conj(cq) * cp) + Tm[:, 0, 1] ** 2 * abs(cp) ** 2
        Pn = Tm[:, 1, 0] ** 2 * abs(cq) ** 2 + 2 * Tm[:, 1, 0] * Tm[:, 1, 1] * np.real(np.conj(cq) * cp) + Tm[:, 1, 1] ** 2 * abs(cp) ** 2
        e = 0.5 * (Pn + (k2m + K.G2 * r) * Qn)
        row[f'{tau_of(s):.0f}'] = float(e.sum() / rho0)
    kk = np.sqrt(k2m); kc = np.sqrt(kcont2)
    res_band = dict(resonant_modes_klat_6p74_7p27=int(((kk >= 6.74) & (kk <= 7.27)).sum()),
                    their_continuum_k_range=[float(kc[(kk >= 6.74) & (kk <= 7.27)].min()) if ((kk >= 6.74) & (kk <= 7.27)).any() else None,
                                             float(kc[(kk >= 6.74) & (kk <= 7.27)].max()) if ((kk >= 6.74) & (kk <= 7.27)).any() else None],
                    max_klat_of_seeded=float(kk.max()), klat_of_axis_k7=float(np.sqrt(K.fd4_k2_axis(7 * L / (2 * np.pi), {'N32': 32, 'N48': 48, 'N64': 64}.get(name, 10 ** 6), L))))
    grid[name] = dict(phi_plus_int_over_rho0_N48=row, **res_band)
    print(name, row, res_band, flush=True)
res['grid_dependence_fixed_orbit_control'] = grid

# ---------- (ii) ensemble mean and variance at tau 800 (N48 symbol) ----------
Mt = np.load(matchf)
lookup = {int(v): ii for ii, v in enumerate(Mt['keys'])}
key = (nvec ** 2).sum(1)
ep = np.array([Mt['evolved_phi'][lookup[int(x)]] for x in key])  # (modes, 2, 2): (q,p) <- 2 complex normals
k2m = k2_grid(48)
vals, ids, so = control(k2m, s_list[-1:])
y = so.y[:, -1]; r = y[:2] @ y[:2]; Tm = y[4:].reshape(len(vals), 2, 2)[ids]
A = np.zeros((len(k2m), 2, 2)); A[:, 0, 0] = 0.5 * (k2m + K.G2 * r); A[:, 1, 1] = 0.5
B = Tm @ ep / np.sqrt(V)                       # per mode c = B z,  E[z z^dag] = I
H = np.conj(B).transpose(0, 2, 1) @ A @ B               # energy of mode = z^dag H z (each of +-n counted once: both are in 'use')
# modes n and -n share the same z (conjugate): the pair energy is 2 z^dag H z  -> count each pair once with factor 2
half = (nvec[:, 0] > 0) | ((nvec[:, 0] == 0) & (nvec[:, 1] > 0)) | ((nvec[:, 0] == 0) & (nvec[:, 1] == 0) & (nvec[:, 2] > 0))
Hp = 2 * H[half]
mean = float(np.real(np.trace(Hp, axis1=1, axis2=2)).sum() / rho0)
var = float(np.real((Hp @ Hp).trace(axis1=1, axis2=2)).sum() / rho0 ** 2)
# realised value with the actual coefficients
Qn = Tm[:, 0, 0] ** 2 * abs(cq) ** 2 + 2 * Tm[:, 0, 0] * Tm[:, 0, 1] * np.real(np.conj(cq) * cp) + Tm[:, 0, 1] ** 2 * abs(cp) ** 2
Pn = Tm[:, 1, 0] ** 2 * abs(cq) ** 2 + 2 * Tm[:, 1, 0] * Tm[:, 1, 1] * np.real(np.conj(cq) * cp) + Tm[:, 1, 1] ** 2 * abs(cp) ** 2
real = float((0.5 * (Pn + (k2m + K.G2 * r) * Qn)).sum() / rho0)
# effective number of independent contributors
neff = mean ** 2 / var
res['ensemble_tau800'] = dict(mean=mean, std=float(np.sqrt(var)), rel_std=float(np.sqrt(var) / mean), realization=real,
                              realization_z_score=float((real - mean) / np.sqrt(var)), n_eff_pairs=float(neff),
                              pairs=int(half.sum()))
print('ensemble', res['ensemble_tau800'], flush=True)

# ---------- (iv) concatenation machinery check ----------
from scipy.optimize import brentq
solA = K.free_orbit_dense(xm, vm, s800 + 1e-9)
f_rv = lambda s: (lambda y: y[0] * y[2] + y[1] * y[3])(solA.sol(s))
ss = np.linspace(0, 1.2, 24001); rv = f_rv(ss); jj = np.flatnonzero((rv[:-1] < 0) & (rv[1:] >= 0))[0]
sp0 = brentq(f_rv, ss[jj], ss[jj + 1], xtol=1e-14)
solP, TP, _, _ = K.one_free_period(E0, Q0)
def r_concat(s):
    s = np.atleast_1d(s); o = np.empty_like(s)
    a = s < sp0; b = ~a
    if a.any(): o[a] = (lambda y: y[0] ** 2 + y[1] ** 2)(solA.sol(s[a]))
    if b.any(): o[b] = (lambda y: y[0] ** 2 + y[1] ** 2)(solP.sol(np.mod(s[b] - sp0, TP)))
    return o
s_f = np.linspace(0, s800, 40001); sL, sR = s_f[:-1], s_f[1:]; sM = 0.5 * (sL + sR)
seed = K.phi_seed_by_k2(rundir / 'initial.bin')
rec = K.propagate(seed['k2'], K.G2 * r_concat(sL), K.G2 * r_concat(sM), K.G2 * r_concat(sR), s_f, [len(s_f) - 1])
e = K.energy_from_T(rec[len(s_f) - 1], seed, K.G2 * r_concat(np.array([s800]))[0])
res['concat_constant_EQ_check'] = dict(value_tau800=float(e.sum() / rho0), target_A=2.0917128, period=TP, first_pericentre_s=sp0)
print('concat check', res['concat_constant_EQ_check'], flush=True)
res['seconds'] = time.time() - t0
out.write_text(json.dumps(res, indent=1) + '\n')
