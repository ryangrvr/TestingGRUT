"""PROTOTYPE (too slow for tau 800; used only for a tau<=20 cross-check of skeptic_tangent.c, see SKEPTIC.md).
st_rerun_skeptic: OWN tangent solver for the preparation-phase response (independent equations and method).

What is ours:
  * equations: physical (cosmic-time) field perturbations  dphi'' + 3H dphi' + (k^2/a^2 + V_rr + xi R P_X) dphi = 0
    for the six canonical reals phi = sqrt2 (Re X, Im X, Re Y, Im Y, Re S, Im S); no conformal canonical variables,
    no Ward projection, R = 6 (Hdot + 2 H^2) taken from the background's Hdot (not Codex's trace formula);
    gauge term dropped (e = 0): its relative size c|g|^2 = e^2 |g|^2 / k^2 is <= 4e-7 at k >= 2 and a <= 56;
  * V_rr: sympy Hessian of OUR potential (rp.potential written out) in the six real fields;
  * seed: central difference of OUR preparation rp.prepare(b, delta=+-1e-5), zero initial velocities;
  * observables: dA = phi_X.dphi_X/|phi_X|^2;  dN = a^3 (dphi1 phi2' + phi1 dphi2' - dphi2 phi1' - phi2 dphi1'),
    normalised by |<N_X>_{64-80}| from OUR C engine on a dense 0.01 grid (bg_own.json);
  * integrator: scipy LSODA (Adams/BDF multistep), not DOP853; output every 0.005 tau; crossing by linear
    interpolation between output points.
What is NOT ours: the 19-variable background right-hand side is Codex's reheating.rhs (Python), co-integrated with the
tangent.  That background was checked separately against our own C engine (bg_own.json: a, H agree to 1e-10).
Usage: python3 -B tangent_own.py <b73_reality_framework_dir> <scratch_dir> <tag> <tau_end> [method] [rtol] [kset]
"""
import sys, json, time, math
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

RF, scratch, tag, tau_end = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], float(sys.argv[4])
method = sys.argv[5] if len(sys.argv) > 5 else 'LSODA'
rtol = float(sys.argv[6]) if len(sys.argv) > 6 else 1e-10
kset = sys.argv[7] if len(sys.argv) > 7 else 'main'
HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3] / 'rp_reheating' / 'code'))
import rp  # our own preparation and potential
sys.path.insert(0, str(RF / 'grut_reheating_preparation_20261026' / 'code'))
import reheating as CX  # Codex background rhs ONLY

beta, decay = (3., .5) if tag == 'quiet3' else (5., 1.)
P = rp.P
KS18 = [.02, .05, .1, .2, .4, .8, 1.2, 1.6, 2., 2.13, 2.4, 2.8, 3.2, 4., 5.6, 5.8, 6., 8.]
if kset == 'main':
    ks = np.array(KS18 + [.005, .01, 2.04, 2.16, 2.18, 2.2])
else:  # fine band scan
    ks = np.round(np.arange(1.80, 2.6001, .01), 4)
nk = len(ks)

# ---- own Hessian of the potential in six canonical real fields (sympy) ----
f = sp.symbols('f0:6', real=True)
x = (f[0] + sp.I * f[1]) / sp.sqrt(2); y = (f[2] + sp.I * f[3]) / sp.sqrt(2); s = (f[4] + sp.I * f[5]) / sp.sqrt(2)
ab2 = lambda w: sp.expand(w * sp.conjugate(w))
V = (P['lx'] * ab2(x) ** 2 + P['M'] ** 2 * ab2(y) + P['ly'] * ab2(y) ** 2 + P['ls'] * (ab2(s) - sp.Rational(P['h'] ** 2) / 2) ** 2
     + 2 * sp.re(sp.expand(P['kap'] * sp.conjugate(y) * x ** 3 + P['mu'] * sp.conjugate(s) * y * y)))
V = sp.expand(V)
hessf = sp.lambdify([f], sp.hessian(V, f), 'numpy')
gradf = sp.lambdify([f], [sp.diff(V, q) for q in f], 'numpy')
# check against rp.potential's Wirtinger force: dV/dphi_re + i dV/dphi_im = sqrt2 * F  (F = dV/dz*)
ztest = np.array([70.3 + 4.1j, -2.0 + .7j, .3 - 28.1j])
phit = np.sqrt(2) * np.array([ztest[0].real, ztest[0].imag, ztest[1].real, ztest[1].imag, ztest[2].real, ztest[2].imag])
_, Ftest = rp.potential(ztest)
g = np.array(gradf(phit)); grad_check = float(np.max(np.abs((g[::2] + 1j * g[1::2]) - np.sqrt(2) * Ftest)) / np.max(np.abs(Ftest)))

# ---- own seed ----
def zprep(delta):
    z, info = rp.prepare(beta, delta=delta)
    assert info['residual'] < 1e-8, info  # root() may report success=False at its 1e-14 tol; residual is what matters
    return z
d = 1e-5
zp, zm = zprep(d), zprep(-d)
dz = (zp - zm) / (2 * d)
seed = np.sqrt(2) * np.array([dz[0].real, dz[0].imag, dz[1].real, dz[1].imag, dz[2].real, dz[2].imag])
v0, _ = CX.prepared(beta)
z0own = zprep(0.)
prep_diff = float(np.max(np.abs(np.array([z0own[0], z0own[1], z0own[2]]) - (v0[2:8:2] + 1j * v0[3:8:2]))))

XI = CX.XI; ETA2 = CX.ETA2
PX = np.diag([1., 1., 0, 0, 0, 0])
k2 = ks ** 2


def rhs(t, yv):
    v = yv[:19]
    dv = CX.rhs(t, v, beta, decay, XI, ETA2)
    a = math.exp(v[15]); H = v[16]; Ric = 6 * (dv[16] + 2 * H * H)
    phi = np.sqrt(2) * v[2:8]
    Mm = np.array(hessf(phi), dtype=float) + XI * Ric * PX
    D = yv[19:19 + 6 * nk].reshape(nk, 6); Dd = yv[19 + 6 * nk:].reshape(nk, 6)
    Ddd = -3 * H * Dd - (k2 / (a * a))[:, None] * D - D @ Mm.T
    return np.concatenate([dv, Dd.ravel(), Ddd.ravel()])


bgown = json.loads((HERE.parents[1] / 'results' / 'bg_own.json').read_text())
norm = abs(bgown[tag]['NX_mean_64_80_dense'])


def observables(t, Y):
    v = Y[:19]; a = np.exp(v[15]); phi = np.sqrt(2) * v[2:4]; phid = np.sqrt(2) * v[8:10]
    D = Y[19:19 + 6 * nk].reshape(nk, 6, -1); Dd = Y[19 + 6 * nk:].reshape(nk, 6, -1)
    r2 = phi[0] ** 2 + phi[1] ** 2
    amp = (phi[0] * D[:, 0] + phi[1] * D[:, 1]) / r2
    dN = a ** 3 * (D[:, 0] * phid[1] + phi[0] * Dd[:, 1] - D[:, 1] * phid[0] - phi[1] * Dd[:, 0])
    NX = a ** 3 * (phi[0] * phid[1] - phi[1] * phid[0])
    return amp.T, (dN / norm).T, a, NX


y0 = np.concatenate([v0, np.tile(seed, nk), np.zeros(6 * nk)])
dt_out = .005
edges = np.unique(np.r_[np.arange(0., tau_end, 20.), tau_end])
T, AMP, CHG, AA, NXs = [], [], [], [], []
t0 = time.time(); nfev = 0; yv = y0
for lo, hi in zip(edges[:-1], edges[1:]):
    te = np.round(np.arange(lo, hi + 1e-9, dt_out), 6)
    sol = solve_ivp(rhs, (lo, hi), yv, method=method, t_eval=te, rtol=rtol, atol=rtol * 1e-3)
    assert sol.success, sol.message
    nfev += sol.nfev; yv = sol.y[:, -1]
    amp, chg, a, nx = observables(sol.t, sol.y)
    keep = slice(None) if lo == 0 else slice(1, None)
    T.append(sol.t[keep]); AMP.append(amp[keep]); CHG.append(chg[keep]); AA.append(a[keep]); NXs.append(nx[keep])
    print(f'{tag} {method} chunk {lo:.0f}-{hi:.0f} nfev {sol.nfev} elapsed {time.time() - t0:.1f}s', flush=True)
T = np.concatenate(T); AMP = np.concatenate(AMP); CHG = np.concatenate(CHG); AA = np.concatenate(AA); NXs = np.concatenate(NXs)
np.savez_compressed(scratch / f'tan_{tag}_{method}_{kset}_{tau_end:g}.npz', t=T, amp=AMP, charge=CHG, a=AA, NX=NXs, ks=ks,
                    seed=seed)


def first_cross(series, level, subset=None):
    idx = np.arange(nk) if subset is None else np.asarray(subset)
    m = np.max(np.abs(series[:, idx]), axis=1); w = np.flatnonzero(m >= level)
    if not len(w): return None
    i = int(w[0]); j = idx[int(np.argmax(np.abs(series[i, idx])))]
    y1, y2 = abs(series[i - 1, j]), abs(series[i, j])
    tc = T[i - 1] + (level - y1) / (y2 - y1) * (T[i] - T[i - 1]) if y2 != y1 else T[i]
    return dict(tau=float(tc), k=float(ks[j]), a=float(np.interp(tc, T, AA)), bracket=[float(T[i - 1]), float(T[i])])


i18 = [int(np.argmin(abs(ks - k))) for k in KS18]
res = dict(tag=tag, method=method, rtol=rtol, kset=kset, ks=ks.tolist(), nfev=int(nfev), seconds=time.time() - t0,
           grad_check_vs_rp_force=grad_check, prep_own_vs_codex_max_abs=prep_diff, seed=seed.tolist(),
           amplitude_per_radian_t0=float(AMP[0, 0]), norm=norm, crossings={})
for grid, sub in [('18k', i18), ('all', None)]:
    for sd in [.03, 1e-6, 1e-9]:
        res['crossings'][f'{grid}_{sd:g}'] = dict(
            amp10=first_cross(sd * AMP, .1, sub), amp1=first_cross(sd * AMP, 1., sub), charge1=first_cross(sd * CHG, 1., sub),
            max_amp=float(sd * np.abs(AMP[:, sub if sub is not None else slice(None)]).max()),
            max_charge=float(sd * np.abs(CHG[:, sub if sub is not None else slice(None)]).max()))
w80 = T <= 80
res['through80_003'] = dict(max_amp_18k=float(.03 * np.abs(AMP[w80][:, i18]).max()),
                            argmax_amp_18k=[float(T[w80][np.unravel_index(np.abs(AMP[w80][:, i18]).argmax(), AMP[w80][:, i18].shape)[0]]),
                                            float(ks[i18][np.unravel_index(np.abs(AMP[w80][:, i18]).argmax(), AMP[w80][:, i18].shape)[1]])],
                            max_amp_all=float(.03 * np.abs(AMP[w80]).max()),
                            max_charge_18k=float(.03 * np.abs(CHG[w80][:, i18]).max()), max_charge_all=float(.03 * np.abs(CHG[w80]).max()))
out = HERE.parents[1] / 'results' / f'tangent_own_{tag}_{method}_{kset}_{tau_end:g}.json'
out.write_text(json.dumps(res, indent=1, default=float) + '\n')
print(json.dumps({k: v for k, v in res.items() if k not in ('seed', 'ks')}, indent=1, default=float))
