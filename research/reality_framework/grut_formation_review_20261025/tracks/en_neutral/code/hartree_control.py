"""en_neutral item 3c: same-seed LINEAR phi driven by MEASURED homogeneous masses from the nonlinear N48 history.
OWN code.  If linear phi on the measured homogeneous (Hartree) mass g2<|chi_X|^2>(s) + 3 lamP <phi^2>(s) already
stays near the nonlinear 0.3%, the suppression is a homogeneous detuning/pump effect; any remaining gap is due to
inhomogeneous scattering.  Controls (all start from the same tau-20 realization, own generator):
  L0 free homogeneous orbit (Codex's control), but integrated on the trace grid with linear interpolation of the
     mass (validates the interpolation scheme against the fine RK4 result 2.0917)
  L2 measured zero-mode pump g2 |chibar|^2(s) only
  L3 measured Hartree g2 <|chi_X|^2>(s) + 3 lamP <phi^2>(s)
  L4 measured zero-mode pump + measured fluctuation shift g2 S(s) WITHOUT the oscillating part of S
     (S replaced by its running mean over ~1 radial period)
Energy = sum_modes 1/2 (|pi|^2 + (k_lat^2 + g2 <|chi|^2>_used) |phi|^2) / rho0, raw (no subtraction), as in Codex.
usage: python3 -I hartree_control.py <published_dir> <rundir> <out_json>
"""
import sys, json, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import en_lib as E

pub, rundir, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
t0 = time.time()
tr = np.genfromtxt(rundir / 'state.tsv', names=True, delimiter='\t')
s, tau = tr['eta'], tr['tau']
rb2 = tr['q0'] ** 2 + tr['q1'] ** 2
S = tr['X2'] - rb2
M = np.load(pub / 'MATCH_L12.npz'); L, N = 12.0, 48
n, coef = E.realization_coefficients(M, L, seed=1, cutoff=8.0)
k2 = E.lattice_k2(n, N, L)
uk, inv = np.unique(np.round(k2, 10), return_inverse=True)
cq, cp = coef[:, 2], coef[:, 5]
cq2, cp2, cqp = abs(cq) ** 2, abs(cp) ** 2, np.real(np.conj(cq) * cp)
rho0 = tr['rho'][0]
# free orbit on the trace grid (own RK4, fine step) for L0
orbit, _ = E.rk4_orbit_and_transfer(M['light_mean_q'][:2], M['light_mean_p'][:2], np.array([0.0]), s, h=1e-3)
r_free = (orbit[:, :2] ** 2).sum(1)
# running mean of S over ~ one radial period (0.85 in s) for L4
dsm = np.median(np.diff(s))
def runmean_s(x, width=0.85):
    out_ = np.empty_like(x); c = np.r_[0, np.cumsum(x * np.r_[np.diff(s), dsm])]
    for i in range(len(x)):
        a = np.searchsorted(s, s[i] - width / 2); b = np.searchsorted(s, s[i] + width / 2)
        out_[i] = (c[b] - c[a]) / max(s[min(b, len(s) - 1)] - s[a], 1e-12) if b > a else x[i]
    return out_
S_smooth = runmean_s(S)
pumps = dict(L0_free_orbit=(E.G2 * r_free, E.G2 * r_free),
             L2_measured_zero_mode=(E.G2 * rb2, E.G2 * rb2),
             L3_measured_hartree=(E.G2 * tr['X2'] + 3 * E.LP * tr['phi2'], E.G2 * tr['X2']),
             L4_zero_mode_plus_smoothed_S=(E.G2 * (rb2 + S_smooth) + 3 * E.LP * tr['phi2'], E.G2 * (rb2 + S_smooth)))
out_idx = np.arange(0, len(s), 10)  # every 0.5 tau
nl = (tr['Kphi'] + tr['Gphi'] + tr['Uphi'] + tr['Umix']) / rho0
res = dict(label='OWN same-seed linear phi on measured homogeneous masses (N48 reproduction trace)', rho0=float(rho0), controls={})
for name, (m2, mint) in pumps.items():
    nk = len(uk); T = np.zeros((nk, 2, 2)); T[:, 0, 0] = 1; T[:, 1, 1] = 1
    def f(m, T):
        d = np.empty_like(T); d[:, 0, :] = T[:, 1, :]; d[:, 1, :] = -(uk + m)[:, None] * T[:, 0, :]; return d
    vals = np.zeros(len(out_idx)); band = np.zeros(len(out_idx)); jo = 0
    bandsel = (np.sqrt(k2) >= 6.74) & (np.sqrt(k2) <= 7.27)
    for i in range(len(s)):
        if jo < len(out_idx) and i == out_idx[jo]:
            t = T[inv]
            Q = t[:, 0, 0] ** 2 * cq2 + 2 * t[:, 0, 0] * t[:, 0, 1] * cqp + t[:, 0, 1] ** 2 * cp2
            P = t[:, 1, 0] ** 2 * cq2 + 2 * t[:, 1, 0] * t[:, 1, 1] * cqp + t[:, 1, 1] ** 2 * cp2
            e = P + (k2 + mint[i]) * Q
            vals[jo] = e.sum() / rho0; band[jo] = e[bandsel].sum() / rho0; jo += 1
        if i == len(s) - 1: break
        h = s[i + 1] - s[i]; a0, a1 = m2[i], m2[i + 1]; am = 0.5 * (a0 + a1)
        a = f(a0, T); b = f(am, T + h / 2 * a); c = f(am, T + h / 2 * b); d = f(a1, T + h * c)
        T = T + h / 6 * (a + 2 * b + 2 * c + d)
    tt = tau[out_idx]
    def at(x, arr): return float(arr[int(np.argmin(abs(tt - x)))])
    res['controls'][name] = dict(at={str(x): dict(total=at(x, vals), band=at(x, band)) for x in (20, 60, 70, 80, 100, 150, 200, 300, 400, 600, 800)},
                                 final=float(vals[-1]), max=float(vals.max()), mean_700_800=float(vals[tt >= 700].mean()),
                                 det_err=float(np.max(abs(np.linalg.det(T) - 1))))
    print(name, 'final %.4g max %.4g mean700-800 %.4g det %.2e' % (vals[-1], vals.max(), vals[tt >= 700].mean(), res['controls'][name]['det_err']),
          {x: round(at(x, vals), 5) for x in (60, 100, 200, 400, 800)}, flush=True)
res['nonlinear_phi_plus_int'] = {str(x): float(nl[int(np.argmin(abs(tau - x)))]) for x in (20, 60, 70, 80, 100, 150, 200, 300, 400, 600, 800)}
res['nonlinear_mean_700_800'] = float(nl[tau >= 700].mean())
res['g2S_smoothed_at'] = {str(x): float(E.G2 * S_smooth[int(np.argmin(abs(tau - x)))]) for x in (40, 60, 70, 80, 100, 150, 200, 400, 800)}
res['seconds'] = time.time() - t0
out.write_text(json.dumps(res, indent=1) + '\n')
print(json.dumps(res['nonlinear_phi_plus_int']), res['nonlinear_mean_700_800'])
print(json.dumps(res['g2S_smoothed_at']))
