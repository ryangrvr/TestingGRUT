"""en_neutral item 3: why does phi stay at ~0.3% in the nonlinear continuation?  OWN analysis code.

Inputs: the dense trace (Delta tau = 0.05) and 2-tau snapshots of a reproduction run of Codex's engine
(run_codex_engine.py; bitwise-equivalent to the published N48 history), Codex's published MATCH_L12.npz.

Diagnostics (all own):
  A. zero-mode orbit: own energy E0 = 1/2|pbar|^2 + lamX/4 |chibar|^4, charge Q0, E0/E_circ(Q0), osculating r_max/r_min;
     fluctuation variance S = <|chi_X|^2> - |chibar|^2; pump modulation of <|chi_X|^2> and of |chibar|^2 per window.
  B. per-window (6 osculating radial periods) pump statistics: modulation depth of <|chi_X|^2> and of |chibar|^2,
     S/<|chi_X|^2>.  (Floquet exponents are computed period by period in floquet_periods.py; a first-pass finite-time
     Lyapunov estimate was dropped because the intra-period amplitude swing biased it to 0.4-0.7 even for stable modes.)
  C. phi spectrum in the nonlinear run vs the same-seed fixed-orbit linear control, shell-binned on the lattice-k
     scale, every 2 tau; the energy in the matched resonance band 6.74 <= k <= 7.27.
usage: python3 -I mechanism.py <published_dir> <rundir> <out_json>
"""
import sys, json, struct, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import en_lib as E

pub, rundir, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
t0 = time.time()
tr = np.genfromtxt(rundir / 'state.tsv', names=True, delimiter='\t')
s, tau = tr['eta'], tr['tau']
rb2 = tr['q0'] ** 2 + tr['q1'] ** 2
E0 = 0.5 * (tr['p0'] ** 2 + tr['p1'] ** 2) + E.LX / 4 * rb2 ** 2
Q0 = tr['q0'] * tr['p1'] - tr['q1'] * tr['p0']
S = tr['X2'] - rb2
rho = tr['rho'][0]
ecc = np.array([E.orbit_turning_points(e, q) for e, q in zip(E0, Q0)])
rmin_osc, rmax_osc = ecc[:, 0], ecc[:, 1]
E0_over_circ = E0 / np.array([E.circular_energy(q) for q in Q0])


NPER = 6
centers = [20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100, 120, 150, 200, 250, 300, 400, 500, 600, 700, 760]
win = []
for tc in centers:
    j = int(np.argmin(abs(tau - tc)))
    Tr, rmn, rmx = E.radial_period(E0[j], Q0[j])
    W = NPER * Tr
    m = (s >= s[j]) & (s <= s[j] + W)
    if s[m][-1] < s[j] + 0.98 * W: break
    sg = s[m]
    # pump modulation depth over the window
    def depth(x): return float((x.max() - x.min()) / (x.max() + x.min()))
    # oscillating part of the pump: rms of (X2 - mean) relative to mean, and of rb2
    def osc(x): return float(np.std(x) / np.mean(x))
    win.append(dict(tau_center=float(tau[j]), tau_end=float(tau[m][-1]), W_s=float(W), radial_period_s=float(Tr),
                    E0_over_rho0=float(E0[j] / rho), E0_over_Ecirc=float(E0_over_circ[j]), rmax_over_rmin=float(rmx / rmn),
                    Q0_over_Q=float(Q0[j] / tr['nQ'][j]),
                    S_over_X2_mean=float(np.mean(S[m] / tr['X2'][m])), g2S_mean=float(E.G2 * np.mean(S[m])),
                    pump_depth_X2=depth(tr['X2'][m]), pump_depth_zero_mode=depth(rb2[m]),
                    pump_rms_X2=osc(tr['X2'][m]), pump_rms_zero_mode=osc(rb2[m])))
    print('win', tc, {k: round(v, 3) for k, v in win[-1].items() if isinstance(v, float)}, flush=True)

# ---- C. phi spectra: nonlinear snapshots vs same-seed fixed-orbit linear control ----
M = np.load(pub / 'MATCH_L12.npz'); L, N = 12.0, 48
a0, C = float(M['a_match']), float(M['C_match'])
snaps = sorted(rundir.glob('state_s*.bin'))
def read_snap(path):
    with open(path, 'rb') as f:
        n_, L_, s_, bits = struct.unpack('<iddi', f.read(24))
        a = np.fromfile(f, dtype='<f4' if bits == 32 else '<f8').reshape(6, n_, n_, n_)
    return s_, a.astype(np.float64)
nf = np.fft.fftfreq(N) * N
NX, NY, NZ = np.meshgrid(nf, nf, nf, indexing='ij')
k2_all = E.lattice_k2(np.stack([NX, NY, NZ], -1), N, L)
klat = np.sqrt(k2_all)
edges = np.r_[0, np.arange(0.25, 12.01, 0.5), 20.0]
bins = np.clip(np.digitize(klat.ravel(), edges) - 1, 0, len(edges) - 2)
band = (klat >= 6.74) & (klat <= 7.27)
# linear control on the same realization: per-mode transfer for unique k2 of band modes
n_r, coef = E.realization_coefficients(M, L, seed=1, cutoff=8.0)
k2r = E.lattice_k2(n_r, N, L); kr = np.sqrt(k2r)
uk, inv = np.unique(np.round(k2r, 10), return_inverse=True)
spec_rows = []
s_snap = []; arrs = []
for path in snaps:
    s_, a = read_snap(path); s_snap.append(s_)
s_snap = np.array(s_snap)
orbit, T = E.rk4_orbit_and_transfer(M['light_mean_q'][:2], M['light_mean_p'][:2], uk, s_snap, h=1e-3)
cq, cp = coef[:, 2], coef[:, 5]
cq2, cp2, cqp = abs(cq) ** 2, abs(cp) ** 2, np.real(np.conj(cq) * cp)
for idx, path in enumerate(snaps):
    s_, a = read_snap(path)
    tau_s = 20 + a0 * s_ + 0.5 * C * s_ * s_
    j = int(np.argmin(abs(s - s_)))
    X2m = tr['X2'][j]
    phk = np.fft.fftn(a[2]) / N ** 3; pik = np.fft.fftn(a[5]) / N ** 3
    mass2 = E.G2 * X2m + 3 * E.LP * tr['phi2'][j]
    e_mode = 0.5 * (abs(pik) ** 2 + (k2_all + mass2) * abs(phk) ** 2)  # Hartree-homogeneous quadratic phi energy per mode
    e_mode[0, 0, 0] = 0.0
    eb_nl = np.bincount(bins, weights=e_mode.ravel(), minlength=len(edges) - 1)
    # linear control per bin (homogeneous orbit mass)
    t = T[idx][inv]; r = orbit[idx, :2] @ orbit[idx, :2]
    Qm = t[:, 0, 0] ** 2 * cq2 + 2 * t[:, 0, 0] * t[:, 0, 1] * cqp + t[:, 0, 1] ** 2 * cp2
    Pm = t[:, 1, 0] ** 2 * cq2 + 2 * t[:, 1, 0] * t[:, 1, 1] * cqp + t[:, 1, 1] ** 2 * cp2
    e_lin = Pm + (k2r + E.G2 * r) * Qm  # pair = 2 x 1/2
    bl = np.clip(np.digitize(kr, edges) - 1, 0, len(edges) - 2)
    eb_lin = np.bincount(bl, weights=e_lin, minlength=len(edges) - 1)
    bandsel = (kr >= 6.74) & (kr <= 7.27)
    # occupation-like measure in the matched band with the homogeneous-Hartree frequency (NL) and orbit-mean frequency
    w2 = k2_all[band] + mass2
    nb_nl = float(np.mean((abs(pik[band]) ** 2 + w2 * abs(phk[band]) ** 2) / (2 * np.sqrt(w2))) * L ** 3)
    spec_rows.append(dict(tau=float(tau_s), nl_bins=(eb_nl / rho).tolist(), lin_bins=(eb_lin / rho).tolist(),
                          nl_band=float(e_mode[band].sum() / rho), lin_band=float(e_lin[bandsel].sum() / rho),
                          nl_band_mean_occupation_hartree=nb_nl,
                          nl_total_quadratic=float(e_mode.sum() / rho), lin_total=float(e_lin.sum() / rho),
                          nl_phi_zero_mode=float(0.5 * (abs(pik[0, 0, 0]) ** 2 + mass2 * abs(phk[0, 0, 0]) ** 2) / rho),
                          nl_frac_k_gt_8=float(e_mode[klat > 8].sum() / e_mode.sum()), nl_frac_k_lt_1p8=float(e_mode[klat < 1.8].sum() / e_mode.sum())))
spec_tau = np.array([r['tau'] for r in spec_rows])
nlb = np.array([r['nl_band'] for r in spec_rows]); lnb = np.array([r['lin_band'] for r in spec_rows])
# orbit-averaged band energy (average over ~ radial period: running mean over 5 snapshots = 10 tau)
def runmean(x, w=5):
    return np.convolve(x, np.ones(w) / w, mode='same')
summary_band = {str(int(tt)): dict(nl=float(runmean(nlb)[np.argmin(abs(spec_tau - tt))]), lin=float(runmean(lnb)[np.argmin(abs(spec_tau - tt))]))
                for tt in (24, 30, 40, 50, 60, 66, 70, 76, 80, 90, 100, 120, 150, 200, 300, 400, 600, 796)}

# late spectrum: where is phi's energy at tau 800 (NL), relative to the linear control
late = spec_rows[-1]
res = dict(label='OWN analysis of a reproduction run of Codex early3d (N48, L12, seed 1); own pump statistics, own spectra and own linear control',
           windows_periods=NPER,
           zero_mode_series={str(int(tt)): dict(E0_over_rho0=float(E0[np.argmin(abs(tau - tt))] / rho),
                                                  E0_over_Ecirc=float(E0_over_circ[np.argmin(abs(tau - tt))]),
                                                  rmax_over_rmin=float(rmax_osc[np.argmin(abs(tau - tt))] / rmin_osc[np.argmin(abs(tau - tt))]),
                                                  S_over_X2=float(S[np.argmin(abs(tau - tt))] / tr['X2'][np.argmin(abs(tau - tt))]))
                             for tt in (20, 40, 50, 55, 60, 65, 70, 80, 100, 150, 200, 300, 400, 600, 800)},
           windows=win, band_energy_runmean10tau=summary_band, spectra_edges=edges.tolist(),
           late_spectrum=dict(tau=late['tau'], nl_bins=late['nl_bins'], lin_bins=late['lin_bins'], nl_frac_k_gt_8=late['nl_frac_k_gt_8'],
                              nl_frac_k_lt_1p8=late['nl_frac_k_lt_1p8']),
           spectra_every_2tau=[dict(tau=r['tau'], nl_band=r['nl_band'], lin_band=r['lin_band'], nl_total_quadratic=r['nl_total_quadratic'],
                                    lin_total=r['lin_total'], nl_band_mean_occupation_hartree=r['nl_band_mean_occupation_hartree'],
                                    nl_frac_k_gt_8=r['nl_frac_k_gt_8'], nl_frac_k_lt_1p8=r['nl_frac_k_lt_1p8']) for r in spec_rows],
           seconds=time.time() - t0)
out.write_text(json.dumps(res, indent=1) + '\n')
np.savez_compressed(out.with_suffix('.npz'), tau=tau[::10], E0=E0[::10], Q0=Q0[::10], S=S[::10], X2=tr['X2'][::10], rb2=rb2[::10],
                    E0_over_circ=E0_over_circ[::10], rmax_over_rmin=(rmax_osc / rmin_osc)[::10], spec_tau=spec_tau,
                    nl_bins=np.array([r['nl_bins'] for r in spec_rows]), lin_bins=np.array([r['lin_bins'] for r in spec_rows]), edges=edges)
print(json.dumps(summary_band, indent=0))
print(json.dumps(res['zero_mode_series'], indent=0))
