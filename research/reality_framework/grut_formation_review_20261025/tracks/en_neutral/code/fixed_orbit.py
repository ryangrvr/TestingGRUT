"""en_neutral item 2: INDEPENDENT reproduction of Codex's fixed-orbit (same-seed) linear phi control, plus the
Floquet exponent of delta-phi on the matched homogeneous reduced orbit, plus the convention comparison.

Own code throughout (en_lib.py).  Data read: Codex's published MATCH_L12.npz (git 7d9c93b) and, for the comparison
only, Codex's published LINEAR_PHI_CONTROL.json and FIELD_CURVES.npz.

usage: python3 -I fixed_orbit.py <published_dir> <out_json> [h]
"""
import sys, json, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import en_lib as E

pub = Path(sys.argv[1]); out = Path(sys.argv[2]); h = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3
t0 = time.time()
M = np.load(pub / 'MATCH_L12.npz')
L, N = 12.0, 48
a0, C = float(M['a_match']), float(M['C_match'])
mq, mp = M['light_mean_q'], M['light_mean_p']

# ---- realization (same documented key) and ensemble covariance ----
n, coef = E.realization_coefficients(M, L, seed=1, cutoff=8.0)
k2lat = E.lattice_k2(n, N, L)
q, p = E.fields_from_coefficients(n, coef, mq, mp, N)
rho0 = E.energy_density(q, p, L)
Qtot = float(np.mean(q[0] * p[1] - q[1] * p[0]))
r0 = float(mq[0] ** 2 + mq[1] ** 2)
E0 = 0.5 * float(mp[0] ** 2 + mp[1] ** 2) + E.LX / 4 * r0 * r0
Q0 = float(mq[0] * mp[1] - mq[1] * mp[0])
# initial phi + interaction (homogeneous |chi|^2 in the interaction, as in the linear control) of the realization
cq, cp = coef[:, 2], coef[:, 5]
def phi_energy_homog(cq2, cp2, cqp, r):
    return float(np.sum(2 * 0.5 * (cp2 + (k2lat + E.G2 * r) * cq2)))
e_phi0_real = phi_energy_homog(abs(cq) ** 2, abs(cp) ** 2, np.real(np.conj(cq) * cp), r0)
# same with local |chi|^2 and lamP term (= the nonlinear diagnostic E_phi + E_mix at tau 20), real-space
dx = L / N
Ephi_nl0 = 0.5 * np.mean(p[2] ** 2) + E.fd4_gradient_energy_density(q[2], dx) + np.mean(E.LP / 4 * q[2] ** 4)
Emix_nl0 = float(np.mean(E.G2 / 2 * q[2] ** 2 * (q[0] ** 2 + q[1] ** 2)))
# ensemble expectation per pair: <|c|^2> = (C_qq) / V with C = ep ep^T
lookup = {int(v): i for i, v in enumerate(M['keys'])}
ix = np.array([lookup[int(x)] for x in (n * n).sum(1)])
ep = M['evolved_phi'][ix]
Cens = ep @ ep.transpose(0, 2, 1) / L ** 3
e_phi0_ens = float(np.sum(2 * 0.5 * (Cens[:, 1, 1] + (k2lat + E.G2 * r0) * Cens[:, 0, 0])))

# ---- fixed homogeneous orbit + per-mode transfer (unique k2) ----
uk, inv = np.unique(np.round(k2lat, 10), return_inverse=True)
tau_out = np.arange(20.0, 800.0 + 1e-9, 0.5)
s_out = 2 * (tau_out - 20) / (a0 + np.sqrt(a0 * a0 + 2 * C * (tau_out - 20)))
orbit, T = E.rk4_orbit_and_transfer(mq[:2], mp[:2], uk, s_out, h=h)
rows = []
cq2, cp2, cqp = abs(cq) ** 2, abs(cp) ** 2, np.real(np.conj(cq) * cp)
for j in range(len(s_out)):
    t = T[j][inv]; r = orbit[j, :2] @ orbit[j, :2]
    Q = t[:, 0, 0] ** 2 * cq2 + 2 * t[:, 0, 0] * t[:, 0, 1] * cqp + t[:, 0, 1] ** 2 * cp2
    P = t[:, 1, 0] ** 2 * cq2 + 2 * t[:, 1, 0] * t[:, 1, 1] * cqp + t[:, 1, 1] ** 2 * cp2
    e_real = float(np.sum(P + (k2lat + E.G2 * r) * Q))  # 2 * 1/2 per pair
    Ce = t @ Cens @ t.transpose(0, 2, 1)
    e_ens = float(np.sum(Ce[:, 1, 1] + (k2lat + E.G2 * r) * Ce[:, 0, 0]))
    # occupation-type (zero-point subtracted) energy of the ensemble, using instantaneous omega
    w2 = k2lat + E.G2 * r; w = np.sqrt(w2)
    nocc = (Ce[:, 1, 1] + w2 * Ce[:, 0, 0]) / (2 * w) * L ** 3 - 0.5  # per mode (one complex coeff)
    e_sub = float(np.sum(2 * w * nocc) / L ** 3)
    rows.append((tau_out[j], s_out[j], e_real / rho0, e_ens / rho0, e_sub / rho0, e_real / E0, r))
rows = np.array(rows)
det_err = float(np.max(abs(np.linalg.det(T[-1]) - 1)))
Eorb = 0.5 * (orbit[:, 2:] ** 2).sum(1) + E.LX / 4 * ((orbit[:, :2] ** 2).sum(1)) ** 2
orbit_E_drift = float(np.max(abs(Eorb / E0 - 1)))

def first(col, lev):
    w = np.flatnonzero(rows[:, col] >= lev); return float(rows[w[0], 0]) if len(w) else None

# per-shell dominant growth at tau 800: which k2 carries the energy
t = T[-1][inv]; r = rows[-1, 6]
Qf = t[:, 0, 0] ** 2 * cq2 + 2 * t[:, 0, 0] * t[:, 0, 1] * cqp + t[:, 0, 1] ** 2 * cp2
Pf = t[:, 1, 0] ** 2 * cq2 + 2 * t[:, 1, 0] * t[:, 1, 1] * cqp + t[:, 1, 1] ** 2 * cp2
ef = Pf + (k2lat + E.G2 * r) * Qf
order = np.argsort(ef)[::-1]
kcont = 2 * np.pi / L * np.sqrt((n * n).sum(1))
top = [dict(n=n[i].tolist(), k_continuum=float(kcont[i]), k_lattice=float(np.sqrt(k2lat[i])), energy_share=float(ef[i] / ef.sum())) for i in order[:8]]
klat_w = float(np.sum(np.sqrt(k2lat) * ef) / ef.sum())

# ---- Floquet exponent on the matched reduced orbit ----
Tr, rmin, rmax = E.radial_period(E0, Q0)
# integrate orbit from pericentre over one radial period on a fine grid (own RK4)
ns = 20000
sg = np.linspace(0, Tr, ns + 1)
x0 = np.array([rmin, 0.0]); v0 = np.array([0.0, Q0 / rmin])
orb1, _ = E.rk4_orbit_and_transfer(x0, v0, np.array([0.0]), sg, h=Tr / ns)
rr = (orb1[:, :2] ** 2).sum(1)
period_check = float(abs(np.sqrt(rr[-1]) / rmin - 1))
kgrid = np.linspace(0.5, 12.0, 1151)
mu, _ = E.floquet_on_rseries(rr, sg, kgrid ** 2)
ipk = int(np.argmax(mu))
# band edges (mu > 0.5 max)
band = kgrid[mu > 0.5 * mu[ipk]]
# convergence: half grid
mu2, _ = E.floquet_on_rseries(rr[::2], sg[::2], kgrid ** 2)
# late fit of the linear control (energy ~ exp(2 mu s)) from tau 300 to 800
late = rows[:, 0] >= 300
mu_fit = float(np.polyfit(rows[late, 1], np.log(rows[late, 3]), 1)[0] / 2)

# ---- Codex comparison ----
cod = json.load(open(pub / 'LINEAR_PHI_CONTROL.json'))
ct = np.array([r['tau'] for r in cod['rows']]); cv = np.array([r['phi_plus_interaction_raw_fraction'] for r in cod['rows']])
cv_on = np.interp(rows[:, 0], ct, cv)
rel = np.abs(rows[:, 2] / cv_on - 1)
fc = np.load(pub / 'FIELD_CURVES.npz')
tN, fN = fc['N48_tau'], fc['N48_phi_plus_interaction_energy_fraction']

res = dict(
    label='INDEPENDENT (own generator, real-space FD4 energy, RK4 orbit/transfer, own Floquet); seed key reproduced from the documented recipe',
    h_rk4=h, N=N, L=L, pairs=int(len(n)), unique_k2=int(len(uk)),
    rho0_own=rho0, rho0_codex=105467.43730088255, rho0_rel_diff=rho0 / 105467.43730088255 - 1,
    nQ_own=Qtot, nQ_codex=8113.477086664833,
    zero_mode=dict(E0=E0, Q0=Q0, E0_over_rho0=E0 / rho0, E0_over_Ecirc=E0 / E.circular_energy(Q0), rmin=rmin, rmax=rmax,
                   rmax_over_rmin=rmax / rmin, radial_period_s=Tr, radial_period_check=period_check),
    initial_tau20=dict(phi_plus_int_homog_realization_frac=e_phi0_real / rho0, phi_plus_int_homog_ensemble_frac=e_phi0_ens / rho0,
                       Ephi_nl_frac=Ephi_nl0 / rho0, Emix_nl_frac=Emix_nl0 / rho0, Ephi_plus_mix_nl_frac=(Ephi_nl0 + Emix_nl0) / rho0,
                       codex_curve_N48_tau20=float(fN[0]), codex_linear_tau20=float(cv[0])),
    linear_control=dict(
        at={str(int(tt)): dict(own_realization=float(rows[np.argmin(abs(rows[:, 0] - tt)), 2]),
                              own_ensemble=float(rows[np.argmin(abs(rows[:, 0] - tt)), 3]),
                              own_ensemble_zero_point_subtracted=float(rows[np.argmin(abs(rows[:, 0] - tt)), 4]),
                              codex=float(cv_on[np.argmin(abs(rows[:, 0] - tt))]),
                              codex_nonlinear_N48=float(fN[np.argmin(abs(tN - tt))]))
            for tt in (20, 40, 60, 70, 80, 100, 134, 200, 300, 400, 500, 600, 650, 700, 750, 800)},
        final_realization=float(rows[-1, 2]), final_ensemble=float(rows[-1, 3]), final_zero_point_subtracted=float(rows[-1, 4]),
        final_realization_over_E0=float(rows[-1, 5]),
        max_realization=float(rows[:, 2].max()),
        first_unit_realization=first(2, 1.0), first_unit_ensemble=first(3, 1.0), first_unit_subtracted=first(4, 1.0),
        first_10pct_realization=first(2, 0.1), first_1pct_realization=first(2, 0.01),
        codex_final=float(cv[-1]), codex_max=float(cod['maximum_fraction']), codex_first_unit=cod['first_sampled_unit_fraction'],
        max_rel_diff_vs_codex=float(rel.max()), median_rel_diff_vs_codex=float(np.median(rel)),
        transfer_det_err=det_err, orbit_energy_drift=orbit_E_drift,
        late_growth_fit_mu_per_s_tau300_800=mu_fit,
        tau800_energy_weighted_k_lattice=klat_w, tau800_top_modes=top),
    floquet_matched_orbit=dict(mu_max=float(mu[ipk]), k_lattice_at_max=float(kgrid[ipk]),
                               band_half_max=[float(band.min()), float(band.max())],
                               mu_max_half_grid=float(mu2.max()),
                               mu_at=dict({f'{kk:.2f}': float(np.interp(kk, kgrid, mu)) for kk in (2, 3, 4, 5, 6, 6.5, 7, 7.5, 8, 9, 10)}),
                               note='k here is sqrt of the FD4 lattice k^2 (equal to continuum k for a continuum Hill equation)'),
    seconds=time.time() - t0)
out.write_text(json.dumps(res, indent=1) + '\n')
np.savez_compressed(out.with_suffix('.npz'), rows=rows, kgrid=kgrid, mu=mu)
print(json.dumps({k: v for k, v in res.items() if k not in ('linear_control',)}, indent=1))
print(json.dumps({k: v for k, v in res['linear_control'].items() if k != 'tau800_top_modes'}, indent=1))
