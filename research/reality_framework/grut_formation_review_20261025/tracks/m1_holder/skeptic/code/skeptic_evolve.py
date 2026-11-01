#!/usr/bin/env python3
"""m1_holder_skeptic (Claude, adversarial reviewer; own code). Third, independent implementation of the homogeneous
reheating system used by Codex f8186fa (prescribed background) and by track m1_holder (rpholder.c, Jordan-frame engine):
plain scipy solve_ivp DOP853 on Codex's PRESCRIBED background, written from the equations in DERIVATIONS.md, not from
either code.  Nothing here runs Codex's evolve.c or the track's rpholder.c.

State y = (I, J, rho_r, ln a, X, Y, Sigma (complex)) with
  I' = J, J' = -(3H + Gamma) J - m_I^2 I, rho_r' = Gamma J^2 - 4 H rho_r, H^2 = (J^2/2 + m_I^2 I^2/2 + rho_r)/3,
  R = 2 m_I^2 I^2 - J^2 (= rho_I - 3 p_I in units 8 pi G = 1),  m_I = 3, Gamma = 0.5  [Q f8186fa DERIVATIONS.md]
  z_k'' + 3 H z_k' + dV/dz_k^* = 0  for X, Y, Sigma with
  V = lx|X|^4 + ly|Y|^4 + M^2|Y|^2 + ls(|S|^2 - h^2/2)^2 + 2 Re[kap Y* X^3 + u(R) S* Y^2] + xiX R |X|^2,  u = mu + lock R.
Charges N_k = -2 a^3 Im(z_k^* z_k'), Gauss N_X + 3 N_Y + 6 N_S = 0.
Outputs: late means of N_X over 640-800 (Codex/track window), 720-800, 800-1000; running 80-wide window means of N_X,
N_Y, N_S to identify the partner of the slow N_X variation; R_min; optional Sigma-phase finite differences.
Usage: skeptic_evolve.py [--end 1000] [--rtol 1e-10] [--phase]    (one core, ~1-2 min per run)
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
H_GEV = 134284273491.30815
LX, LY, LS = 1e-3, 1.0, 1.0
KAP = 0.0034586847211042523
M, HV, MU = 4028528204739.244 / H_GEV, 5371370939652.326 / H_GEV, 961402368973.181 / H_GEV
XIX = -5.0 / 6.0
MI, GAM = 3.0, 0.5
LOCK = (10j - MU) / 12.0
# initial states [Q]: original = d4747cb preparation; holder = f8186fa holder_model.json (both checked by the track and by
# skeptic_point.py to 1e-15; using the published digits here keeps this run independent of our own root solves)
Z_ORIG = np.array([76.39978266076902 + 0j, -2.111509854283515 - 0.9393290404896636j, 0 - 28.293142125868524j])
Z_HOLD = np.array([87.19556213287578 + 0j, -5.75424796328559 + 0j, -5.2146407257463994e-15 - 28.387181944361622j])
CODEX = dict(ref=-8120.739695645271, ref_half=-8118.6145409703995, holder=-8267.43703207837, holder_half=-8267.128654142336,
             ratio=1.0180645288398746, Rmin=-1.0632716653512906, phase_fd=[0.5543843710942298, 0.5543843215267685])
TRACK = dict(ref=-8120.740082447965, holder=-8267.437810464611, ratio=1.0180645761996148, ratio_720_800=1.018292946395092,
             ratio_800_1000=1.0196124009638576, holder_800_1000=-8277.723936923412, ref_800_1000=-8118.500647008936, Rmin=-1.0636,
             phase_fd=[0.5543819905788894, 0.5543819412672235])


def make_rhs(lock):
    def rhs(t, y):
        I, J, rad, lna = y[0], y[1], y[2], y[3]
        X = complex(y[4], y[5]); Y = complex(y[6], y[7]); S = complex(y[8], y[9])
        vX = complex(y[10], y[11]); vY = complex(y[12], y[13]); vS = complex(y[14], y[15])
        H = np.sqrt((0.5 * J * J + 0.5 * MI * MI * I * I + rad) / 3.0)
        R = 2 * MI * MI * I * I - J * J
        u = MU + lock * R
        x2 = (X * X.conjugate()).real; y2 = (Y * Y.conjugate()).real; s2 = (S * S.conjugate()).real
        gX = 2 * LX * x2 * X + 3 * KAP * Y * X.conjugate() ** 2 + XIX * R * X
        gY = (M * M + 2 * LY * y2) * Y + KAP * X ** 3 + 2 * u.conjugate() * S * Y.conjugate()
        gS = 2 * LS * (s2 - 0.5 * HV * HV) * S + u * Y * Y
        aX = -3 * H * vX - gX; aY = -3 * H * vY - gY; aS = -3 * H * vS - gS
        return np.array([J, -(3 * H + GAM) * J - MI * MI * I, GAM * J * J - 4 * H * rad, H,
                         vX.real, vX.imag, vY.real, vY.imag, vS.real, vS.imag,
                         aX.real, aX.imag, aY.real, aY.imag, aS.real, aS.imag])
    return rhs


def initial(z, v=None):
    y = np.zeros(16)
    y[0] = np.sqrt(6.0) / 3.0                      # 4.5 I^2 = 3  ->  H = 1, R = 12
    y[4:10] = np.c_[z.real, z.imag].ravel()
    if v is not None:
        y[10:16] = np.c_[v.real, v.imag].ravel()
    return y


def charges(Y):
    a3 = np.exp(3 * Y[3])
    N = []
    for k in range(3):
        z = Y[4 + 2 * k] + 1j * Y[5 + 2 * k]; v = Y[10 + 2 * k] + 1j * Y[11 + 2 * k]
        N.append(-2 * a3 * (np.conj(z) * v).imag)
    return np.array(N)


def window_mean(t, f, t0, t1):
    m = (t >= t0 - 1e-9) & (t <= t1 + 1e-9)
    if m.sum() < 2:
        return float('nan')
    return float(np.trapezoid(f[m], t[m]) / (t[m][-1] - t[m][0]))


def run(tag, lock, z, end=1000.0, rtol=1e-10, atol=1e-12, dt=0.05, v=None):
    t0 = time.time()
    t_eval = np.arange(0.0, end + dt / 2, dt)
    sol = solve_ivp(make_rhs(lock), (0.0, end), initial(z, v), method='DOP853', rtol=rtol, atol=atol, t_eval=t_eval)
    if not sol.success:
        raise RuntimeError(sol.message)
    t = sol.t; Y = sol.y
    N = charges(Y)                                  # 3 x nt
    gauss = N[0] + 3 * N[1] + 6 * N[2]
    scale = np.maximum(np.abs(N[0]) + 3 * np.abs(N[1]) + 6 * np.abs(N[2]), 1.0)
    R = 2 * MI * MI * Y[0] ** 2 - Y[1] ** 2
    res = dict(tag=tag, lock=[lock.real, lock.imag], end=end, rtol=rtol, nfev=int(sol.nfev), seconds=time.time() - t0,
               NX_640_800=window_mean(t, N[0], 640, 800), NX_720_800=window_mean(t, N[0], 720, 800),
               NY_640_800=window_mean(t, N[1], 640, 800), NS_640_800=window_mean(t, N[2], 640, 800),
               max_abs_gauss_rel=float(np.max(np.abs(gauss) / scale)), final_a=float(np.exp(Y[3, -1])),
               R_range_sampled=[float(R.min()), float(R.max())], sample_dt=dt)
    if end >= 1000:
        res.update(NX_800_1000=window_mean(t, N[0], 800, 1000), NY_800_1000=window_mean(t, N[1], 800, 1000), NS_800_1000=window_mean(t, N[2], 800, 1000),
                   NX_900_1000=window_mean(t, N[0], 900, 1000))
        # running 80-wide window means from 400 to end-80, step 20 -> slow modulation of N_X and its partner
        starts = np.arange(400.0, end - 80 + 1e-9, 20.0)
        run80 = [dict(t0=float(s), NX=window_mean(t, N[0], s, s + 80), NY=window_mean(t, N[1], s, s + 80), NS=window_mean(t, N[2], s, s + 80)) for s in starts]
        res['running_80_means'] = run80
        nx = np.array([r['NX'] for r in run80]); ny = np.array([r['NY'] for r in run80]); ns = np.array([r['NS'] for r in run80])
        res['slow_modulation'] = dict(NX_min=float(nx.min()), NX_max=float(nx.max()), NX_peak_to_peak_rel=float((nx.max() - nx.min()) / abs(nx.mean())),
                                      dNY_times3_range=float(3 * (ny.max() - ny.min())), dNS_times6_range=float(6 * (ns.max() - ns.min())),
                                      corr_NX_3NY=float(np.corrcoef(nx, 3 * ny)[0, 1]), corr_NX_6NS=float(np.corrcoef(nx, 6 * ns)[0, 1]))
        # instantaneous (unaveraged) N_X, N_Y, N_S samples late: amplitude of fast oscillation
        m = t >= 900
        res['late_instantaneous'] = dict(NX_std=float(N[0][m].std()), NY_std=float(N[1][m].std()), NS_std=float(N[2][m].std()),
                                         NY_mean=float(N[1][m].mean()), NS_mean=float(N[2][m].mean()),
                                         absX_max=float(np.max(np.hypot(Y[4][m], Y[5][m]))), absY_max=float(np.max(np.hypot(Y[6][m], Y[7][m]))),
                                         S_over_vev_mean=float(np.mean(np.hypot(Y[8][m], Y[9][m])) / (HV / np.sqrt(2))))
    # coarse N_X(t) record for the notes
    keep = [0.5, 1, 2, 3, 5, 10, 20, 50, 100, 200, 400, 640, 800]
    if end >= 1000:
        keep += [900, 1000]
    res['NX_at_tau'] = {str(k): float(np.interp(k, t, N[0])) for k in keep if k <= end}
    return res, t, N


def r_min(lock, z):
    """fine minimum of R(t) on 0 <= t <= 3 from a dense-output integration."""
    sol = solve_ivp(make_rhs(lock), (0.0, 3.0), initial(z), method='DOP853', rtol=1e-11, atol=1e-13, dense_output=True)
    Rf = lambda tt: 2 * MI * MI * sol.sol(tt)[0] ** 2 - sol.sol(tt)[1] ** 2
    tt = np.linspace(0, 3, 30001); Rv = Rf(tt); i = int(np.argmin(Rv))
    m = minimize_scalar(Rf, bounds=(tt[max(i - 2, 0)], tt[min(i + 2, len(tt) - 1)]), method='bounded', options=dict(xatol=1e-10))
    return float(m.fun), float(m.x)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--end', type=float, default=1000.0)
    p.add_argument('--rtol', type=float, default=1e-10)
    p.add_argument('--phase', action='store_true', help='also run the four Sigma-phase perturbed holder runs to 800')
    p.add_argument('--out', default='skeptic_evolve.json')
    a = p.parse_args()
    out = dict(method='scipy solve_ivp DOP853 on the prescribed background of f8186fa DERIVATIONS.md; own rhs; trapezoid window means on dt = 0.05 samples')
    Rm, tm = r_min(LOCK, Z_HOLD)
    out['R_min_holder'] = dict(R_min=Rm, at_tau=tm, codex_sampled=CODEX['Rmin'], track=TRACK['Rmin'])
    ref, t, Nref = run('original', 0j, Z_ORIG, end=a.end, rtol=a.rtol)
    hol, t, Nhol = run('holder', LOCK, Z_HOLD, end=a.end, rtol=a.rtol)
    out['original'] = ref; out['holder'] = hol
    ratio = hol['NX_640_800'] / ref['NX_640_800']
    cmp = dict(NX_original_own=ref['NX_640_800'], NX_original_codex=CODEX['ref'], NX_original_track=TRACK['ref'],
               original_rel_diff_codex=ref['NX_640_800'] / CODEX['ref'] - 1, original_rel_diff_track=ref['NX_640_800'] / TRACK['ref'] - 1,
               NX_holder_own=hol['NX_640_800'], NX_holder_codex=CODEX['holder'], NX_holder_track=TRACK['holder'],
               holder_rel_diff_codex=hol['NX_640_800'] / CODEX['holder'] - 1, holder_rel_diff_track=hol['NX_640_800'] / TRACK['holder'] - 1,
               ratio_640_800_own=ratio, ratio_codex=CODEX['ratio'], ratio_track=TRACK['ratio'],
               ratio_720_800_own=hol['NX_720_800'] / ref['NX_720_800'], ratio_720_800_track=TRACK['ratio_720_800'],
               H_ratio_own=ratio ** (-2 / 3), H_new_GeV_own=H_GEV * ratio ** (-2 / 3))
    if a.end >= 1000:
        cmp.update(ratio_800_1000_own=hol['NX_800_1000'] / ref['NX_800_1000'], ratio_800_1000_track=TRACK['ratio_800_1000'],
                   ratio_900_1000_own=hol['NX_900_1000'] / ref['NX_900_1000'],
                   NX_holder_800_1000_own=hol['NX_800_1000'], NX_holder_800_1000_track=TRACK['holder_800_1000'],
                   NX_original_800_1000_own=ref['NX_800_1000'], NX_original_800_1000_track=TRACK['ref_800_1000'])
        # ratio from running windows: spread
        rr = np.array([h['NX'] for h in hol['running_80_means']]) / np.array([r['NX'] for r in ref['running_80_means']])
        cmp['ratio_running_80_windows'] = dict(t0=[h['t0'] for h in hol['running_80_means']], ratio=rr.tolist(), min=float(rr.min()), max=float(rr.max()),
                                               H_ratio_spread=[float(rr.max() ** (-2 / 3)), float(rr.min() ** (-2 / 3))])
    out['comparison'] = cmp
    if a.phase:
        fd = {}
        for step in (1e-3, 5e-4):
            vals = []
            for sgn in (-1, 1):
                zz = Z_HOLD.copy(); zz[2] *= np.exp(1j * sgn * step)
                r, _, _ = run(f'phase{sgn:+d}_{step:g}', LOCK, zz, end=800.0, rtol=a.rtol)
                vals.append(r['NX_640_800'])
            fd[str(step)] = float((np.log(abs(vals[1])) - np.log(abs(vals[0]))) / (2 * step))
        out['phase_sensitivity'] = dict(own=fd, codex_fd=CODEX['phase_fd'], track_fd=TRACK['phase_fd'])
    json.dump(out, open(RES / a.out, 'w'), indent=1, default=float)
    print(json.dumps({k: out[k] for k in out if k not in ('original', 'holder')}, indent=1, default=float))
    for k in ('original', 'holder'):
        print(k, {kk: vv for kk, vv in out[k].items() if kk not in ('running_80_means',)})


if __name__ == '__main__':
    main()
