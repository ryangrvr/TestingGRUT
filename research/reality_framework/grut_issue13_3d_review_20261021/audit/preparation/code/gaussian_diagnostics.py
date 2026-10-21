"""Initial-condition diagnostics + quasi-linear predictions for Codex's dense Gaussian preparations.

Usage: python3 -I gaussian_diagnostics.py <codex_raw_dir> <out_json> [run ...]
Regenerates Codex's noise exactly as gaussian_run.noise() does (PCG64, rfftn, exp(-k^2/2), zero k=0,
unit realization rms) and checks it against the recorded means/covariances before use.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
import numpy as np
from scipy.special import erfinv

sys.path.insert(0, str(Path(__file__).resolve().parent))
import quasilinear as ql  # noqa: E402

TGRID = np.round(np.arange(0, 48.0001, 0.1), 10)
RUNS = {'gaussian_s101': (192, 0.4, 101), 'gaussian_s202': (192, 0.4, 202), 'gaussian_large_s101': (288, 0.4, 101)}


def codex_noise(N, dx, seed, shift_S=0.0):
    rng = np.random.default_rng(seed)
    k = 2 * np.pi * np.fft.fftfreq(N, dx)
    kr = 2 * np.pi * np.fft.rfftfreq(N, dx)
    k2 = k[:, None, None] ** 2 + k[None, :, None] ** 2 + kr[None, None, :] ** 2
    filt = np.exp(-.5 * k2)
    q = np.empty((3, N, N, N))
    for f in range(3):
        F = np.fft.rfftn(rng.standard_normal((N, N, N)))
        F *= filt
        F.flat[0] = 0
        x = np.fft.irfftn(F, s=(N, N, N), axes=(0, 1, 2))
        x -= x.mean()
        x /= x.std()
        q[f] = x
    q[2] += shift_S
    return q


def lattice_k2_rfft(N, dx):
    n = np.fft.fftfreq(N) * N
    nr = np.fft.rfftfreq(N) * N
    s = lambda m: 4 * np.sin(np.pi * m / N) ** 2 / dx ** 2  # noqa: E731
    k2 = s(n)[:, None, None] + s(n)[None, :, None] + s(nr)[None, None, :]
    kphys2 = ((2 * np.pi / (N * dx)) ** 2) * (n[:, None, None] ** 2 + n[None, :, None] ** 2 + nr[None, None, :] ** 2)
    w = np.full(k2.shape, 2.0)
    w[:, :, 0] = 1
    if N % 2 == 0:
        w[:, :, -1] = 1
    return k2, kphys2, w


def analyze(q, N, dx, kappa=2.0, rec=None, label='', field_marks=(0.3, 0.5, 1.0), chunk=1500):
    F = kappa ** (1 / 3)
    k2, kp2, w = lattice_k2_rfft(N, dx)
    Qu, Qw, QS = (np.fft.rfftn(q[i]) for i in range(3))
    norm = 1.0 / N ** 6
    zero_S = float(q[2].mean())  # box mean of S noise (non-zero only for a forced shift)
    QS.flat[0] = 0  # zero mode handled separately
    pS = w * np.abs(QS) ** 2 * norm
    cuS = w * np.real(Qu * np.conj(QS)) * norm
    cwS = w * np.real(Qw * np.conj(QS)) * norm
    # bin in lattice k2: exact unique values below 1.0, width-0.002 bins above (power above k2=30 < 1e-13)
    flat = k2.ravel()
    low = flat < 1.0
    uniq = np.unique(np.round(flat[low], 12))
    edges_hi = np.arange(1.0, 30.0 + 1e-9, 0.002)
    centers = np.concatenate([uniq, 0.5 * (edges_hi[:-1] + edges_hi[1:])])
    idx = np.empty(flat.shape, dtype=np.int64)
    idx[low] = np.searchsorted(uniq, np.round(flat[low], 12))
    hi = ~low
    idx[hi] = len(uniq) + np.clip(np.searchsorted(edges_hi, flat[hi], side='right') - 1, 0, len(edges_hi) - 2)
    keep = flat <= 30.0
    nb = len(centers)
    PS_b = np.bincount(idx[keep], weights=pS.ravel()[keep], minlength=nb)
    CU_b = np.bincount(idx[keep], weights=cuS.ravel()[keep], minlength=nb)
    CW_b = np.bincount(idx[keep], weights=cwS.ravel()[keep], minlength=nb)
    used = np.flatnonzero((PS_b > 0) | (CU_b != 0) | (CW_b != 0))
    cen = centers[used]
    G = np.zeros((len(TGRID), nb))
    J = np.zeros(len(TGRID))
    Ju = np.zeros(len(TGRID))
    var = np.zeros(len(TGRID))
    for s in range(0, len(used), chunk):
        sl = used[s:s + chunk]
        g, _ = ql.s_growth(centers[sl], kappa, TGRID)
        Au, Aw = ql.x_response(centers[sl], kappa, TGRID)
        G[:, sl] = g
        J += -np.sum(ql.EPS_X * F * ql.EPS_S * g * (Au * CU_b[sl] + Aw * CW_b[sl]), axis=1)
        Ju += -np.sum(ql.EPS_X * F * ql.EPS_S * g * Au * CU_b[sl], axis=1)
        var += np.sum(g ** 2 * PS_b[sl], axis=1)
    a = ql.a_of(TGRID)
    sig_f = ql.EPS_S * np.sqrt(var) / a
    g0, _ = ql.s_growth(np.array([0.0]), kappa, TGRID)
    S0_forced = ql.EPS_S * zero_S * g0[:, 0] / a  # homogeneous growth of any supplied zero mode
    # bare cubic -<chi_f^3> from the realization at tau = 0,1,...,30
    gi = np.interp(k2, cen, G[0, used])  # placeholder shape check
    tc = np.arange(0, 30.0001, 1.0)
    cube = []
    for t in tc:
        j = int(round(t / 0.1))
        gi = np.interp(k2, cen, G[j, used])
        f = np.fft.irfftn(QS * gi, s=(N, N, N), axes=(0, 1, 2)) * ql.EPS_S
        cube.append(float(np.mean(f ** 3)))
    J3 = -np.interp(TGRID, tc, cube, right=0.0)
    chi0 = ql.zero_mode(TGRID, J, kappa)
    chi0_c = ql.zero_mode(TGRID, J + J3, kappa)
    S0 = chi0 / a + S0_forced
    S0c = chi0_c / a + S0_forced
    out = {'label': label, 'N': N, 'dx': dx, 'L': N * dx, 'kappa': kappa, 'zero_mode_shift_S': ql.EPS_S * zero_S,
           'C_uS': float(np.sum(cuS)), 'C_wS': float(np.sum(cwS)), 'var_S_noise': float(np.sum(pS)),
           'lowest_shell_k': float(2 * np.pi / (N * dx)), 'a_cross_lowest_shell': float(np.sqrt(kappa ** (2 / 3) + (2 * np.pi / (N * dx)) ** 2))}
    kphys = np.sqrt(kp2)
    for kc in (0.1, 0.25, 0.35, 0.5):
        out[f'S_power_frac_k_lt_{kc}'] = float(np.sum(pS[kphys < kc]) / np.sum(pS))
    out['S_power_frac_lowest_shell'] = float(np.sum(pS[np.isclose(kphys, 2 * np.pi / (N * dx))]) / np.sum(pS))
    comp = []
    if rec is not None:
        for i in range(len(rec)):
            t = rec['t'][i]
            if t > 30.5:
                break
            j = int(round(t / 0.1))
            comp.append({'a': float(rec['a'][i]), 'mean_S_recorded': float(rec['mean_S'][i]), 'S0_pred_cross': float(S0[j]),
                         'S0_pred_cross_plus_bare_cubic': float(S0c[j]), 'sigma_f_pred': float(sig_f[j]),
                         'max_abs_S_recorded': float(max(-rec['min_S'][i], rec['max_S'][i]))})
    out['comparison'] = comp
    tot = np.sqrt(sig_f ** 2 + S0c ** 2)
    marks = {}
    for thr in field_marks:
        j = int(np.argmax(tot >= thr)) if np.any(tot >= thr) else len(TGRID) - 1
        gi = np.interp(k2, cen, G[j, used])
        f = np.fft.irfftn(QS * gi, s=(N, N, N), axes=(0, 1, 2)) * ql.EPS_S / a[j]
        full = f + S0c[j]
        full_nocubic = f + S0[j]
        # ensemble statistical spread of <sign S> for this spectrum (Gaussian field, arcsine law)
        Pth = np.exp(-kp2) * gi ** 2
        Pth.flat[0] = 0
        C = np.fft.irfftn(Pth, s=(N, N, N), axes=(0, 1, 2))
        rho = np.clip(C / C.flat[0], -1, 1)
        sd_m = float(np.sqrt(np.mean(2 / np.pi * np.arcsin(rho))))
        wpow = pS * gi ** 2
        marks[f'rms_{thr}'] = {'a': float(a[j]), 'S0': float(S0c[j]), 'S0_cross_only': float(S0[j]), 'sigma_f': float(sig_f[j]),
                               'S0_over_sigma_f': float(S0c[j] / sig_f[j]),
                               'mean_sign_fluct_only': float(np.mean(np.sign(f))),
                               'mean_sign_with_zero_mode': float(np.mean(np.sign(full))),
                               'pplus_with_zero_mode': float(np.mean(full > 0)),
                               'mean_sign_zero_mode_without_bare_cubic': float(np.mean(np.sign(full_nocubic))),
                               'ensemble_sd_of_mean_sign': sd_m,
                               'growth_weighted_power_frac_k_lt_0.25': float(np.sum(wpow[kphys < 0.25]) / np.sum(wpow)),
                               'growth_weighted_power_frac_lowest_shell': float(np.sum(wpow[np.isclose(kphys, 2 * np.pi / (N * dx))]) / np.sum(wpow)),
                               'effective_k_peak': float(np.sqrt(np.sum(wpow * kp2) / np.sum(wpow))),
                               'bare_cubic_over_cross_source': float(abs(J3[j]) / max(abs(J[j]), 1e-300)),
                               'erfinv_m': float(erfinv(np.clip(np.mean(np.sign(full)), -0.999999, 0.999999)))}
    out['predictions'] = marks
    sel = slice(0, len(TGRID), 10)
    out['series'] = {'a': a[sel].round(4).tolist(), 'S0_cross': (chi0 / a)[sel].tolist(), 'S0_total': S0c[sel].tolist(),
                     'sigma_f': sig_f[sel].tolist(), 'J_cross': J[sel].tolist(), 'J_bare_cubic': J3[sel].tolist()}
    return out


def main():
    raw, outp = Path(sys.argv[1]), Path(sys.argv[2])
    names = sys.argv[3:] or list(RUNS)
    res = json.loads(outp.read_text()) if outp.exists() else {}
    for name in names:
        N, dx, seed = RUNS[name]
        q = codex_noise(N, dx, seed)
        rec_noise = json.loads((raw / f'gaussian_N{N}_s{seed}.noise.json').read_text())
        cov = [[float(np.mean(q[i] * q[j])) for j in range(3)] for i in range(3)]
        dev = float(np.max(np.abs(np.array(cov) - np.array(rec_noise['covariance']))))
        rec = np.genfromtxt(raw / f'{name}.csv', delimiter=',', names=True)
        r = analyze(q, N, dx, rec=rec, label=name)
        r['regeneration_check'] = {'max_abs_covariance_deviation_vs_codex_record': dev, 'numpy': np.__version__,
                                   'codex_numpy': rec_noise['numpy_version']}
        res[name] = r
        p = r['predictions']['rms_0.5']
        print(name, 'covdev=%.2e' % dev, 'C_uS=%.4f C_wS=%.4f' % (r['C_uS'], r['C_wS']),
              'a=%.3f S0=%.3g sig=%.3g m_lin=%.4f sd_m=%.4f' % (p['a'], p['S0'], p['sigma_f'], p['mean_sign_with_zero_mode'], p['ensemble_sd_of_mean_sign']), flush=True)
        outp.write_text(json.dumps(res, indent=1))
        del q


if __name__ == '__main__':
    main()
