"""Initial-condition diagnostics + quasi-linear zero-mode prediction for Codex's sparse 64-mode seeds.

Usage: python3 -I sparse_diagnostics.py <raw_dir> <out_json>
<raw_dir> holds Codex's extracted raw/ (mode tables and csv series), read only.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import quasilinear as ql  # noqa: E402

RAW = Path(sys.argv[1]) if __name__ == '__main__' else None
OUT = Path(sys.argv[2]) if __name__ == '__main__' else None

RUNS = {  # run -> (mode table, N, dx, kappa) ; dt/grid variants share main_s101's polynomial
    'main_s101': ('modes_L76.8_s101.txt', 192, 0.4, 2.0),
    'main_s202': ('modes_L76.8_s202.txt', 192, 0.4, 2.0),
    'main_s303': ('modes_L76.8_s303.txt', 192, 0.4, 2.0),
    'kappa4_s101': ('modes_L76.8_s101.txt', 192, 0.4, 4.0),
    'large_s101': ('modes_L115.2_s101.txt', 288, 0.4, 2.0),
    'grid256_s101': ('modes_L76.8_s101.txt', 256, 0.3, 2.0),
    'grid384_s101': ('modes_L76.8_s101.txt', 384, 0.2, 2.0),
    'dt_half_s101': ('modes_L76.8_s101.txt', 192, 0.4, 2.0),
}
TGRID = np.round(np.arange(0, 48.0001, 0.05), 10)


def load(name):
    a = np.loadtxt(RAW / name)
    return a[::, :3].astype(int), a[:, 3:]


def triads(n):
    """Signed-mode triples with s1 n1 + s2 n2 + s3 n3 = 0 (ordered)."""
    signed = [(m, s) for m in range(len(n)) for s in (1, -1)]
    look = {tuple(s * n[m]): (m, s) for m, s in signed}
    out = []
    for m1, s1 in signed:
        for m2, s2 in signed:
            key = tuple(-(s1 * n[m1] + s2 * n[m2]))
            if key in look:
                m3, s3 = look[key]
                out.append((m1, s1, m2, s2, m3, s3))
    return out


def cube_mean(z, tri):
    """<f^3> for f = sum Re(z_m e^{i theta_m}) = 1/2 sum (z e^{i th} + z* e^{-i th})."""
    tot = 0.0
    for m1, s1, m2, s2, m3, s3 in tri:
        p = 1.0 + 0j
        for m, s in ((m1, s1), (m2, s2), (m3, s3)):
            p *= z[..., m] if s == 1 else np.conj(z[..., m])
        tot = tot + p
    return np.real(tot) / 8.0


def field_on_grid(n, cos_c, sin_c, N, sub):
    """Evaluate sum c cos(theta)+s sin(theta), theta = 2 pi n.(i+.5)/N, on every sub-th point."""
    idx = (np.arange(0, N, sub) + 0.5) / N
    f = np.zeros((len(idx),) * 3)
    for m in range(len(n)):
        ex = np.exp(2j * np.pi * n[m, 0] * idx)
        ey = np.exp(2j * np.pi * n[m, 1] * idx)
        ez = np.exp(2j * np.pi * n[m, 2] * idx)
        ph = ex[:, None, None] * ey[None, :, None] * ez[None, None, :]
        f += cos_c[m] * ph.real + sin_c[m] * ph.imag
    return f


def analyze(run, table, N, dx, kappa):
    n, c = load(table)
    rec = np.genfromtxt(RAW / f'{run}.csv', delimiter=',', names=True)
    return analyze_arrays(run, table, n, c, N, dx, kappa, rec)


def analyze_arrays(run, table, n, c, N, dx, kappa, rec=None):
    L = N * dx
    kphys = np.linalg.norm(n, axis=1) * 2 * np.pi / L
    k2 = ql.khat2_int(n, N, dx)
    F = kappa ** (1 / 3)
    cu, su, cw, sw, cS, sS = c.T
    pS = (cS ** 2 + sS ** 2) / 2
    CuS = (cu * cS + su * sS) / 2
    CwS = (cw * cS + sw * sS) / 2
    order = np.argsort(kphys)
    g, _ = ql.s_growth(k2, kappa, TGRID)
    Au, Aw = ql.x_response(k2, kappa, TGRID)
    # J(t) = -<dF2 chi_S> ; dF2_m = eps_X F (Au cu + Aw cw) [cos part], same for sin; chi_S,m = eps_S g c_S
    J = -np.sum(ql.EPS_X * F * ql.EPS_S * g * (Au * CuS + Aw * CwS), axis=1)
    J_u = -np.sum(ql.EPS_X * F * ql.EPS_S * g * Au * CuS, axis=1)
    J_w = -np.sum(ql.EPS_X * F * ql.EPS_S * g * Aw * CwS, axis=1)
    chi0 = ql.zero_mode(TGRID, J, kappa)
    chi0_u = ql.zero_mode(TGRID, J_u, kappa)
    chi0_w = ql.zero_mode(TGRID, J_w, kappa)
    a = ql.a_of(TGRID)
    sig_f = ql.EPS_S * np.sqrt(np.sum(g ** 2 * pS, axis=1)) / a  # physical rms of fluctuations
    S0 = chi0 / a
    # cubic source from triads of the linear fluctuation field (bare -<chi_f^3>), for magnitude comparison
    tri = triads(n)
    z = ql.EPS_S * g * (cS - 1j * sS)[None, :]
    cub = cube_mean(z, tri) if tri else np.zeros_like(TGRID)
    J3 = -cub
    chi0_c = ql.zero_mode(TGRID, J + J3, kappa)
    # comparison with Codex's recorded nonlinear series (read only)
    comp = []
    for i in range(len(rec) if rec is not None else 0):
        t = rec['t'][i]
        if t > 30.5:
            break
        j = int(round(t / 0.05))
        comp.append({'a': float(rec['a'][i]), 'mean_S_recorded': float(rec['mean_S'][i]),
                     'S0_pred_cross': float(S0[j]), 'S0_pred_cross_plus_bare_cubic': float(chi0_c[j] / a[j]),
                     'sigma_f_pred': float(sig_f[j]), 'max_abs_S_recorded': float(max(-rec['min_S'][i], rec['max_S'][i]))})
    # a_nl: total predicted rms reaches 0.3 (0.2 and 0.5 as sensitivity)
    tot = np.sqrt(sig_f ** 2 + S0 ** 2)
    out = {'run': run, 'table': table, 'N': N, 'dx': dx, 'L': L, 'kappa': kappa,
           'k_min': float(kphys[order[0]]), 'k_lowest6': kphys[order[:6]].round(4).tolist(),
           'S_power_lowest6': pS[order[:6]].round(5).tolist(),
           'S_power_k_lt_0p35': float(pS[kphys < 0.35].sum()), 'S_power_k_lt_0p5': float(pS[kphys < 0.5].sum()),
           'n_modes_k_lt_0p5': int((kphys < 0.5).sum()),
           'a_cross_k0': float(np.sqrt(kappa ** (2 / 3))), 'a_cross_kmin': float(np.sqrt(kappa ** (2 / 3) + k2[order[0]])),
           'C_uS': float(CuS.sum()), 'C_wS': float(CwS.sum()),
           'J0_codex_convention': float(-2 * F * F * ql.EPS_S * ql.EPS_X * CuS.sum()),
           'n_triads_ordered': len(tri), 'comparison': comp}
    marks = {}
    for thr in (0.2, 0.3, 0.5, 1.0):
        j = int(np.argmax(tot >= thr)) if np.any(tot >= thr) else len(TGRID) - 1
        aj = float(a[j])
        lam0 = float(ql.lam_exponent(0.0, kappa, aj))
        lamk = float(ql.lam_exponent(k2[order[0]], kappa, aj))
        cosc = ql.EPS_S * g[j] * cS / aj
        sinc = ql.EPS_S * g[j] * sS / aj
        f = field_on_grid(n, cosc, sinc, N, max(1, N // 96))
        full = f + S0[j]
        marks[f'rms_{thr}'] = {
            'a': aj, 'S0': float(S0[j]), 'sigma_f': float(sig_f[j]), 'S0_over_sigma_f': float(S0[j] / sig_f[j]),
            'S0_from_u_only': float(chi0_u[j] / aj), 'S0_from_w_only': float(chi0_w[j] / aj),
            'S0_cross_plus_bare_cubic': float(chi0_c[j] / aj),
            'bare_cubic_source_over_cross_source': float(abs(J3[j]) / max(abs(J[j]), 1e-300)),
            'WKB_gap_efolds_k0_minus_kmin': lam0 - lamk,
            'mean_sign_fluct_only': float(np.mean(np.sign(f))),
            'mean_sign_with_zero_mode': float(np.mean(np.sign(full))),
            'pplus_with_zero_mode': float(np.mean(full > 0)),
            'mean_S_fluct_only_grid': float(f.mean()),
            'growth_weighted_power_lowest_mode_fraction': float(g[j, order[0]] ** 2 * pS[order[0]] / np.sum(g[j] ** 2 * pS)),
            'growth_weighted_C_uS': float(np.sum(g[j] * CuS) / np.sqrt(np.sum(g[j] ** 2 * pS))),
            'growth_weighted_C_wS': float(np.sum(g[j] * CwS) / np.sqrt(np.sum(g[j] ** 2 * pS))),
        }
    out['predictions'] = marks
    # time series (subsampled) for plotting / audit
    sel = slice(0, len(TGRID), 20)
    out['series'] = {'a': a[sel].round(4).tolist(), 'S0_cross': S0[sel].tolist(), 'S0_u_only': (chi0_u / a)[sel].tolist(),
                     'S0_w_only': (chi0_w / a)[sel].tolist(), 'sigma_f': sig_f[sel].tolist(), 'J_cross': J[sel].tolist(),
                     'J_bare_cubic': J3[sel].tolist()}
    return out


def main():
    res = {}
    for run, (table, N, dx, kappa) in RUNS.items():
        res[run] = analyze(run, table, N, dx, kappa)
        p = res[run]['predictions']['rms_0.3']
        print(run, 'kmin=%.3f C_uS=%.3f C_wS=%.3f' % (res[run]['k_min'], res[run]['C_uS'], res[run]['C_wS']),
              'a_nl=%.3f S0=%.3g sig=%.3g ratio=%.3g m_pred=%.3f gap=%.2f cub/cross=%.2g' % (
                  p['a'], p['S0'], p['sigma_f'], p['S0_over_sigma_f'], p['mean_sign_with_zero_mode'],
                  p['WKB_gap_efolds_k0_minus_kmin'], p['bare_cubic_source_over_cross_source']), flush=True)
    OUT.write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
