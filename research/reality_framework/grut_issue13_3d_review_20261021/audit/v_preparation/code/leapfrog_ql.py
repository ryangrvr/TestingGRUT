"""Independent quasi-linear check of the sparse-preparation zero mode (verifier of track 'preparation').

Written from expanding.cpp directly, NOT from the audit's quasilinear.py. Differences from the audit:
  * every amplitude is advanced with the solver's own KDK leapfrog at the run's dt (default 0.06), using the lattice
    eigenvalue khat^2 (exact for a plane wave under the 7-point Laplacian), so the linear map is the solver's map;
  * the homogeneous X background is advanced with the same leapfrog (not the analytic orbit);
  * the S zero mode chi0 is advanced with the box-averaged kick (a^2 - |X_bg|^2) chi0 - <dF^2 chi_S>, where
    dF^2 = 2(U du + W dw) per mode, cos and sin parts evolved separately ( <.> = (cos.cos + sin.sin)/2 ).
Second-order in the preparation amplitudes; S^2-driven X response and cubic S terms are omitted [E].
Preparation inputs [P]: du = 1e-3 F f_u, dw = 1e-3 F f_w, S = 1e-2 f_S; p_u = omega q_w, p_w = -omega q_u,
p_S = rate q_S (expanding.cpp initial()).

Usage: python3 -I leapfrog_ql.py <codex_raw_dir> <out_json>
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
import numpy as np

EPS_X, EPS_S, RATE = 1e-3, 1e-2, 1 / 40.0


def evolve(n, c, N, dx, kappa, dt=0.06, T=48.0, xscale=(1.0, 1.0), extra=None):
    """Return dict of time series. c columns: cu su cw sw cS sS. xscale scales (u, w) X perturbations."""
    k2 = np.sum(4 * np.sin(np.pi * np.asarray(n, float) / N) ** 2, axis=1) / dx ** 2
    F = kappa ** (1 / 3); x = F * F; om = np.sqrt(0.1 * x)
    cu, su, cw, sw, cS, sS = [np.asarray(v, float) for v in np.asarray(c).T]
    cu, su = cu * xscale[0], su * xscale[0]
    cw, sw = cw * xscale[1], sw * xscale[1]
    # state per mode and part (cos=0, sin=1): du, dw, chi; momenta
    du = np.stack([EPS_X * F * cu, EPS_X * F * su]); dw = np.stack([EPS_X * F * cw, EPS_X * F * sw])
    ch = np.stack([EPS_S * cS, EPS_S * sS])
    U, W = F, 0.0
    pU, pW = om * W, -om * U
    pdu, pdw = om * dw, -om * du
    pch = RATE * ch
    # separate u-only and w-only sourced zero modes (linear in X amplitudes, so they add)
    c0 = np.zeros(3); p0 = np.zeros(3)  # total, u-part, w-part

    # to split u and w contributions we evolve two extra copies of the X perturbations
    duU, dwU = np.stack([EPS_X * F * cu, EPS_X * F * su]), np.zeros_like(du)
    duW, dwW = np.zeros_like(du), np.stack([EPS_X * F * cw, EPS_X * F * sw])
    pduU, pdwU = om * dwU, -om * duU
    pduW, pdwW = om * dwW, -om * duW
    steps = int(round(T / dt))
    out = {'t': [0.0], 'a': [1.0], 'S0': [0.0], 'S0_u': [0.0], 'S0_w': [0.0], 'sigma_f': [EPS_S * np.sqrt(np.sum(cS ** 2 + sS ** 2) / 2)],
           'J': [-0.5 * np.sum(2 * (U * du + W * dw) * ch)], 'g_k0': [1.0]}
    g0 = np.array([1.0]); pg0 = np.array([RATE])  # an unsourced k=0 S mode for the head-start comparison

    def kick(h, a):
        nonlocal pU, pW, pdu, pdw, pch, pduU, pdwU, pduW, pdwW, p0, pg0
        f_bg = U * U + W * W
        pU_ = -(0.1 * f_bg) * U; pW_ = -(0.1 * f_bg) * W
        res = []
        for (xu, xw) in ((du, dw), (duU, dwU), (duW, dwW)):
            proj = U * xu + W * xw
            res.append((-k2 * xu - 0.1 * (f_bg * xu + 2 * U * proj), -k2 * xw - 0.1 * (f_bg * xw + 2 * W * proj), 2 * proj))
        fch = -(k2 + f_bg - a * a) * ch
        J = np.array([-0.5 * np.sum(r[2] * ch) for r in res])
        f0 = (a * a - f_bg) * c0 + J
        pU = pU + h * pU_; pW = pW + h * pW_
        pdu = pdu + h * res[0][0]; pdw = pdw + h * res[0][1]
        pduU = pduU + h * res[1][0]; pdwU = pdwU + h * res[1][1]
        pduW = pduW + h * res[2][0]; pdwW = pdwW + h * res[2][1]
        pch = pch + h * fch
        p0 = p0 + h * f0
        pg0 = pg0 + h * (a * a - f_bg) * g0
        return J[0]

    for it in range(1, steps + 1):
        t0, t1 = (it - 1) * dt, it * dt
        a0, a1 = 1 + RATE * t0, 1 + RATE * t1
        kick(0.5 * dt, a0)
        U = U + dt * pU; W = W + dt * pW
        du = du + dt * pdu; dw = dw + dt * pdw
        duU = duU + dt * pduU; dwU = dwU + dt * pdwU
        duW = duW + dt * pduW; dwW = dwW + dt * pdwW
        ch = ch + dt * pch
        c0 = c0 + dt * p0
        g0 = g0 + dt * pg0
        J = kick(0.5 * dt, a1)
        out['t'].append(t1); out['a'].append(a1)
        out['S0'].append(c0[0] / a1); out['S0_u'].append(c0[1] / a1); out['S0_w'].append(c0[2] / a1)
        out['sigma_f'].append(np.sqrt(np.sum(ch ** 2) / 2) / a1)
        out['J'].append(J); out['g_k0'].append(float(g0[0]))
    out = {k: np.array(v) for k, v in out.items()}
    out['_final_ch'] = ch
    return out, k2


def amplitudes_at(n, c, N, dx, kappa, a_target, dt=0.06, xscale=(1.0, 1.0)):
    """Linear S amplitudes (cos, sin) per mode at a_target and the zero mode, for field evaluation."""
    T = round((a_target - 1) * 40 / dt) * dt
    o, k2 = evolve(n, c, N, dx, kappa, dt=dt, T=T, xscale=xscale)
    return o, k2


def field_sign_stats(n, ch, S0, a, N, sub=1):
    """Evaluate S_lin = (sum ch_c cos + ch_s sin)/a on the solver grid (every sub-th point) plus S0."""
    idx = (np.arange(0, N, sub) + 0.5) / N
    f = np.zeros((len(idx),) * 3)
    for m in range(len(n)):
        ex = np.exp(2j * np.pi * n[m][0] * idx); ey = np.exp(2j * np.pi * n[m][1] * idx); ez = np.exp(2j * np.pi * n[m][2] * idx)
        ph = (ex[:, None] * ey[None, :])[:, :, None] * ez[None, None, :]
        f += (ch[0, m] * ph.real + ch[1, m] * ph.imag) / a
    full = f + S0
    return {'mean_sign_fluct': float(np.mean(np.sign(f))), 'mean_sign_total': float(np.mean(np.sign(full))),
            'm_05_total': float((np.mean(full > 0.5) - np.mean(full < -0.5)) / max(np.mean(np.abs(full) > 0.5), 1e-300)),
            'rms_fluct_grid': float(np.sqrt(np.mean(f * f)))}


RUNS = {
    'main_s101': ('modes_L76.8_s101.txt', 192, 0.4, 2.0),
    'main_s202': ('modes_L76.8_s202.txt', 192, 0.4, 2.0),
    'main_s303': ('modes_L76.8_s303.txt', 192, 0.4, 2.0),
    'kappa4_s101': ('modes_L76.8_s101.txt', 192, 0.4, 4.0),
    'large_s101': ('modes_L115.2_s101.txt', 288, 0.4, 2.0),
}


def main():
    raw, outp = Path(sys.argv[1]), Path(sys.argv[2])
    res = {}
    for run, (tab, N, dx, kap) in RUNS.items():
        t = np.loadtxt(raw / tab); n = t[:, :3].astype(int); c = t[:, 3:]
        o, k2 = evolve(n, c, N, dx, kap)
        rec = np.genfromtxt(raw / f'{run}.csv', delimiter=',', names=True)
        comp = []
        for i in range(len(rec)):
            j = int(round(rec['t'][i] / 0.06))
            if rec['a'][i] > 1.80:
                break
            comp.append({'a': float(rec['a'][i]), 'rec': float(rec['mean_S'][i]), 'pred': float(o['S0'][j]),
                         'ratio': float(rec['mean_S'][i] / o['S0'][j]) if o['S0'][j] != 0 else None,
                         'pred_u': float(o['S0_u'][j]), 'pred_w': float(o['S0_w'][j]), 'sigma_f': float(o['sigma_f'][j])})
        tot = np.sqrt(o['sigma_f'] ** 2 + o['S0'] ** 2)
        j5 = int(np.argmax(tot >= 0.5))
        a5 = float(o['a'][j5])
        # head start: ln(g_k0) - ln(|g_kmin|) evaluated with the leapfrog amplitudes
        kmin_i = int(np.argmin(k2))
        o2, _ = evolve(n[[kmin_i]], np.array([[0, 0, 0, 0, 1.0, 0.0]]), N, dx, kap, T=o['t'][j5])
        g_kmin = o2['sigma_f'][-1] * o2['a'][-1] / (EPS_S * np.sqrt(0.5))
        head = float(np.log(o['g_k0'][j5] / g_kmin))
        # field sign statistics at rms 0.5 using amplitudes at that time
        oo, _ = evolve(n, c, N, dx, kap, T=o['t'][j5])
        st = field_sign_stats(n, oo['_final_ch'], oo['S0'][-1], oo['a'][-1], N, sub=2)
        res[run] = {'N': N, 'L': N * dx, 'kappa': kap, 'k_min_phys': float(np.min(np.linalg.norm(n, axis=1)) * 2 * np.pi / (N * dx)),
                    'C_uS': float(np.sum(c[:, 0] * c[:, 4] + c[:, 1] * c[:, 5]) / 2), 'C_wS': float(np.sum(c[:, 2] * c[:, 4] + c[:, 3] * c[:, 5]) / 2),
                    'J0': float(o['J'][0]), 'a_rms05': a5, 'S0_rms05': float(o['S0'][j5]), 'S0_u_rms05': float(o['S0_u'][j5]),
                    'S0_w_rms05': float(o['S0_w'][j5]), 'sigma_f_rms05': float(o['sigma_f'][j5]),
                    'S0_over_sigma_f_rms05': float(o['S0'][j5] / o['sigma_f'][j5]), 'head_start_efolds_k0_vs_kmin': head,
                    'sign_stats_rms05': st, 'comparison': comp}
        good = [x for x in comp if x['a'] >= 1.12 and x['a'] <= 1.601 and abs(x['rec']) > 1e-6]
        dev = max(abs(x['ratio'] - 1) for x in good)
        res[run]['max_abs_rel_dev_1.12_1.60'] = dev
        print(run, 'C_uS=%.4f C_wS=%.4f J0=%.4e a5=%.3f S0=%.4f (u %.4f w %.4f) sig=%.4f ratio=%.3f head=%.2f m_sign=%.3f maxdev(1.12-1.60)=%.4f' % (
            res[run]['C_uS'], res[run]['C_wS'], res[run]['J0'], a5, res[run]['S0_rms05'], res[run]['S0_u_rms05'], res[run]['S0_w_rms05'],
            res[run]['sigma_f_rms05'], res[run]['S0_over_sigma_f_rms05'], head, st['mean_sign_total'], dev), flush=True)
        for x in comp:
            if 1.1 < x['a'] < 1.67:
                print('   a=%.3f rec=% .4e pred=% .4e ratio=%.4f' % (x['a'], x['rec'], x['pred'], x['ratio'] if x['ratio'] else np.nan))
    outp.write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
