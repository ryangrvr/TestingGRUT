#!/usr/bin/env python3
"""Inflationary stationary point, full 6x6 mass spectrum and phase dictionary (own method).

Method differs from Codex point.py: minimization over the six canonical real components from random
complex starts (no aligned-phase ansatz), Newton polish, full canonical Hessian (radial and angular
blocks emerge from it), an analytic large-field approximation, and the gauge-invariant reduction.
Outputs results/point_full.json.
"""
import json, sys
from math import gamma
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import iv

sys.path.insert(0, str(Path(__file__).resolve().parent))
from model import build, radial_V, GEV, COUP

OUT = Path(__file__).resolve().parents[1] / 'results'
SMAX = 9.1e-6   # Planck 2018 uncorrelated CDM isocurvature (beta_iso<0.038), as used in affleck_dine [Q]
Q = np.array([1., 3., 6.])


def find_min(m, nstart=48, seed=7):
    rng = np.random.default_rng(seed)
    best = None; rows = []
    for k in range(nstart):
        R = rng.uniform(0.1, 40, 3) * np.array([1, .2, .05])
        ph = rng.uniform(-np.pi, np.pi, 3)
        z0 = np.ravel([[np.sqrt(2) * r * np.cos(p), np.sqrt(2) * r * np.sin(p)] for r, p in zip(R, ph)])
        r = minimize(m['V'], z0, jac=m['grad'], method='BFGS', options=dict(gtol=1e-11, maxiter=20000))
        z = r.x
        for _ in range(30):  # Newton polish, least squares handles the gauge zero mode
            dz = np.linalg.lstsq(m['hess'](z), -m['grad'](z), rcond=1e-12)[0]
            z = z + dz
            if np.linalg.norm(dz) < 1e-14 * max(1, np.linalg.norm(z)):
                break
        rows.append((m['V'](z), z))
        if best is None or m['V'](z) < best[0] - 1e-9:
            best = (m['V'](z), z)
    vals = np.array([r[0] for r in rows])
    return best, dict(n_starts=nstart, n_reaching_global=int(np.sum(np.abs(vals - best[0]) < 1e-7 * abs(best[0]))),
                      distinct_values=sorted(set(np.round(vals, 6).tolist()))[:6])


def to_gauge_rep(z):
    """Gauge rotate so that X is real positive; return complex fields."""
    X = (z[0] + 1j * z[1]) / np.sqrt(2); Y = (z[2] + 1j * z[3]) / np.sqrt(2); S = (z[4] + 1j * z[5]) / np.sqrt(2)
    a = -np.angle(X)
    X, Y, S = X * np.exp(1j * a), Y * np.exp(3j * a), S * np.exp(6j * a)
    return X, Y, S


def analyse(m, label):
    (vmin, z), starts = find_min(m)
    X, Y, S = to_gauge_rep(z)
    R = np.abs([X, Y, S])
    Hs = m['hess'](z)
    # Rotate to a real representative so real/imag blocks separate.
    zr = np.array([np.sqrt(2) * X.real, np.sqrt(2) * X.imag, np.sqrt(2) * Y.real, np.sqrt(2) * Y.imag,
                   np.sqrt(2) * S.real, np.sqrt(2) * S.imag])
    Hr = m['hess'](zr)
    ev, evec = np.linalg.eigh(Hr)
    re_idx, im_idx = [0, 2, 4], [1, 3, 5]
    radial = np.linalg.eigvalsh(Hr[np.ix_(re_idx, re_idx)])
    angH = Hr[np.ix_(im_idx, im_idx)]
    ang_ev, ang_vec = np.linalg.eigh(angH)
    offblock = float(np.abs(Hr[np.ix_(re_idx, im_idx)]).max())
    # phase variables: field z_i = z0_i exp(-i theta_i), canonical imag comp c_i = -sqrt2 z0_i theta_i
    z0 = np.array([X.real, Y.real, S.real])
    dtheta_dc = -1 / (np.sqrt(2) * z0)          # theta_i per unit canonical imaginary amplitude
    W = np.array([[3., -1., 0.], [0., 2., -1.]])  # psi1 = 3thX - thY, psi2 = 2thY - thS
    wdelta = np.array([6., 0., -1.])              # delta = 6thX - thS = 2psi1 + psi2
    modes = []
    for k in range(3):
        v = ang_vec[:, k]
        th = dtheta_dc * v
        modes.append(dict(mass2_over_H2=float(ang_ev[k]), mass_over_H=float(np.sqrt(max(ang_ev[k], 0))),
                          dtheta_per_canonical=th.tolist(),
                          dpsi1=float(W[0] @ th), dpsi2=float(W[1] @ th), ddelta=float(wdelta @ th),
                          gauge_overlap=float(abs(v @ (np.sqrt(2) * np.abs(z0) * Q * np.sign(-z0)))
                                              / np.linalg.norm(np.sqrt(2) * np.abs(z0) * Q))))
    # Gauge-invariant reduction: K = diag(2R^2); G = (W K^-1 W^T)^-1
    K = np.diag(2 * R**2)
    G = np.linalg.inv(W @ np.linalg.inv(K) @ W.T)
    c = m['coup']
    A = 2 * c['kap'] * R[1] * R[0]**3
    B = 2 * m['mu'] * R[2] * R[1]**2
    Hpsi = np.diag([A, B])
    red = np.sort(np.real(np.linalg.eigvals(np.linalg.inv(G) @ Hpsi)))
    # light mode: canonical amplitude -> psi2, delta
    light = min([md for md in modes if md['mass2_over_H2'] > 1e-8], key=lambda d: d['mass2_over_H2'])
    heavy = max(modes, key=lambda d: d['mass2_over_H2'])
    m2 = light['mass2_over_H2']; nu = np.sqrt(9 / 4 - m2)
    Cnu = 2**(nu - 1.5) * gamma(nu) / gamma(1.5)
    iso = {}
    for N in [50, 60]:
        dchi = Cnu / (2 * np.pi) * np.exp(-(1.5 - nu) * N)
        dd = abs(light['ddelta']) * dchi
        iso[f'N{N}'] = dict(delta_delta_rms_rad=dd, S_over_cot_delta0=dd, ratio_to_Smax_at_cot1=dd / SMAX)
    # stochastic equilibrium of psi2 in its valley V=-B cos psi2 (Starobinsky-Yokoyama)
    kpar = 8 * np.pi**2 * B / 3
    eq = dict(exponent_8pi2B_over_3H4=kpar, mean_cos_psi2=float(iv(1, kpar) / iv(0, kpar)),
              prob_abs_psi2_gt_pi_over_2=float(_prob_gt(kpar, np.pi / 2)),
              relaxation_efolds_3H2_over_m2=3 / m2)
    # radial light mode -> amplitude isocurvature
    rv, rvec = np.linalg.eigh(Hr[np.ix_(re_idx, re_idx)])
    k0 = int(np.argmin(rv)); vr = rvec[:, k0]
    dlnRX = abs(vr[0]) / (np.sqrt(2) * R[0])
    nur = np.sqrt(9 / 4 - rv[k0]); Cr = 2**(nur - 1.5) * gamma(nur) / gamma(1.5)
    rad = {}
    for N in [50, 60]:
        d = dlnRX * Cr / (2 * np.pi) * np.exp(-(1.5 - nur) * N)
        rad[f'N{N}'] = dict(dlnRX_rms=d, S_if_dlnY_dlnRX_5=5 * d, ratio_to_Smax=5 * d / SMAX)
    # vector mass at the point
    mA = 1e-7 * np.sqrt(np.sum(2 * Q**2 * R**2))
    # analytic large-field approximation (derived in NOTES)
    lam_eff = c['lam_X'] - 1.5 * c['kap']**(4 / 3) * (2 * c['lam_Y'])**(-1 / 3)
    xa = np.sqrt(m['c_H'] / (2 * lam_eff)) if m['c_Y'] == 0 and m['c_S'] == 0 else np.nan
    ya = (c['kap'] / (2 * c['lam_Y']))**(1 / 3) * xa
    sa = (m['mu'] * ya**2 / (2 * c['lam_S']))**(1 / 3)
    out = dict(label=label, c_H=m['c_H'], c_Y=m['c_Y'], c_S=m['c_S'], starts=starts,
               V_min_over_H4=vmin, moduli_over_H=R.tolist(), moduli_GeV=(R * m['Hinf']).tolist(),
               phases_rel_X_rad=[float(np.angle(Y)), float(np.angle(S))],
               grad_norm=float(np.linalg.norm(m['grad'](zr))),
               spectrum_all_mass2_over_H2=ev.tolist(), radial_mass2_over_H2=radial.tolist(),
               angular_mass2_over_H2=ang_ev.tolist(), real_imag_offblock_max=offblock,
               angular_modes=modes, reduced_gauge_invariant_mass2=red.tolist(),
               G_psi_metric=G.tolist(), A_over_H4=A, B_over_H4=B,
               light_mode=light, heavy_mode=heavy,
               isocurvature_light_phase=iso, stochastic_equilibrium=eq, radial_isocurvature=rad,
               vector_mass_over_H=mA,
               analytic=dict(lam_eff=lam_eff, x=xa, y=ya, s=sa),
               Higgs_amplitude_over_vacuum=float(R[2] / (m['h'] / np.sqrt(2))),
               validity=dict(M_Y_over_H=m['M'], kap_X2_over_MY2_at_min=c['kap'] * R[0]**2 / m['M']**2,
                             kap_X2_over_MY2_single_field=c['kap'] * (m['c_H'] / (2 * c['lam_X'])) / m['M']**2,
                             lamY_Y2_over_MY2=c['lam_Y'] * R[1]**2 / m['M']**2,
                             Y_quadratic_elimination_over_actual=float(c['kap'] * R[0]**3 / m['M']**2 / R[1])))
    return out


def _prob_gt(k, a):
    t = np.linspace(-np.pi, np.pi, 200001); w = np.exp(k * np.cos(t))
    return np.trapezoid(w * (np.abs(t) > a), t) / np.trapezoid(w, t)


def fd_check(m, z, eps=1e-5):
    g0 = m['hess'](z); fd = np.zeros_like(g0)
    for i in range(6):
        e = np.zeros(6); e[i] = eps
        fd[:, i] = (m['grad'](z + e) - m['grad'](z - e)) / (2 * eps)
    return float(np.abs(fd - g0).max() / np.abs(g0).max())


def main():
    res = {}
    base = build()
    res['baseline'] = analyse(base, 'stated point, c_H=1, H=4.7e13 GeV')
    R = res['baseline']['moduli_over_H']
    zr = np.array([np.sqrt(2) * R[0], 0, -np.sqrt(2) * R[1], 0, -np.sqrt(2) * R[2], 0])
    res['baseline']['hessian_fd_rel_err'] = fd_check(base, zr)
    res['baseline']['radial_V_check'] = radial_V(R, base) - res['baseline']['V_min_over_H4']
    # one-field comparison quantities (affleck_dine): delta theta = H/(2 pi sqrt2 R), S = 6 dtheta
    R1 = np.sqrt(1 / (2 * COUP['lam_X']))
    res['one_field_reference'] = dict(R_over_H=R1, S_over_cot=6 / (2 * np.pi * np.sqrt(2) * R1))
    # sensitivity to tagged extra Hubble masses on Y and Sigma (NOT in 7737487)
    sens = []
    for cY in [-1., 0., 1.]:
        for cS in [-1., 0., 1.]:
            if cY == 0 and cS == 0:
                continue
            a = analyse(build(c_Y=cY, c_S=cS), f'c_Y={cY}, c_S={cS}')
            sens.append(dict(c_Y=cY, c_S=cS, moduli_over_H=a['moduli_over_H'],
                             angular_mass_over_H=[md['mass_over_H'] for md in a['angular_modes']],
                             light_ddelta=a['light_mode']['ddelta'], light_dpsi1=a['light_mode']['dpsi1'],
                             S_over_cot_N50=a['isocurvature_light_phase']['N50']['S_over_cot_delta0'],
                             radial_mass2=a['radial_mass2_over_H2']))
    res['sensitivity_extra_Hubble_masses'] = sens
    # H dependence at c_H=1 (the stated point fixes couplings; H is the review's supplied maximum)
    hs = []
    for H in [1e12, 1e13, 4.7e13]:
        a = analyse(build(Hinf=H), f'H={H}')
        hs.append(dict(H_GeV=H, moduli_over_H=a['moduli_over_H'], angular_mass_over_H=[md['mass_over_H'] for md in a['angular_modes']],
                       radial_mass2=a['radial_mass2_over_H2'], M_over_H=a['validity']['M_Y_over_H'],
                       S_over_cot_N50=a['isocurvature_light_phase']['N50']['S_over_cot_delta0']))
    res['H_scan'] = hs
    OUT.mkdir(exist_ok=True)
    (OUT / 'point_full.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    b = res['baseline']
    print('moduli GeV', b['moduli_GeV'], 'V', b['V_min_over_H4'], 'starts', b['starts'])
    print('radial m2', b['radial_mass2_over_H2'], 'angular m2', b['angular_mass2_over_H2'])
    print('reduced', b['reduced_gauge_invariant_mass2'], 'offblock', b['real_imag_offblock_max'], 'fd', b['hessian_fd_rel_err'])
    print('light', {k: b['light_mode'][k] for k in ['mass_over_H', 'dpsi1', 'dpsi2', 'ddelta', 'gauge_overlap']})
    print('heavy', {k: b['heavy_mode'][k] for k in ['mass_over_H', 'dpsi1', 'dpsi2', 'ddelta', 'gauge_overlap']})
    print('iso', b['isocurvature_light_phase'], 'eq', b['stochastic_equilibrium'])
    print('radial iso', b['radial_isocurvature'])
    print('analytic', b['analytic'], 'Higgs/vac', b['Higgs_amplitude_over_vacuum'], 'validity', b['validity'], 'mA/H', b['vector_mass_over_H'])
    print('one-field', res['one_field_reference'])
    for s in sens:
        print(s['c_Y'], s['c_S'], np.round(s['moduli_over_H'], 4), np.round(s['angular_mass_over_H'], 4), round(s['S_over_cot_N50'], 5))
    for s in hs:
        print(s)


if __name__ == '__main__':
    main()
