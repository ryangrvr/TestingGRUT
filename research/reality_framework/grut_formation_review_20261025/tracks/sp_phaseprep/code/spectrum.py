#!/usr/bin/env python3
"""sp_phaseprep: phase spectrum, light source-phase mode and its decay constant at H_i (own code).

Reuses our gp_matching model (code/model.py: sympy-exact canonical 6x6 Hessian, random complex starts) at
H_i = 1.1211167304108865e13 GeV, with 4.7e13 GeV as a control against the gp_matching values (0.378 H, f_eff 0.377 H).

Definitions (H = 1 units):
  - canonical real components z_i = (x_i1 + i x_i2)/sqrt2; the Hessian is the mass-squared matrix.
  - gauge-invariant phases (gp_matching convention z = R e^{-i theta}): psi1 = 3thX - thY (kappa), psi2 = 2thY - thS
    (mu source phase), delta6 = 6thX - thS = 2 psi1 + psi2. Codex's supplied mismatch (Sigma -> -|Sigma| e^{i delta},
    X, Y fixed) shifts psi2 and delta6 by the same amount, with psi1 = 0.
  - f_eff = 1 / |d psi2 / d phi_c| along the light angular eigenvector (phi_c canonical), i.e. psi2 = phi_c / f_eff.
  - per-e-fold de Sitter kick of psi2: H / (2 pi f_eff) [D, standard].
  - coarse-grained variance of a free massive mode (Bunch-Davies), modes with k/(aH) < x_c at the end of inflation:
      sigma^2(x_c) = int_0^{x_c} dx/x P(x),  P(x) = (H^2/8pi) x^3 e^{-pi Im nu} |H^(1)_nu(x)|^2,  nu = sqrt(9/4 - m^2/H^2)
    (exact mode functions; for x_c << 1 and real nu this is (H/2pi)^2 C_nu^2 x_c^(3-2nu)/(3-2nu),
     C_nu = 2^(nu-3/2) Gamma(nu)/Gamma(3/2); the x_c -> 1, nu -> 3/2 limit is Starobinsky-Yokoyama 3H^4/(8 pi^2 m^2)).

Usage: python3 spectrum.py ../results/spectrum.json
"""
import json
import sys
from math import gamma
from pathlib import Path
import numpy as np
import mpmath as mp

TR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(TR / 'gp_matching' / 'code'))
from model import build  # noqa: E402  (our gp_matching code)
from point_full import analyse  # noqa: E402

H_STAR = 1.1211167304108865e13


def P_mode(x, m2):
    """Dimensionless power P/H^2 of a free scalar of mass^2 m2 (H^2 units) at k/(aH) = x (exact BD)."""
    nu = mp.sqrt(mp.mpf(9) / 4 - m2)
    h = mp.hankel1(nu, x)
    pref = mp.exp(-mp.pi * mp.im(nu))
    return float(mp.re(x**3 * pref * abs(h)**2 / (8 * mp.pi)))


def sigma2_coarse(m2, xc):
    """int_0^xc dx/x P(x) (H^2 units). Small-x part analytic for real nu, numerical quadrature above."""
    nu2 = 9 / 4 - m2
    if nu2 > 0:
        nu = np.sqrt(nu2)
        Cn = 2**(nu - 1.5) * gamma(nu) / gamma(1.5)
        x0 = min(1e-4, xc)
        low = Cn**2 / (4 * np.pi**2) * x0**(3 - 2 * nu) / (3 - 2 * nu)
    else:
        x0 = min(1e-4, xc)
        # heavy: P ~ x^3 at small x; integrate numerically from 1e-8
        low = float(mp.quad(lambda t: P_mode(t, m2) / t, [1e-8, x0]))
    if xc <= x0:
        return low
    hi = float(mp.quad(lambda t: P_mode(t, m2) / t, mp.linspace(x0, xc, 6)))
    return low + hi


def summarize(label, H):
    m = build(Hinf=H)
    a = analyse(m, label)
    light = a['light_mode']; heavy = a['heavy_mode']
    f_eff = 1 / abs(light['dpsi2'])
    R = a['moduli_over_H']
    m2 = light['mass2_over_H2']
    out = dict(H_GeV=H, moduli_over_H=R, V_min=a['V_min_over_H4'],
               radial_mass_over_H=np.sqrt(np.array(a['radial_mass2_over_H2'])).tolist(),
               angular_mass_over_H=[md['mass_over_H'] for md in a['angular_modes']],
               light=dict(mass_over_H=light['mass_over_H'], mass2=m2, dpsi1=light['dpsi1'], dpsi2=light['dpsi2'],
                          ddelta6=light['ddelta'], dtheta=light['dtheta_per_canonical'], f_eff_over_H=f_eff,
                          f_eff_over_sqrt2_absSigma=f_eff / (np.sqrt(2) * R[2]),
                          kick_per_efold_rad=1 / (2 * np.pi * f_eff)),
               heavy=dict(mass_over_H=heavy['mass_over_H'], dpsi1=heavy['dpsi1'], dpsi2=heavy['dpsi2']),
               A_over_H4=a['A_over_H4'], B_over_H4=a['B_over_H4'], G_psi=a['G_psi_metric'],
               reduced_gauge_invariant_mass2=a['reduced_gauge_invariant_mass2'],
               starts=a['starts'], vector_mass_over_H=a['vector_mass_over_H'])
    # Sigma-quartic-with-linear-source relations [D]: m_t^2 = 2 lam_S |S|^2 (h<<|S|), m_r^2 = 3 m_t^2
    hS = 1e11 / H
    out['sigma_quartic_relations'] = dict(mt2_pred=2 * (R[2]**2 - hS**2 / 2), mr2_pred=2 * (3 * R[2]**2 - hS**2 / 2))
    # stochastic equilibrium of psi2 at fixed radii, P ~ exp(kappa cos psi2), kappa = 8 pi^2 B / 3
    from scipy.special import iv
    kap = 8 * np.pi**2 * a['B_over_H4'] / 3
    t = np.linspace(-np.pi, np.pi, 400001); w = np.exp(kap * (np.cos(t) - 1))
    Z = np.trapezoid(w, t)
    out['fixed_radius_equilibrium'] = dict(kappa=kap, mean_cos=float(iv(1, kap) / iv(0, kap)),
                                           rms_psi2=float(np.sqrt(np.trapezoid(w * t**2, t) / Z)),
                                           gaussian_rms=float(1 / np.sqrt(kap)),
                                           relax_efolds_3_over_m2=3 / m2)
    # coarse-grained exact-BD spread of the light mode (and of each radial / heavy mode) for several cut-offs
    cuts = [1.0, 1 / 40, 1 / 100, 1 / 500]
    out['coarse_grained'] = {}
    for name, mm2, fe in [('light_phase', m2, f_eff)]:
        rows = []
        for xc in cuts:
            s2 = sigma2_coarse(mm2, xc)
            rows.append(dict(x_c=xc, sigma_canonical_over_H=float(np.sqrt(s2)), sigma_psi2_rad=float(np.sqrt(s2) / fe)))
        out['coarse_grained'][name] = rows
    out['coarse_grained']['SY_light_3_over_8pi2m2'] = dict(
        sigma_canonical_over_H=float(np.sqrt(3 / (8 * np.pi**2 * m2))), sigma_psi2_rad=float(np.sqrt(3 / (8 * np.pi**2 * m2)) / f_eff))
    # CMB-scale amplitude of the light mode (free-field decay), N_* = 50, 60
    nu = np.sqrt(9 / 4 - m2); Cn = 2**(nu - 1.5) * gamma(nu) / gamma(1.5)
    out['cmb_amplitude'] = {f'N{N}': dict(dpsi2_rms_per_lnk=float(Cn / (2 * np.pi * f_eff) * np.exp(-(1.5 - nu) * N)),
                                          n_iso_minus_1=float(3 - 2 * nu)) for N in (50, 55, 60)}
    return out


def mass_scan(f_eff, cuts=(1 / 40, 1 / 100, 1 / 500), targets=(0.03, 0.08)):
    """Phase spread vs a (supplied) phase mass m_t, at fixed f_eff: what mass makes the spread <= target."""
    ms = np.round(np.concatenate([np.linspace(0.3, 1.45, 24), np.linspace(1.5, 3.0, 7)]), 4)
    rows = []
    for mt in ms:
        r = dict(m_over_H=float(mt))
        for xc in cuts:
            r[f'sigma_psi_xc_{1/xc:.0f}'] = float(np.sqrt(sigma2_coarse(mt**2, xc)) / f_eff)
        r['sigma_psi_xc_1'] = float(np.sqrt(sigma2_coarse(mt**2, 1.0)) / f_eff)
        rows.append(r)
    need = {}
    for xc in list(cuts) + [1.0]:
        key = f'sigma_psi_xc_{1/xc:.0f}'
        for tg in targets:
            ok = [r['m_over_H'] for r in rows if r[key] <= tg]
            need[f'{key}<={tg}'] = min(ok) if ok else None
    return dict(f_eff_over_H=f_eff, rows=rows, minimal_mass=need)


def main(out):
    res = dict(H_star=summarize('H_i = 1.1211e13 GeV, c_H = 1', H_STAR),
               control_47e13=summarize('control 4.7e13 GeV', 4.7e13))
    res['mass_scan'] = mass_scan(res['H_star']['light']['f_eff_over_H'])
    Path(out).write_text(json.dumps(res, indent=1, default=float) + '\n')
    for k in ('H_star', 'control_47e13'):
        r = res[k]
        print(k, 'moduli', np.round(r['moduli_over_H'], 6), 'rad', np.round(r['radial_mass_over_H'], 4),
              'ang', np.round(r['angular_mass_over_H'], 5))
        print('  light', r['light'])
        print('  B', r['B_over_H4'], 'eq', r['fixed_radius_equilibrium'], 'quartic rel', r['sigma_quartic_relations'])
        print('  coarse', r['coarse_grained'])
        print('  cmb', r['cmb_amplitude'])
    print(res['mass_scan']['minimal_mass'])


if __name__ == '__main__':
    main(sys.argv[1])
