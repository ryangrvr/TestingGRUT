#!/usr/bin/env python3
"""Every annihilation / number-changing channel of a symmetric X/Xbar population at the light point (own code).

Kinematic thresholds, cross sections and rate/H through the epochs for:
 2->2 into lighter states (s s, s, D D, A', SM through the s-h mixing), D D -> X Xbar,
 the 4->2 cannibal processes (threshold amplitudes from amp42.py), the core-catalysed conversion of bb211d1
 (quoted sigma v), and object/anti-object collisions. Background: n_Q = Y_Q s(T_vis).
usage: python3 channels.py <path to profile_Q200000_lambda0.03_refined.csv>
"""
import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))  # own code dir (python -I safe)
from common import *

OUT = Path(__file__).resolve().parents[1] / 'results'
ME, MMU, ALPHA = 0.00051099895, 0.1056583755, 1 / 137.035999084


def unit_width_low(M, photon_coefficient=10.):
    """SM-Higgs-like width below the two-pion threshold: e, mu and photons (portal-debris convention)."""
    lep = sum(ml**2 * M / (8 * np.pi * VH**2) * max(1 - 4 * ml**2 / M**2, 0)**1.5 for ml in (ME, MMU))
    return lep + photon_coefficient**2 * ALPHA**2 * M**3 / (256 * np.pi**3 * VH**2)


def main(profile):
    amp = json.loads((OUT / 'amplitudes.json').read_text())
    A = amp['threshold_4to2']['A_GeV_minus8']
    lam3 = 2 * G**2 * V
    res = {}
    # thresholds (cold pair at rest has sqrt(s) = 2 m_X)
    th = dict(
        s_s=dict(sqrt_s_needed=2 * MS, open_cold=bool(2 * MX >= 2 * MS), KE_per_particle_needed=MS - MX),
        s_single=dict(sqrt_s_needed=MS, open_cold=False, KE_per_particle_needed=MS / 2 - MX),
        D_D=dict(sqrt_s_needed=2 * MD, open_cold=bool(MX >= MD), KE_per_particle_needed=MD - MX,
                 note='D is Z2-odd and stable; D D -> X Xbar is open, so D only stores dark energy'),
        Aprime=dict(m_Aprime=M_AP, ratio_to_2mX=M_AP / (2 * MX), open_cold=False),
        muons=dict(sqrt_s_needed=2 * MMU, open_cold=bool(2 * MX >= 2 * MMU)),
        three_to_two=dict(XXXbar_to_X_s=bool(3 * MX >= MX + MS), XXXbar_to_X_D='forbidden by D -> -D',
                          note='charge conservation requires X X Xbar -> X + neutral; no neutral state lighter than 2 m_X'),
        four_to_one_s='forbidden: invariant mass >= 4 m_X != m_s',
        four_to_two=dict(XXXbarXbar_to_XXbar=True, XXXXbar_to_XX=True, XXXbarXbar_to_DD=bool(4 * MX >= 2 * MD),
                         XXXbarXbar_to_ss=bool(4 * MX >= 2 * MS)),
        on_cores=dict(
            conversion_Xbar_core_to_X='open: E_out = m_X + 2 omega_Q = %.4f GeV' % (MX + 2 * OMEGA_REF),
            Xbar_absorption_with_s=bool(MX + OMEGA_REF >= MS), Xbar_absorption_with_s_margin_GeV=MX + OMEGA_REF - MS,
            X_capture_with_s=bool(MX - OMEGA_REF >= MS),
            X_capture_two_body='closed: no hidden quantum lighter than m_X - omega_Q = %.4f GeV' % (MX - OMEGA_REF),
            X_capture_three_body_XX_core='open: E_out = 2 m_X - omega_Q (rate not computed)',
            inverse_conversion_threshold_E='E_X >= m_X + 2 omega_Q (closes after any redshift of the daughter)'),
        bound_states=dict(alpha_eff=lam3**2 / (16 * np.pi * MX**2), criterion_alpha_mX_over_ms=lam3**2 / (16 * np.pi * MX**2) * MX / MS,
                          note='Yukawa bound states need alpha m_X/m_s >~ 1.7; none, and a bound pair could not decay in the hidden sector'))
    res['thresholds'] = th
    # portal annihilation through s-h mixing
    theta = KAPPA * V * VH / (MH**2 - MS**2)
    sqs = 2 * MX
    gam_unit = unit_width_low(sqs)
    sv_portal = lam3**2 * theta**2 * gam_unit / (MX * (sqs**2 - MS**2)**2)
    gam_s_vis = theta**2 * unit_width_low(MS)
    res['portal'] = dict(theta=theta, unit_width_at_2mX=gam_unit, sigma_v_GeV2=sv_portal, sigma_v_cm3s=sv_portal * GEV2_CM3S,
                         s_visible_width_GeV=gam_s_vis,
                         note='mu mu closed for cold pairs; on-shell s -> SM branch ~1e-27 (portal-debris) so the hidden sector cannot drain energy')
    # D D -> X Xbar
    Mdd = (2 * HD**2 * V) * lam3 / (4 * MD**2 - MS**2)
    beta = np.sqrt(1 - MX**2 / MD**2)
    sv_dd = Mdd**2 * beta / (32 * np.pi * MD**2)
    res['DD_to_XXbar'] = dict(M=Mdd, beta=beta, sigma_v_GeV2=sv_dd, sigma_v_cm3s=sv_dd * GEV2_CM3S)
    # object geometry
    a = np.loadtxt(profile, delimiter=',', skiprows=1)
    r, F = a[:, 0], a[:, 1]
    w = r * r * F * F
    c = np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(r)); c = np.concatenate([[0], c]) / c[-1]
    r_half = np.interp(0.5, c, r) / MX; r90 = np.interp(0.9, c, r) / MX
    M_obj = E_PER_FREE_MASS * Q_REF * MX
    res['object'] = dict(r_half_GeV_inv=r_half, r90_GeV_inv=r90, mass_GeV=M_obj)
    # conversion (bb211d1, quoted)
    sv_conv_cold = 1.13039e-16 / GEV2_CM3S
    sv_conv_hot = 8.1225e-15 / GEV2_CM3S
    res['conversion_quoted'] = dict(sigma_v_cold_GeV2=sv_conv_cold, sigma_v_daughter_GeV2=sv_conv_hot,
                                    source='bb211d1 REPORT.md section 2 (cold Q=200000 ground profile, weak-wave)')
    # epochs
    Ts = [86., 40., 20., 10., 7., 5., 3., 1., 0.3, 0.1, 0.03, 0.01, 3e-3, 1e-3, 1e-4, 1e-5, 2.5e-6]
    rows = []
    for T in Ts:
        nQ = YQ * entropy(T)
        H = hubble(T)
        row = dict(T_vis_GeV=T, n_Q_GeV3=nQ, H_GeV=H, nQ_over_mX3=nQ / MX**3)
        # 4->2 per Xbar
        row['G42_asym_over_H'] = A['b'] * nQ**3 / 6 / H
        for p in (0.1, 1., 10., 100.):
            npl, nmi = (1 + p) * nQ, p * nQ
            g = (A['a'] * npl**2 * nmi / 4 + A['b'] * npl**3 / 6 + A['b'] * npl * nmi**2 / 6) / H
            row[f'G42_p{p:g}_over_H'] = g
        # portal per Xbar (symmetric p=10 background)
        row['Gportal_p10_over_H'] = sv_portal * 11 * nQ / H
        # conversion on Q=200000 cores holding all net charge
        ncore = nQ / Q_REF
        row['Gconv_cold_over_H'] = ncore * sv_conv_cold / H
        # elastic X Xbar (kinetic equilibrium) at v ~ 0.3: sigma = |M|^2/(64 pi^2 s) * 4pi
        sig_el = amp['threshold_elastic']['XXbar']**2 / (16 * np.pi * 4 * MX**2)
        row['Gelastic_v0.3_over_H'] = nQ * sig_el * 0.3 / H
        # object/anti-object collisions: velocity needed for Gamma = H with sigma = pi (2 r90)^2
        row['v_needed_obj_collision'] = H / (ncore * np.pi * (2 * r90)**2)
        rows.append(row)
    res['epochs'] = rows
    # temperatures where Gamma = H
    def T_cross(key):
        from scipy.optimize import brentq
        def f(lT):
            T = np.exp(lT); nQ = YQ * entropy(T); H = hubble(T)
            if key == 'asym':
                return np.log(A['b'] * nQ**3 / 6 / H)
            if key.startswith('p'):
                p = float(key[1:]); npl, nmi = (1 + p) * nQ, p * nQ
                return np.log((A['a'] * npl**2 * nmi / 4 + A['b'] * npl**3 / 6 + A['b'] * npl * nmi**2 / 6) / H)
            if key == 'conv':
                return np.log(nQ / Q_REF * sv_conv_cold / H)
        return float(np.exp(brentq(f, np.log(1e-8), np.log(1e3))))
    res['Gamma_equals_H_T_vis_GeV'] = {k: T_cross(k) for k in ('asym', 'p0.1', 'p1', 'p10', 'p100', 'conv')}
    (OUT / 'channels.json').write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: v for k, v in res.items() if k != 'epochs'}, indent=1, default=float))
    for r_ in rows:
        print({k: f'{v:.3g}' for k, v in r_.items()})


if __name__ == '__main__':
    main(sys.argv[1])
