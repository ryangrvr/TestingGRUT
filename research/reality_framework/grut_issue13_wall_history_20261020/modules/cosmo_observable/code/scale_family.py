"""Part 4: scale dependence of the persistent-network consequence at fixed couplings (lambda_phi = 0.01, g_X = 0.1),
sigma = (2 sqrt(2 lambda_phi)/3) v^3. Writes results/scale_family.json.

v_max(dom): the largest v for which a persistent network in the scaling regime (rho_w = A sigma / t, backreaction
included) has not reached rho_w = rho_others by T = T_CMB. Variants: rho_Lambda = 0 or Planck-like; A in
{0.643 (LMS matter era), 0.7, 0.8, 0.9}; freezing (rho_w ~ a^-1) after the onset of acceleration.
"""
import json
import math
import os
import time
import numpy as np
from scipy.optimize import brentq
import core as C
import observables as O

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')

def log_ratio_today(lnsig, A, rho_L, freeze):
    sig = math.exp(lnsig)
    if freeze:
        h = C.history(sig, A=A, n_post=1.0, f_switch=1e30, rho_L=rho_L, switch_on='all', T_start=10.0,
                      freeze_at_acceleration=True)
    else:
        h = C.history(sig, A=A, n_post=1.0, f_switch=1e30, rho_L=rho_L, switch_on='all', T_start=10.0)
    rw = C.wall_density_at(h, C.T0) if h['switched'] else h['rho_w_end']
    other = C.RHO_R0 + C.RHO_M0 + rho_L
    return math.log(rw / other), h

def sigma_dom_today(A, rho_L, freeze=False):
    f = lambda ls: log_ratio_today(ls, A, rho_L, freeze)[0]
    ls = brentq(f, math.log(1e-7), math.log(1e-2), xtol=1e-10)
    return math.exp(ls)

def main():
    t0 = time.time()
    out = {'module': 'cosmo_observable', 'part': '4 scale family', 'lambda_phi': C.LAMBDA_PHI, 'g_X': C.G_X}
    # ---- v_max(dom) -----------------------------------------------------------------------------------------------
    vd = {}
    for A in (C.A_MATTER_LMS, C.A_LOW, C.A_CENTRAL, C.A_HIGH):
        for rl_name, rl in (('rhoL=planck_like', C.RHO_L_PLANCK), ('rhoL=0', 0.0)):
            if rl_name == 'rhoL=0' and A != C.A_CENTRAL:
                continue
            s = sigma_dom_today(A, rl)
            v = C.v_from_sigma(s)
            vd['A=%.3f_%s' % (A, rl_name)] = {'sigma_max_GeV3': s, 'cuberoot_MeV': s**(1/3.) * 1e3,
                                              'v_max_MeV': v * 1e3, 'm_X_max_MeV': C.G_X * v * 1e3,
                                              'r_max_product_bound_63.93': 63.93 / (C.G_X * v)}
    s = sigma_dom_today(C.A_CENTRAL, C.RHO_L_PLANCK, freeze=True)
    v = C.v_from_sigma(s)
    vd['A=0.800_rhoL=planck_like_freeze_after_acceleration_onset'] = {
        'sigma_max_GeV3': s, 'cuberoot_MeV': s**(1/3.) * 1e3, 'v_max_MeV': v * 1e3, 'm_X_max_MeV': C.G_X * v * 1e3,
        'r_max_product_bound_63.93': 63.93 / (C.G_X * v)}
    # analytic comparison: subdominant scaling today, rho_w0 = rho_c0  -> sigma = 3 M_P^2 H0^2 t0 / A
    H0 = C.kmsmpc_to_gev(67.36)
    t0p = C.AGE_PLANCK_GYR * C.GYR_S / C.HBAR_GEV_S
    vd['analytic_rho_w0_eq_rho_c0_A0.8'] = {'sigma_GeV3': 3 * C.M_P**2 * H0**2 * t0p / C.A_CENTRAL,
                                            'note': 'no backreaction, LCDM t0 = 13.797 Gyr, H0 = 67.36; order-of-magnitude check'}
    out['v_max_dom'] = vd
    vals = [e['v_max_MeV'] for k, e in vd.items() if 'v_max_MeV' in e]
    out['v_max_dom_range_MeV'] = [min(vals), max(vals)]

    # ---- v_max(CMB) from the reconciled LMS bound (subdominant-network bound) ---------------------------------------
    with open(os.path.join(RES, 'normalization.json')) as f:
        nrm = json.load(f)
    out['v_max_CMB'] = {k: {kk: e[kk] for kk in ('sigma_GeV3', 'v_max_MeV', 'm_X_max_MeV_gX0.1', 'r_max_from_mX_r<63.93GeV')}
                        for k, e in nrm['bounds_table'].items()}

    # ---- scan over v: T_dom and R_H for n = 1, 2 (rho_L = 0) --------------------------------------------------------
    scan = []
    for lv in np.linspace(math.log10(3e-3), math.log10(3e4), 36):
        v = 10**lv
        sig = C.sigma_tree(v)
        row = {'v_GeV': v, 'sigma_GeV3': sig, 'm_X_GeV': C.G_X * v}
        for n in (1.0, 2.0):
            Tst = max(100.0, 10 * v) if v > 10 else 100.0
            h = C.history(sig, A=C.A_CENTRAL, n_post=n, f_switch=1.0, rho_L=0.0, T_start=Tst, switch_on='all')
            s = O.summarize(h, with_accel=False)
            row['n=%d' % n] = {'T_dom_GeV': s['T_switch_GeV'], 'z_dom': s['z_switch'],
                               'H0_pred_kmsMpc': s['H0_pred_kmsMpc'], 'R_H_vs70': s['R_H']['70'],
                               'Omega_w_T_rec': s['fractions']['T=%.4g_GeV' % O.T_REC]['rho_w/rho_tot'],
                               'DeltaNeff_equiv_T1MeV': s['fractions']['T=0.001_GeV']['DeltaNeff_equiv']}
        scan.append(row)
    out['scan_rhoL0_A0.8_fswitch1'] = scan
    out['scan_note'] = ('T_dom here is where rho_w = rho_rad + rho_m (switch_on=all): scaling is assumed until walls '
                        'dominate the total, then rho_w ~ a^-n. For v >~ 1 GeV this is in the radiation era and agrees '
                        'with rho_w = rho_rad to <1e-3 at v = 1 TeV.')

    # ---- characteristic v: domination at T = 1 MeV (before BBN), at T_rec, at T_eq-like ------------------------------
    def T_dom_of_v(lv):
        v = 10**lv
        h = C.history(C.sigma_tree(v), A=C.A_CENTRAL, n_post=1.0, f_switch=1.0, rho_L=0.0,
                      T_start=max(100.0, 10 * v), switch_on='all')
        return h['T_switch_GeV'] if h['switched'] else 0.0
    char = {}
    for name, Tt in (('T_dom=1MeV', 1e-3), ('T_dom=0.1MeV', 1e-4), ('T_dom=T_rec', O.T_REC), ('T_dom=T_CMB_today_rhoL0', C.T0 * 1.000001)):
        lv = brentq(lambda x: math.log(max(T_dom_of_v(x), C.T0 * 0.5) / Tt), -3.0, 5.0, xtol=1e-9)
        char[name] = {'v_GeV': 10**lv, 'sigma_GeV3': C.sigma_tree(10**lv)}
    # v at which the persistent network alone (rho_L = 0) would give H0 = 70 km/s/Mpc (R_H = 1), n = 1 and n = 2
    def RH_of_v(lv, n):
        v = 10**lv
        h = C.history(C.sigma_tree(v), A=C.A_CENTRAL, n_post=n, f_switch=1.0, rho_L=0.0, T_start=100.0, switch_on='all')
        return O.summarize(h, with_accel=False)['R_H']['70']
    for n in (1.0, 2.0):
        lv = brentq(lambda x: math.log(RH_of_v(x, n)), math.log10(0.03), math.log10(0.5), xtol=1e-9)
        char['R_H=1_vs70_rhoL0_n=%d' % n] = {'v_GeV': 10**lv, 'sigma_GeV3': C.sigma_tree(10**lv),
                                             'sigma_over_reconciled_LMS_bound': C.sigma_tree(10**lv) / 8.409393636408877e-10}
    out['characteristic_v'] = char
    out['T_dom_scaling_note'] = 'radiation era, no backreaction: T_dom ~ sigma^(1/2) g_*^(-1/4) ~ v^(3/2)'
    out['runtime_s'] = time.time() - t0
    with open(os.path.join(RES, 'scale_family.json'), 'w') as f:
        json.dump(out, f, indent=1)
    return out

if __name__ == '__main__':
    o = main()
    print(json.dumps({k: o[k] for k in ('v_max_dom', 'v_max_dom_range_MeV', 'characteristic_v', 'runtime_s')}, indent=1))
    for r in o['scan_rhoL0_A0.8_fswitch1']:
        print('v=%.3g GeV  Tdom=%s  R_H(n=1)=%.3g  R_H(n=2)=%.3g  Om_w(rec)=%.3g' % (
            r['v_GeV'], r['n=1']['T_dom_GeV'], r['n=1']['R_H_vs70'], r['n=2']['R_H_vs70'], r['n=1']['Omega_w_T_rec']))
