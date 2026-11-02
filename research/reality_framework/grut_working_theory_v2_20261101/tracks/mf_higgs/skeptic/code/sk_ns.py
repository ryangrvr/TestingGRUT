#!/usr/bin/env python3
"""mf_higgs SKEPTIC: the spectator neutrino N_s (not in the bench; mf_hybrid skeptic N8), its loop-induced Higgs mass
during inflation, as a function of the seesaw m-tilde_s that the track fixed at 0.05 eV.  [D/E, model choices]
  delta m_h^2(s) = -(y_nus^2/8 pi^2) M(s)^2 [ln(M(s)^2/mu^2) - 1],  M(s) = M_soft + y_s s,  y_nus^2 = mt_s M_soft/v^2.
  * the largest mt_s that keeps delta m_h^2 >= -0.66 H_i^2 (the inflationary tolerance at xi_EW = 0.1, M_t = 172.86,
    tree constant-H HM; reproduced in sk_infl) under 'vanish at mu = M_soft';
  * the s-odd piece (2 M_soft y_s s): the CP-branch difference of the inflationary Higgs mass;
  * the smallest mt_s that makes N_s decay before it dominates (N_s from s decay at a ~ 1.8e3 with 1.3% of rho,
    matter-like afterwards: Gamma_Ns > H at a_dom = a_dec / 0.013, H ~ K H_i/a^2) -- the seesaw is then not forced;
  * Sigma m_nu with five singlets at 0.05 eV each (real Casas-Ibarra R): 0.25 eV.
Output: results/sk_ns.json
"""
import sys
sys.dont_write_bytecode = True
import math
import sk_common as C

H, VEW = 1.30002426e11, 174.1
ys, sstar = 4.6e-3, 0.003 * C.MPL / 1.30002426e11
Ms_inf = ys * sstar
TOL = 0.66
out = dict(M_s_inflation_over_Hi=Ms_inf, tolerance_Hi2=TOL, rows={})
for Msoft in (0.3, 1.0, 1.5):
    def dm2(mt, mu=Msoft, M=Ms_inf + Msoft):
        y2 = mt * 1e-9 * Msoft * H / VEW ** 2
        c = y2 * M ** 2 / (8 * math.pi ** 2)
        return -c * (math.log(M ** 2 / mu ** 2) - 1), c
    v05, c05 = dm2(0.05)
    mt_max = 0.05 * TOL / (-v05)
    vp, _ = dm2(0.05, M=Ms_inf + Msoft)
    vm, _ = dm2(0.05, M=Ms_inf - Msoft)
    # decay before N_s domination: Gamma = y^2 M_soft/(8 pi) (both L H and conj.), need >= H(a_dom)
    a_dec, f_dec, K = 1.8e3, 0.013, 2.05
    a_dom = a_dec / f_dec
    H_dom = K / a_dom ** 2
    y2_min = 8 * math.pi * H_dom / Msoft
    mt_min = y2_min * VEW ** 2 / (Msoft * H) * 1e9
    out['rows'][f'M_soft={Msoft}'] = dict(dm2_at_0p05eV_vanish_at_Msoft=v05, coeff_c=c05,
                                         mt_s_max_for_tolerance_eV=mt_max,
                                         CP_branch_difference_Hi2=abs(vp - vm),
                                         mt_s_min_to_decay_before_domination_eV=mt_min)
    print(Msoft, out['rows'][f'M_soft={Msoft}'])
out['Sigma_m_nu_5_singlets_0p05eV'] = 0.25
out['bounds_quoted'] = 'Planck 2018 + BAO: Sigma m_nu < 0.12 eV (1807.06209); DESI DR2 + CMB < 0.064 eV (2503.14744) [Q via mf_hybrid skeptic]'
out['dark_radiation_alternative'] = 'no seesaw: N_s massless/stable -> Delta N_eff ~ 0.04 (mf_hybrid skeptic), inside Planck N_eff = 2.99 +- 0.17'
C.jdump(out, C.RES / 'sk_ns.json')
