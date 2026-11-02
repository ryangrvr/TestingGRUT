"""Smaller checks of mf_hybrid [D/E/C]: gravitational X production, wall timescales, neutrino-mass sum of the seesaw
singlets, spectator isocurvature with mf_ideal's p, the spectator's one-loop mass from N_s, the induced holder portal
leak, the light-neutrino-mass sensitivity of the Pauli-blocked history.
Usage: python3 -B sk_misc.py -> ../results/sk_misc.json
"""
import math, json, os
import numpy as np
import sk_common as C

mP = C.mP
out = dict(label='D/E/C (skeptic)', note=__doc__)
d, P = C.load_bench('HFa')

# 1. gravitational X pairs from R = Rbar + R1 cos(2 m_chi t), V = (1/2) xi R1 cos(w t) x^2 per real d.o.f.
#    energy rate = A^2 w/(64 pi) per real d.o.f. with w = 2 m_chi (one oscillation quantum -> one pair of energy w)
A = (5 / 6.) * 3 * 0.85 * 3
w = 2 * P['m_chi']
rate_mine = 2 * A ** 2 * w / (64 * math.pi)          # complex X: 2 real d.o.f.
rate_track = d and 101.06494311084468
out['grav_X'] = dict(A_Hi2=A, rate_mine_Hi5=rate_mine, rate_track_Hi5=rate_track, ratio_mine_track=rate_mine / rate_track,
                     fraction_of_M1r_EX=rate_mine / 1.887e5,
                     note='D: Gamma(phi->xx) = sigma^2/(32 pi m) for -sigma phi x^2/2 gives dE/dt = A^2 w/(64 pi) per real dof; '
                          'the track used w = m_chi (x2 low). Still ~1e-3 of E_X')

# 2. domain walls of the real waterfall if they form (lambda/4 (chi^2-v^2)^2, tadpole bias c)
lam, v, c = P['lam'], P['v'], P['cb']
sig_w = (2 * math.sqrt(2) / 3) * math.sqrt(lam) * v ** 3
dV = 2 * c * v
H_dom = sig_w / (3 * mP ** 2)          # rho_w ~ sigma_w H (scaling, A ~ 1) equals 3 mP^2 H^2
H_col = dV / sig_w                     # pressure dV ~ sigma_w H
out['walls_if_formed'] = dict(tension_Hi3=sig_w, bias_dV_Hi4=dV, dV_over_V0=dV / P['V0'], H_domination=H_dom,
                              t_domination_rad=1 / (2 * H_dom), H_collapse=H_col, t_collapse_rad=1 / (2 * H_col),
                              walls_dominate_before_collapse=H_dom > H_col,
                              note='E: scaling-network estimates; if the Pauli-limited waterfall fragments and walls '
                                   'percolate, they dominate at t ~ 1/(2 H_dom) before the bias removes them')

# 3. light-neutrino masses from the seesaw singlets (Casas-Ibarra: sum_j mtilde_j = sum_i m_i sum_j |R_ji|^2 >= sum m_i)
mt = 0.05
nN = 4
out['neutrino_mass_sum'] = dict(mtilde_each_eV=mt, n_singlets=nN, sum_mtilde_eV=nN * mt,
                                sum_m_nu_if_R_real_eV=nN * mt, planck2018_bao_95=0.12, desi_dr2_cmb_95=0.0642,
                                required_sum_j_absR2_for_sum_mnu_0p06=nN * mt / 0.06,
                                note='D: with real R the light masses sum to sum_j mtilde_j = 0.20 eV (0.25 eV with a fifth, '
                                     'seesaw-coupled N_s). Planck 2018+BAO (1807.06209) < 0.12 eV, DESI DR2+CMB '
                                     '(2503.14744) < 0.064 eV (retrospective). Either a complex R with '
                                     'sum_j |R_ji|^2 ~ 3 (cancellations in m_nu, large CP phases) or smaller mtilde '
                                     '(slower N decays) is needed')

# 4. isocurvature: beta = (H/2 pi s_*)^2 p^2 / A_s; Planck 2018 axion-I (uncorrelated, n_II = 1) beta < 0.038
amp = 1 / (2 * math.pi * P['s_star'])
b1 = amp ** 2 / C.AS
p_max = math.sqrt(0.038 / (1 - 0.038) / b1)
D_ideal = [4.42, 4.91, 5.07]              # mf_ideal isocurvature.json, charge-only and total D on fast histories
out['isocurvature'] = dict(beta_per_p2=b1, p_max_planck_0p038=p_max,
                           beta_at_mf_ideal_D=[b1 * x * x for x in D_ideal],
                           s_star_needed_MP_at_mf_ideal_D=[P['s_star'] / mP * x / p_max for x in D_ideal],
                           note='C/Q: Planck 2018 X (1807.06211) beta_iso < 0.038; with mf_ideal p = 4.4-5.1 the bench s_* '
                                '= 3e-3 M_P gives beta = 0.075-0.098 (retrospective conflict); p must come from mf_chain')

# 5. spectator one-loop mass from N_s (mass y_s s, 2 Majorana d.o.f.): V'' = -(y^4 s^2/(32 pi^2))(12 L - 4)
ys, s = 4.6e-3, P['s_star']
base = ys ** 4 * s ** 2 / (32 * math.pi ** 2)
out['spectator_Ns_loop'] = dict(y_s=ys, ms2_over_Hi2_at_L=[dict(L=L, ms2=-base * (12 * L - 4)) for L in (-1, 0, 1. / 3, 1, 2)],
                                drift_lnS_per_efold_for_ms2_0p03=0.03 / 3,
                                note='D: the N_s loop gives |m_s^2| = 0.004 (12L - 4) H^2 at s_*, L = ln(y^2 s^2/mu^2); '
                                     '<= 0.01 H^2 only for |L - 1/3| <~ 0.2 (a renormalisation choice). N_s needs a '
                                     'soft Majorana mass (< 1.5 H_i) and a seesaw Yukawa to reach the SM; otherwise it is '
                                     'dark radiation (Delta N_eff ~ 0.04 at f_s = 1.3%)')

# 6. induced holder portal A_s s|Sigma|^2 ~ c 2|g_h mu| (Y loop), s -> Sigma_radial mixing -> A_L A_L
gh, mu, h, lamS, ms = 2.1887e-4, 7.159456, 40.0, 0.5, 3.0
mrho2 = 2 * lamS * h * h
rows = []
for cc, lab in ((1 / (16 * math.pi ** 2), '1/16pi^2'), (1 / (4 * math.pi ** 2), '1/4pi^2')):
    Aind = cc * 2 * gh * mu
    G = Aind ** 2 * ms ** 3 / (32 * math.pi * (mrho2 - ms * ms) ** 2)
    rows.append(dict(coefficient=lab, A_sSigma_Hi=Aind, Gamma_leak_Hi=G, BR_vs_Gamma_s=G / 1.2629e-6))
out['induced_holder_portal'] = dict(rows=rows, mf_ideal_R_leak='0.36-5.8 x (1.3e-6 H_i/Gamma_s) [Q mf_ideal NOTES 10.6]',
                                    note='E: zero only if the induced portal vanishes at the stated scale (M_Y benchmark); '
                                         'the track says "no hidden heat" without this channel')

# 7. Gaussian wrong-vacuum probability for the linear ratios found in sk_waterfall
try:
    wf = json.load(open(os.path.join(C.OUT, 'sk_waterfall_HFa.json')))
    out['walls_linear'] = {k: dict(ratio_min=v['ratio_amp_min'], t=v['t_ratio_amp_min'], t_ratio_lt2=v.get('t_ratio_lt2'),
                                   frad_t03=v['frad_t03'])
                           for k, v in wf.items() if isinstance(v, dict) and 'ratio_amp_min' in v}
except Exception as e:  # noqa
    out['walls_linear'] = str(e)

# 8. Pauli-blocked history at smaller mtilde (sum m_nu = 0.06 eV with real R: mtilde = 0.015 eV each)
import sk_pauli as SP
SP.MNU = 0.015e-9
s1, _ = SP.run(blocking=True, Gth=0.0)
s2, _ = SP.run(blocking=False)
out['pauli_mtilde_0p015'] = dict(pauli_Gth0_pertfree=s1, unblocked=s2)
print(json.dumps(out, indent=1, default=float)[:6000])
C.dump('sk_misc.json', out)
