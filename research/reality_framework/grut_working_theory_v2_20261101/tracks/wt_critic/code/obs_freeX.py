"""wt_critic (Claude, own code): what the free-X relic of M1r shows an observer.  Arithmetic [D] on supplied or
track numbers; literature thresholds are [R].

- single-X nucleon cross section through phi-h mixing theta (trilinear A phi|H|^2 or b phi^2|H|^2) and through
  lambda_HX |X|^2|H|^2:  sigma = g^2 y_N^2 mu^2/(pi m_phi^4) with y_N = theta f_N m_N / v_H (complex scalar X,
  X X phi vertex 2 g m_X from g^2 phi^2 |X|^2), and sigma = lambda^2 (f_N m_N/m_h^2)^2 mu^2/(4 pi m_X^2);
- X X -> X X self-interaction sigma/m from REVIEW_M1's sigma = 2.16627561e-7/m_X^2;
- fugacity of the all-gas relic z = exp(5/2 - S/Q) (non-relativistic Boltzmann limit) and the condensation
  temperature T_s where mu_gas = m + T ln z reaches the bulk Omega_b = 0.775865;
- free-streaming mass for the quoted comoving lambda_fs.
"""
import json, math, sys, os
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'obs_freeX.json')
GEV2_TO_CM2 = 0.3893794e-27
GEV_TO_G = 1.78266192e-24
mX, g, mN, fN, vH, mh = 0.0924089637375362, 0.1, 0.938919, 0.30, 246.22, 125.25
mphi2 = 6 * mX ** 2
mu = mX * mN / (mX + mN)
yN_per_theta = fN * mN / vH
sig_theta2 = g ** 2 * yN_per_theta ** 2 * mu ** 2 / (math.pi * mphi2 ** 2) * GEV2_TO_CM2     # cm^2 per theta^2
sig_lHX2 = (fN * mN / mh ** 2) ** 2 * mu ** 2 / (4 * math.pi * mX ** 2) * GEV2_TO_CM2       # cm^2 per lambda^2
theta_A = lambda A: A * vH / mh ** 2
theta_b = lambda b: 2 * b * (10 * mX) * vH / mh ** 2
cases = {
    'trilinear_A_2.7e-10GeV (fx_observables XLZD value for Q=2e5 objects)': sig_theta2 * theta_A(2.7e-10) ** 2,
    'b_2.2e-10 (fx_observables XLZD value for Q=2e5 objects)': sig_theta2 * theta_b(2.2e-10) ** 2,
    'lambda_HX_9e-11': sig_lHX2 * 9e-11 ** 2,
    'lambda_HX_1e-3 (would thermalise the sectors)': sig_lHX2 * 1e-6,
    'theta_1e-3': sig_theta2 * 1e-6,
}
sidm = 2.16627561e-7 / mX ** 3 * GEV2_TO_CM2 / GEV_TO_G
rho_m = 0.3111 * 2.775e11 * 0.6766 ** 2      # M_sun / Mpc^3 (Planck-like) [R]
out = dict(sigma_Xn_per_theta2_cm2=sig_theta2, sigma_Xn_per_lambdaHX2_cm2=sig_lHX2,
           sigma_Xn_cm2=cases, theta_per_A_GeVinv=vH / mh ** 2, theta_per_b=2 * 10 * mX * vH / mh ** 2,
           sidm_sigma_over_m_cm2_per_g=sidm)
fug = {}
for SQ in (3.47, 4.68, 4.94, 5.13):
    z = math.exp(2.5 - SQ)
    fug['%.2f' % SQ] = dict(z=z, T_s_over_mX=(1 - 0.775865) / (-math.log(z)))
out['fugacity_all_gas'] = fug
out['M_fs_Msun'] = {('%.1f pc' % l): 4 * math.pi / 3 * rho_m * (l * 1e-6) ** 3 for l in (0.3, 0.5)}
out['closed_channels_at_rest'] = dict(XXbar_to_phi=dict(need=math.sqrt(6), have=2.0), XXbar_to_phiphi=dict(need=2 * math.sqrt(6), have=2.0),
                                      XXbarX_to_Xphi=dict(need=1 + math.sqrt(6), have=3.0))
json.dump(out, open(OUT, 'w'), indent=1)
print(json.dumps(out, indent=1))
