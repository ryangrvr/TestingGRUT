"""Verifier v_vacua: collect the verdict table from results/*.json -> ../results/verdicts.json"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, '..', 'results')
v1 = json.load(open(os.path.join(R, 'v1_qsector.json')))
v2 = json.load(open(os.path.join(R, 'v2_nucleation.json')))
v3 = json.load(open(os.path.join(R, 'v3b_sm_check.json')))
B = v3['B_buttazzo_matching_smdr_running']
C = v3['C_smdr_matching_smdr_running']
D = v3['D_tension_heldout']
F = v3['F_decomposition_eff']
row1d = {(r['Kprime'], r['x']): r for r in v1['1d_numeric_KV']}
walls = {w['m_w']: w for w in v1['1f_linear_wall_sech2']}
lc = v2['LCDM'][0]
rt = v2['rates']['rows']['1e-577']
pot = B['potential_PDG_central']

V = []


def add(n, claim, track, mine, verdict, note=''):
    V.append({'#': n, 'claim (track section)': claim, 'track': track, 'verifier': mine, 'verdict': verdict,
              'note': note})


add(1, 'E1/E2 static equations, signs (1.1)', 'eps\' - 2K\'Lam/K = mu; Lam = V + eps - mu q',
    'same; G = -3H^2 g checked: %s' % v1['1a_dS']['G_plus_3H2_g_is_zero'], 'CONFIRMED')
add(2, 'a2 = 2K\'^2/(K^2 eps\'\') = 4/(3 K m_q^2), a3, dq/dx, dK/K (1.3)', 'as stated',
    'a2 %s, a3 %s, dK/K = a2 %s' % (v1['1b_series']['a2_equals_track'], v1['1b_series']['a3_equals_track'],
                                    v1['1b_series']['dK/K equals a2']), 'CONFIRMED',
    'H^2 = Lam/(3K(q)) = x/(3K0) + O(x^3): no x^2 term in the expansion rate (as in equations 5.2)')
add(3, 'exact slope dLam/dV (1.3)', 'K^2 D/(K^2 D - 4 Lam K\'^2)',
    'equal on shell (after E1): %s' % v1['1c_slope']['equals_track_formula_on_shell'], 'CONFIRMED',
    'holds only on the static family (uses E1); the off-shell implicit derivative has an extra 2(mu-eps\')KK\' term')
add(4, 'numeric static roots, KV eps (1.3)', '0.010190 (K\'=1, x=0.01); 0.0100177 (K\'=0.3)',
    '%.6f; %.7f' % (row1d[(1.0, 0.01)]['Lam_exact'], row1d[(0.3, 0.01)]['Lam_exact']), 'CONFIRMED')
add(5, 'a2 x at Planck depth; "<0.5% down to |V| ~ 1e71" (0, 1.3)', '4e-3; <0.5% to 1e71',
    '%.2e at 7.5e70; %.2e at 1e71; 0.5%% reached at %.2e' % (
        v1['1g_sizes']['a2 x at x=7.5e70 (track: ~4e-3)'],
        v1['1g_sizes']['a2 x at x=1e71 (track says <0.5% down to 1e71)'], v1['1g_sizes']['x at which a2 x = 0.5%']),
    'CORRECTED', '0.57% at 1e71 (trivial); conclusion unchanged')
add(6, 'q-coupled matter: Lam_i = V_i - eps*(mu - theta_i); difference formula (1.4)', 'as stated',
    'sympy: residual %s; first-order term %s' % (v1['1e_coupled']['Lam_i - [V_i - eps*(mu - theta_i)]'],
                                                 v1['1e_coupled']['first-order term check (should be 0)']), 'CONFIRMED')
add(7, 'junction: [Kx]=0, [n.dK]=S/3=-sigma, -K[Kx_ab]=S_ab-h_ab S/3 (1.5)', 'as stated',
    'distributional parts of G_ab and R recomputed: %s' % v1['1f_junction']['distributional_parts_confirmed'],
    'CONFIRMED')
add(8, 'linear wall: kink -0.994/-0.997 sigma; tracking deviation 1.8% (m w=10), 0.12% (m w=30) (1.5)',
    '-0.994, -0.997; 1.8%; 0.12%',
    '%.4f, %.4f; %.2f%%; %.3f%% (analytic 2/(mw)^2 = 0.222%%)' % (
        walls[0.01]['kink_corrected_exp(m zo)'], walls[0.1]['kink_corrected_exp(m zo)'],
        100 * walls[10.0]['max_rel_dev_from_local_tracking(|z|<2w)'],
        100 * walls[30.0]['max_rel_dev_from_local_tracking(|z|<2w)']),
    'CORRECTED', 'the m w = 30 value is 0.22%, not 0.12% (track numerics); scaling 2/(m w)^2 confirmed')
add(9, 'V4: Minkowski pi t^4/3; MD 81pi/1485 t^4; LCDM 0.1230/H0^4, eta0 = 3.18/H0, t0 = 13.79 Gyr (2.1)',
    '0.1230; 3.18; 13.79', '%.4f; %.3f; %.2f' % (lc['V4_H0^4'], lc['eta0_H0'], lc['t0_Gyr']), 'CONFIRMED',
    'Buttazzo\'s 0.15 corresponds to Omega_m ~ 0.27 (0.143); the track flags the difference correctly')
add(10, 'volume criterion Gamma/H^4 < 9/(4pi); dS slope (4pi/3)/H^3 (2.2)', '0.716',
    '%s = %.4f' % (v2['volume_criterion']['threshold'], v2['volume_criterion']['numeric']), 'CONFIRMED')
add(11, 'rate table 1e-577 Gyr^-1 Gpc^-3 (2.3)', '-741.3; -574.8; -573.6; 1e574 Gyr; 1e145 Gyr',
    '%.1f; %.1f; %.1f; 1e%.0f; 1e%.0f' % (rt['log10_Gamma_GeV4'], rt['log10_P_past_cone'], rt['log10_Gamma/H_L^4'],
                                          rt['log10_t_dS_worldline_Gyr'], rt['log10_t_Minkowski_Gyr']), 'CONFIRMED',
    'rates quoted from CMS 1707.09301 abstract (checked at abstract level)')
add(12, 'AdS crunch t_c = pi sqrt(3K/|Lam|) (2.4)', '3e-41 s; 2e-38 s; 55 Gyr',
    '%.2e s; %.2e s; %.1f Gyr' % (v2['AdS_crunch']['times']['7.5e70']['t_c_s'],
                                  v2['AdS_crunch']['times']['2e65']['t_c_s'],
                                  v2['AdS_crunch']['times']['rho_Lambda']['t_c_Gyr']), 'CONFIRMED')
cd = {r['dMt']: r for r in v3['E_CDL']}
add(13, 'CDL quench moves the stability boundary by < 1e-3 GeV (2.5)', 'eps/thr = 5934 at +0.001 GeV',
    'eps/thr = %.0f at +0.001; crossing eps = thr at dMt ~ 3e-6 GeV (%.2f at 1e-6, %.2f at 1e-5)' % (
        cd[0.001]['eps/thr(pos)'], cd[1e-06]['eps/thr(pos)'], cd[1e-05]['eps/thr(pos)']), 'CONFIRMED',
    'conclusion holds; ratios quoted at >= 0.01 GeV depend on the tension definition (full-path sigma gives 0.17 '
    'at 0.05 GeV) because thin wall fails there')
add(14, 'track 1- to 3-loop beta functions and gamma_h (transcription from 1307.3536 App. A/B)', 'own transcription',
    'equal to SMDR to 1e-13 (1-2 loop), <2e-5 rel. (3-loop, rounded coefficients), 4-loop QCD 1.5e-6', 'CONFIRMED')
bp = B['ref_point_at_MPl_3loop']
add(15, 'validation at 1307.3536 inputs: lambda(M_Pl), g_i, y_t (5.3)', '-0.01437; 0.61537 0.50546 0.48731 0.38248',
    '%.5f; %.5f %.5f %.5f %.5f' % (bp['lam'], bp['g1'], bp['g2'], bp['g3'], bp['yt']), 'CONFIRMED')
add(16, 'M_t,crit at PDG 2025 point, msbar criterion (5.4)', '170.94, mu* = 4.7e17',
    'Buttazzo matching + SMDR running: %.3f (mu* %.2e); SMDR full matching: %.3f (3-loop), %.3f (loop order 5)' % (
        B['PDG2025']['Mt_crit'], B['PDG2025']['mu_star'], C['Mtcrit_PDG2025']['loopOrder3'][0],
        C['Mtcrit_PDG2025']['loopOrder5'][0]),
    'CORRECTED', 'pipeline reproduced to 2 MeV; with current matching (4-loop QCD y_t, correct g3) the line is '
                 '171.03 +- 0.11, = Hiller et al. 171.04 (the track listed 171.04 but did not resolve the 0.1 GeV gap)')
add(17, 'slopes 0.495 (M_h), 0.353 per 0.0009 (alpha_s) (5.4)', '0.495; 0.353',
    'B: %.3f; %.3f. S: %.3f; %.3f' % (B['PDG2025']['dMtcrit/dMh'], B['PDG2025']['dMtcrit per 0.0009 as'],
                                      C['dMtcrit/dMh'], C['dMtcrit per 0.0009 as']), 'CONFIRMED')
add(18, 'theory shifts: dlam +0.072, dyt -0.087, 2-3 loop -0.109; 4-loop "several times smaller" (5.4)',
    '+0.072, -0.087, -0.109', '%+.3f, %+.3f, %+.3f; loop order 4/5 minus 3: %+.4f/%+.4f' % (
        B['shifts']['dlam +0.0003'], B['shifts']['dyt +0.0005'], B['shifts']['RGE 2-loop minus 3-loop'],
        B['shifts']['RGE loopOrder4 minus 3'], B['shifts']['RGE loopOrder5 minus 3']), 'CONFIRMED')
add(19, 'matching error budget 0.113 GeV covers the matching (5.4, 7)', '0.113 (Buttazzo +-0.0003, +-0.0005)',
    'SMDR - Buttazzo couplings at Mt: y_t %.5f, g3 %.4f, lambda %+.5f -> M_t,crit %+.3f (y_t %+.3f, g3 %+.3f)' % (
        C['couplings_at_Mt: SMDR minus Buttazzo-fit']['(172.56, 125.2, 0.118)']['yt'],
        C['couplings_at_Mt: SMDR minus Buttazzo-fit']['(172.56, 125.2, 0.118)']['g3'],
        C['couplings_at_Mt: SMDR minus Buttazzo-fit']['(172.56, 125.2, 0.118)']['k'],
        D['matching difference S - B (GeV)'], F['B with SMDR yt'] - F['B (all Buttazzo)'],
        F['B with SMDR g3'] - F['B (all Buttazzo)']),
    'CORRECTED', 'net offset 0.088 GeV is inside 0.113, but the y_t shift (1.9x its quoted error, missing 4-loop QCD '
                 'that 1307.3536 itself mentions) and a g3(Mt) 0.0017 above the 4-loop QCD value are outside the budget')
add(20, 'inverse: M_h,crit = 128.54 at M_t = 172.56 (5.4)', '128.54',
    'B: %.2f; S: %.2f' % (B['Mh_crit_at_Mt=172.56'], C['Mh_crit_at_Mt=172.56']), 'CORRECTED', '128.36 with SMDR matching')
add(21, 'FN double condition (M_t, M_h) at mu* = 1.22e19 / 2.435e18 / 1e18 (5.4)',
    '174.20/132.03; 172.61/128.64; 171.71/126.77',
    '%.2f/%.2f; %.2f/%.2f; %.2f/%.2f' % tuple(B['FN_double']['1.2209e19'] + B['FN_double']['2.435e18']
                                             + B['FN_double']['1e18']), 'CONFIRMED',
    'with Buttazzo matching (same caveat on matching as #16)')
add(22, 'potential at PDG central: V(2.435e18) = -7.5e70 GeV^4 = 3e117 rho_L; lambda = 0 at 4.2e10 GeV (5.4, 3.1)',
    '-7.48e70; 3e117; 4.2e10',
    '%.3g (1-loop lambda_eff); %.2g; %.2g' % (pot['2.435e+18']['V_eff_GeV4'], pot['|V(2.435e18)|/rho_Lambda'],
                                               pot['lambda zero crossing GeV']), 'CORRECTED',
    'depth confirmed; the lambda = 0 crossing is 4.7e10 (track picks a coarse grid point; 6.2e9 not 5.7e9 at the '
    '1307.3536 point)')
tb = D['B (Buttazzo matching)']
ts = D['S (SMDR matching), same error budget']
add(23, 'retrospective tension z = 3.3 (direct), 2.8 (+0.3 MC), 1.8 (xsec) (5.5)', '3.34; 2.84; 1.84',
    'B: %.2f; %.2f; %.2f. S: %.2f; %.2f; %.2f' % (
        tb['direct 172.56+-0.31 | exp+param+matching']['z'], tb['direct 172.56+-0.31 | +criterion+0.3 MC']['z'],
        tb['xsec 172.4+-0.7 | exp+param+matching']['z'], ts['direct 172.56+-0.31 | exp+param+matching']['z'],
        ts['direct 172.56+-0.31 | +criterion+0.3 MC']['z'], ts['xsec 172.4+-0.7 | exp+param+matching']['z']),
    'CORRECTED', '3.1 / 2.7 / 1.7 sigma with current matching; still a retrospective tension, not a decisive failure')
add(24, 'held out: FCC-ee sigma 0.13 GeV, z ~ 12; HL-LHC z ~ 5.6 (5.6)', '0.1325; 12.25; 5.57',
    'B: %.4f; %.2f; %.2f. S: %.2f; %.2f' % (D['heldout FCC-ee (track inputs), B']['sigma'],
                                            D['heldout FCC-ee (track inputs), B']['z(direct)'],
                                            D['heldout HL-LHC (0.2, 0.00045, 0.02), B']['z(direct)'],
                                            D['heldout FCC-ee, S']['z(direct)'], D['heldout HL-LHC, S']['z(direct)']),
    'CORRECTED', 'arithmetic confirmed; with current matching 11.6 and 5.3; future precisions themselves NOT CHECKED')
add(25, 'eff criterion: M_t,crit 171.00, eff - msbar = +0.064 (5.3, 5.4)', '171.00; +0.064',
    '1-loop lambda_eff: %.3f; %+.3f' % (F['eff_1loop_B'], F['eff_1loop_minus_msbar_B']), 'CONFIRMED',
    'to 0.013 GeV with a 1-loop lambda_eff; the 2-loop piece was not independently coded')
add(26, 'cross-check with literature: Bednyakov 171.44 at (125.7, 0.1185); 1307.3536 own formula 171.12 vs quoted 171.53',
    '171.38; inconsistency 0.4 GeV in source',
    'B: %.2f; S: %.2f (loop 5 %.2f); formula: %.3f' % (
        B['Mtcrit_Bednyakov_inputs(125.7,0.1185) msbar'][0], C['Mtcrit_Bednyakov(125.7,0.1185) loop3/loop5'][0][0],
        C['Mtcrit_Bednyakov(125.7,0.1185) loop3/loop5'][1][0], B['Buttazzo linear formula at Mh=125.15']),
    'CONFIRMED', 'SMDR matching reproduces Bednyakov to 0.03 and Hiller (171.04; alpha_s 0.1215; M_h 127.97) to '
                 '0.01-0.05')
add(27, 'PDG 2025 inputs (5.1)', '125.20+-0.11; 172.56+-0.31; 172.4+-0.7; 0.1180+-0.0009; 80.3692+-0.0133',
    'read from pdg-2025-v0.2.3.sqlite and QCD review eq. (9.25): identical', 'CONFIRMED')
add(28, 'readings (a)-(d) and the fixed-point argument (3)', 'see track 3.5',
    'logic re-derived; (c)-worldline -> (a); (d) not implied by the action; Volovik counting', 'CONFIRMED',
    '(b) "only if Gamma = 0" is an interpretive [A] step (it imposes (c)\'s end-point condition on (b)); under (a) '
    'with uncalibrated QV1, matching rho_L,obs > 0 needs exact criticality (splitting = rho_L), i.e. the equality')
add(29, 'novelty and credit (4, 6.4)', 'FN, Degrassi/Buttazzo/Bednyakov, Klinkhamer, AHHKM, KV, Volovik',
    'missing: Shaposhnikov-Wetterich 0912.0208; Holthausen-Lim-Lindner 1112.2415; Bezrukov et al. 1205.2893; '
    'Kawai-Okada 1110.2303 / Hamada-Kawai-Kawana 1509.05955 (far-future selection of Lambda, MPP); '
    'Martin-Robertson 1907.02500 / Alam et al. (matching behind Hiller 171.04)', 'CORRECTED',
    'credit additions only; the bookkeeping claim stands')
json.dump(V, open(os.path.join(R, 'verdicts.json'), 'w'), indent=1)
for v in V:
    print(v['#'], v['verdict'], '|', v['claim (track section)'], '|', v['track'], '|', v['verifier'])
