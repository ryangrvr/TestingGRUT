#!/usr/bin/env python3
"""m1_ledger_skeptic (Claude, own code, algebra only, no engine runs): adversarial recomputation of the m1_ledger
track's numbers from Codex's raw packet inputs, read through git show.  Everything here is mpmath (50 digits) or
plain arithmetic on quoted [Q] numbers; no Codex code is imported, no reheating engine is run.

Checks: (1) frozen-test coefficients from (v, g, M, Q, D); (2) the old target decomposition (object + released pairs);
(3) Codex's own yield formula re-evaluated from its row numbers (g* = 106.75) -- a convention check that does not use
any engine; (4) the H_i^{3/2} exponent from Codex's scale_0.5 / scale_2.0 control rows (factor 4 in scale, not used
by m1_ledger); (5) the lambda chain; (6) the clock/window ambiguity of lambda from Codex's charge_clock entries;
(7) the D = 4.4e-10 GeV target against the Planck-2018 value; (8) A from the M1 heavy parameters and the hw_gates
formula alpha6 = kappa^2 mu <Sigma>/(M_Y^4 - 4 mu^2 <Sigma>^2); (9) the Gamma t0 chain from Codex's quoted channel
rates (old and new); (10) the pulse rows of m1_ledger re-read with the charge-conserving (eta > 1) bookkeeping;
(11) the cold-completion and 10%-kinetics rows.
Output: ../results/skeptic_algebra.json
"""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps = 50
OUT = Path(__file__).resolve().parents[1] / 'results'
R = {}

# ----------------------------------------------------------------------------------------------------- [Q] inputs
v = mp.mpf('0.924089637375362')                 # d4747cb CANDIDATE_M1.json vPhi_GeV
gX = mp.mpf('0.1')
mX_codex = mp.mpf('0.09240896373753621')        # CANDIDATE_M1.json mX_GeV
D = mp.mpf('4.4e-10')                           # GeV per unit entropy, d4747cb/1fcc6b7
Q = 200000
M_new = mp.mpf('17324.832346421026')            # DEFORMATION_RESULTS primary mass_GeV
M_old = mp.mpf('17064.537454385838')            # 1fcc6b7 ABUNDANCE_LEDGER object_mass_GeV
YQ_old = mp.mpf('4.8294009076988175e-9')
late_mass_old = mp.mpf('18221.721841256607')    # d71cc01 FIXED_YIELD_RESULTS Q200000 late_mass_per_object_GeV
late_X_old = mp.mpf('12521.802854958796')       # late_X_particles_per_object
scalars_old = mp.mpf('6261.211326289298')       # source_scalar_number_per_object
late_D_old = mp.mpf('0.6197976197995704')
codex = dict(YQ_new=mp.mpf('5.079414232725842e-09'), c_eta=mp.mpf('0.06678046736338095'), c_p=mp.mpf('2.133560934726762'),
             eta_min=mp.mpf('0.8502556152297386'), p_max=mp.mpf('0.004686999952631147'),
             old_fraction=mp.mpf('0.9364942349053573'), old_mult=mp.mpf('1.0678122328228328'),
             old_pair_max=mp.mpf('0.0049296413205678555'), old_req_pairs_eta1=mp.mpf('0.031306064370482646'),
             old_core_frac_no_pairs=mp.mpf('0.1835032174604583'))

# (1) frozen test ---------------------------------------------------------------------------------------------
mX = gX * v
eps = M_new / Q
r = eps / mX
c_eta = 1 / r - 1
c_p = 2 / r
YQ_new = D / eps
R['frozen'] = dict(
    mX_GeV=mX, mX_matches_codex=abs(mX / mX_codex - 1), eps_GeV=eps, eps_over_mX=r,
    YQ_new=YQ_new, YQ_new_rel_codex=YQ_new / codex['YQ_new'] - 1,
    c_eta=c_eta, c_eta_rel_codex=c_eta / codex['c_eta'] - 1, c_p=c_p, c_p_rel_codex=c_p / codex['c_p'] - 1,
    eta_min=1 - mp.mpf('0.01') / c_eta, p_max=mp.mpf('0.01') / c_p,
    eta_min_rel_codex=(1 - mp.mpf('0.01') / c_eta) / codex['eta_min'] - 1, p_max_rel_codex=(mp.mpf('0.01') / c_p) / codex['p_max'] - 1,
    identity_check='(rho/s)/D = eta + (1-eta)/r + 2p/r = 1 + (1/r - 1)(1-eta) + (2/r) p  [exact]',
    m1_ledger_values=dict(c_eta=0.06678047, c_p=2.13356093, eta_min=0.8502556, p_max=0.0046870, YQ=5.079414233e-9))

# (2) old target ----------------------------------------------------------------------------------------------
m_eff_old = D / YQ_old
eps_old = M_old / Q
R['old_target'] = dict(
    D_over_YQ_old_GeV=m_eff_old, over_mX=m_eff_old / mX,
    d71cc01_late_mass_per_charge_GeV=late_mass_old / Q, identity_YQ_old_times_late_mass_over_D=YQ_old * late_mass_old / Q / D - 1,
    decomposition_object_plus_pairs_GeV=M_old + late_X_old * mX, remainder_D_particles_GeV=late_mass_old - M_old - late_X_old * mX,
    implied_mD_GeV=(late_mass_old - M_old - late_X_old * mX) / late_D_old,
    released_pairs_per_charge=scalars_old / Q, codex_required_pairs_eta1=codex['old_req_pairs_eta1'],
    old_object_fraction=YQ_old * eps_old / D, rel_codex=YQ_old * eps_old / D / codex['old_fraction'] - 1,
    old_multiplier=D / (YQ_old * eps_old), rel_codex_mult=D / (YQ_old * eps_old) / codex['old_mult'] - 1,
    kin_f_coefficient=2 * mX / m_eff_old, old_pair_max_1pct=mp.mpf('0.01') / (2 * mX / m_eff_old), rel_codex_pair_max=(mp.mpf('0.01') / (2 * mX / m_eff_old)) / codex['old_pair_max'] - 1,
    old_core_fraction_no_pairs=(1 - m_eff_old / mX) / (1 - eps_old / mX), rel_codex_core_frac=((1 - m_eff_old / mX) / (1 - eps_old / mX)) / codex['old_core_frac_no_pairs'] - 1,
    reading='D/Y_Q,old = 0.98593 m_X = object 0.923316 m_X + 2 x 0.031306 m_X released pairs (+1e-6 D particles): the old '
            'target was object + scalar-release debris, calibrated through v (d71cc01 scale_status), as m1_ledger says')

# (3) Codex's yield formula from its own row numbers (no engine) -----------------------------------------------
MRED = mp.mpf('2.435e18'); GSTAR = mp.mpf('106.75')
rows = {  # MATCHED_RESULTS / DEFORMED_MATCHED_RESULTS [Q]
    'b3_g0.5_baseline': dict(NX='-8120.740082505459', a='56.3606178466993', T='7702499481901.287', H='133349619398.33566', Y='5.026475516056739e-09', scale='1.0'),
    'b5_g1_baseline': dict(NX='-5210.315329485499', a='48.59925395539165', T='7700753422824.31', H='133349619398.33566', Y='5.033450222058231e-09', scale='1.0'),
    'scale_0.5': dict(NX='-8120.74008406811', a='56.36061782457464', T='5446489617877.453', H='66674809699.16783', Y='1.7771274617767915e-09', scale='0.5'),
    'scale_2.0': dict(NX='-8120.740086950753', a='56.36061785535858', T='10892979229798.16', H='266699238796.67133', Y='1.4217019699288572e-08', scale='2.0'),
    'b3_g0.5_matched_oldobj': dict(NX='-8120.740082531376', a='56.360617837163446', T='7768548315222.136', H='135646367122.97003', Y='5.156893366462146e-09', scale='1.0172235041614452'),
    'M1_b3_g0.5_matched': dict(NX='-8120.740082517851', a='56.36061784655097', T='7729445947446.704', H='134284273491.30815', Y='5.079414232732486e-09', scale='1.0070090495735164'),
    'M1_b5_g1_matched': dict(NX='-5210.315329489032', a='48.5992539548102', T='7724122783045.604', H='134160195609.41566', Y='5.079414232730545e-09', scale='1.006078579112091'),
}
yf = {}
for k, d in rows.items():
    NX, a, T, H, Y = (mp.mpf(d[x]) for x in ('NX', 'a', 'T', 'H', 'Y'))
    s = 2 * mp.pi ** 2 / 45 * GSTAR * T ** 3
    Yrec = abs(NX) * H ** 3 / (a ** 3 * s)
    yf[k] = dict(Y_from_row_numbers=Yrec, rel_to_row_Y=Yrec / Y - 1, H_over_H0=H / mp.mpf('133349619398.33566'))
R['codex_yield_formula_from_rows'] = yf
# exponents from Codex's scale controls (factor 4 in scale) and the matched rows
def ex(num, den, key):
    rn, rd = rows[num], rows[den]
    lam = mp.mpf(rn['scale']) / mp.mpf(rd['scale'])
    return dict(lambda_ratio=lam, exponent=mp.log(mp.mpf(rn[key]) / mp.mpf(rd[key])) / mp.log(lam))
R['exponent_from_codex_rows'] = dict(
    Y_scale2_over_scale05=ex('scale_2.0', 'scale_0.5', 'Y'), T_scale2_over_scale05=ex('scale_2.0', 'scale_0.5', 'T'),
    Y_scale2_over_base=ex('scale_2.0', 'b3_g0.5_baseline', 'Y'), Y_scale05_over_base=ex('scale_0.5', 'b3_g0.5_baseline', 'Y'),
    NX_spread_over_factor4=dict(s05=mp.mpf(rows['scale_0.5']['NX']) / mp.mpf(rows['b3_g0.5_baseline']['NX']) - 1,
                                s2=mp.mpf(rows['scale_2.0']['NX']) / mp.mpf(rows['b3_g0.5_baseline']['NX']) - 1),
    a_spread_over_factor4=dict(s05=mp.mpf(rows['scale_0.5']['a']) / mp.mpf(rows['b3_g0.5_baseline']['a']) - 1,
                               s2=mp.mpf(rows['scale_2.0']['a']) / mp.mpf(rows['b3_g0.5_baseline']['a']) - 1),
    Y_M1_b3=ex('M1_b3_g0.5_matched', 'b3_g0.5_baseline', 'Y'), Y_M1_b5=ex('M1_b5_g1_matched', 'b5_g1_baseline', 'Y'),
    Y_oldobj_b3=ex('b3_g0.5_matched_oldobj', 'b3_g0.5_baseline', 'Y'),
    eta2_effect_expected='(lambda^2-1) x 3e-15 ~ 1e-14 relative; the observed 1e-10 spread in N_X and a over a factor 4 is integrator noise (rtol 2e-11)')

# (5) lambda chain ---------------------------------------------------------------------------------------------
Yb3 = mp.mpf(rows['b3_g0.5_baseline']['Y']); Yb5 = mp.mpf(rows['b5_g1_baseline']['Y'])
tgt_old = D / eps_old; tgt_new = D / eps
R['lambda_chain'] = dict(
    old_point_yield_over_YQ_old=Yb3 / YQ_old, b5=Yb5 / YQ_old,
    lambda_M1_b3=(tgt_new / Yb3) ** (mp.mpf(2) / 3), codex=mp.mpf('1.0070090495735164'),
    lambda_M1_b5=(tgt_new / Yb5) ** (mp.mpf(2) / 3), codex_b5=mp.mpf('1.006078579112091'),
    lambda_oldobj_b3=(tgt_old / Yb3) ** (mp.mpf(2) / 3), codex_oldobj_b3=mp.mpf('1.0172235041614452'),
    lambda_oldobj_b5=(tgt_old / Yb5) ** (mp.mpf(2) / 3), codex_oldobj_b5=mp.mpf('1.0162835956038305'),
    identity_1p0518=dict(lhs=tgt_new / YQ_old, rhs=(Yb3 / YQ_old) * mp.mpf('1.0070090495735164') ** mp.mpf(1.5)),
    brief_naive_scale=(tgt_new / YQ_old) ** (mp.mpf(2) / 3),
    H_M1_GeV=mp.mpf('133349619398.33566') * (tgt_new / Yb3) ** (mp.mpf(2) / 3), codex_H=mp.mpf('134284273491.30815'),
    note='m1_ledger quotes the (5,1) old-object repair as 1.016277 (+1.628%); Codex row b5_g1_matched has 1.0162835956 (+1.628%): same to 5 digits')

# (6) clock/window ambiguity of lambda from Codex's own charge_clock entries ---------------------------------------
clock = {'20': '-8117.881574752486', '50': '-8104.618318826644', '80': '-8147.109179920771', '200': '-8133.914222541242', '800': '-8116.072835674296'}
NXw = abs(mp.mpf(rows['b3_g0.5_baseline']['NX']))
amb = {}
for t, n in clock.items():
    ratio = abs(mp.mpf(n)) / NXw                      # Y scales with N_X at fixed everything else
    amb['tau_' + t] = dict(NX=mp.mpf(n), Y_over_window=ratio, lambda_if_this_clock=(tgt_new / (Yb3 * ratio)) ** (mp.mpf(2) / 3))
lams = [x['lambda_if_this_clock'] for x in amb.values()]
R['clock_ambiguity'] = dict(entries=amb, lambda_window_640_800=(tgt_new / Yb3) ** (mp.mpf(2) / 3),
                            lambda_min=min(lams), lambda_max=max(lams),
                            codex_instantaneous_tau80_over_tau800=abs(mp.mpf(clock['80'])) / abs(mp.mpf(clock['800'])),
                            m1_ledger_window_ratio_tau80_over_tau800=mp.mpf('1.0044324235464697'),
                            our_engine_NX_window_64_80=mp.mpf('8156.734641992738'), our_engine_NX_window_640_800=mp.mpf('8120.740082437151'),
                            remark='m1_ledger pairs its window ratio 1.00443 with Codex instantaneous clocks (-8147.1/-8116.1 = 1.00382): '
                                   'different estimators of the same drift; the window ratio is the right one for the quoted yields')

# (7) the target D itself --------------------------------------------------------------------------------------
Och2 = mp.mpf('0.1200'); dOch2 = mp.mpf('0.0012')     # [R] Planck 2018 VI, TT,TE,EE+lowE+lensing
rhoc_h2 = mp.mpf('1.05368e-5')                       # [R] PDG: critical density / h^2, GeV cm^-3
s0 = mp.mpf('2891.2')                                # [R] PDG: present entropy density, cm^-3 (T0 = 2.7255 K)
D_planck = Och2 * rhoc_h2 / s0
R['target_D'] = dict(D_codex=D, D_planck2018=D_planck, codex_over_planck=D / D_planck - 1, planck_fractional_uncertainty=dOch2 / Och2,
                     lambda_shift_if_D_planck=(D_planck / D) ** (mp.mpf(2) / 3) - 1, lambda_band_from_Och2=((1 + dOch2 / Och2) ** (mp.mpf(2) / 3) - 1),
                     remark='the +0.7009% calibration step is the same size as the rounding of the target (0.61%) and smaller than '
                            'its +-1% observational band; m1_ledger did not raise this; [R] inputs recalled, not fetched')

# (8) A from the M1 heavy parameters (hw_gates formula alpha6 = kappa^2 mu <Sigma>/(M_Y^4 - 4 mu^2 <Sigma>^2), <Sigma> = h/sqrt2)
kap = mp.mpf('0.0034586847211042523'); mu = mp.mpf('961402368973.181'); h = mp.mpf('5371370939652.326'); MY = mp.mpf('4028528204739.244')
lam = mp.mpf('1.0070090495735164'); H0 = mp.mpf('133349619398.33566')
sig = h / mp.sqrt(2)
A_M1 = kap ** 2 * mu * sig / (MY ** 4 - 4 * mu ** 2 * sig ** 2)
r2 = 2 * mu ** 2 * h ** 2 / MY ** 4
A_old_full = mp.mpf('2.108859793121751e-31')        # 1fcc6b7 COUPLING_RELATIONS.json
A_new_codex = mp.mpf('2.079605512546138e-31')
B_codex = mp.mpf('9.242691166871723e-31')
R['A_coupling'] = dict(A_from_M1_parameters=A_M1, codex_A=A_new_codex, rel=A_M1 / A_new_codex - 1,
                       r_squared=r2, r=mp.sqrt(r2), two_A_over_B=2 * A_M1 / B_codex,
                       A_old_full=A_old_full, A_old_over_lambda2=A_old_full / lam ** 2, rel_to_codex_new=A_old_full / lam ** 2 / A_new_codex - 1,
                       m1_ledger_used_A_old=mp.mpf('2.10886e-31'), m1_ledger_A_new=mp.mpf('2.0796057165545546e-31'),
                       heavy_ratios=dict(MY_over_H=MY / (lam * H0), h_over_H=h / (lam * H0), mu_over_H=mu / (lam * H0)),
                       scaling='kappa fixed, mu,h,M_Y ~ lambda: A ~ lambda^-2 exactly since r^2 is invariant')

# (9) Gamma t0 chain from Codex's quoted channel rates (old: e160356/8986d34; new: d4747cb COLD_LOSS_RESULTS) -------
AGE = mp.mpf('4.3549488e17')
coh_old = mp.mpf('3.91958726382564e-50'); born_old = mp.mpf('1.0399701293581432e-46')
coh_new = mp.mpf('6.099550997404576e-51'); born_new = mp.mpf('9.471343360338484e-47')
c6_old = A_old_full * MRED ** 2; c6_new = A_new_codex * MRED ** 2
G_old = (coh_old + born_old) * c6_old ** 2 * AGE
G_new = (coh_new + born_new) * c6_new ** 2 * AGE
R['gamma_t0'] = dict(c6_old=c6_old, c6_new=c6_new, codex_c6_new=mp.mpf('1233044.8995126376'),
                     Gamma_t0_old=G_old, codex_old_e160356=mp.mpf('7.08367857848652e-17'), rel_old=G_old / mp.mpf('7.08367857848652e-17') - 1,
                     m1_ledger_old=mp.mpf('7.083679968300075e-17'), m1_ledger_old_rel=mp.mpf('7.083679968300075e-17') / mp.mpf('7.08367857848652e-17') - 1,
                     Gamma_t0_new=G_new, codex_new=mp.mpf('6.271629359413795e-17'), rel_new=G_new / mp.mpf('6.271629359413795e-17') - 1,
                     ratio_new_over_old=G_new / G_old, profile_factor=(coh_new + born_new) / (coh_old + born_old), coupling_factor=(c6_new / c6_old) ** 2,
                     lambda_minus4=lam ** -4, codex_relative_included_change_at_fixed_profile=mp.mpf('0.9724482663579905'),
                     coherent_share_new=coh_new / (coh_new + born_new), coherent_share_old=coh_old / (coh_old + born_old),
                     age_check_13p8Gyr_julian=mp.mpf('13.8e9') * mp.mpf('365.25') * 86400)

# (10) pulse rows re-read with charge conservation (eta = 1 + (Nbar - Nplus)/Q0, p = Nbar/Q0) ---------------------
def term(eta, p):
    return c_eta * (1 - eta) + c_p * p
pulses = {'5%': dict(xbar='0.002088', nplus=None), '10%': dict(xbar='0.008448', nplus='94.859'), '20%': dict(xbar='0.03447', nplus=None), '40%': dict(xbar='0.1459', nplus=None)}
pr = {}
for k, d in pulses.items():
    p = mp.mpf(d['xbar'])
    npl = mp.mpf(d['nplus']) / Q if d['nplus'] else mp.mpf(0)
    eta = 1 + p - npl
    pr[k] = dict(p=p, Nplus_over_Q0=npl, eta_charge_conserving=eta, term_m1_ledger_eta1=term(1, p), term_charge_conserving=term(eta, p),
                 ratio=term(eta, p) / term(1, p), coefficient_on_Nbar=c_p - c_eta, coefficient_on_Nplus=c_eta)
R['pulse_rows'] = dict(rows=pr, remark='m1_ledger reads the emitted Xbar as p with eta = 1; charge conservation makes the core charge rise, '
                                        'eta = 1 + (Nbar - Nplus)/Q0 (ed73791 remark), so the frozen term is (c_p - c_eta) Nbar/Q0 + c_eta Nplus/Q0 = '
                                        '2.0668 Nbar/Q0 + 0.0668 Nplus/Q0: 3% below m1_ledger rows; verdict unchanged (5% passes, 10% fails x1.75)')

# (11) cold completions and the 10% kinetics row ---------------------------------------------------------------------
cc = [('5%', '199563.6967569675', '436.3032430325056', '17289.28104310131', '0.00027515575036596296'),
      ('10%', '198210.13106309355', '1789.8689369064627', '17178.96535550169', '0.0011274534942455894'),
      ('20%', '192174.09397143824', '7825.906028561768', '16686.59469287748', '0.00490314775719769'),
      ('40%', '159038.49838783825', '40961.50161216174', '13969.584034711766', '0.024817648841077844')]
ccr = {}
for k, Qf, Nf, Ef, ex_c in cc:
    Qf, Nf, Ef, ex_c = map(mp.mpf, (Qf, Nf, Ef, ex_c))
    eta = Qf / Q
    excess = (Ef + Nf * mX) / M_new - 1                  # mass of (core + free X at rest) over the reference mass
    ccr[k] = dict(eta=eta, Qf_plus_Nf_minus_Q=Qf + Nf - Q, excess_from_masses=excess, codex_excess=ex_c, rel=excess / ex_c - 1,
                  frozen_term=term(eta, 0), frozen_over_actual=term(eta, 0) / excess, EQ_over_mX_final=Ef / Qf / mX)
R['cold_completion'] = ccr
Qc = mp.mpf('188371.170'); nbar = mp.mpf('0.559'); nplus = mp.mpf('11629.389'); Ekin = mp.mpf('2526.289'); Egas = mp.mpf('14156.238')
eta_k = Qc / Q; p_k = nbar / Q
R['kinetics_10pct'] = dict(eta=eta_k, p=p_k, charge_closure=Qc + nplus - nbar - Q, free_X_from_ledger=(1 - eta_k + p_k) * Q, codex_free_X=nplus,
                           frozen_term=term(eta_k, p_k), m1_ledger_frozen=mp.mpf('0.003888856814259084'),
                           gas_kinetic_over_reference_mass=Ekin * mX / M_new, codex_quoted=mp.mpf('0.013475'),
                           gas_rest_plus_kinetic_check=(nplus + nbar) + Ekin, codex_gas_energy=Egas,
                           codex_resolved_excess=mp.mpf('0.007266'), total_if_heat_stays=mp.mpf('0.007266') + Ekin * mX / M_new)

# (12) kin_solve translation ----------------------------------------------------------------------------------------
conv = c_p / (2 * mX / m_eff_old)
R['kin_translation'] = dict(new_term_per_old_f=conv, equals_YQ_ratio=YQ_new / YQ_old, rel=conv / (YQ_new / YQ_old) - 1,
                            rows={name: dict(f=mp.mpf(f), p=mp.mpf(f) / (2 * mX / m_eff_old), term=conv * mp.mpf(f)) for name, f in
                                  (('R4.9665', '0.0045855'), ('R5', '0.0047584'), ('R7', '0.0288444'), ('R41.7', '1.4921315'))},
                            crossing_R=dict(old=5 + mp.log(mp.mpf('0.01') / mp.mpf('0.0047584')) / (mp.log(mp.mpf('0.0288444') / mp.mpf('0.0047584')) / 2),
                                            new=5 + mp.log(mp.mpf('0.01') / conv / mp.mpf('0.0047584')) / (mp.log(mp.mpf('0.0288444') / mp.mpf('0.0047584')) / 2)))

def conv_json(o):
    if isinstance(o, dict):
        return {k: conv_json(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [conv_json(v) for v in o]
    if isinstance(o, mp.mpf):
        return float(o)
    return o
OUT.mkdir(parents=True, exist_ok=True)
(OUT / 'skeptic_algebra.json').write_text(json.dumps(dict(method=__doc__, results=conv_json(R)), indent=1))
print(json.dumps(conv_json(R), indent=1))
