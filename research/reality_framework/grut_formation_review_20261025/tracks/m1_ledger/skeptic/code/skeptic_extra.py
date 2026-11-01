#!/usr/bin/env python3
"""m1_ledger_skeptic, second pass (Claude, own code, algebra only; no engine, no Codex code run).

Extra adversarial checks on tracks/m1_ledger that the first skeptic pass did not make:
 (1) the basis of the H_i^{3/2} exponent: r = rho_rad/(M_red^2 H_i^2) rebuilt from Codex's (T_vis, H) row numbers
     must be the same number on every row of a history (T ~ (H_i M_red)^{1/2} is then exact, not assumed);
 (2) the history dependence of the calibration: lambda(3,0.5) vs lambda(5,1);
 (3) the sensitivity of lambda to the target D (dlnlambda = 2/3 dlnD) and the clock band, combined;
 (4) the pulse rows of m1_ledger re-read with the ACTUAL N_plus and N_phi of bf8ab3b SOURCE_RESULTS.json for all four
     pulses (the first pass had N_plus only for the 10% pulse and set it to zero elsewhere);
 (5) the interpolation dependence of the kin_solve crossing R (log-linear vs linear between the R = 5 and R = 7 points);
 (6) Codex's late_half_window_relative_charge_change translated into a lambda shift.
Inputs are [Q] numbers read through git show (d4747cb MATCHED/DEFORMED_MATCHED_RESULTS.json, CONSEQUENCES.json;
bf8ab3b SOURCE_RESULTS.json; REVIEW_KINETIC.md).  Output: ../results/skeptic_extra.json
"""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps = 40
OUT = Path(__file__).resolve().parents[1] / 'results'
R = {}
MRED = mp.mpf('2.435e18'); GSTAR = mp.mpf('106.75')
D = mp.mpf('4.4e-10'); Q = 200000
v = mp.mpf('0.924089637375362'); mX = mp.mpf('0.1') * v
eps = mp.mpf('17324.832346421026') / Q
c_eta = mX / eps - 1; c_p = 2 * mX / eps

# (1) r = rho_rad/(M_red^2 H^2) from (T_vis, H) on each Codex row  ----------------------------------------------
rows = {
    'b3_g0.5_baseline': ('7702499481901.287', '133349619398.33566', '5.026475516056739e-09'),
    'b3_g0.5_matched_oldobj': ('7768548315222.136', '135646367122.97003', '5.156893366462146e-09'),
    'scale_0.5': ('5446489617877.453', '66674809699.16783', '1.7771274617767915e-09'),
    'scale_2.0': ('10892979229798.16', '266699238796.67133', '1.4217019699288572e-08'),
    'M1_b3_g0.5_matched': ('7729445947446.704', '134284273491.30815', '5.079414232732486e-09'),
    'b5_g1_baseline': ('7700753422824.31', '133349619398.33566', '5.033450222058231e-09'),
    'M1_b5_g1_matched': ('7724122783045.604', '134160195609.41566', '5.079414232730545e-09'),
}
rr = {}
for k, (T, H, Y) in rows.items():
    T, H = mp.mpf(T), mp.mpf(H)
    rho = mp.pi ** 2 / 30 * GSTAR * T ** 4
    rr[k] = dict(r_rad=rho / (MRED ** 2 * H ** 2), T_over_sqrt_HM=T / mp.sqrt(H * MRED))
b3 = [rr[k]['r_rad'] for k in rows if k.startswith(('b3', 'scale', 'M1_b3'))]
b5 = [rr[k]['r_rad'] for k in rows if 'b5' in k]
R['radiation_constant_from_rows'] = dict(rows=rr,
    b3_spread=(max(b3) - min(b3)) / min(b3), b5_spread=(max(b5) - min(b5)) / min(b5),
    remark='r is the same on every (3,0.5) row (baseline, scale 0.5, 2.0, old-object match, M1 match) to the spread quoted, '
           'so T_vis ~ (H_i M_red)^{1/2} exactly and Y ~ H_i^{3/2} follows from N_X, a unchanged; the engine run adds nothing')

# (2) history dependence ---------------------------------------------------------------------------------------
lam_b3 = mp.mpf('1.0070090495735164'); lam_b5 = mp.mpf('1.006078579112091')
R['history_dependence'] = dict(lambda_b3=lam_b3, lambda_b5=lam_b5, relative=lam_b5 / lam_b3 - 1,
                               remark='the two supplied reheating histories differ by -0.092% in lambda; the quoted +0.7009% is the (3,0.5) choice')

# (3) D sensitivity and the combined band ------------------------------------------------------------------------
Yb3 = mp.mpf('5.026475516056739e-09')
D_planck = mp.mpf('0.1200') * mp.mpf('1.05368e-5') / mp.mpf('2891.2')
clock = {'20': '-8117.881574752486', '50': '-8104.618318826644', '80': '-8147.109179920771', '200': '-8133.914222541242', '800': '-8116.072835674296'}
NXw = mp.mpf('8120.740082505459')
def lam(Dv, nx=NXw):
    return ((Dv / eps) / (Yb3 * nx / NXw)) ** (mp.mpf(2) / 3)
lams = {}
for Dname, Dv in (('codex_4.4e-10', D), ('planck_central', D_planck), ('planck_minus_1pct', D_planck * mp.mpf('0.99')), ('planck_plus_1pct', D_planck * mp.mpf('1.01'))):
    for t, n in clock.items():
        lams[Dname + '_tau' + t] = lam(Dv, abs(mp.mpf(n)))
    lams[Dname + '_window'] = lam(Dv)
R['lambda_band'] = dict(values=lams, min=min(lams.values()), max=max(lams.values()),
                        dlnlambda_per_dlnD=mp.mpf(2) / 3, D_planck_GeV=D_planck,
                        remark='every value is a legal reading of the same calculation (target rounding, clock choice); the spread is '
                               '%.2f%% to %.2f%% against the quoted +0.7009%%' % (100 * (min(lams.values()) - 1), 100 * (max(lams.values()) - 1)))

# (4) pulse rows with the actual N_plus and N_phi (bf8ab3b SOURCE_RESULTS.json) ------------------------------------
src = {'5%': ('417.6076896560908', '16.976768462493535', '0.859392456960631'),
       '10%': ('1689.5063628733974', '89.35622473151756', '5.503174650773861'),
       '20%': ('6893.863229804531', '824.6722008107597', '53.68529897323827'),
       '40%': ('29185.758436050142', '10596.704346059894', '589.5194150258485')}
m1_rows = {'5%': 0.0044548752317094785, '10%': 0.018024322776571684, '20%': 0.07354384542003148, '40%': 0.3112865403766345}
first_pass = {'5%': 0.004315437615854739, '10%': 0.017491835030053957, '20%': 0.07124192271001574, '40%': 0.3015432701883173}
pr = {}
for k, (nm, npl, nphi) in src.items():
    nm, npl, nphi = mp.mpf(nm) / Q, mp.mpf(npl) / Q, mp.mpf(nphi) / Q
    eta = 1 + nm - npl                       # core charge rises by N_minus - N_plus
    t_bar = c_eta * (1 - eta) + c_p * nm     # p = N_minus/Q0
    t_phi = c_eta * (1 - eta) + c_p * (nm + nphi)   # each phi -> X Xbar adds one pair (eta unchanged)
    pr[k] = dict(Nminus_over_Q0=nm, Nplus_over_Q0=npl, Nphi_over_Q0=nphi, eta=eta, free_X_over_Q0=1 - eta + nm,
                 term_m1_ledger=m1_rows[k], term_first_skeptic_pass=first_pass[k], term_with_actual_Nplus=t_bar, term_with_Nplus_and_phi_pairs=t_phi,
                 ratio_to_m1_ledger=t_phi / m1_rows[k], factor_over_1pct=t_phi / mp.mpf('0.01'))
R['pulse_rows_actual_Nplus_Nphi'] = dict(rows=pr,
    remark='charge-conserving bookkeeping with the measured N_plus (not zero) and the emitted phi counted as future pairs: the first skeptic '
           'pass underestimated the 20%/40% rows by 0.4%/1.2% (N_plus set to 0); net effect on m1_ledger rows -3.0/-2.9/-2.3/+0.2%; verdicts unchanged')

# (5) kin_solve crossing R: interpolation dependence ----------------------------------------------------------------
f5, f7 = mp.mpf('0.0047584'), mp.mpf('0.0288444'); conv = c_p / (2 * mX / (D / mp.mpf('4.8294009076988175e-9')))
allow_old, allow_new = mp.mpf('0.01'), mp.mpf('0.01') / conv
def cross(allow, mode):
    if mode == 'log':
        return 5 + mp.log(allow / f5) / (mp.log(f7 / f5) / 2)
    return 5 + 2 * (allow - f5) / (f7 - f5)
R['kin_crossing_interpolation'] = dict(log_linear=dict(old=cross(allow_old, 'log'), new=cross(allow_new, 'log')),
                                       linear=dict(old=cross(allow_old, 'lin'), new=cross(allow_new, 'lin')),
                                       conv_new_term_per_f=conv, allowed_f_new=allow_new,
                                       remark='with only the R = 5 and R = 7 points, the crossing is 5.82 (log-linear, m1_ledger and REVIEW_KINETIC) or '
                                              '5.44 (linear); the M1 shift (-0.056 log-linear, -0.041 linear) is robust in sign, the absolute R is not')

# (6) half-window drift -> lambda ----------------------------------------------------------------------------------
hw = mp.mpf('0.0002617013414302205')
R['half_window_to_lambda'] = dict(codex_late_half_window_relative_charge_change=hw, lambda_shift=(1 + hw) ** (mp.mpf(2) / 3) - 1,
                                  remark='the drift inside the averaging window moves lambda by 1.7e-4, i.e. the 8th quoted digit of lambda is noise at the 4th')

def conv_json(o):
    if isinstance(o, dict):
        return {k: conv_json(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [conv_json(v) for v in o]
    if isinstance(o, mp.mpf):
        return float(o)
    return o
OUT.mkdir(parents=True, exist_ok=True)
(OUT / 'skeptic_extra.json').write_text(json.dumps(dict(method=__doc__, results=conv_json(R)), indent=1))
out = conv_json(R)
print('r spread b3 %.2e b5 %.2e' % (out['radiation_constant_from_rows']['b3_spread'], out['radiation_constant_from_rows']['b5_spread']))
print('r values', {k: '%.12f' % v['r_rad'] for k, v in out['radiation_constant_from_rows']['rows'].items()})
print('history dependence', out['history_dependence']['relative'])
print('lambda band', out['lambda_band']['min'], out['lambda_band']['max'], out['lambda_band']['remark'])
for k, r in out['pulse_rows_actual_Nplus_Nphi']['rows'].items():
    print(k, 'm1 %.5f first %.5f actualNplus %.5f withphi %.5f ratio %.4f x%.2f' % (r['term_m1_ledger'], r['term_first_skeptic_pass'], r['term_with_actual_Nplus'], r['term_with_Nplus_and_phi_pairs'], r['ratio_to_m1_ledger'], r['factor_over_1pct']))
print('crossing', out['kin_crossing_interpolation']['log_linear'], out['kin_crossing_interpolation']['linear'])
print('half window -> lambda', out['half_window_to_lambda']['lambda_shift'])
