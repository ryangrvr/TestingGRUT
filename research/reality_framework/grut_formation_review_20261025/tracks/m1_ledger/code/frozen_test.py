#!/usr/bin/env python3
"""m1_ledger (Claude, own code): the frozen mass test, its coefficients, and the (eta, p) translation of the
computed populations.  Items 1 and 3 of the brief.  Pure algebra; no Codex code imported or run.

Quoted inputs [Q] are marked with their source commit.  Output: ../results/frozen_test.json
"""
import json
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parents[1] / 'results'

# ---- quoted inputs ------------------------------------------------------------------------------------------------
v = 0.924089637375362                 # [Q d4747cb CANDIDATE_M1.json]
gX = 0.1
mX = gX * v                           # 0.0924089637 GeV
D = 4.4e-10                           # GeV per unit visible entropy, the fixed target [Q 1fcc6b7/d4747cb]
Q = 200000
M_new = 17324.832346421026            # GeV [Q d4747cb DEFORMATION_RESULTS primary]
eps_new = M_new / Q                   # E(Q)/Q, GeV
M_old = 17064.537454385838            # GeV [Q 1fcc6b7 ABUNDANCE_LEDGER]
eps_old = M_old / Q
YQ_old = 4.8294009076988175e-9        # [Q 1fcc6b7; rp.py; the supplied primordial yield of the heavy spot point]
YQ_new_codex = 5.079414232725842e-09  # [Q d4747cb]
c_eta_codex, c_p_codex = 0.0667805, 2.13356
eta_min_codex, p_max_codex = 0.85026, 0.004687

# ---- item 1: derivation --------------------------------------------------------------------------------------------
# rho/s = Y_Q [eta eps_Q + (1-eta) m_X + 2 m_X p]; with Y_Q eps_Q = D:
# (rho/s)/D = eta + (1-eta) m_X/eps_Q + 2 p m_X/eps_Q = 1 + (1-eta)(m_X/eps_Q - 1) + (2 m_X/eps_Q) p
r_new = eps_new / mX
c_eta = mX / eps_new - 1              # = (1 - eps/m)/(eps/m)
c_p = 2 * mX / eps_new
YQ_new = D / eps_new
tol = 0.01
eta_min = 1 - tol / c_eta
p_max = tol / c_p

r_old = eps_old / mX
m_eff_old = D / YQ_old                # the mass per net charge the OLD target assumed
old_objects_fraction = YQ_old * eps_old / D
YQ_mult_old_repair = D / (YQ_old * eps_old)   # 1.0678: multiplier the old object needed
c_p_old_convention = 2 * mX / m_eff_old       # kin_solve's f = 2.0286 P/n_Q

item1 = dict(
    mX_GeV=mX, D_GeV=D, eps_new_GeV=eps_new, eps_new_over_mX=r_new, codex_eps_over_mX=0.93739999,
    YQ_new=YQ_new, YQ_new_codex=YQ_new_codex, YQ_new_rel_diff=YQ_new / YQ_new_codex - 1,
    YQ_new_times_eps_new_GeV=YQ_new_codex * eps_new,
    coefficient_unbound_charge=c_eta, codex_coefficient_unbound_charge=c_eta_codex,
    coefficient_pairs=c_p, codex_coefficient_pairs=c_p_codex,
    eta_min_at_p0=eta_min, codex_eta_min=eta_min_codex, p_max_at_eta1=p_max, codex_p_max=p_max_codex,
    identity='(rho/s)/D = 1 + (mX/eps-1)(1-eta) + (2mX/eps) p, exact given Y_Q eps = D and the three-reservoir ledger',
    charge_constraint='N_plus/Q_net = 1 - eta + p >= 0, so eta > 1 is legal when p > 0 (ed73791 remark confirmed)',
    old_target=dict(YQ_old=YQ_old, eps_old_GeV=eps_old, eps_old_over_mX=r_old,
                    mass_per_net_charge_assumed_GeV=m_eff_old, mass_per_net_charge_assumed_over_mX=m_eff_old / mX,
                    note='the old target was NOT an object E/Q: D/Y_Q,old = 0.98592 m_X is a v-calibrated supplied mass '
                         'per net charge (d71cc01 FIXED_YIELD_RESULTS: late_mass_per_object 18221.72 GeV vs object 17064.54)',
                    old_objects_deliver_fraction_of_D=old_objects_fraction, codex_quoted=0.936494,
                    YQ_multiplier_for_old_objects=YQ_mult_old_repair, codex_quoted_multiplier=1.0678122,
                    kin_solve_f_coefficient=c_p_old_convention, codex_kin_f_coefficient=2.0286),
)

# ---- item 3: translate computed populations -------------------------------------------------------------------------
def term(eta, p):
    return c_eta * (1 - eta) + c_p * p

def row(name, label, eta, p, source, note=''):
    t = term(eta, p)
    return dict(scenario=name, label=label, eta=eta, p=p, term=t, mass_over_D=1 + t, passes_1pct=bool(t <= tol),
                factor_over_allowance=t / tol, source=source, note=note)

rows = []
# (a) pt_relax: Xbar per net charge BEFORE annihilation (p if nothing annihilates)
for p in (0.069, 0.110):
    rows.append(row('pt_relax equilibrium Xbar/charge at the pilot energy (no annihilation)', '[E]', 1.0, p,
                    'tracks/pt_relax NOTES: 0.069-0.110',
                    'p here = Xbar per net charge if none annihilate; survival fraction needed to pass at eta=1: %.4f' % (p_max / p)))
for p in (1.9, 2.9):
    rows.append(row('pt_relax equilibrium with e=1e-7 vector heat thermalised (no annihilation)', '[E]', 1.0, p,
                    'tracks/pt_relax NOTES: 1.9-2.9', 'survival fraction needed: %.5f' % (p_max / p)))
rows.append(row('classical prethermal plateau 1.75 Xbar/charge (no annihilation)', '[C]', 1.0, 1.75,
                'tracks/pt_relax / Codex 5f5513b', 'relaxes classically at ~1e-3 per conformal unit vs ~1e11-1e12 units available'))
# (b) kin_solve case L: f = 2.0286 P/n_Q is the pair MASS over the old target budget per net charge; p = f/2.0286
conv = c_p / c_p_old_convention       # new term = conv * f_old
kin = [('e=1e-6 band gas-only, R=4.9665', 0.0045855), ('diagnostic R=5 gas-only', 0.0047584),
       ('diagnostic R=7 gas-only', 0.0288444), ('e=1e-7 band gas-only, R=41.7', 1.4921315)]
for name, f in kin:
    p = f / c_p_old_convention
    rows.append(row('kin_solve ' + name, '[C]', 1.0, p, 'REVIEW_KINETIC.md reproduction table',
                    'f=%.5f (old convention, 2.0286 P/n_Q); new term = %.5f f; room for (1-eta) <= %.4f' % (f, conv, max(0.0, (tol - c_p * p) / c_eta))))
# crossing in R: log-linear between R=5 and R=7
f5, f7 = 0.0047584, 0.0288444
slope = np.log(f7 / f5) / 2
R_cross_old = 5 + np.log(0.01 / f5) / slope
f_new_allow = tol / conv
R_cross_new = 5 + np.log(f_new_allow / f5) / slope
# (c) Codex's own eta-type numbers
rows.append(row('cold dilute masks hold 30% of charge, read as eta (p=0)', '[C Codex ed73791]', 0.30, 0.0, 'ed73791 REPORT',
                'Codex itself says a mask fraction is not final eta'))
rows.append(row('cold dilute masks hold 32% of charge, read as eta (p=0)', '[C Codex ed73791]', 0.32, 0.0, 'ed73791 REPORT'))
rows.append(row('dense primary masks 49.4-51.7% (p=0), lower value', '[C Codex ed73791]', 0.494, 0.0, 'ed73791 REPORT'))
rows.append(row('hot elliptic control 0.6% in dense masks (64^3), p=0', '[C Codex ed73791]', 0.006, 0.0, 'ed73791 REPORT', 'not spatially converged; 96^3 gives 3.45%'))
rows.append(row('hot elliptic control 3.45% (96^3), p=0', '[C Codex ed73791]', 0.0345, 0.0, 'ed73791 REPORT'))
for amp, xbar in (('5%', 0.002088), ('10%', 0.008448), ('20%', 0.03447), ('40%', 0.1459)):
    rows.append(row('pulse %s: measured outgoing Xbar/Q0 as raw leftover p (eta=1, no annihilation)' % amp, '[C Codex bf8ab3b]', 1.0, xbar,
                    'bf8ab3b REPORT (direct Xbar/Q0; the receiver net charge 0.0020/0.0080/0.0304/0.0949 is a lower bound)'))
# Codex's conditional cold completion: all Xbar converted, free X stay free, E(Q_final)/Q_final re-solved
cc = [('5%', 0.9978184837848374, 0.9375216239335298, 0.00027515575036596296),
      ('10%', 0.9910506553154678, 0.9379011211496998, 0.0011274534942455894),
      ('20%', 0.9608704698571912, 0.9396341123676687, 0.00490314775719769),
      ('40%', 0.7951924919391913, 0.9505328004426522, 0.024817648841077844)]
for amp, eta, eq, codex_excess in cc:
    frozen = term(eta, 0.0)
    actual = eta * eq / r_new + (1 - eta) / r_new - 1       # same ledger with the re-solved E/Q at the smaller Q
    rows.append(dict(scenario='pulse %s conditional cold completion (all Xbar converted, p=0)' % amp, label='[C Codex bf8ab3b]',
                     eta=eta, p=0.0, term=frozen, mass_over_D=1 + frozen, passes_1pct=bool(frozen <= tol),
                     factor_over_allowance=frozen / tol, source='bf8ab3b COLD_COMPLETION_RESULTS.json',
                     codex_excess_with_resolved_EQ=codex_excess, our_excess_with_resolved_EQ=actual,
                     note='frozen coefficient (fixed eps_Q) understates the excess by %.2fx because E/Q rises at smaller Q' % (actual / frozen)))
# 10% pulse kinetics (c481530)
eta_k, nbar_k = 188371.170 / Q, 0.559 / Q
rows.append(dict(scenario='10% pulse kinetics tau=1000: Q_core/Q0=0.94186, 0.559 Xbar per core', label='[C Codex c481530]',
                 eta=eta_k, p=nbar_k, term=term(eta_k, nbar_k), mass_over_D=1 + term(eta_k, nbar_k),
                 passes_1pct=bool(term(eta_k, nbar_k) <= tol), factor_over_allowance=term(eta_k, nbar_k) / tol,
                 source='c481530 REPORT', codex_stationary_excess=0.007266, codex_gas_kinetic_energy_fraction=0.013475,
                 note='frozen term %.5f; Codex 0.7266%% with re-solved E/Q; gas still holds 1.3475%% of the reference mass as kinetic energy '
                      '(2.07%% total if it does not leave); fixed-profile experiment, not a formation ensemble' % term(eta_k, nbar_k)))

passing = [r['scenario'] for r in rows if r['passes_1pct']]
item3 = dict(
    conversion_new_term_per_old_f=conv,
    kin_solve_crossing_R_old_convention=float(R_cross_old), kin_solve_crossing_R_new_convention=float(R_cross_new),
    f_old_allowed_by_new_test=f_new_allow, rows=rows, scenarios_passing=passing,
    what_decides_eta='a formation ensemble: expanding 3D evolution from the actual (matched) release spectrum with the changed '
                     'scalar cubic, giving the charge-weighted object spectrum y_Q and the core/tail membership; mask fractions '
                     '(30-52%) and the fixed-profile pulses are not that ensemble',
    what_decides_p='chemical relaxation of the Xbar inventory (number-changing 2<->4 and condensate-assisted channels, pt_relax) '
                   'followed by annihilation/conversion with the ACTUAL heat R (vector heat ~ 1/e) and object capture; '
                   'kin_solve (forced-thermal start, gas only) and c481530 (fixed profile, one pulse) are partial pieces',
)

OUT.mkdir(exist_ok=True)
(OUT / 'frozen_test.json').write_text(json.dumps(dict(method=__doc__, item1=item1, item3=item3), indent=1, default=float))
print(json.dumps(item1, indent=1, default=float))
print('conv', conv, 'R_cross old/new', R_cross_old, R_cross_new, 'f allowed', f_new_allow)
for r in rows:
    print('%-95s eta %.5f p %.6f term %.5f %s x%.1f' % (r['scenario'][:95], r['eta'], r['p'], r['term'], 'PASS' if r['passes_1pct'] else 'fail', r['factor_over_allowance']))
