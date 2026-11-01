"""fx_formation item 4 (Claude, 1 Nov 2026): verdict table from family.json, eq_history.json, eq_saha.json, kinetics.json.

Mass per net charge of the relic (m_X units):  M = eta <E/Q> + (1 - eta) * 1 + 2 p   (free X at rest mass, pairs 2 m_X)
Reference (the calibration's object): 0.9374 m_X.  Frozen test (REVIEW_M1, reference form): 0.0667805 (1 - eta) + 2.13356 p.
"""
import json, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
FXC = os.path.join(HERE, '..', '..', 'fx_chain', 'results', 'heat_onset.json')
sys.path.insert(0, HERE)

REF = 0.9373999905262296
fam = json.load(open(os.path.join(RES, 'family.json')))
EH = json.load(open(os.path.join(RES, 'eq_history.json')))
KI = json.load(open(os.path.join(RES, 'kinetics.json')))
f3 = fam['fit3']


def eq_of_Q(Q):
    return f3['Omega_b'] + f3['A'] * Q ** (-1 / 3) + f3['B'] * Q ** (-2 / 3) + f3['C'] / Q


out = {}
out['E_over_Q_of_mean_Q'] = {('%.0e' % Q): eq_of_Q(Q) for Q in (2e5, 3e5, 1e8, 1e12, 1e13, 1e16, 1e20)}
out['E_over_Q_of_mean_Q']['2e+05_solver'] = REF
# Y_h(R) interpolation and R_crit for the two release choices
HO = json.load(open(FXC))
Ycrit = EH['Y_crit']['Y_crit']
crit = {}
for rel in ('heavy_release', 'pilot_state'):
    pts = sorted([(c['R'], c['energy_per_charge_conf']) for c in HO['chain'] if c['phi_in_equilibrium'] and c['label'].startswith(rel)])
    lab = {c['label']: c for c in HO['chain'] if c['phi_in_equilibrium'] and c['label'].startswith(rel)}
    Rs, Ys = [], []
    for L, c in lab.items():
        o = c['onset_untilted']
        Rs.append(c['R']); Ys.append(o['nQ_over_mX3'] / (2 * math.pi ** 2 / 45 * 3 * o['T_h_over_mX'] ** 3))
    o = np.argsort(Rs); Rs = np.array(Rs)[o]; Ys = np.array(Ys)[o]
    # log Y linear in log(1 + R/R1) with R1 = 1 (smooth monotone interpolation)
    x = np.log1p(Rs); y = np.log(Ys)
    xs = np.linspace(x[0], x[-1], 4001); ys = np.interp(xs, x, y)
    j = int(np.argmin(np.abs(ys - math.log(Ycrit))))
    Rc = float(np.expm1(xs[j]))
    # e_crit: R_tot(e) = R_vec(1e-6) (1e-6/e) + R_med(e); R_vec(1e-6) band = 4.0538, R_med 0.913 (1e-6) .. 1.144 (1e-7)
    Rv, Rm = 4.053801019129129, 0.9127346451110923
    e_c = 1e-6 * Rv / max(Rc - Rm, 1e-9) if Rc > Rm else None
    crit[rel] = dict(R=list(map(float, Rs)), Y_h=list(map(float, Ys)), R_crit=Rc, e_crit_band=e_c)
out['survival_threshold'] = dict(Y_crit_bulk=Ycrit, S_per_charge_max=1 / Ycrit, Y_crit_finite_Q=EH['Y_crit_finite_Q'], by_release=crit)
# case table
cases = []
for key, pgas in (('R0_heavy|Y_onset', 0.0), ('med_e1e6_heavy|Y_onset', 0.0), ('tot_e1e6_band_heavy|Y_onset', 0.0045855 / 2.028545),
                  ('tot_e1e6_band_pilot|Y_onset', 0.0045855 / 2.028545), ('tot_e1e7_band_heavy|Y_onset', 1.4921315 / 2.028545)):
    H = EH['histories'][key]
    survive = H['T_ev'] is None
    row = dict(case=key, Y_h=H['Y_h'], S_per_charge=1 / H['Y_h'], T_ev=H['T_ev'], T_s=H['T_s'], survive=survive,
               eta_min_equilibrium=H['eta_min_below_T6'])
    if survive:
        k = KI['histories'][key]['capture_freezeout']
        lo = k['s1_Q1e+14']['one_minus_eta_final']; hi = k['s0.01_Q1e+20']['one_minus_eta_final']
        mid = k['s1_Q1e+20']['one_minus_eta_final']
        g = KI['histories'][key]['growth']
        Qlo = min(v[-1]['Q_coag'] for v in g.values() if v); Qhi = max(v[-1]['Q_ostwald'] for v in g.values() if v)
        p = 0.0
        etas = (1 - hi, 1 - mid, 1 - lo)
        EQ = (eq_of_Q(Qlo), eq_of_Q(Qhi))
        Ms = [e * q + (1 - e) + 2 * p for e in etas for q in EQ]
        row.update(eta_final_range=[1 - hi, 1 - lo], eta_final_mid=1 - mid, mean_Q_range=[Qlo, Qhi], EQ_range=list(EQ), p=p,
                   M_per_charge_range=[min(Ms), max(Ms)], dM_over_ref_range=[min(Ms) / REF - 1, max(Ms) / REF - 1],
                   frozen_test_term_reference_form=[0.0667805 * lo, 0.0667805 * hi])
    else:
        p = pgas
        M = 1 + 2 * p
        row.update(eta_final_range=[0.0, 0.0], note='condensed phase evaporates between T_ev and T_s; re-nucleation barrier/T >= 8e4; '
                   'no BEC (fugacity 0.08-0.3); relic = free X + gas-only pairs (kin_solve f)', p=p, M_per_charge_range=[M, M],
                   dM_over_ref_range=[M / REF - 1, M / REF - 1], frozen_test_term_reference_form=[0.0667805 + 2.13356 * p] * 2)
    cases.append(row)
out['cases'] = cases
# spinodal passage and other numbers
sp = KI['spinodal_passage']
out['spinodal'] = dict(n_sp=sp['n_sp'], c1=sp['c1'], c2=sp['c2'], T_vis_GeV=sp['T_vis_GeV'], H_over_mX=sp['H_over_mX'],
                       Q_frag_range=[min(c['Q_half_cube'] for c in sp['cases']), max(c['Q_cube'] for c in sp['cases'])],
                       Gamma_over_H_at_growth=[c['Gamma_over_H'] for c in sp['cases']],
                       sudden_quench_reference=dict(n=88.0, Q_cube=2.96e5, note='fastest cube at F = 1 (m1_vacuum), not the physical regime'))
nuc = KI['nucleation_from_vapour']
out['renucleation_min_barrier_over_T'] = min(r['barrier_over_T'] for r in nuc if 'barrier_over_T' in r)
out['vapour_fugacity_range'] = [min(r['fugacity'] for r in nuc if 'fugacity' in r), max(r['fugacity'] for r in nuc if 'fugacity' in r)]
json.dump(out, open(os.path.join(RES, 'summary.json'), 'w'), indent=1, default=float)
print('Y_crit = %.4f (S/Q < %.3f); finite-Q: %s' % (Ycrit, 1 / Ycrit, {k: round(v['Y_crit'], 4) for k, v in EH['Y_crit_finite_Q'].items()}))
for rel, c in crit.items():
    print('%s: R_crit = %.2f, e_crit(band) = %s' % (rel, c['R_crit'], '%.2e' % c['e_crit_band'] if c['e_crit_band'] else '-'))
for r in cases:
    print('%-30s Y_h=%.3f S/Q=%.2f survive=%s T_ev=%s T_s=%s eta=%s <Q>=%s E/Q=%s M=%s dM/ref=%s test=%s' % (
        r['case'], r['Y_h'], r['S_per_charge'], r['survive'], r['T_ev'], r['T_s'],
        ['%.3f' % x for x in r['eta_final_range']], ['%.1e' % x for x in r.get('mean_Q_range', [])], ['%.5f' % x for x in r.get('EQ_range', [])],
        ['%.4f' % x for x in r['M_per_charge_range']], ['%+.4f' % x for x in r['dM_over_ref_range']], ['%.4f' % x for x in r['frozen_test_term_reference_form']]))
print('E/Q at <Q>:', {k: round(v, 5) for k, v in out['E_over_Q_of_mean_Q'].items()})
print('spinodal:', out['spinodal'])
print('renucleation min barrier/T: %.2e; fugacity %s' % (out['renucleation_min_barrier_over_T'], out['vapour_fugacity_range']))
