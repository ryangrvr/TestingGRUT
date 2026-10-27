#!/usr/bin/env python3
"""st_clocks_skeptic: collect per-claim skeptic verdicts from the own results files into results/skeptic_verdicts.json."""
import json
from pathlib import Path
R = Path(__file__).resolve().parent.parent / 'results'
w = json.loads((R / 'widths_skeptic.json').read_text())
b = json.loads((R / 'background_skeptic.json').read_text())['histories']
o = json.loads((R / 'ordering_skeptic.json').read_text())['histories']
p = json.loads((R / 'population_skeptic.json').read_text())
c = json.loads((R / 'criterion_skeptic.json').read_text())
s = w['scalar']['theta_minus_pi_over_2']
unit = lambda tag: [r for r in b[tag]['sources']['YL']['at_crossings'] if r['kind'] == 'amplitude_unity' and r['seed'] == 0.03][0]
V = {
 'W1': dict(refuted=False, own=s['Gamma_L_GeV'], rel_vs_codex=w['rel_vs_codex']['Gamma_L'], method='real-field vertex sum (sympy), recursive phase space'),
 'W2': dict(refuted=False, own=[w['vector']['e1e-07']['Gamma_GeV'], w['vector']['e1e-06']['Gamma_GeV']],
            rel_vs_codex=[w['rel_vs_codex']['Gamma_A_e1e_7'], w['rel_vs_codex']['Gamma_A_e1e_6']], method='explicit polarisation sum + angular quadrature'),
 'W3': dict(refuted=False, own=w['vector_median_ratio']),
 'M1': dict(refuted=False, own=[w['medians']['YL'], w['medians']['A_e1e_6'], w['medians']['A_e1e_7']]),
 'C1': dict(refuted=False, own=[b['quiet3']['C_own'], b['quiet5']['C_own']], rel_vs_st=[b['quiet3']['C_rel_vs_st_clocks'], b['quiet5']['C_rel_vs_st_clocks']],
            rel_vs_codex=[b['quiet3']['C_rel_vs_codex'], b['quiet5']['C_rel_vs_codex']], tau_star=[b['quiet3']['tau_star_own'], b['quiet5']['tau_star_own']]),
 'A1': dict(refuted=False, median_a_YL=[b[t]['sources']['YL']['number_quantiles'][2]['a_own'] for t in ('quiet3', 'quiet5')]),
 'F1': dict(refuted=False, own=[unit(t)['number_fraction_own'] for t in ('quiet3', 'quiet5')]),
 'F2': dict(refuted=False, own=[unit(t)['energy_fraction_own'] for t in ('quiet3', 'quiet5')],
            rel_vs_codex=[unit(t)['energy_rel_vs_codex'] for t in ('quiet3', 'quiet5')],
            codex_data_pchip=[unit(t)['codex_data_pchip_over_total_own'] for t in ('quiet3', 'quiet5')],
            codex_data_spline=[unit(t)['codex_data_spline_over_total_own'] for t in ('quiet3', 'quiet5')]),
 'F3': dict(refuted=False, late_limit_rel=[b[t]['sources']['YL']['late_limit_rel'] for t in ('quiet3', 'quiet5')]),
 'Q1': dict(refuted=False, Q_2p13=[o[t]['Q_cube']['k=2.13'] for t in ('quiet3', 'quiet5')], Q_at_crossing_k_quiet5=o['quiet5']['Q_cube'].get('k=2')),
 'T1': dict(refuted=False, scope_overstated=True,
            orbit_criterion_0p03=[c[t]['seeds']['seed_0.03']['orbit_amplitude_criterion_first_sample'] for t in ('quiet3', 'quiet5')],
            orbit_measure_at_codex_samples_0p03=[c[t]['seeds']['seed_0.03']['orbit_amplitude_measure_at_codex_sample'] for t in ('quiet3', 'quiet5')],
            orbit_measure_at_codex_samples_1e6=[c[t]['seeds']['seed_1e-06']['orbit_amplitude_measure_at_codex_sample'] for t in ('quiet3', 'quiet5')],
            orbit_criterion_1e6_running_max_800=[c[t]['seeds']['seed_1e-06']['orbit_amplitude_running_max_at_800'] for t in ('quiet3', 'quiet5')],
            seed_1e9_codex_criterion_E_tau=[o[t]['seed_1e9']['extrapolation']['fit_500_800']['tau_cross'] for t in ('quiet3', 'quiet5')],
            cond_assisted_by_unit_from_0_1_5=[[o[t]['condensate_assisted']['from_%g_to_%g' % (lo, o[t]['unit_sample'])]['codex_data'] for lo in (0, 1, 5)] for t in ('quiet3', 'quiet5')]),
 'N1': dict(refuted=False, NYL=[p[t]['NYL_late'] for t in ('quiet3', 'quiet5')], R_Y=[p[t]['R_Y_homogeneous'] for t in ('quiet3', 'quiet5')]),
}
(R / 'skeptic_verdicts.json').write_text(json.dumps(V, indent=2) + '\n')
print(json.dumps(V, indent=1))
