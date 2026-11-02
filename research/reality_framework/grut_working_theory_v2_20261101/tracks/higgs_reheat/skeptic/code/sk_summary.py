"""Collects the skeptic's headline numbers into ../results/sk_summary.json (no new computation)."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, '..', 'results')
ld = lambda f: json.load(open(os.path.join(R, f)))
m = ld('sk_modes.json')
mc = ld('modes_canon_benchmark.json')
edges = {}
for f in ('sk_edges_0_1.json', 'sk_edges_2.json', 'sk_edges_3.json', 'sk_edges_4.json'):
    edges.update(ld(f)['lower_edges_N=60'])
band = {}
for f in ('sk_edges_band_0_1.json', 'sk_edges_band_2_3.json'):
    band.update(ld(f))
inf = ld('sk_infl.json')
th = ld('sk_thermal.json')
bb = ld('sk_bubbles.json')
dl = ld('track_code_plasma_deadline.json')
bench = m['benchmark_xiEW=0.1']
crit = m['criterion_sensitivity_xiEW=0.1']
S = dict(
    benchmark_D60={k: dict(track_norm_own_code=v['track cH=1.0']['D'], canonical_own_code=v['canon cH=1.0']['D'],
                           canonical_track_code_copy=mc[k]['D60'],
                           canonical_cH_half_2_track_code_copy=[mc[k]['D60_cH05'], mc[k]['D60_cH2']],
                           canonical_band_over_criteria=[min(vv for kk, vv in crit[k].items()
                                                             if kk.startswith('canon') and not isinstance(vv, dict)),
                                                         max(vv for kk, vv in crit[k].items()
                                                             if kk.startswith('canon') and not isinstance(vv, dict))])
                   for k, v in bench.items()},
    lower_edge_xi_EW_N60={k: dict(track_norm=v['track']['xi_EW_lower'], canonical=v['canon']['xi_EW_lower'],
                                  canonical_band=[band[k][x]['xi_EW_lower'] for x in ('f=0.5', 'gauss R=1/k_c', 'f=2')]
                                  if k in band else None) for k, v in edges.items()},
    plasma_deadline_D_of_tau_th=dict(own_code_sharp={k: {n: crit[k][f'{n} D(tau_th)'] for n in ('track', 'canon')}
                                                     for k in ('Mt=172.87 H=1.29e+11', 'Mt=172.87 H=1.51e+11')},
                                     track_code_tanh_ramp_worst_corner=dl),
    curvature_mass={k: {h: dict(track_one_loop=v[h]['m2_track_one_loop'], tree_h0_units=v[h]['m2_tree_in_h0_units_12xi_e2Gamma'],
                                genuine_loop=v[h]['m2_genuine_one_loop_shift'], canonical=v[h]['m2_canonical_with_loop'],
                                var_exact_canonical=v[h]['var_superhorizon_canonical_phi'],
                                var_exact_h0=v[h]['var_superhorizon_in_h0_units'], var_track=v[h]['var_track_value_check_m2_0937'])
                            for h in v} for k, v in inf.items()},
    xi_running={k: {h: [v[h]['xi_factor_indep_1loop_Buttazzo'], v[h]['xi_factor_track']] for h in v} for k, v in inf.items()},
    HM_tree={k: {h: [v[h]['HM_tree_indep']['B'], v[h]['xi_EW_min_tree_indep_N60']] for h in v} for k, v in inf.items()},
    thermal_S3_over_T={k: {T: [r['track']['S3_over_T'], r['canon']['S3_over_T']] for T, r in v.items()}
                       for k, v in th.items() if k.startswith('Mt')},
    thermal_bubbles_log10=bb,
    history_sensitivity={k: v['canon']['D'] for k, v in m['history_sensitivity'].items()})
json.dump(S, open(os.path.join(R, 'sk_summary.json'), 'w'), indent=1, default=float)
print(json.dumps(S, indent=1, default=float)[:6000])
