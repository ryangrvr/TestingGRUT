#!/usr/bin/env python3
"""rp_critic: collect this track's own numbers next to Codex's and the three tracks' (quoted) into results/summary.json."""
import json, glob, os, sys
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
out = dict(track='rp_critic', packet='e160356', note='own code; Codex and track values quoted [Q]')
codex_six = {'b2_g0.03': (-0.19127890, 0.3326972, 0.0035042871), 'b2_g0.1': (-0.19876909, 0.5595552, 0.0035064833),
             'b2_g0.3': (-0.21868353, 0.7842706, 0.0035116238), 'b10_g0.03': (-0.18147886, 0.1847296, 1.7397750),
             'b10_g0.1': (-0.18152782, 0.3250488, 1.7071335), 'b10_g0.3': (-0.17971470, 0.5181728, 1.6142281)}
rp_line = {'b3_g0.5': 1.0454, 'b5_g1': 1.0474, 'b10_g3': 1.0194, 'b3_g0.9': 1.1482, 'b12_g0.1_NYL': 2.22325}
reh = {}
for f in sorted([f for f in glob.glob(os.path.join(R, 'rc_b*.json')) if os.path.basename(f) != 'rc_boundary.json']):
    d = json.load(open(f)); m = d['late_mean']
    key = os.path.basename(f)[3:-5]
    reh[key] = dict(b=d['b'], gamma=d['gamma'], tau_end=d['tau_end'], eps=m['eps'], NX=m['NX'], NYL=m['NYL'], NYH=m['NYH'],
                    yield_over_YQ=abs(m['Y_over_YQ']), dilution=m['dilution'], charge_ratio=d['charge_ratio_to_instant'],
                    a_f=d['final']['a'], T_f_GeV=d['final']['T_GeV'], timing=d['timing'], checks=d['checks'],
                    windows=d.get('windows'), prep_residual=d['prep_residual'], R0=d['R0'])
out['reheating_own'] = reh
out['reheating_codex_six_quoted'] = {k: dict(eps=v[0], yield_over_YQ=v[1], NYL=v[2]) for k, v in codex_six.items()}
out['reheating_rp_line_quoted'] = rp_line
out['boundary_own'] = json.load(open(os.path.join(R, 'rc_boundary.json')))
out['source_own'] = json.load(open(os.path.join(R, 'rc_source.json')))
out['arith_own'] = json.load(open(os.path.join(R, 'rc_arith.json')))
json.dump(out, open(os.path.join(R, 'summary.json'), 'w'), indent=1)
for k, v in reh.items():
    print(f"{k:16s} tau_end {v['tau_end']:6g} eps {v['eps']:.5f} NX {v['NX']:.2f} NYL {v['NYL']:.5g} Y/YQ {v['yield_over_YQ']:.5f} D {v['dilution']:.3f} a_f {v['a_f']:.3f}")
