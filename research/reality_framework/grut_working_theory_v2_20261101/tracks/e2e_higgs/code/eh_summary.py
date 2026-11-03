#!/usr/bin/env python3
"""e2e_higgs: collect the headline numbers of post.json, xiscan.json, inflation.json, control_tend.json
-> results/summary.json.  Model choices throughout; [C] on the argued criteria [E]."""
import sys
sys.dont_write_bytecode = True
import json
import eh_common as E

P = json.load(open(E.RES / 'post.json'))
X = json.load(open(E.RES / 'xiscan.json'))
I = json.load(open(E.RES / 'inflation.json'))
T = json.load(open(E.RES / 'control_tend.json'))
r = lambda x, n=3: round(x, n) if isinstance(x, float) else x
S = dict(note='e2e_higgs headline numbers; xi_H = 0.1 read as xi(M_t); PDG M_t = 172.57 +- 0.29, M_h = 125.20 +- 0.11, '
              'alpha_s = 0.1180 +- 0.0009 [R]; every choice a model choice; nothing validated or excluded.',
         post_grav={}, post_media={}, deadlines={}, R_fluct={}, xi_lower_edge_post_grav={}, inflation={})
for k, o in P.items():
    S['post_grav'][k] = dict(Mt_crit=r(o['Mt_crit_grav']), band=[r(x) for x in o['Mt_crit_band']], N=r(o['N_central']),
                             D_172p57=o['D_grav_coarse']['172.57'], D_172p86=o['D_grav_coarse']['172.86'],
                             Mt_crit_by_N={a: r(b) for a, b in o['Mt_crit_by_N'].items()},
                             alpha_s={a: dict(Mt_crit=r(b['Mt_crit']), D_central=b['D_central_Mt']) for a, b in o.get('alpha_s', {}).items()},
                             Mh={a: r(b) for a, b in o.get('Mh', {}).items()})
    S['post_media'][k] = {a: dict(D_max_174to174p6=b['D_max_174_to_174p6'], safe_factor=b['smallest_safe_factor_at_174p6'])
                          for a, b in o['media'].items()}
    S['deadlines'][k] = o['deadline']
    S['R_fluct'][k] = {a: dict(Mt_crit=r(b['Mt_crit']), D_172p86=b['D_172p86']) for a, b in o['R_rms_variants'].items()}
for k, o in X.items():
    S['xi_lower_edge_post_grav'][k] = {a: r(b, 4) for a, b in o['xi_EW_lower_edge_grav'].items()}
for k, o in I.items():
    S['inflation'][k] = dict(N0=r(o['N0']), caps={a: r(b) for a, b in o['Mt_cap_xi0.1'].items()},
                             xi_min_const_tree={m: r(o['rows'][m]['tree const xi_min'], 4) for m in ('172.57', '172.86')},
                             xi_min_const_loop={m: r(o['rows'][m]['loop const xi_min'], 4) for m in ('172.57', '172.86')},
                             xi_min_inst_exact={m: r(o['rows'][m]['tree inst_exact xi_min'], 4) for m in ('172.57', '172.86')})
S['control_t_end'] = T
S['media_low_xi'] = 'xi_EW = 0.05 runs (1.3e11): D_grav(172.57) = +45 (KIN, LAT_P0); with N+SM media D(174.6) = -inf / -6554'
E.jdump(S, E.RES / 'summary.json')
print(json.dumps(S['post_grav']['KIN@1.3e+11'], indent=0)[:800])
