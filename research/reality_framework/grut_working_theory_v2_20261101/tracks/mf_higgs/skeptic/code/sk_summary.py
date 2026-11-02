#!/usr/bin/env python3
"""mf_higgs SKEPTIC: collect the verdicts and the recomputed numbers -> results/sk_summary.json"""
import sys
sys.dont_write_bytecode = True
import json

import sk_common as C

R = {n: json.load(open(C.RES / f'sk_{n}.json')) for n in ('infl', 'postD', 'plasma', 'medium', 'g6', 'windows', 'ns', 'upper')}
inf, pd = R['infl'], R['postD']
inst = inf['instantaneous']
summary = dict(
    track='mf_higgs', skeptic='mf_higgs_skeptic (Claude, multi-field pass; GitHub issue #13 context)',
    labels='[C] own code (skeptic/code), [D] derived, [E] estimate/argued criterion, [A] argued, [Q] quoted, [R] literature',
    model_choice_note='Every benchmark, history, criterion, closure and renormalisation condition is a model choice; nothing here '
                      'validates or excludes anything; data comparisons are retrospective.',
    checks={
        'H1': dict(refuted=False, recomputed={'xi_min tree A (172.57/172.86/174.6)': [inf['rows'][m]['A xi_min'] for m in ('172.57', '172.86', '174.6')],
                                              'xi_min flat 2-loop B': [inf['rows'][m]['B xi_min'] for m in ('172.57', '172.86', '174.6')],
                                              'xi_min 1-loop dS C': [inf['rows'][m]['C xi_min'] for m in ('172.57', '172.86', '174.6')],
                                              'cap const H': inf['caps_const_H']}),
        'H2': dict(refuted=True, recomputed=dict(caps_with_endpoint=inst['caps_with_endpoint'], xi_min_with_endpoint=inst['xi_min_with_endpoint'],
                                                 rows=inst['rows'], gaussian_tail_at_eps1=inf['gaussian_tail_own_modes'])),
        'H6_H7': dict(refuted=False, recomputed={k: {kk: v.get(kk) for kk in ('Mt_crit', 'Mt_crit_band', 'Mt_crit_N60', 'Mt_crit als=0.1171',
                                                                            'Mt_crit als=0.1189', 'Mt_crit Mh=125.09', 'Mt_crit Mh=125.31',
                                                                            'D_at', 'deadlines', 'D(172.57, als=0.1171)')}
                                                 for k, v in pd.items() if k in ('EXP', 'MC_TR', 'MC_P0@1.3e11', 'MC_P0@2.4e11')},
                      controls=pd.get('controls_EXP')),
        'H5_H8': dict(refuted=False, recomputed=dict(plasma=R['plasma'], medium_factor=R['medium'])),
        'H4': dict(refuted=False, recomputed=R['windows']),
        'H9': dict(refuted=False, recomputed=dict(upper=R['upper'], lower_edge_grav=pd.get('lower_edge_grav'))),
        'H10': dict(refuted=False, recomputed=dict(virial=R['g6']['virial'], floquet=R['g6']['floquet'])),
        'H3': dict(refuted=False, recomputed=R['ns']),
    })
C.jdump(summary, C.RES / 'sk_summary.json')
print('written', C.RES / 'sk_summary.json')
