#!/usr/bin/env python3
"""Skeptic of inf_higgs: collect the headline numbers of skeptic/results/*.json into sk_summary.json."""
import sys
sys.dont_write_bytecode = True
import json

import sk_common as C


def L(n):
    return json.load(open(C.RES / n))


def interp_zero(pts):
    pts = sorted(pts)
    for (a, da), (b, db) in zip(pts[:-1], pts[1:]):
        if da < 0 <= db:
            return a + (0 - da) * (b - a) / (db - da)
    return None


def main():
    S = {}
    b = L('sk_bench.json')
    S['regression_D60'] = {k: dict(D_172_57=v['D60']['172.57'], D_172_86=v['D60']['172.86'], Mt_crit_60=v['Mt_crit_60'],
                                   Mt_crit_Nobs=v['Mt_crit_Nobs'], tau_c_172_86=v['tau_c60']['172.86']) for k, v in b.items()}
    S['track_values_for_comparison'] = {
        'HT6t|MD': [-365, -23.8, 172.917, 172.936], 'T1|MD': [-371, -23.9, 172.917, 172.936],
        'E1|MD': [-371, -24.0, 172.917, 172.936], 'HT6t|G05': [-312, -6.2, 172.875, 172.884],
        'M1r|G05': [-310, -5.1, 172.874, None]}
    a = L('sk_aterm.json')
    S['Aterm_gravitational_only'] = {k: dict(D_172_57=v['D60']['172.57'], D_172_86=v['D60']['172.86'], Mt_crit=v['Mt_crit_60'],
                                             dm2_S1_over_H2=v['m2_extra_at_S1'], kappa_A=v['kappa_A_median'])
                                     for k, v in a.items()}
    ht = [(0.0, b['HT6t|MD|iH0']['D60']['172.86'])] + [(v['A_TeV'], v['D60']['172.86']) for k, v in a.items()
                                                       if v['name'] == 'HT6t' and v['A_TeV'] >= 0]
    t1 = [(0.0, b['T1|MD|iH0']['D60']['172.86'])] + [(-v['A_TeV'], v['D60']['172.86']) for k, v in a.items()
                                                     if v['name'] == 'T1' and v['A_TeV'] < 0]
    S['Aterm_D_bound_at_172_86_TeV'] = dict(HT6t=interp_zero(ht), T1=interp_zero(t1),
                                            note='destabilising sign; gravitational only; D60 = 0 at M_t = 172.86')
    S['Aterm_with_S1_plasma'] = {k: dict(Mt_crit=v['Mt_crit'], D60=v['D60']) for k, v in L('sk_aterm_plasma.json').items()}
    p = L('sk_plasmaD.json')
    S['plasmaD_crit_extended_grid'] = {f'{n}|{r}': dict(track_grid=v['hmult=20']['Mt_crit'], extended=v['hmult=1000']['Mt_crit'],
                                                         D_174_6_extended=v['hmult=1000']['D60'][-1])
                                       for n, o in p.items() for r, v in o.items()}
    S['inflation_HN'] = L('sk_infl.json')
    S['induced_Rven_mass'] = L('sk_induced.json')
    S['sensitivities'] = L('sk_sens.json')
    S['inputs'] = L('sk_inputs.json')
    S['misc'] = L('sk_misc.json')
    S['repro'] = L('sk_repro.json')
    S['xcheck'] = L('sk_xcheck.json')
    C.jdump(S, C.RES / 'sk_summary.json')
    print(json.dumps(S['Aterm_D_bound_at_172_86_TeV']))
    print({k: (v['track_grid'], v['extended']) for k, v in S['plasmaD_crit_extended_grid'].items() if v['track_grid'] != v['extended']})


if __name__ == '__main__':
    main()
