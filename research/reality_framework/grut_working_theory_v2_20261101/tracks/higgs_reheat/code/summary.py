"""Collect the headline numbers of track higgs_reheat into ../results/summary.json.  Own code; reads only our results."""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
L = lambda f: json.load(open(os.path.join(RES, f)))
inf, bg, md, ed, th, ct = L('inflation.json'), L('background.json'), L('modes.json'), L('edges.json'), L('thermal.json'), \
    L('controls.json')

S = {'labels': '[C] computed here; [E] estimate/criterion choice; [A] argued; [Q] quoted (abstract level)'}

# 1. SM potential
S['1_SM_potential [C]'] = {f'Mt={Mt}': {k: inf[f'Mt={Mt}']['H_free-X 1.29e11']['flat'][k] for k in ('Lambda_I', 'h_max')}
                           for Mt in (172.33, 172.6, 172.87)}
S['1_xi_running [C]'] = {'xi(M_t)=0.1 gives xi(H_i)': inf['Mt=172.6']['H_free-X 1.29e11']['xi_at_mu'],
                         'factor (xi_H - 1/6)/(xi_EW - 1/6)': ed['Mt=172.6 H=1.29e+11']['xi_H_over_xi_EW_running_factor']}

# 2. inflation
rows = {}
for Mt in (172.33, 172.6, 172.87):
    for hl in ('H_free-X 1.29e11', 'H_Q-matter 1.51e11'):
        r = inf[f'Mt={Mt}'][hl]
        b = r['xiEW=0.1']
        rows[f'Mt={Mt} {hl}'] = dict(
            h_max_over_H=r['h_max_over_H'],
            m2_curv_over_H2_one_loop=b['m2_curv_over_H2_one_loop'],
            dS_barrier_h_b_over_H=b['cR=1.0,loop=True']['h_b_over_H'],
            B_HM_one_loop=b['cR=1.0,loop=True']['B_HM'], B_HM_tree_RGI=b['cR=1.0,loop=False']['B_HM'],
            log10_bubbles_N60_one_loop=b['cR=1.0,loop=True']['log10_N_bubbles']['60'],
            xi_EW_min_N60_one_loop=r['xi_min']['60']['xi(M_t), one loop dS'],
            xi_EW_min_N60_tree=r['xi_min']['60']['xi(M_t), tree RGI'],
            xi_EW_min_N50_one_loop=r['xi_min']['50']['xi(M_t), one loop dS'],
            xi_H_min_N60_one_loop=r['xi_min']['60']['xi(H_i), one loop dS'],
            sigma2_per_component_light=r['variance_units_H2']['stochastic_light_3H^4/(8pi^2 m^2)'],
            sigma2_per_component_dS_superhorizon=r['variance_units_H2']['dS_superhorizon_eps=1.0'],
            end_of_inflation_log10_patches_above_flat_hmax_N60=r['tails_at_end_of_inflation'][
                'Fokker-Planck full V | flat h_max']['log10_expected_patches']['60'])
S['2_inflation [C]'] = rows

# 3a background
S['3a_background [C]'] = dict(R0=bg['R0'], first_R_negative_tau=bg['first_R_negative_tau'],
                              R_min_first=bg['R_minima'][0], R_min_second=bg['R_minima'][1],
                              radiation_equality_tau=bg['radiation_equality_tau'],
                              engine_vs_scipy=bg['engine_vs_scipy_max_rel_dev'],
                              T_max=[(t['H_GeV'], t['T_max_GeV'], t['tau_Tmax']) for t in bg['temperatures']],
                              T_at_tau_0p05=[t['T_at_tau']['0.05'] for t in bg['temperatures']])


# 3b gravitational-only window
def mf_cross(rows):
    xs = sorted((float(k), 2 * v['max_sigma2'] ** 0.5 / v['h_max_over_H']) for k, v in rows.items())
    for (x1, y1), (x2, y2) in zip(xs, xs[1:]):
        if x1 > 1 and y1 < 1 <= y2:
            return float(np.exp(np.interp(0.0, [np.log(y1), np.log(y2)], [np.log(x1), np.log(x2)])))
    return None


win = {}
for k, e in ed.items():
    fac = e['xi_H_over_xi_EW_running_factor']
    gk = 'grav_only ' + k.replace('H=', 'H=')
    mf = mf_cross(md[gk]) if gk in md else None
    win[k] = dict(benchmark_xi_H=e['benchmark']['xi_H'], benchmark_D60=e['benchmark']['D']['60']['D'],
                  benchmark_D60_cH05=e['benchmark']['D_cH']['cH=0.5']['60']['D'],
                  benchmark_D60_cH2=e['benchmark']['D_cH']['cH=2.0']['60']['D'],
                  benchmark_tau_c=e['benchmark']['D']['60'].get('tau_c'),
                  N60=e['N=60'], N50=e['N=50'],
                  upper_mean_field_2sigma_eq_hmax_xi_H=mf,
                  upper_mean_field_xi_EW=(1 / 6 + (mf - 1 / 6) / fac if mf else None),
                  upper_cH_bracket_xi_H=(e['N=60 cH=2.0 (from scan)']['smallest_unsafe_xi_H_above_1'],
                                         e['N=60 cH=0.5 (from scan)']['smallest_unsafe_xi_H_above_1']))
S['3b_gravitational_only_window [C][E]'] = win

# 3c plasma
pl = {}
for pk in ('plasma Mt=172.6 H=1.29e+11', 'plasma Mt=172.87 H=1.51e+11'):
    rows = md[pk]
    pl[pk] = {kk: dict(D60=v['danger']['best']['60']['D'], runaway=v['runaway'],
                       maxsigL2_ratio=v['max_sigma2_L'] / v['sigma2_L_0']) for kk, v in rows.items()
              if kk.startswith('xi=0.088') or kk.startswith('xi=100.0') or kk.startswith('xi=0.05')}
S['3c_plasma [C][A]'] = pl
S['3c_prompt_large_xi [C]'] = md.get('prompt plasma, large xi, Mt=172.87 H=1.51e11, to tau=2.5')

# 3d thermal
S['3d_thermal [C][E]'] = {k: dict(S3_over_T={T: b['S3_over_T'] for T, b in v['bounce'].items()},
                                  kappa_T={T: b['kappa_T'] for T, b in v['bounce'].items()},
                                  barrier_h_over_T={T: b['barrier_h_over_T'] for T, b in v['bounce'].items()},
                                  log10_bubbles=v['log10_bubbles_in_past_light_cone']) for k, v in th.items()}
S['controls [C]'] = {k: {kk: vv['D60'] for kk, vv in v.items()} for k, v in ct.items()}


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, float) and not np.isfinite(o):
        return str(o)
    return o


S = clean(S)
json.dump(S, open(os.path.join(RES, 'summary.json'), 'w'), indent=1, default=float)
print(json.dumps(S, indent=1, default=float)[:12000])
