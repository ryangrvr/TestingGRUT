#!/usr/bin/env python3
"""Assemble evidence, convergence differences and inherited observable outputs."""
from pathlib import Path
import json, math
import numpy as np

ROOT=Path(__file__).resolve().parent
PARENT=ROOT.parent/'grut_candidate_1_20261018/code'


def main():
    case_names=['control','small','stress','time','coarse','space','domain']
    runs={name:json.loads((ROOT/f'evolution_{name}.json').read_text()) for name in case_names}
    convergence={}
    for name in ['time','coarse','space','domain']:
        a=runs['stress']['trajectory']; b=runs[name]['trajectory']
        assert [x['tau'] for x in a]==[x['tau'] for x in b]
        convergence[name]={
            'max_core_Q_fraction_r8_difference':max(abs(x['core_Q_fraction']['8']-y['core_Q_fraction']['8']) for x,y in zip(a,b)),
            'max_core_Q_fraction_r12_difference':max(abs(x['core_Q_fraction']['12']-y['core_Q_fraction']['12']) for x,y in zip(a,b)),
            'max_profile_error_difference':max(abs(x['phase_aligned_profile_error']-y['phase_aligned_profile_error']) for x,y in zip(a,b)),
            'max_central_S_difference':max(abs(x['central_S']-y['central_S']) for x,y in zip(a,b)),
            'max_quadrupole_norm_difference':max(abs(x['quadrupole_norm']-y['quadrupole_norm']) for x,y in zip(a,b)),
        }
    old=json.loads((PARENT/'results.json').read_text())
    probe=json.loads((PARENT/'reporter_results.json').read_text())
    direct=[]
    for p in probe['phase_results']:
        k=p['k_over_gv'];delta=p['phase_delta_s_wave']
        direct.append({'momentum_over_gv':k,'delta0':delta,
                       'S0_phase':2*delta,'s_wave_partial_cross_section_times_gv_squared':4*math.pi*math.sin(delta)**2/k**2,
                       'scope':'inherited fixed-background classical elastic partial wave; not total plane-wave cross section'})
    omega0=(.1)**.25
    summaries=[]
    for name in case_names:
        run=runs[name]; s=run['summary']; row=run['trajectory'][-1]
        summaries.append({**s,'final_core_Q_fraction_r8':row['core_Q_fraction']['8'],
                           'final_core_Q_fraction_r12':row['core_Q_fraction']['12'],
                           'final_profile_error':row['phase_aligned_profile_error'],
                           'final_quadrupole_norm':row['quadrupole_norm'],
                           'max_absolute_E_drift':s['initial_E_grid']*s['max_E_relative_drift'],
                           'max_local_global_gradient_energy_difference':max(abs(x['spectral_local_gradient_total_difference']) for x in run['trajectory']),
                           'free_escape_energy_gap':399.4384624223012-s['initial_E_grid'],
                           'conditional_nonfree_charge_lower_bound':(399.4384624223012-s['initial_E_grid'])/(1-omega0)})
    branches=json.loads((ROOT/'branch_results.json').read_text())
    margins=[x['delta_Ebar'] for k in ['all_resolved_two_lump_pairs','all_resolved_lump_plus_free','all_resolved_branch_combinations'] for x in branches[k]]
    result={'scope':'same-law leading classical scalar stability follow-up',
            'evolution_summaries':summaries,'matched_run_differences':convergence,
            'direct_scattering_observables':direct,
            'inherited_weak_gravity':old['checks'],
            'minimum_resolved_fragmentation_excess':min(margins),
            'all_declared_evolution_validity_checks_pass':all(x['numerical_validity_pass'] for x in summaries),
            'all_declared_operating_diagnostics_pass':all(x['operating_diagnostic_pass'] for x in summaries),
            'caveats':['finite time and selected preparations','periodic Fourier discretization',
                       'sampled/found branch envelope only','no quantum lifetime or full Einstein evolution',
                       'reporter source/detector and empirical relevance unconstructed']}
    (ROOT/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_valid':result['all_declared_evolution_validity_checks_pass'],
                      'all_core_tests_pass':result['all_declared_operating_diagnostics_pass'],
                      'minimum_resolved_fragmentation_excess':result['minimum_resolved_fragmentation_excess'],
                      'matched_run_differences':convergence},indent=2))


if __name__=='__main__':main()
