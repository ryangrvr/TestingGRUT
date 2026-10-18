#!/usr/bin/env python3
"""Evidence checks for the bounded same-law stability follow-up.

Numerical acceptance and internal consistency are not empirical validation,
global minimization, a long lifetime, or a quantum/gravitational calculation.
"""
from pathlib import Path
import hashlib,json
import numpy as np

ROOT=Path(__file__).resolve().parent


def main():
    branch=json.loads((ROOT/'branch_results.json').read_text())
    refined=json.loads((ROOT/'branch_thresholds_refined.json').read_text())
    independent=json.loads((ROOT/'independent_comparison.json').read_text())
    newton=json.loads((ROOT/'independent_refinement_results.json').read_text())
    analysis=json.loads((ROOT/'analysis.json').read_text())
    prescription=json.loads((ROOT/'evolution_prescription.json').read_text())
    checks=[]
    def check(name,passed,evidence):checks.append({'name':name,'passed':bool(passed),'evidence':evidence})
    check('same_original_background',prescription['background_sha256']==hashlib.sha256((ROOT.parent/'grut_candidate_1_20261018/code/results_profile.npz').read_bytes()).hexdigest(),prescription['background_sha256'])
    check('evolution_source_matches_execution',prescription['source_sha256']==hashlib.sha256((ROOT/'evolve.py').read_bytes()).hexdigest(),prescription['source_sha256'])
    check('all_112_scalar_bvp_attempts_accepted',len(branch['attempts'])==112 and all(x['accepted'] and x['status']==0 for x in branch['attempts']),len(branch['attempts']))
    for name in ['turning','binding']:
        record=refined[name]
        check(name+'_root_converged',record['root_converged'],record['root_converged'])
        check(name+'_domain_energy_convergence',abs(record['Ebar_domain_difference'])<1e-6,record['Ebar_domain_difference'])
        check(name+'_domain_charge_convergence',abs(record['Qbar_domain_difference'])<1e-6,record['Qbar_domain_difference'])
    check('all_resolved_fragmentation_margins_positive',analysis['minimum_resolved_fragmentation_excess']>0,analysis['minimum_resolved_fragmentation_excess'])
    radial=branch['branch_radial_stability']
    check('radial_branch_negative_control',all(x['lowest_eigenvalues'][0]>0 for x in radial if x['Omega']==.95) and all(x['lowest_eigenvalues'][0]<0 for x in radial if x['Omega']==.97),[{k:x[k] for k in ['Omega','N','lowest_eigenvalues']} for x in radial])
    check('independent_gradient_preflight',independent['gradient_preflight_max_relative_discrepancy']<1e-6,independent['gradient_preflight_max_relative_discrepancy'])
    for x in independent['comparison']:
        a=abs(x['meshes']['800']['Ebar_error_vs_BVP']);b=abs(x['meshes']['1600']['Ebar_error_vs_BVP'])
        check('independent_energy_mesh_convergence_'+str(x['fraction']),3.5<a/b<4.5,{'coarse_error':a,'fine_error':b,'ratio':a/b})
    for x in newton['runs']:
        end=x['final']
        check('independent_FV_stationarity_'+str(x['fraction']),end['full_gradient_inf']<1e-8,end['full_gradient_inf'])
    for x in analysis['evolution_summaries']:
        name=x['case']['name']
        check(name+'_global_charge_conservation',x['max_Q_relative_drift']<1e-9,x['max_Q_relative_drift'])
        check(name+'_Hamiltonian_energy_error',x['max_E_relative_drift']<2e-3,x['max_E_relative_drift'])
        check(name+'_sampled_core_retention',x['min_core_Q_fraction_r8']>.9,x['min_core_Q_fraction_r8'])
        check(name+'_same_charge_below_free_escape_energy',x['free_escape_energy_gap']>0,x['free_escape_energy_gap'])
        check(name+'_consistent_regional_gradient_energy',x['max_local_global_gradient_energy_difference']<1e-5,x['max_local_global_gradient_energy_difference'])
    indexed={x['case']['name']:x for x in analysis['evolution_summaries']}
    ratio=indexed['stress']['max_E_relative_drift']/indexed['time']['max_E_relative_drift']
    check('time_refinement_reduces_energy_error',ratio>3,ratio)
    for name in ['time','space','domain']:
        error=analysis['matched_run_differences'][name]
        check(name+'_core_retention_comparison',error['max_core_Q_fraction_r8_difference']<.002,error['max_core_Q_fraction_r8_difference'])
        check(name+'_profile_comparison',error['max_profile_error_difference']<.005,error['max_profile_error_difference'])
    report={'scope':'bounded leading classical numerical evidence',
            'passed':all(x['passed'] for x in checks),'number_of_checks':len(checks),'checks':checks,
            'declaration':'Protocol conservation/retention gates were frozen before trajectories; extra convergence/cross-method checks were consolidated after execution, not statistical preregistration.',
            'open':['all field configurations and partitions','arbitrarily long evolution','quantum lifetime','physical apparatus','empirical calibration','nonlinear Einstein evolution']}
    (ROOT/'VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'passed':report['passed'],'checks':len(checks),'failed':[x['name'] for x in checks if not x['passed']]},indent=2))
    return 0 if report['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
