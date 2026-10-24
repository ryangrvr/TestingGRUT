#!/usr/bin/env python3
import hashlib,json
import numpy as np
from streaming import ROOT,PARENT,C,H0,OR,OM,travel,mass_accounting


def main():
    d=json.loads((ROOT/'STREAMING_RESULTS.json').read_text());runs=d['runs'];checks=[]
    def check(name,ok,value):checks.append(dict(name=name,passed=bool(ok),evidence=value))
    errors=[]
    for af,u,end in [(1e-11,.707,OR/OM),(1e-6,2.532,OR/OM),(1e-7,20.,1e-3)]:
        q=af*u;exact=C*q/(H0*np.sqrt(OR))*(np.arcsinh(end/q)-np.arcsinh(af/q))
        errors.append(abs(travel(af,u,end,True)/exact-1))
    check('Collisionless quadrature agrees with independent radiation-era closed form',max(errors)<1e-10,errors)
    # Stable kinetic-energy arithmetic must retain momenta below E-1 roundoff.
    p=1e-12;a=1e-7;af=1.;q=a*p
    expected=q/af;k=(q/af)**2/(np.sqrt(1+(q/af)**2)+1)
    recovered=np.sqrt(k*(k+2))
    check('Comoving-momentum receiver preserves extremely cold momentum without subtracting rest energy',abs(recovered/expected-1)<1e-14,float(recovered/expected-1))
    ne=max(r['maximum_number_absolute_error'] for r in runs);qe=max(r['maximum_charge_absolute_error'] for r in runs)
    ee=max(abs(s['energy_ledger_relative_error']) for r in runs for s in r['snapshots'])
    check('All production runs conserve number, charge and energy including core exchange and expansion',ne<1e-10 and qe<1e-10 and ee<1e-8,dict(number=ne,charge=qe,energy=ee))
    spread={}
    for label in ['T20MeV','T1keV','z1e6']:
        group=sorted([r for r in runs if r['label']==label],key=lambda r:r['nodes'])
        vals=[next(s for s in r['snapshots'] if s['name']=='z5')['converted_post_escape_mean_Mpc'] for r in group]
        spread[label]=max(vals)/min(vals)-1
    check('Receiver horizon is stable across all three grid/time/source-bin refinements',max(spread.values())<.002,spread)
    hashes={r['rate_table_sha256'] for r in runs};expected=hashlib.sha256((PARENT/'RATE_RESULTS.json').read_bytes()).hexdigest()
    check('Every production run uses the same pinned ground-object rate table',hashes=={expected},len(hashes))
    first=json.loads((ROOT/'FIRST_CONVERSION_RESULTS.json').read_text())['rows']
    ge=max(r['optical_quadrature_32_vs_64_max_relative'] for r in first)
    check('Independent first-event optical-depth quadrature converges from order 32 to 64',ge<1e-9,ge)
    m=d['mass'];e=m['omega_over_m'];E0=np.sqrt(1.5)
    check('Conversion energy and charge use the chemical potential, with consistent daughter mass accounting',
        abs(m['cold_ejected_K_over_m']-2*e)<1e-14 and abs(m['initial_daughter_ejected_K_over_m']-(E0-1+2*e))<1e-14
        and abs(m['remaining_rest_matter_ratio_without_retuning']+m['released_energy_fraction_of_original_calibrated_rest_matter']-1)<1e-14,m)
    # The physical result must explicitly exclude the failed path estimator.
    rejection=json.loads((ROOT/'REJECTED_PATH_DIAGNOSTIC.json').read_text())
    bounded=all(not s['active_grid_path_estimates_certified'] and
       s['converted_total_traveled_path_mean_upper_Mpc']>=s['converted_total_traveled_path_mean_lower_Mpc']
       for r in runs for s in r['snapshots'])
    check('Rejected active-grid path estimator is fenced from physical conclusions',not rejection['passed'] and bounded,rejection['active_grid_path_relative_error'])
    out=dict(scope='Same-author verification of a conditional imposed-source transport calculation, not independent review or an observational likelihood.',
        checks=checks,passed=sum(r['passed'] for r in checks),total=len(checks),rejected_diagnostic=rejection,
        open_gates=['Actual formation spectrum and cosmological production history','Evolving object profiles and additional interactions','Elastic-angle transport and a mixed-species perturbation calculation','Lyman-alpha/BBN likelihood and independent replication'])
    (ROOT/'VERIFICATION_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for r in checks:print('PASS' if r['passed'] else 'FAIL',r['name'],r['evidence'])
    assert all(r['passed'] for r in checks)


if __name__=='__main__':main()
