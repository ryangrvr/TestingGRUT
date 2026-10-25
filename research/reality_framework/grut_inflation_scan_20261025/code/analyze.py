"""Publish compact numerical evidence and the conditional Y-source budget."""
import csv,json,hashlib,shutil
from pathlib import Path
import numpy as np
from build import WORK
from model import YQ,minimum,parameters
from run import ROOT,write

MX=.09240896373753621
DM=.44e-9


def raw(name,legacy):
    p=WORK/name
    if not p.exists():
        WORK.mkdir(parents=True,exist_ok=True)
        shutil.move(str(ROOT/legacy),str(p))
    if (ROOT/legacy).exists():
        assert (ROOT/legacy).read_bytes()==p.read_bytes()
        (ROOT/legacy).unlink()
    return json.loads(p.read_text()),hashlib.sha256(p.read_bytes()).hexdigest()


def forced_matrix(a,u,pi,p):
    A=a*a*p['M']**2;B=2*a*p['mu']*u[2];J=1e-3*u[0]**3
    matrix=np.array([[A+B.real,B.imag],[B.imag,A-B.real]])
    yr=np.linalg.solve(matrix,-np.array([J.real,J.imag]));yp=yr[0]+1j*yr[1]
    Ap=2*a*p['M']**2;Bp=2*p['mu']*(u[2]+a*pi[2]);Jp=3e-3*u[0]**2*pi[0]
    derivative=np.array([[Ap+Bp.real,Bp.imag],[Bp.imag,Ap-Bp.real]])
    vr=np.linalg.solve(matrix,-np.array([Jp.real,Jp.imag])-derivative@yr)
    return yp,vr[0]+1j*vr[1]


def summarize_ensemble(e,eid,members):
    fields=['H_GeV','cH','delta','epsX','NX','NY','NS','Neff','EX','EY','ES','JY_approx',
            'YX_over_YQ','late_charge_drift_fraction','max_Gauss_over_initial_Ncirc',
            'max_Hamiltonian_work_error_over_initial_energy']+[f'final_state_{i}' for i in range(12)]
    for i,r in enumerate(e['rows']):
        members.append(dict(ensemble_id=eid,member=i,**{k:r[k] for k in fields}))
    return {k:v for k,v in e.items() if k!='rows'}|dict(ensemble_id=eid)


def main():
    v,vhash=raw('candidate_validation_raw.json','CANDIDATE_VALIDATION.json')
    initial,ihash=raw('initial_ensembles_raw.json','INITIAL_ENSEMBLES.json')
    members=[];cs=[];iis=[];high=[];pair_budgets=[];ychecks=[]
    for i,e in enumerate(initial['rows']):iis.append(summarize_ensemble(e,f'initial_{i}',members))
    for i,e in enumerate(v['candidate_ensembles']):cs.append(summarize_ensemble(e,f'candidate_{i}',members))
    pair_limit=.01*DM/(2*MX*YQ)
    gammaY=3*1e-6*1e12/(256*np.pi**3)
    for c in v['high_fraction_controls']:
        i=c['candidate_index'];large=summarize_ensemble(c['larger_seed_ensemble'],f'candidate_{i}_sigma003',members)
        ctrls=c['tighter_and_longer_controls'];shorts=[]
        for r in ctrls:
            # Full window moments and canonical endpoint retained; sampled histories unnecessary here.
            shorts.append({k:v for k,v in r.items() if k!='history'})
        radial={k:v for k,v in c['conditional_radial_release'].items() if k!='history'}
        high.append({k:v for k,v in c.items() if k not in ['larger_seed_ensemble','tighter_and_longer_controls','conditional_radial_release']}
                    |dict(larger_seed_ensemble=large,tighter_and_longer_controls=shorts,conditional_radial_release=radial))
        r=ctrls[0];w=r['windows']['600-1000'];old=r['windows']['300-600']
        H=r['H_GeV'];adec=np.sqrt(H/gammaY);eydec=w['JY_approx']*1e12/H*adec
        pairs=1.5*w['JY_approx']/abs(w['NX'])
        state=np.array(r['final_canonical_state']);u=state[:6:2]+1j*state[1:6:2];pi=state[6::2]+1j*state[7::2]
        a=r['a_end'];p=parameters(H);A=a*a*p['M']**2;B=2*a*p['mu']*u[2]
        yp,ypd=forced_matrix(a,u,pi,p);d=u[1]-yp;vel=pi[1]-ypd-d/a
        ey=abs(vel)**2+A*abs(d)**2+np.real(np.conj(B)*d*d)
        step=1e-5
        fd=(forced_matrix(a+step,u+step*pi,pi,p)[0]-forced_matrix(a-step,u-step*pi,pi,p)[0])/(2*step)
        direct=r['history'][-1]['EY_exc_quadratic']
        ychecks.append(dict(candidate_index=i,matrix_vs_C_energy_relative_error=float(abs(ey-direct)/direct),
                            forced_derivative_finite_difference_relative_error=float(abs(fd-ypd)/max(abs(ypd),1e-30))))
        pair_budgets.append(dict(candidate_index=i,H_GeV=H,cH=r['cH'],delta=r['delta'],
            raw_epsX=w['epsX'],YX_over_YQ=w['YX_over_YQ'],JY_approx=w['JY_approx'],
            JY_window_relative_shift=abs(w['JY_approx']-old['JY_approx'])/abs(w['JY_approx']),
            final_Y_quadratic_mass_matrix_positive=bool(A>abs(B)),
            final_Y_selfinteraction_to_quadratic_force_ratio=float(2*abs(u[1])**2/A),
            final_Y_mass_splitting_relative_parameter=float(abs(B)/A),
            a_at_assumed_Y_decay=adec,Y_excitation_energy_over_X_at_assumed_decay=eydec/w['EX'],
            extra_Y_decay_pairs_per_abs_X_charge=pairs,one_percent_pair_allowance_per_net_charge=pair_limit,
            ratio_to_one_percent_pair_allowance=pairs/pair_limit,
            extra_DM_fraction_if_yield_diluted_to_YQ_and_Y_unthermalized=2*MX*YQ*pairs/DM,
            energy_equivalent_quartic_fraction_after_Y_injection=abs(w['epsX'])/(1+eydec/w['EX'])**.75,
            energy_proxy_scope='Only a comoving-energy ratio mapped to the quartic circular-charge scale; NOT a solved mixed-plasma remainder or thermal relic.',
            assumptions='Adiabatic conserved late Y number extrapolated to perturbative Y->3X decay; zero mean Y charge; no thermalization/number changing; visible RD. Small Y mass splitting omitted. Not an evolved decay history.'))
    with (ROOT/'ENSEMBLE_MEMBERS.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(members[0]),lineterminator="\n");w.writeheader();w.writerows(members)
    write(ROOT/'ENSEMBLE_SUMMARY.json',dict(initial_ensembles=iis,candidate_ensembles=cs,
        high_fraction_controls=high,member_count=len(members),member_data='ENSEMBLE_MEMBERS.csv',
        raw_candidate_output_sha256=vhash,raw_initial_output_sha256=ihash,
        scope=v['scope'],SEM_scope='Finite random-sample error only. Not inflationary-model, spatial, numerical or relic uncertainty.'))
    write(ROOT/'Y_SOURCE_RESULTS.json',dict(pair_limit=pair_limit,Gamma_Y_3body_GeV=gammaY,rows=pair_budgets,
                                          independent_real_matrix_and_derivative_checks=ychecks))
    points=[]
    grid=json.loads((ROOT/'GRID_RESULTS.json').read_text())
    for g in grid['groups']:
        p=g['initial_point'];fg=p['vector_mass_over_H']/1e-7
        points.append(dict(H_GeV=g['H_GeV'],cH=g['cH'],gauge_Higgs_norm_over_H=fg,
            vector_mass_over_H=p['vector_mass_over_H'],
            assumed_BD_Goldstone_sqrt_power_over_H=1/(2*np.pi),
            assumed_BD_gauge_angle_sqrt_power=1/(2*np.pi*fg),
            assumed_superhorizon_AL_power_over_H2_per_k_over_H2=1/(4*np.pi**2*p['vector_mass_over_H']**2)))
    write(ROOT/'VECTOR_INITIALIZATION.json',dict(rows=points,
        assumption='Stationary inflationary scalar background, light vector, Bunch-Davies state and superhorizon Goldstone/Proca limit. Gauge-angle normalization is gauge dependent; no post-release mode mixing/occupation evolved.',
        after_release='Use full evolving Higgs norm and coupled physical scalar/vector fluctuations. These initialization data do not supply a final vector abundance.'))
    coarse=list(csv.DictReader((ROOT/'SCAN_RESULTS.csv').open()))
    refine=list(csv.DictReader((ROOT/'REFINEMENT_RESULTS.csv').open()))
    local=list(csv.DictReader((ROOT/'LOCAL_WINDOW_RESULTS.csv').open()))
    cfg=json.loads((ROOT/'REFINED_CANDIDATES.json').read_text())
    nbase=len(coarse)+len(refine)+cfg['n_shifted']+len(local)
    nensembles=len(members)
    # 4 tighter/longer releases + 1 radial preparation for each of three high-fraction points.
    count=nbase+nensembles+15
    control=json.loads((ROOT/'CONTROL_RESULTS.json').read_text())
    best=cs[0];larger=high[0]['larger_seed_ensemble'];pb=pair_budgets[0]
    finding=dict(production_release_count=count,production_grid_release_count=nbase,
        phase_ensemble_member_count=nensembles,additional_high_point_control_releases=15,
        control_script_releases_not_included_in_production_count=8,
        refined_nominal_candidates=23,nominal_candidates_above_019=3,
        baseline_high_scale_mean=iis[0]['eps_net'],baseline_high_scale_SEM=iis[0]['eps_net_sem'],
        baseline_high_scale_local_rms=iis[0]['local_eps_rms'],
        best_supplied_phase_H_GeV=best['H_GeV'],best_cH=best['cH'],best_delta=best['delta'],
        best_sigma001_epsilon=best['eps_net'],best_sigma001_SEM=best['eps_net_sem'],
        best_sigma003_epsilon=larger['eps_net'],best_sigma003_SEM=larger['eps_net_sem'],
        best_Y_source=pb,
        no_window_claim_supported=False,
        a_complete_create_and_protect_cosmology_constructed=False,
        baseline_action_has_a_phase_selector=False,
        final_relic_population_computed=False,
        outcome='A reproducible conditional supplied-phase release clears the X-only regularity/yield/0.19 tests. It does not clear the complete leftover or preparation problem. Its Y source exceeds the unthermalized pair allowance; a real decay/thermalization calculation is still required.',
        max_grid_Gauss_error=max(float(r['max_Gauss_over_initial_Ncirc']) for r in coarse+refine+local),
        max_grid_work_error=max(float(r['max_Hamiltonian_work_error_over_initial_energy']) for r in coarse+refine+local),
        controls=control)
    write(ROOT/'FINDINGS.json',finding)
    print('FINDINGS',count,'production releases;',nensembles,'ensemble members; best phase epsilon',best['eps_net'],larger['eps_net'],
          'Y budget ratio',pb['ratio_to_one_percent_pair_allowance'])


if __name__=='__main__':main()
