"""Non-independent controls; report discretization limits as well as ledgers."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import quad
from expansion import grid, deposit, free_stream_control, clock, OMEGA, ROOT, initial
from hybrid import profiled_operator, Hybrid
import transport


def file(tag):return json.loads((ROOT/(tag+'.json')).read_text())
def sample(d,a):return min(d['snapshots'],key=lambda r:abs(r['a_over_a0']-a))
def taurow(d,t):return min(d['snapshots'],key=lambda r:abs(r['tau']-t))


def exchange_roundtrips(n,omega,cycles):
    e=grid(0,8,'translated',n)
    v=np.zeros_like(e);v[np.argmin(abs(e-1.4))]=1
    start=float(e@v)
    for k in range(cycles):
        for shift in (2*omega,-2*omega):
            i=np.flatnonzero(v>0)
            if np.min(e[i]+shift)<e[0]-1e-11 or np.max(e[i]+shift)>e[-1]+1e-11:
                raise ValueError('two-state test lost support')
            v=deposit(e,e[i]+shift,v[i])
    return {'shift_bins':n,'Omega':omega,'cycles':cycles,'number_error':float(v.sum()-1),
            'mean_energy_error':float(e@v-start),
            'spurious_variance_mX2':float(v@((e-start)**2))}


def main():
    c=clock();eps=c['H0_mX2_over_ncore0']; checks=[]
    add=lambda name,passed,**kw:checks.append({'name':name,'pass':bool(passed),**kw})
    a=10
    integral=quad(lambda x:np.exp(-x)/eps,0,np.log(a),epsabs=1e-8)[0]
    exact=(1-1/a)/eps
    add('FRW_optical_time_integral',abs(integral/exact-1)<1e-12,numerical=integral,analytic=exact)
    fam=json.loads((ROOT/'FAMILY.json').read_text())
    fe=max(abs(x['energy_interpolation_error_mX']) for x in fam['holdouts'])
    fo=max(abs(x['omega_interpolation_error']) for x in fam['holdouts'])
    add('held_out_stationary_family',fe<.001 and fo<1e-6,max_energy_error_mX=fe,max_omega_error=fo)
    primary=file('expanding96');parent=json.loads((ROOT.parent/'grut_m1_kinetics_20261028'/'translated_full_96.json').read_text())
    u=taurow(primary,1000);v=parent['snapshots'][-1]
    diff={k:abs(u[k]/v[k]-1) for k in ('Nplus','Nminus','gas_kinetic_mX')}
    add('early_expanding_versus_static',max(diff.values())<.03,relative_differences=diff,
        note='Differences include real redshift and 8 versus 12 angular nodes; close agreement is not exact equivalence.')
    # An independently called supplied-kernel path must reduce exactly to the
    # inherited core operator at the reference profile and frequency.
    e=grid(0,8,'translated',32)
    z=np.load(ROOT.parent/'grut_m1_kinetics_20261028'/'CORE_KERNELS.npz')
    alt,alt_omit,alt_error=profiled_operator(e,OMEGA,z['k'],z['W'])
    old,old_omit,_,_=transport.core_operator(e)
    diff_operator=float(np.max(abs(alt-old))); diff_omit=float(np.max(abs(alt_omit-old_omit)))
    add('supplied_kernel_matches_inherited_operator',diff_operator<1e-10 and diff_omit<1e-10,
        max_matrix_difference=diff_operator,max_omission_difference=diff_omit)
    # Directional derivative tests the Q-dependent collision Jacobian and
    # redshift-work row against a finite difference of the independently
    # evaluated nonlinear RHS. Passive omitted/defect rows are excluded.
    model=Hybrid(e,'full',8,eps,.00025,True);source,_=initial(e);state=np.r_[source,0.,0.,0.]
    rng=np.random.default_rng(7261);direction=rng.normal(size=len(state));direction[-3:]=0
    step=.01;numeric=(model.rhs(0,state+step*direction)-model.rhs(0,state-step*direction))/(2*step)
    analytic=model.jac(0,state)@direction;dim=2*len(e)
    jerr=float(np.linalg.norm(numeric[:dim]-analytic[:dim])/np.linalg.norm(numeric[:dim]))
    werr=float(abs(numeric[-2]-analytic[-2]))
    add('profiled_jacobian_directional_difference',jerr<1e-7 and werr<1e-6,population_relative_error=jerr,work_row_absolute_error=werr)
    qnodes=list(range(186000,204001,3000));mat_errors=[]
    for q in qnodes:
        a=np.load(ROOT/f'core_Q{q}.npz'); mat_errors.append(np.max(a['checks'][:,:2]))
    add('profiled_kernel_unitarity_and_reciprocity',max(mat_errors)<1e-9,max_error=float(max(mat_errors)),charge_nodes=qnodes)
    holdouts=[]
    for tau in (0,1,10,100,1000):
        row=file('KERNEL_HOLDOUT_'+str(tau)); D=np.array(row['rate_sums']['direct']); P=np.array(row['rate_sums']['interpolated'])
        holdouts.append({'tau':row['tau'],'Q':row['q'],'coverage_ranked_by_interpolated_rate':row['coverage_by_interpolated_event_rate'],
                         'relative_total_rate_error':float(P.sum()/D.sum()-1),'sampled_bins':row['count_selected']})
    add('selected_population_weighted_kernel_holdouts',min(r['coverage_ranked_by_interpolated_rate'] for r in holdouts)>.98 and
        max(abs(r['relative_total_rate_error']) for r in holdouts)<.01,rows=holdouts,
        interpretation='Direct stationary scattering at sampled populated nodes, not a bound on unvisited momenta or on rate cancellation.')
    runs=['expanding64_replay','expanding96','ringdown005early_expanding64',
          'ringdown01_expanding64','ringdown005_expanding64','ringdown005_expanding96']
    for tag in runs:
        d=file(tag);ch=d['checks']
        add('fixed_ledger_'+tag,ch['energy_plus_work_max_absolute_drift_mX']<1e-5 and ch['minimum_saved_population']>-1e-7,
            drift_mX=ch['energy_plus_work_max_absolute_drift_mX'],minimum_population=ch['minimum_saved_population'])
    for coarse,fine in [('expanding64_replay','expanding96'),('ringdown005_expanding64','ringdown005_expanding96')]:
        a=sample(file(coarse),10);b=sample(file(fine),10)
        dd={k:abs(a[k]/b[k]-1) for k in ('Nplus','Qcore','gas_kinetic_mX')}
        add('fixed_refinement_'+fine,max(dd.values())<.01,relative_difference=dd,
            interpretation='selected integrated observables only, not a bound on the distribution or physical model')
    free={}
    for tag in ('expanding32','expanding64_replay','expanding96','ringdown005_expanding96'):
        d=file(tag); z=np.load(ROOT/(tag+'.npz')); i=int(np.argmin(abs(z['x']-np.log(2))))
        free[tag]=free_stream_control(z['E'],z['Y'][:-2,i],[5.])[0]
    add('free_streaming_energy_refinement',abs(free['expanding96']['relative_kinetic_error'])<.015,
        two_to_ten_control=free,threshold_is_engineering_only=True)
    for tag in ('hybrid32','hybrid64','hybrid96','hybrid_h005_64','hybrid_h005_96','hybrid64_fineomega',
                'profiled32','profiled64','profiled96','profiled_h005_64','profiled_h005_96',
                'profiled64_to10','profiled_h005_64_to10','profiled96_to10','profiled_h005_96_to10',
                'profiled64_fineomega','profiled64_tight','profiled64_angles12',
                'profiled_h0025_32','profiled_h0025_64','profiled_h0025_96',
                'profiled_h0025early_64','profiled_h0025_64_to10','profiled_h0025_96_to10'):
        path=ROOT/(tag+'.json')
        if not path.exists():continue
        d=file(tag);last=d['snapshots'][-1];ch=d['checks']
        add('hybrid_ledger_'+tag,ch['max_energy_minus_generator_defect_drift_mX']<.001 and
            abs(last['interpolated_generator_energy_defect_mX'])<.1 and ch['minimum_saved_population']>-1e-7,
            defect_mX=last['interpolated_generator_energy_defect_mX'],
            drift_after_defect_mX=ch['max_energy_minus_generator_defect_drift_mX'],minimum_population=ch['minimum_saved_population'],
            interpretation='energy defect is explicitly measured; it is not a bound on projection diffusion or profile error')
    for coarse,fine in [('profiled64','profiled96'),('profiled_h005_64','profiled_h005_96'),
                        ('profiled_h0025_64','profiled_h0025_96')]:
        if not (ROOT/(coarse+'.json')).exists() or not (ROOT/(fine+'.json')).exists():continue
        a=file(coarse)['snapshots'][-1];b=file(fine)['snapshots'][-1]
        delta={k:abs(a[k]-b[k]) for k in ('Nplus','Qcore','gas_kinetic_mX','stationary_cold_excess_fraction')}
        add('profiled_resolution_'+fine,delta['stationary_cold_excess_fraction']<.0001,absolute_differences=delta,
            interpretation='Selected scalar outputs only; moving-grid exchange remains numerically diffusive.')
    for control in ('profiled64_fineomega','profiled64_tight','profiled64_angles12'):
        if not (ROOT/(control+'.json')).exists():continue
        a=file('profiled64')['snapshots'][-1];b=file(control)['snapshots'][-1]
        delta={k:abs(a[k]-b[k]) for k in ('Nplus','gas_kinetic_mX','stationary_cold_excess_fraction')}
        add('profiled_control_'+control,delta['stationary_cold_excess_fraction']<.0001,absolute_differences=delta)
    cycles=[exchange_roundtrips(n,om,10) for n in (32,64,96) for om in (OMEGA,.88435)]
    baseline=max(x['spurious_variance_mX2'] for x in cycles if x['Omega']==OMEGA)
    add('aligned_vs_off_grid_exchange_control',baseline<1e-13 and all(x['spurious_variance_mX2']>1e-7 for x in cycles if x['Omega']!=OMEGA),
        rows=cycles,interpretation='Off-grid round trips demonstrably diffuse despite conserving number and mean energy; finer hybrid results are still approximate.')
    # Source is a compacted signed finite-window spectrum from another local
    # worktree. Check exact copied bytes and compare their energy moments.
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'sources').glob('*'))}
    d=file('ringdown005_expanding96');s=d['snapshots'][0]
    metadata=json.loads((ROOT/'sources'/'ringdown_h005.json').read_text())['signed_spectrum_diagnostics']
    expected=sum(metadata[w][species]['frequency_times_number_energy_mX'] for w in ('0_1000','1000_10000') for species in ('plus','minus','neutral'))
    err=s['gas_energy_mX']-expected
    add('radial_source_energy_moments',abs(err)<.1,loaded_gas_energy_mX=s['gas_energy_mX'],
        signed_energy_mX=expected,difference_mX=err,source_is_not_positivity_certified=True,sha256=hashes)
    if (ROOT/'profiled_h0025_96.json').exists():
        d=file('profiled_h0025_96');s=d['snapshots'][0]
        meta=json.loads((ROOT/'sources'/'ringdown_h0025.json').read_text())['signed_spectrum_diagnostics']
        expected=sum(meta[w][species]['frequency_times_number_energy_mX'] for w in ('0_1000','1000_10000') for species in ('plus','minus','neutral'))
        err=s['gas_energy_mX']-expected
        add('finest_radial_source_energy_moments',abs(err)<.1,loaded_gas_energy_mX=s['gas_energy_mX'],
            signed_energy_mX=expected,difference_mX=err,source_is_not_positivity_certified=True)
    h=[]
    for label in ('h01','h005','h0025'):
        row=json.loads((ROOT/'sources'/f'ringdown_{label}.json').read_text())
        h.append({'spacing':row['configuration']['h'],'late_probe_energy_mX':row['late_1000_to_10000']['probe_energy_mX'],
                  'late_core_excess_mX':row['checkpoints'][-1]['core20_excess_over_same_grid_stationary_mX']})
    ratio=(h[2]['late_probe_energy_mX']-h[1]['late_probe_energy_mX'])/h[2]['late_probe_energy_mX']
    add('radial_resolution_trend',0<ratio<.05,rows=h,last_relative_change=ratio,
        interpretation='Improving resolution trend, not proof of continuum convergence or of a quantum particle source.')
    out={'same_author':True,'independent_replication':False,'checks':checks,
         'all_engineering_checks_pass':all(c['pass'] for c in checks),
         'limitations':'Tolerances validate selected controls, not missing physics, local inhomogeneity, formation, late excitation or the final 1% mass criterion.'}
    (ROOT/'VERIFICATION.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2))
    if not out['all_engineering_checks_pass']:raise SystemExit(1)

if __name__=='__main__':main()
