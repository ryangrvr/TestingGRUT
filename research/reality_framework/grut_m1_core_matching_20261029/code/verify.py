"""Same-author checks of the stored matching, restart and transport results."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicHermiteSpline

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT.parent/'grut_m1_local_cloud_20261029'
PARENT=ROOT.parent/'grut_m1_expansion_clock_20261028'
def read(p):return json.loads(p.read_text())

def main():
    tests=[]
    def add(name,ok,**detail):tests.append({'name':name,'pass':bool(ok),**detail})
    for tag,tol in [('h005',.04),('h0025',.01)]:
        d=read(ROOT/f'late_{tag}_analysis.json');c=d['checks'];start=max(abs(x) for x in d['restart_differences_in_unscaled_engine_units'].values())
        add('restart_state_'+tag,start<1e-12,max_abs_difference=start)
        add('late_energy_ledger_'+tag,c['segment_energy_plus_absorber_max_abs_drift_mX']<tol,
            max_drift_mX=c['segment_energy_plus_absorber_max_abs_drift_mX'])
        add('late_charge_ledger_'+tag,c['segment_charge_plus_absorber_max_abs_drift']<1e-6,
            max_drift=c['segment_charge_plus_absorber_max_abs_drift'])
        anti=[x['minus']['positive_projection_number'] for x in d['late_window_sources']]
        neutral=[x['neutral']['positive_projection_number'] for x in d['late_window_sources']]
        add('late_antiparticle_neutral_window_stability_'+tag,max(anti)/min(anti)<1.01 and max(neutral)/min(neutral)<1.01,
            direct_antiparticles=anti,neutral_parents=neutral,
            limitation='Window stability is not a quantum particle identification or a bound on future emission.')
    a=read(ROOT/'MATCHING_SURFACES.json')['field_shell_40_to_80']
    add('shell_energy_flux_consistency',abs(a['energy_residual_boundary_matching_mX'])<.005,**{k:a[k] for k in ('delta_energy_mX','in_minus_out_integrated_probe_energy_mX','energy_residual_boundary_matching_mX')})
    add('shell_charge_flux_consistency',abs(a['charge_residual_boundary_matching'])<.003,**{k:a[k] for k in ('delta_charge','in_minus_out_integrated_probe_charge','charge_residual_boundary_matching')})
    family=read(PARENT/'FAMILY.json')['nodes'];mass=CubicHermiteSpline([r['Q'] for r in family],[r['energy_mX'] for r in family],[r['Omega'] for r in family],extrapolate=False)
    for tag in ('h005_64','h005_96','h0025_64','h0025_96'):
        d=read(ROOT/f'late_kin_{tag}.json');c=d['checks']
        err=d['initial_pulse_energy_mX']-float(mass(d['initial_Q_core']))-d['initial_gas_energy_mX']-d['held_unassigned_recomputed_mX']
        add('source_remainder_identity_'+tag,abs(err)<1e-8,error_mX=err)
        add('transport_ledger_'+tag,c['max_energy_plus_work_minus_generator_defect_drift_mX']<2e-5 and c['minimum_saved_population']>-1e-8,**c)
    s=read(ROOT/'SCENARIO_AUDIT.json');parts=s['energy_decompositions_at_scale_factor_10']
    for k,v in parts.items():
        error=v['cold_excess_percent']+v['gas_kinetic_percent']+v['held_remainder_percent']-v['full_percent']
        add('mass_decomposition_'+k,abs(error)<1e-9,error_percent=error)
    a=read(ROOT/'late_kin_h0025_64.json')['full_modeled_final_excess_percent_with_remainder']
    b=read(ROOT/'late_kin_h0025_96.json')['full_modeled_final_excess_percent_with_remainder']
    add('late_transport_energy_node_comparison',abs(a-b)<.002,percent_64=a,percent_96=b)
    a=read(ROOT/'ANGULAR_AUDIT.json')
    archived=[r for r in a['3D_runs'] if 'matches_archived_initial_sha256' in r]
    matched=all(r['matches_archived_initial_sha256'] for r in archived)
    diff=max(abs(k) for r in archived for row in r['common_time_comparison_to_archived_trace'] for k in (row['relative_energy_difference'],row['relative_charge_difference']))
    add('3d_initial_and_archived_trace',matched and diff<1e-12,max_common_time_relative_difference=diff,
        limitation='Same numerical code and same initial data; not independent replication.')
    er=max(r['max_relative_energy_drift'] for r in a['3D_runs']);qr=max(r['max_relative_charge_drift'] for r in a['3D_runs'])
    add('3d_conservation',er<1e-6 and qr<2e-12,max_relative_energy_drift=er,max_relative_charge_drift=qr)
    ortho=max(x['orthonormality_max_error'] for x in a['multipole_controls'].values())
    add('angular_quadrature_orthonormality',ortho<1e-12,max_error=ortho)
    one=a['multipole_controls']['multipoles96']['r20_difference'];two=a['multipole_controls']['multipoles96_double']['r20_difference']
    err=max(abs(x['field_power_fraction_by_l'][2]-y['field_power_fraction_by_l'][2]) for x,y in zip(one,two))
    add('angular_snapshot_precision_at_r20',err<2e-5,max_l2_fraction_difference=err)
    low=read(ROOT/'multipoles96_t20.json')['result']['rows'];high=read(ROOT/'multipoles96_double_t20.json')['result']['rows']
    select=lambda rows:next(r for r in rows if r['radius_mX_inverse']==24 and r['case']=='difference_to_control' and r['channel']=='neutral_real')
    f32=sum(select(low)['field_power_by_l']);f64=sum(select(high)['field_power_by_l'])
    add('weak_neutral_difference_requires_double_snapshot',f32==0 and f64>0,power_float32=f32,power_float64=f64,
        limitation='This precision contrast does not establish physical signal accuracy at such a small amplitude.')
    fine=read(ROOT/'mixing_cont_nv8001.json');coarse=read(ROOT/'mixing_cont_nv4001.json')
    nerr=max(abs(x) for r in fine['rows'] if r['mX_t']<=1e7 for x in r['relative_number_error'])
    verr=max(abs(x/y-1) for a,b in zip(fine['rows'][:5],coarse['rows']) for x,y in zip(a['density_cv'],b['density_cv']))
    add('continuous_mixing_number_and_radius_resolution',nerr<2e-5 and verr<.003,max_number_error=nerr,max_cv_relative_difference=verr)
    asym=fine['asymptotic_check']['density_cv'];num=fine['rows'][-1]['density_cv'];rel=max(abs(x/y-1) for x,y in zip(num,asym))
    add('continuous_mixing_asymptotic_control',rel<.001,relative_difference=rel,finite=num,asymptotic=asym)
    discrete=read(ROOT/'mixing_discrete_sweep.json');previous=read(OLD/'shot_noise_dr2p5.json')
    oldcv=previous['rows'][2];newcv=discrete['rows'][0]
    err=max(abs(newcv['plus']['density_cv']/oldcv['poisson_cv_plus']-1),abs(newcv['minus']['density_cv']/oldcv['poisson_cv_minus']-1))
    add('discrete_shell_sweep_vs_parent_spatial_grid',err<.003,max_relative_difference=err)
    mono=fine['mono_shell_control'];a=mono['radius_inner'];b=mono['radius_outer'];width=b-a;nc=s['ncore_over_mX3']
    integral=quad(lambda r:1/r**2,a,b,epsabs=1e-24,epsrel=1e-10)[0]
    independent=np.sqrt(integral/(4*np.pi*nc*width**2));err=abs(independent/mono['exact_poisson_cv']-1)
    add('mono_speed_shell_vs_direct_radial_quad',err<1e-9,relative_difference=err)
    add('clock_conversion',abs(s['unit_optical_time_mX_t']*nc-1)<1e-12 and abs(s['tau_at_mX_t_1e7']-nc*1e7)<1e-18,
        optical_time_one_mX_t=s['unit_optical_time_mX_t'])
    result={'same_author':True,'independent_replication':False,'checks':tests,'all_pass':all(t['pass'] for t in tests),
            'scope':'Checks code consistency and selected numerical controls, not the physical angular distribution, quantum source positivity, mixing, future emission, formation, or a cosmological 1% verdict.'}
    (ROOT/'VERIFICATION.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2))
    if not result['all_pass']:raise SystemExit(1)

if __name__=='__main__':main()
