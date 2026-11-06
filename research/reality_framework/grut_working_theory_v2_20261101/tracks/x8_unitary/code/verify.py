#!/usr/bin/env python3
"""Verify the scoped X8 physics, replay, and failed acceptance gates."""
import argparse
import hashlib
import json
import math
from pathlib import Path

TRACK=Path(__file__).resolve().parents[1]
RESULTS=TRACK/'results'
EVIDENCE=('controls','source','paired_source','preparation','transport','memory')


def read(name,folder=RESULTS):
    return json.loads((folder/f'{name}.json').read_text())


def close(a,b):
    if isinstance(a,dict):
        return all(k in b and close(v,b[k]) for k,v in a.items())
    if isinstance(a,list):
        return len(a)==len(b) and all(close(x,y) for x,y in zip(a,b,strict=True))
    if isinstance(a,float):
        return math.isclose(a,b,rel_tol=1e-8,abs_tol=2e-11)
    return a==b


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path)
    replay=p.add_mutually_exclusive_group()
    replay.add_argument('--replay',type=Path)
    replay.add_argument('--spot-replay',type=Path)
    args=p.parse_args()
    controls,first,paired,prepared,transport,memory=(read(name) for name in EVIDENCE)
    checks={}

    def check(name,value):
        checks[name]=bool(value)

    inputs=read('inputs')
    for path,digest in inputs['sha256'].items():
        check(f'input {path}',hashlib.sha256((TRACK.parent/path).read_bytes()).hexdigest()==digest)

    s=controls['static']
    check('unmatched contact term eliminated against exact spectrum',max(
        r['absolute_error'] for r in s['contact_checks'])<1e-13 and
          min(abs(r['old_contact_excess']) for r in s['contact_checks'])>1e-3)
    response=s['periodic_response']
    check('all-mode kinetic remainder has a first-order tail',all(
        abs(r['asymptotic_order']-1)<0.02 for r in response))
    paired_orders=[r['paired_asymptotic_order'] for r in response]
    check('periodic reference converges independently',max(
        r['reference_refinement_error'] for r in response)<1e-7)
    check('paired remainder improves the quadratic cutoff order',min(paired_orders)>1.5 and
          sum(paired_orders)/len(paired_orders)>1.85)
    check('paired continuum reference accuracy',max(
        r['ladder'][-1]['paired_relative_error'] for r in response)<2e-6)
    check('paired finite-amplitude static refinement',abs(
        s['finite_amplitude'][-2]['paired']/s['finite_amplitude'][-1]['paired']-1)<1e-3)
    check('homogeneous helicity solver recovered',max(controls['homogeneous'].values())<1e-11)
    spatial=controls['spatial_ode']
    check('independent matrix exponential ODE and second-order evolution',spatial['comparisons'][-1]['max_abs_error']<1e-6
          and spatial['comparisons'][0]['max_abs_error']>10*spatial['comparisons'][-1]['max_abs_error'])
    check('unitary matter split reverses',spatial['reversibility_max_abs_error']<1e-10)
    check('force from the same mass energy',max(controls[k]['force_abs_error'] for k in ('action','paired_action'))<1e-9)
    check('pressure from the matched action',max(controls[k]['pressure_relative_error']
          for k in ('action','paired_action'))<1e-7)
    check('mass magnitude preserved exactly',controls['action']['mass_squared_error']<1e-13)
    check('static and massless controls',max(abs(v) for c in controls['vacuum'].values() for v in c.values())<1e-11)

    check('original instantaneous cutoff miss preserved',not first['matched_comparisons']['wide']['gate_pass']
          and first['matched_comparisons']['wide']['instantaneous_peak_scaled_error']>0.02)
    check('original wide budget miss preserved',first['cases']['wide']['max_pressure_work_error']>1e-5)
    check('paired order-two instantaneous miss preserved',not paired['matched_comparisons']['paired_wide']['gate_pass']
          and paired['matched_comparisons']['paired_wide']['instantaneous_peak_scaled_error']>0.02)
    for name,data in (('paired',paired),('third-order',prepared)):
        cases=data['cases'].values()
        check(f'{name} pressure work',max(r['max_pressure_work_error'] for r in cases)<1e-5)
        check(f'{name} raw and matched fermion work',max(
            max(r['max_fermion_work_error'],r['max_fermion_work_error_matched']) for r in cases)<1e-5)
        check(f'{name} total momentum',max(r['max_comoving_momentum_error'] for r in cases)<1e-9)
        check(f'{name} occupied projectors',max(r['orthonormality_error'] for r in cases)<1e-10)
        check(f'{name} gradient and directional stress',all(
            r['max_scalar_gradient_fraction']>1e-3 and r['max_fermion_anisotropic_stress']>1e-5 for r in cases))
    pairs=prepared['matched_comparisons']
    check('third-order both cutoff ladders pass',all(pairs[k]['gate_pass'] for k in ('wide','cutoff_later')))
    check('third-order independent discretizations pass',all(
        pairs[k]['gate_pass'] for k in ('wide_half','longitudinal','transverse','radial','collocation','bandwidth')))
    check('preparation choice shifts work beyond tolerance',all(not pairs[k]['gate_pass']
          and pairs[k]['instantaneous_peak_scaled_error']>0.04 for k in ('later','wide_later')))

    stationary=transport['stationary_controls']
    check('rotating rates from independent two-body angular integral',stationary['width_max_abs_error']<1e-13)
    check('constant-rotation spin update positive and matched',stationary['bloch_max_abs_error']<1e-13
          and stationary['max_bloch_norm']<=1)
    check('zero-drive total time dilation',stationary['zero_drive_max_abs_error']<1e-13)
    check('drive reversal and branch closure',stationary['helicity_drive_reversal_error']<1e-13
          and stationary['anti_aligned_threshold_matches'])
    t=transport['cases']
    check('same exact H55 checkpoint and restart',transport['checkpoint_sha256']==memory['checkpoint_sha256']
          and max(r['timestep_restart_error'] for r in t.values())<1e-12)
    check('rotating collision substep energy and charge',all(r['collision_substep']['energy_error']<1e-8
          and r['collision_substep']['charge_error']<1e-11 for r in t.values()))
    check('live expansion pressure budget and Pauli bounds',all(r['budget_max']<1e-6 and
          r['bloch_norm_max']<1+1e-9 and r['pauli_min']>=-1e-10 and r['pauli_max']<=1+1e-10 for r in t.values()))
    check('daughter source independently integrated',all(r['radiation_source_integral_error']<1e-3 and
          r['final']['radiation_sm_fraction']>0 for r in t.values()))
    check('no qualifying H55 splice',all(not r['gate_pass'] and not r['canonical_scalar_free_projection']['gate_pass']
          and r['sampled_passing_points']==0 and r['best_window_max_gap']>0.8 for r in t.values()))
    check('timestep rejection robust',transport['timestep_max_abs_differences']['scalar_work']<0.01)
    check('local collision approximation misses threshold diagnostic',all(
        0.19<r['validity']['all']['integrated']<0.25 for r in t.values()))

    benchmark=memory['stationary']
    check('two-time stationary positive channel matches',benchmark['probability_relative_error']<1e-8
          and benchmark['energy_relative_error']<1e-3)
    check('two-time stationary negative branch leakage negligible',benchmark['negative_branch_probability']<1e-15)
    check('sampled memory calculation remains perturbative',memory['maximum_exact_probability']<0.01
          and memory['archive_positivity_roundoff_adjustment']<1e-10)
    check('two-time numerical controls',all(memory[k]<limit for k,limit in
          (('conformal_timestep_refinement_max_relative_error',1e-4),
           ('spectral_spacing_refinement_max_relative_error',1e-4),
           ('background_timestep_refinement_max_relative_error',2e-3),
           ('local_quadrature_refinement_max_relative_error',1e-3),
           ('uv_refinement_max_relative_error',1e-8))))
    rows=[row for window in memory['windows'] for row in window['modes']]
    check('finite-time discrepancy reproduced on frozen controls',len(rows)==32 and all(
          r['frozen_energy']/r['frozen_local_energy']>2. for r in rows)
          and max(abs(r['exact_energy']/r['frozen_energy']-1) for r in rows)<0.11)
    check('finite exposure dependence preserved',max(
          abs(r['exact_to_local_energy']-r['alternative_envelope_energy_ratio']) for r in rows)>1)
    check('no physical relic or empirical acceptance promoted',all(not r['scientific_acceptance']
          for r in (first,paired,prepared,transport,memory)))

    if args.replay:
        for name in EVIDENCE:
            check(f'fresh {name} replay',close(read(name),read(name,args.replay)))
    if args.spot_replay:
        for name in ('controls','transport','memory'):
            check(f'fresh {name} replay',close(read(name),read(name,args.spot_replay)))
        for name,data in (('paired_reference',paired),('uv_reference',prepared)):
            check(f'fresh {name} replay',close(data['cases'][name],read(name,args.spot_replay)))

    result=dict(passed=sum(checks.values()),total=len(checks),checks=checks,
                fixed_preparation_synthetic_source_gate_pass=
                checks['third-order both cutoff ladders pass'] and checks['third-order independent discretizations pass']
                and checks['third-order pressure work'] and checks['third-order total momentum'],
                preparation_independence_pass=all(pairs[k]['gate_pass'] for k in ('later','wide_later')),
                H55_spatial_source_certified=False,H55_transport_gate_pass=any(r['gate_pass'] for r in t.values()),
                scientific_acceptance=False,inherited_strict_accepted=0,inherited_strict_tested=135)
    if args.output:
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(f'{result["passed"]}/{result["total"]} evidence checks passed')
    failed=[name for name,value in checks.items() if not value]
    if failed:
        raise SystemExit('\n'.join(failed))


if __name__=='__main__':
    main()
