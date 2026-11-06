#!/usr/bin/env python3
"""Verify X7 evidence while retaining the distinct source and interface decisions."""
import argparse
import hashlib
import json
import math
from pathlib import Path

TRACK = Path(__file__).resolve().parents[1]
RESULTS = TRACK/'results'
EVIDENCE = ('static','controls','evolution','radial','overlap','pulse')
SPOT_CASES = ('gauss_refined','gauss_cutoff','gauss_heavier')


def read(name,folder=RESULTS):
    return json.loads((folder/f'{name}.json').read_text())


def close(a,b):
    if isinstance(a,dict):
        # The first archived uniform runs predate explicit quadrature metadata.
        return all(k in b and close(v,b[k]) for k,v in a.items())
    if isinstance(a,list):
        return len(a)==len(b) and all(close(x,y) for x,y in zip(a,b,strict=True))
    if isinstance(a,float):
        return math.isclose(a,b,rel_tol=1e-8,abs_tol=2e-11)
    return a==b


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    replay=parser.add_mutually_exclusive_group()
    replay.add_argument('--replay',type=Path)
    replay.add_argument('--spot-replay',type=Path)
    args=parser.parse_args()
    static,controls,evolution,radial,overlap,pulse=(read(name) for name in EVIDENCE)
    checks={}

    def check(name,criterion):
        checks[name]=bool(criterion)

    inputs=read('inputs')
    for path,digest in inputs['sha256'].items():
        check(f'input {path}',hashlib.sha256((TRACK.parent/path).read_bytes()).hexdigest()==digest)
    check('signed mass moments cancel',max(map(abs,static['cancellation_residuals']))<1e-13)
    check('signed edge falls at cubic rate',all(
        abs(row['signed_order']-3)<0.04 for rows in static['edge_orders_96_to_192'].values() for row in rows))
    check('physical-only edge grows linearly',all(
        abs(row['unsigned_order']+1)<0.03 for rows in static['edge_orders_96_to_192'].values() for row in rows))
    check('finite matrix static curvature',max(row['relative_error'] for row in static['curvature_checks'])<1e-5)
    check('independent finite matrix spectral curvature',max(
        row['absolute_error'] for row in static['spectral_curvature_checks'])<1e-13)
    refs=static['references']
    check('independent periodic reference stable',all(abs(r['periodic']/r['refined']-1)<1e-8 for r in refs))
    check('large periodic volume approaches bulk',all(abs(r['large_periodic']/r['bulk']-1)<1e-8 for r in refs))
    check('static harmonics converge against continuum',max(
        h['continuum_relative_error'] for r in static['cylinders']
        if r['longitudinal']==192 and r['spacing']==0.125 for h in r['harmonics'])<0.003)
    check('static quadrature limitation retained',max(
        h['continuum_relative_error'] for r in static['cylinders']
        if r['longitudinal']==192 and r['spacing']==0.5 for h in r['harmonics'])>0.03)
    check('static physical limit approached as auxiliaries get heavier',all(
        abs(r['highest_auxiliary_mass']-r['physical'])<abs(r['heavier']-r['physical'])
        <abs(r['refined']-r['physical']) for r in refs))
    check('massless response vanishes',not any(static['massless'].values()) and static['massless_continuum']==0)
    subthreshold=next(r for r in static['frequency_response'] if r['a']==1. and r['harmonic']==1
                      and r['frequency_squared']>0.)
    full=subthreshold['unbounded_auxiliary_difference']
    last=subthreshold['cylinders'][-1]
    check('finite auxiliary response persists without the cutoff',abs(full)>1e-3 and
          abs(last['total_difference']/full-1)<2e-4 and abs(last['edge'])<2e-7)
    check('homogeneous mode occupations converge at second order',
          controls['helicity'][-1]['max_occupation_error']<1e-6 and
          controls['helicity'][0]['max_occupation_error']>10*controls['helicity'][-1]['max_occupation_error'])
    check('homogeneous signed stress agrees',max(controls['helicity'][-1]['relative_errors'].values())<1e-5
          and controls['helicity'][-1]['field_spread']<1e-10)
    check('mass variation gives reciprocal force',controls['action']['force_error']<1e-9)
    check('volume variation gives pressure',controls['action']['pressure_relative_error']<1e-7)
    check('fixed-scale split reverses',controls['action']['reversibility_max_error']<1e-10)
    check('static and massless vacuum controls',max(abs(v) for d in controls['vacuum'].values()
          for v in d.values())<1e-11)
    all_cases=[*evolution['cases'].values(),*radial['cases'].values()]
    check('coupled energy pressure budget',max(r['max_pressure_work_error'] for r in all_cases)<1e-5)
    check('raw separate fermion work',max(r['max_fermion_work_error'] for r in all_cases)<1e-5)
    check('matched separate fermion work',max(r['max_fermion_work_error_matched'] for r in all_cases)<1e-5)
    check('translation momentum budget',max(r['max_comoving_momentum_error'] for r in all_cases)<1e-9)
    check('occupied orbitals remain orthonormal',max(max(r['orthonormality_error_by_sector'])
          for r in all_cases)<1e-10)
    check('stress independently cross-checked',max(abs(v) for d in evolution['stress_crosscheck'].values()
          for v in d.values())<1e-12)
    check('spatial wave samples directional stress',all(r['max_fermion_anisotropic_stress']>1e-5
          and r['max_scalar_gradient_fraction']>1e-3 for r in all_cases))
    base,refined=evolution['comparisons'],radial['comparisons']
    check('initial cutoff and independent spacing controls',all(row['gate_pass']
          for k,row in base.items() if k not in ('transverse_bins','heavier')))
    check('initial radial miss preserved',not base['transverse_bins']['gate_pass']
          and base['transverse_bins']['instantaneous_peak_scaled_error']>0.01)
    check('independent fine radial and cutoff controls',all(row['gate_pass']
          for k,row in refined.items() if k!='gauss_heavier'))
    check('auxiliary mass sensitivity persists',all(not r['gate_pass']
          and r['instantaneous_peak_scaled_error']>0.1 for r in (base['heavier'],refined['gauss_heavier'])))
    check('compact pulse leaves small auxiliary excitation',all(r['gate_pass']
          and r['weighted_end_energy'][0]>0 for r in pulse['cases']))
    check('compact pulse physical energy timestep agreement',abs(
        pulse['cases'][0]['weighted_end_energy'][0]/pulse['cases'][1]['weighted_end_energy'][0]-1)<1e-5)
    old=overlap['reference_gaps']
    cases=overlap['cases']
    check('restart matches archived background and stress',max(
        r['restart_max_relative_difference'] for r in cases.values())<1e-12)
    check('restart reproduces archived projection gaps',all(
        abs(r['initial'][k]-v)<1e-12 for r in cases.values() for k,v in old.items()))
    check('live restart budget and positivity',all(
        r['budget_max']<1e-6 and r['bloch_norm_max']<1+1e-9 and
        r['pauli_min']>=-1e-10 and r['pauli_max']<=1+1e-10 for r in cases.values()))
    check('live overlap samples the entire registered interval',all(
        r['sample_spacing']<2.5e-4 and abs(r['final']['t']-r['initial']['t']-0.12)<1e-12
        for r in cases.values()))
    check('H55 projection fails in both live models and both timestep runs',all(
        not r['gate_pass'] and r['sampled_passing_points']==0 and r['best_window_max_gap']>0.8
        for r in cases.values()))
    check('canonical scalar allocation also fails',all(
        not r['canonical_scalar_free_projection']['gate_pass'] and
        r['canonical_scalar_free_projection']['sampled_passing_points']==0 for r in cases.values()))
    check('vacuum damping transfers energy to visible radiation',all(
        (r['final']['radiation_sm_fraction']>0)==name.startswith('vacuum_decay')
        for name,r in cases.items()))
    check('timestep sensitivity of rejection is small compared with window gap',max(
        d['scalar_work'] for d in overlap['timestep_max_abs_differences'].values())<0.01)
    moments=overlap['archive_moment_check']
    check('exact coherence moment energy identity',max(moments['archive_reconstruction_errors'].values())<1e-12
          and moments['max_continuity_identity_error']<1e-12 and moments['max_bloch_ball_excess']<1e-9)
    check('phase reflection changes work at fixed gas moments',abs(moments['phase_reflected_work_over_rho']
          +old['scalar_work'])<1e-12 and moments['phase_reflection_derivative_separation_over_rho']>0.2)
    frozen=overlap['frozen_tail']
    check('conditional frozen pressure bound misses before half parent energy',
          frozen['energy_half_before_phase_independent_pressure_gate'] and
          frozen['pressure_envelope_at_half_parent_energy']>1e-3)
    check('earlier terminal refinement miss retained',not overlap['inherited_terminal_refinement_pass'])
    check('no scientific acceptance promoted',all(
        not r['scientific_acceptance'] for r in (evolution,radial,overlap)))
    if args.replay:
        for name in EVIDENCE:
            check(f'fresh {name} replay',close(read(name),read(name,args.replay)))
    if args.spot_replay:
        for name in ('static','controls','overlap'):
            check(f'fresh {name} replay',close(read(name),read(name,args.spot_replay)))
        for name in SPOT_CASES:
            check(f'fresh {name} replay',close(radial['cases'][name],read(name,args.spot_replay)))
    spatial_pass=all(row['gate_pass'] for k,row in base.items() if k!='transverse_bins') and all(
        row['gate_pass'] for row in refined.values())
    result=dict(passed=sum(checks.values()),total=len(checks),checks=checks,
                fixed_auxiliary_discretization_gate_pass=checks['initial cutoff and independent spacing controls']
                and checks['independent fine radial and cutoff controls'],
                auxiliary_mass_independence_pass=base['heavier']['gate_pass'] and refined['gauss_heavier']['gate_pass'],
                spatial_source_gate_pass=spatial_pass, H55_transport_gate_pass=any(r['gate_pass'] for r in cases.values()),
                scientific_acceptance=False, inherited_strict_accepted=0,inherited_strict_tested=135)
    if args.output:
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(f'{result["passed"]}/{result["total"]} evidence checks passed')
    failed=[name for name,passed in checks.items() if not passed]
    if failed:
        raise SystemExit('\n'.join(failed))


if __name__=='__main__':
    main()
