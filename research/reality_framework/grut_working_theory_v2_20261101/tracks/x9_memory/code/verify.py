#!/usr/bin/env python3
"""Verify X9 matched preparations, resolved responses, and scope."""
import argparse
import hashlib
import json
import math
from pathlib import Path

TRACK=Path(__file__).resolve().parents[1]
RESULTS=TRACK/'results'
ARCHIVES=('production','envelope','bath','bath_refined','bath_dense','live')


def read(name,folder=RESULTS):
    return json.loads((folder/f'{name}.json').read_text())


def close(a,b):
    if isinstance(a,dict):
        # The two first bath runs recorded a file digest while the unused
        # later part of the history was still being written.
        return all(k in b and close(v,b[k]) for k,v in a.items() if k!='histories_sha256')
    if isinstance(a,list):
        return len(a)==len(b) and all(close(x,y) for x,y in zip(a,b,strict=True))
    if isinstance(a,float):
        return math.isclose(a,b,rel_tol=2e-8,abs_tol=2e-10)
    return a==b


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path)
    replay=p.add_mutually_exclusive_group()
    replay.add_argument('--replay',type=Path)
    replay.add_argument('--spot-replay',type=Path)
    args=p.parse_args()
    production,envelope,first,refined,bath,live=(read(n) for n in ARCHIVES)
    checks={}

    def check(name,criterion):
        checks[name]=bool(criterion)

    for path,digest in read('inputs')['sha256'].items():
        check(f'input {path}',hashlib.sha256((TRACK.parent/path).read_bytes()).hexdigest()==digest)
    for name,row in production['cases'].items():
        match=row['change_at_copy']
        keys=('t','a','H','rho','pressure','sigma','chi','s','dsigma','dchi','ds')
        matched=(*keys,*(f'{kind}_{j}' for j in range(1,5) for kind in ('n','rho','p')))
        check(f'{name} matches initial macroscopic data',max(abs(match[k]) for k in matched)<1e-12)
        check(f'{name} energy and Pauli',row['max_budget_error']<1e-8 and row['max_bloch_norm']<1+1e-9
              and row['pauli_min']>=-1e-10 and row['pauli_max']<=1+1e-10)
        if name.endswith('_reflect'):
            check(f'{name} changes initial work sign',row['normalized_work_at_copy_before']>0
                  and row['normalized_work_at_copy_after']<0)
    pvalues=list(production['comparisons'].values())
    largest_counts=max(abs(row['final_window_relative_difference'][k])
                       for row in pvalues for k in ('n_1','n_2','n_3','n_4','a3_n_1','a3_n_2','a3_n_3','a3_n_4'))
    comoving=max(abs(row['final_window_relative_difference']['comoving_fermion_number']) for row in pvalues)
    check('largest reflected species count difference stays subpercent',largest_counts<0.0012
          and all(row['final_window_within_one_percent'] for row in pvalues))
    check('total comoving count stable in reflected tests',comoving<0.00012)
    check('early work response survives radial and timestep variation',all(
        2.28<row['early_max_work_difference_over_rho']<2.30 and
        0.00037<row['early_max_fermion_fraction_difference']<0.0004 for row in pvalues))
    check('late baseline resolution distinguished from paired effect',max(
        abs(row['n_4']) for row in production['base_final_window_resolution_differences'].values())>largest_counts)
    e=[row['total_relative_bounds'] for row in envelope['cases'].values()]
    check('admissible fixed-history window envelope on all grids',max(row['mean_number'] for row in e)<0.0004
          and max(row['mean_energy'] for row in e)<0.00006)
    check('admissible fixed-history endpoint envelope on all grids',max(row['endpoint_number'] for row in e)<0.0007)
    check('orthogonal response propagation',max(
        row['tangent_error'] for case in envelope['cases'].values()
        for row in case['per_species_absolute_bounds'])<1e-10)

    check('same early scalar trace pinned at all bath stages',
          first['source_trace_sha256']==refined['source_trace_sha256']==bath['source_trace_sha256'])
    pole=bath['stationary_poles']
    check('dispersive real mass reference',max(row['independent_principal_error'] for row in pole)<2e-9)
    check('stationary pole width matches Born scale',max(row['born_width_relative_error'] for row in pole)<0.004)
    c=[row for runs in bath['cases'].values() for row in runs]
    check('positive spectral packet preparations',max(
        max(row[k] for k in ('positive_projection_error','orthogonality','bath_packet_sterile_leakage',
                             'sterile_covariance_difference')) for row in c)<1e-10)
    check('identical initial sterile and bath marginals on measured quantities',max(
        max(abs(v) for v in row['initial_differences'].values())/max(1.,row['mean_initial_energy'])
        for row in c)<1e-11)
    check('driven star avoids taper at instantaneous branches',max(
        row['maximum_drive_branch_cutoff_ratio'] for row in c)<0.5)
    check('driven bath interval below discrete recurrence time',all(
        row['conformal_duration']<row['recurrence_time'] for row in c))
    check('Padé energy plus work identity and unitarity',max(row['max_energy_work_error'] for row in c)<1e-9
          and max(row['max_norm_error'] for row in c)<1e-9)
    check('independent continuous torque converges',max(
        row['continuous_torque_difference'] for row in bath['cases']['timestep_fine'])<2e-5
          and max(row['continuous_torque_difference'] for row in bath['cases']['reference'])>8e-5)
    check('dense output refinement',max(
        row['decimated_window_absolute_differences']['bath_energy']/row['mean_initial_energy'] for row in c)<1e-7)
    b=bath['final_window_refinement_differences']
    check('dense frequency and timestep comparisons',all(
        max(b[k].values())<0.005 for k in ('frequency','timestep','timestep_fine')))
    check('early uniform-grid and sparse-output misses retained',
          first['final_window_refinement_differences']['frequency']['bath_energy']>0.1
          and refined['final_window_refinement_differences']['timestep']['bath_energy']>0.1)
    check('remaining bath cutoff sensitivity retained',b['cutoff']['bath_energy']>0.01)
    for i in (1,2):
        ref=bath['cases']['reference'][i]['final_window_differences']['bath_energy']
        check(f'robust nonzero driven daughter energy contrast mode {i}',abs(ref)>10 and all(
            abs(r[i]['final_window_differences']['bath_energy']/ref-1)<0.005 for r in bath['cases'].values()))

    micro=live['microscopic']
    check('split propagator against dense exponential',micro['split'][-1]['error']<1e-6
          and micro['split'][0]['error']>10*micro['split'][-1]['error'])
    check('Padé propagator against dense exponential',micro['pade']['error']<1e-8
          and micro['pade']['max_energy_work_error']<1e-10)
    check('leading bath amplitude matches direct Born integral',micro['leading_born_relative_error']<1e-5)
    st=live['stationary']
    check('conjugate stationary states have same spectral weights',st['exact_stationary_spectral_weight_difference']<1e-13)
    check('stationary bulk count loses memory before recurrence',all(
        row['recurrence_time']>120 and row['max_early_sterile_count_difference']>0.5
        and row['max_late_sterile_count_difference']<2e-6
        and row['max_late_total_positive_bath_count_difference']<2e-6 for row in st['runs']))
    outgoing=[row['mean_late_chiral_bath_count_difference'] for row in st['runs']]
    check('stationary daughter channel retains phase memory',min(outgoing)>0.10
          and max(outgoing)-min(outgoing)<1e-4)
    feedback=live['live']['runs']
    check('live reciprocal scalar and bath energy and charge',all(
        row['max_energy_error']<1e-11 and row['max_charge_error']<1e-10
        and row['max_pauli_projector_error']<1e-10 for r in feedback.values() for row in r['cases']))
    check('live matched initial inputs',max(abs(v) for r in feedback.values()
          for v in r['differences']['initial'].values())<1e-11)
    check('live scalar and bath retain different energies',all(
        abs(r['differences']['terminal']['scalar_energy'])>0.09
        and abs(r['differences']['terminal']['bath_energy'])>0.45 for r in feedback.values()))
    check('live feedback phase differences stable under timestep',all(
        abs(feedback['base']['differences']['terminal'][k]/feedback['half']['differences']['terminal'][k]-1)<0.01
        for k in ('scalar_energy','bath_energy','sterile_count')))
    check('no relic abundance acceptance promoted',all(
        not read(name)['scientific_acceptance'] for name in ARCHIVES))

    if args.replay:
        for name in ARCHIVES:
            check(f'fresh {name} replay',close(read(name),read(name,args.replay)))
    if args.spot_replay:
        for name in ('bath_dense','live'):
            check(f'fresh {name} replay',close(read(name),read(name,args.spot_replay)))
        spot=read('production',args.spot_replay)
        check('fresh coarse full-history pair',close(production['comparisons']['coarse'],spot['comparisons']['coarse'])
              and all(close(production['cases'][k],v) for k,v in spot['cases'].items()))

    report=dict(passed=sum(checks.values()),total=len(checks),checks=checks,
                H55_reflection_bulk_window_pass=all(row['final_window_within_one_percent'] for row in pvalues),
                H55_full_bath_cutoff_independence_pass=b['cutoff']['bath_energy']<0.01,
                largest_reflected_species_count_difference=largest_counts,
                largest_reflected_comoving_total_count_difference=comoving,
                largest_fixed_history_window_count_envelope=max(row['mean_number'] for row in e),
                physical_spatial_H55_source_certified=False,H55_particle_handoff_certified=False,
                inherited_strict_accepted=0,inherited_strict_tested=135,scientific_acceptance=False)
    if args.output:
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(f'{report["passed"]}/{report["total"]} evidence checks passed')
    fail=[k for k,v in checks.items() if not v]
    if fail:
        raise SystemExit('\n'.join(fail))


if __name__=='__main__':
    main()
