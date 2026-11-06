#!/usr/bin/env python3
"""Verify both X6 gates and preserve the blocked H55 scientific verdict."""
import argparse
import hashlib
import json
import math
from pathlib import Path

TRACK = Path(__file__).resolve().parents[1]
TRACKS = TRACK.parent
RES = TRACK/'results'
REPO = TRACK.parents[4]
PARENT = '7a4c5aa9caf534d806564e89dced4cd66dda28c9'


def read(name,directory=RES):
    return json.loads((directory/name).read_text())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--replay',type=Path)
    args=parser.parse_args()
    checks={}

    def check(name,value):
        checks[name]=bool(value)

    ctl=read('controls.json')
    old=read('source_initial.json')
    source=read('source.json')
    follow=read('source_followup.json')
    edge=read('edge.json')
    transport=read('transport.json')
    check('initial source experiment retained',len(old['cases'])==6)
    check('separate source refinements',len(source['cases'])==9 and len(follow['cases'])==4)
    histories=[*source['cases'].values(),*follow['cases'].values()]
    check('pressure-work total budget',max(s['max_pressure_work_error'] for s in histories)<1e-5)
    check('bare fermion pressure-work identity',max(s['max_separate_fermion_work_error'] for s in histories)<1e-5)
    check('matched fermion pressure-work identity',max(s['max_matched_fermion_work_error'] for s in histories)<1e-5)
    check('longitudinal momentum conservation',max(s['max_comoving_momentum_error'] for s in histories)<1e-10)
    check('all occupied orbitals orthonormal',max(s['orthonormality_error'] for s in histories)<1e-10)
    check('directional pressures reproduce isotropic variation',all(
        abs((s['final']['P_x']+2*s['final']['P_perp'])/3-s['final']['P'])<1e-10 for s in histories))
    check('finite source includes expansion and gradients',all(
        s['final']['a']>1 and s['maximum_gradient_fraction']>0 for s in histories))
    check('homogeneous helicity occupations',ctl['helicity'][-1]['mode_occupation_max_absolute']<1e-5)
    check('homogeneous helicity refinement',ctl['helicity'][-1]['mode_occupation_max_absolute']
          <ctl['helicity'][0]['mode_occupation_max_absolute']/10)
    check('independent spatial matrix ODE',ctl['spatial']['ode'][-1]['max_absolute_error']<1e-5)
    check('fixed-scale reversibility',ctl['spatial']['inversion_max_absolute']<1e-10)
    check('static uniform vacuum',abs(ctl['vacuum']['static_vacuum']['rho_F'])<1e-11)
    check('massless decoupling',abs(ctl['vacuum']['massless_wave']['rho_F'])<1e-11
          and ctl['vacuum']['massless_wave']['Jv']==0)
    check('spatial counterterm finite difference',ctl['matching']['static_curvature_relative_error']<1e-5)
    check('mass vertex yields scalar force',ctl['matching']['field_variation_absolute_error']<1e-8)
    check('isotropic action variation yields pressure',ctl['matching']['pressure_variation_relative_error']<1e-8)
    bare=source['refinement_comparisons']
    matched=source['matched_refinement_comparisons']
    check('original raw cutoff miss retained',not old['refinement_comparisons']['wide_cutoff']['gate_pass'])
    check('spatial quadrature and bandwidth resolved',all(matched[k]['gate_pass'] for k in (
        'spatial_quadrature','scalar_bandwidth')))
    check('independent timestep controls at both cutoffs',all(matched[k]['gate_pass'] for k in (
        'half_step','wide_half_step')))
    check('fine transverse quadrature resolves both work metrics',
          follow['comparisons']['transverse_bins__transverse_64_first']['matched']['gate_pass'])
    check('matched accumulated work stable against wide cutoff',matched['wide_cutoff']['integrated_transfer_error']<0.01)
    check('matched instantaneous work fails wide cutoff',matched['wide_cutoff']['instantaneous_peak_scaled_error']>=0.01)
    check('larger fine cutoff still misses',not follow['comparisons']['wide_half_step__wider_first']['matched']['gate_pass'])
    check('second-order rotating-vacuum improvement',all(
        row['second_error']<row['first_error'] for row in follow['rotating_vacuum']))
    check('second-order preparation still misses cutoff criterion',
          not follow['comparisons']['reference_second__wide_second']['matched']['gate_pass'])
    check('edge dominates reference static matching',edge['cylinders'][0]['harmonics'][0]['edge_fraction']>0.5)
    check('edge decomposition checked by all six static probes',all(
        h['static_curvature_relative_error']<2e-5 for c in edge['cylinders'][:2] for h in c['harmonics']))
    ratios=edge['largest_cylinder_edge_ratios']
    check('sharp-edge nonlocal harmonic scaling',abs(ratios[1]-0.5)<0.02 and abs(ratios[2]-1/3)<0.02)

    pulse=transport['pulse']
    check('compact pulse energy and momentum source ledger',pulse[-1]['work_relative_error']<2e-5
          and pulse[-1]['momentum_relative_error']<2e-5 and abs(pulse[-1]['momentum'])>0.01*pulse[-1]['energy'])
    check('pulse refinement',pulse[-1]['work_relative_error']<pulse[0]['work_relative_error']/3.5
          and pulse[-1]['momentum_relative_error']<pulse[0]['momentum_relative_error']/3.5)
    d=transport['controlled_elastic']
    plain=transport['controlled_no_elastic']
    check('free tail closes energy and momentum',max(d['max_energy_relative_error'],d['max_momentum_relative_error'],
          plain['max_energy_relative_error'],plain['max_momentum_relative_error'])<1e-12)
    check('daughter directional stress gives massless trace',
          abs(d['final_daughter_pressure_x']+2*d['final_daughter_pressure_perp']-d['final_daughter_energy'])<1e-12)
    check('early cut fails pressure with exact energy',not d['trials'][0]['accepted'] and
          d['trials'][0]['maximum_downstream_relative_errors']['energy_total']<1e-12 and
          max(d['trials'][0]['maximum_downstream_relative_errors'][p] for p in ('P_x','P_perp'))>1e-6)
    check('three predefined moving boundaries pass',len(d['trials'])==4 and all(
        t['accepted'] and t['surviving_parent_energy_fraction']>0.5 and
        max(t['maximum_downstream_relative_errors'].values())<1e-6 for t in d['trials'][1:]))
    check('moving boundaries give same downstream transport',d['max_difference_between_accepted_boundaries']<1e-12)
    check('no-elastic pressure bound fails through half energy',
          plain['boundary'] is None and plain['pressure_envelope_at_parent_half_time']>1e-6)
    check('full overlap covariance Pauli admissible',all(
        row['covariance_eigenvalue_range'][0]>-1e-10 and row['covariance_eigenvalue_range'][1]<1+1e-10
        for model in (d,plain) for row in model['full_covariance_controls']))
    check('full matrix overlap pressure reproduces mode formula',max(
        row['pressure_reference_relative_error'] for model in (d,plain) for row in model['full_covariance_controls'])<1e-11)
    check('local spatial coherence scope retained',not transport['local_transport_adequacy_claimed'] and
          d['full_covariance_controls'][-1]['remaining_same_band_spatial_coherence_relative_norm']>0)
    h55=transport['inherited_H55_interface']
    check('H55 energy pressure and work fail',not h55['tested_boundary_pass']
          and abs(h55['dressed_energy_gap'])>1e-6 and abs(h55['free_dephasing_pressure_gap'])>1e-6
          and abs(h55['normalized_fermion_energy_derivative_jump'])>1e-3)
    check('H55 terminal refinement miss retained',not h55['terminal_refinement_pass'])
    check('test interactions never applied to H55',not h55['applied_control_tail_to_H55'])
    check('scientific acceptance remains false',not (source['scientific_acceptance'] or
          follow['scientific_acceptance'] or transport['scientific_acceptance']))

    inputs=read('inputs.json')
    check('X5 parent pinned',inputs['parent']==PARENT)
    check('inherited inputs unchanged',all(
        hashlib.sha256((REPO/path).read_bytes()).hexdigest()==digest for path,digest in inputs['sha256'].items()))
    scan=json.loads((TRACKS/'x2_strip/results/blind_verdict.json').read_text())
    check('original strict result zero of 135',not scan['strict_M2_X_keys'] and all(
        x['grid_points']==135 for x in scan['nodes'].values()))
    if args.replay:
        fresh=read('source.json',args.replay)
        check('fresh source replay all work summaries',all(
            math.isclose(s['integrated_transfer_matched'],fresh['cases'][k]['integrated_transfer_matched'],
                         rel_tol=1e-8,abs_tol=1e-12) for k,s in source['cases'].items()))
        final=read('source_followup.json',args.replay)
        check('fresh followup replay all work summaries',all(
            math.isclose(s['integrated_transfer_matched'],final['cases'][k]['integrated_transfer_matched'],
                         rel_tol=1e-8,abs_tol=1e-12) for k,s in follow['cases'].items()))
        end=read('transport.json',args.replay)['controlled_elastic']
        check('fresh transport boundary and downstream energy',end['boundary']==d['boundary'] and
              math.isclose(end['final_daughter_energy'],d['final_daughter_energy'],rel_tol=1e-10))
        benchmark=read('controls.json',args.replay)
        check('fresh independent ODE comparisons',benchmark['helicity'][-1]['mode_occupation_max_absolute']<1e-5
              and benchmark['spatial']['ode'][-1]['max_absolute_error']<1e-5)
        replay_edge=read('edge.json',args.replay)
        check('fresh sharp-edge decomposition',all(math.isclose(x,y,rel_tol=1e-10) for x,y in zip(
            replay_edge['largest_cylinder_edge_ratios'],ratios,strict=True)))
    result=dict(passed=sum(checks.values()),total=len(checks),checks=checks,
                finite_source_cutoff_gate_pass=False,controlled_averaged_transport_gate_pass=True,
                inherited_H55_interface_gate_pass=False,scientific_acceptance=False,strict_accepted=0,strict_tested=135)
    if args.output:
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(f'{result["passed"]}/{result["total"]} evidence checks passed')
    if not all(checks.values()):
        raise SystemExit(', '.join(name for name,value in checks.items() if not value))


if __name__=='__main__':
    main()
