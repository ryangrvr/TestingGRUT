#!/usr/bin/env python3
"""Check the X11 evidence while preserving the failed physical gates."""
import argparse
import hashlib
import json
import math
from pathlib import Path

TRACK=Path(__file__).resolve().parents[1]
RESULTS=TRACK/'results'
TRACKS=TRACK.parent


def read(folder,name):
    return json.loads((folder/f'{name}.json').read_text())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--results',type=Path,default=RESULTS)
    parser.add_argument('--axial-replay',type=Path)
    args=parser.parse_args()
    production=read(args.results,'production')
    continuation=read(args.results,'continuation')
    law=read(args.results,'axial_law')
    cases=production['cases']
    checks={}

    def check(name,condition):
        checks[name]=bool(condition)

    for name,sha in (production['input_sha256']|law['input_sha256']).items():
        check(f'pinned input {name}',hashlib.sha256((TRACKS/name).read_bytes()).hexdigest()==sha)

    fail=cases['cutoff']
    check('wide cutoff fails at the original step',not fail['completed'] and
          -0.02<fail['last_recorded_time']<0 and 'did not converge' in fail['reason'])
    for name,case in cases.items():
        if name=='cutoff':
            continue
        check(f'{name} total energy ledger',case['max_total_energy_ledger_error']<1e-5)
        check(f'{name} occupation positivity and Bloch norm',case['occupation_range'][0]>-1e-10 and
              case['occupation_range'][1]<1+1e-10 and case['max_bloch_norm_squared_error']<1e-8)
        if case['moment_projection']:
            check(f'{name} native projection',all(row['native_moment_relative_error']<1e-10
                                                  for row in case['moment_projection']))
        if 'max_reciprocal_work_ledger_error' in case:
            measured=production['work_from_quantum_stress'][name]['integrated_quantum_stress_transfer']
            check(f'{name} reciprocal work closes on the integration steps',
                  case['max_reciprocal_work_ledger_error']<1e-8 and
                  math.isclose(measured,case['integrated_comoving_scalar_work'],rel_tol=2e-5,abs_tol=1e-9))
    check('direct X5 homogeneous replay agrees',production['independent_homogeneous']['max_error']<1e-10)
    check('full projected late source reaches the archived diagnostic time',
          math.isclose(cases['late']['terminal']['t'],2.9980937,abs_tol=1e-10))
    check('later phase-sensitive work and pressure are still present',
          abs(cases['late']['terminal']['signed_work_over_density'])>1e-2 and
          abs(cases['late']['terminal']['dephased_pressure_jump'])>1e-6)
    delta=production['comparisons']
    check('bulk timestep stable while radial count remains unresolved',
          abs(delta['timestep']['terminal_relative_changes']['N_fraction'])<.001 and
          abs(delta['radial']['terminal_relative_changes']['N_fraction'])>.01)
    check('cutoff comparison controls timestep independently',
          abs(cases['half_step']['timestep']-cases['cutoff_half_step']['timestep'])<1e-15 and
          abs(delta['cutoff']['terminal_relative_changes']['N_fraction'])<.002 and
          abs(delta['cutoff']['terminal_signed_work_difference_over_baseline_magnitude'])>.5)
    check('small spatial cubes omit the bulk of seed gradient power',
          cases['reference']['seed_fractions']['chi']['gradient']<
          cases['spatial']['seed_fractions']['chi']['gradient']<1e-3 and
          cases['scalar_fine']['seed_fractions']['chi']['gradient']>.05)
    check('waterfall spans multiple cells of the archived scalar lattice',
          production['travel_in_archived_scalar_cells']>10 and
          production['spatially_controlled_H55_production'] is False)

    stages=continuation['cases']
    kinetic, relic, baryon=(stages[k] for k in ('kinetic','abundance','baryon'))
    for name,row in kinetic.items():
        splice=row['initial_splice']
        check(f'{name} native number and energy carry through the kinetic input',
              splice['native_moment_max_relative_error']<1e-10 and abs(splice['energy_gap'])<1e-6 and
              row['max_kinetic_exchange_error']<1e-8)
        check(f'{name} pressure or work prevents a physical splice',
              max(abs(splice['complete_pressure_jump_over_density']),
                  splice['max_abs_sampled_residual_work_last_window'])>1e-6 and
              splice['passed'] is False)
    check('refined kinetic clock reduces pressure-work quadrature error',
          kinetic['at055']['max_kinetic_budget_error']<1e-5 and
          kinetic['at055_half']['max_kinetic_budget_error']<kinetic['at055']['max_kinetic_budget_error']/3 and
          kinetic['legacy_at055']['max_kinetic_budget_error']>1e-4)
    check('entropy result is stable under halving the kinetic timestep',
          abs(kinetic['at055_half']['visible_entropy_comoving_Hi3']/
              kinetic['at055']['visible_entropy_comoving_Hi3']-1)<.01)
    for name,row in relic.items():
        pop=row['population']
        check(f'{name} release clock matches its scale-factor history',
              row['release']['H_consistency']<1e-7 and row['release']['charge_balance']<1e-6)
        check(f'{name} relic receives the calculated charge with fixed scale',
              row['H55_scale_GeV']==5.5e10 and row['holder_c0']==1.55 and
              row['holder_strength_implementation']=='scaled coupling at fixed archived spectator amplitude' and
              not pop['Y_Q_calibrated'] and math.isclose(row['Y_charge'],
              abs(row['release']['N_X'])/row['s_comoving_visible_Hi3'],rel_tol=1e-12))
        check(f'{name} abundance comparison preserves the miss',
              0<row['conditional_Omega_h2']<row['target_Omega_h2'] and
              not row['source_handoff_passed'] and not row['accepted_abundance'])
    check('direct same-source kinetic timestep abundance comparison',
          abs(continuation['kinetic_step_relative_abundance_difference'])<.01)
    check('three spatial-model boundary calculations are present',
          set(continuation['boundary_relative_abundance_difference'])=={'at035','late'})
    check('late boundary and radial refinement change the conditional yield above tolerance',
          abs(continuation['boundary_relative_abundance_difference']['late'])>.01 and
          abs(relic['radial055']['conditional_Omega_h2']/relic['at055']['conditional_Omega_h2']-1)>.01)

    literal=law['literal_positive_parity']
    check('real positive-parity law fails the mass-sum benchmark for any real shape',
          literal['sum_light_masses_eV']>=literal['mass_sum_trace_bound_eV']>.19 and
          literal['mass_sum_trace_bound_eV']>law['Planck_BAO_mass_sum_limit_eV'])
    check('literal one-loop benchmark also exceeds the mass-sum comparison',
          sum(literal['repository_loop_estimate']['sv_tree_plus_loop_eV'])>law['Planck_BAO_mass_sum_limit_eV'])
    check('common phase and phase conjugation controls',max(law['tests'].values())<1e-11)
    inherited=law['inherited_parity_later_widths']
    inv=inherited['normalized_CP_invariants']
    check('axial law suppresses invariants inside the waterfall triplet',
          max(abs(inv[i][j]) for i in range(1,4) for j in range(1,4))<1e-12 and
          min(abs(inv[0][j]) for j in range(1,4))>.24)
    check('parity extension approaches the reference mass sum only at the later width',
          .058<inherited['sum_light_masses_eV']<.060 and
          law['inherited_parity_native_widths']['sum_light_masses_eV']>.09)
    check('fixed phases miss the conditional baryon yield across washout controls',
          all(abs(row['ratio_to_observed'])>10 for row in baryon['clock_later']['washout_cases'].values()))
    check('strict inherited verdict and scientific scope',
          continuation['strict_candidate_ledger']=='0/135' and not continuation['accepted_abundance'] and
          not continuation['accepted_GRUT_law'] and not law['accepted_GRUT_law'])
    if args.axial_replay:
        check('fresh exact axial matching replay',json.loads(args.axial_replay.read_text())==law)

    out=dict(passed=sum(checks.values()),total=len(checks),checks=checks,
             accepted_abundance=False,accepted_GRUT_law=False,strict_candidate_ledger='0/135')
    if args.output:
        args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(f"X11 evidence checks {out['passed']}/{out['total']}")
    for name,passed in checks.items():
        if not passed:
            print('FAIL',name)
    raise SystemExit(out['passed']!=out['total'])


if __name__=='__main__':
    main()
