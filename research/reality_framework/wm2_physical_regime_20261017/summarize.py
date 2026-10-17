"""Root-weighted finite outcome and independent analytic-claim ledger."""
from pathlib import Path
from fractions import Fraction as Q
import json
import sys
sys.path.append(str(Path(__file__).resolve().parent.parent/'wm2_schmidt_receiver_20261017'))
import bounds as E

BASE=Path(__file__).resolve().parent


def run():
    rows=[]
    for i in range(2):
        p=BASE/f'results/S9_N5_root{i}.json'
        rows.append(json.loads(p.read_text()) if p.exists() else None)
    all_present=all(r is not None for r in rows)
    complete=all_present and all(r['status']=='COMPLETE' for r in rows)
    intervals=None
    if complete:
        intervals={name:[E.OLD.rational_bounds(sum(w*Q.from_float(r['physical_probability_intervals'][name][j])
            for w,r in zip([Q(3,5),Q(2,5)],rows)))[j] for j in range(2)] for name in ['macro','full']}
    gate=bool(complete and all(r['finite_instance_gate'] for r in rows) and intervals['macro'][0]>=.5)
    A=dict(status=('PASS' if gate else 'FAIL') if all_present else 'PENDING',gate_pass=gate,
        original_root_weights=['3/5','2/5'],event_intervals=intervals,
        components=[None if r is None else dict(status=r['status'],counts=r['nonroot_type_counts'],
            target_time=100,completed_step_time=r['completed_physical_steps_time'],
            whole_step_state=r['whole_step_state'],stop_reason=r['stop_reason'],
            physical_probability_intervals=r['physical_probability_intervals']) for r in rows])
    present=[r for r in rows if r is not None]
    if present:
        A['resources']=dict(sequential_wall_seconds=sum(r['whole_worker_wall_seconds'] for r in present),
            sequential_CPU_seconds=sum(r['whole_worker_cpu_seconds'] for r in present),
            maximum_worker_peak_MiB=max(r['peak_resident_MiB'] for r in present),worker_limit_seconds=1200,
            equal_error_or_budget_speed_claim=False)
    B=dict(status='BUILDER_DERIVED_NOT_EXTERNALLY_REVIEWED',uniform_independent_storage_disorder=True,
        N0=2**28,N0_is_sufficient_not_minimal=True,formation_lower='1/2',
        frozen_protocol='uniform time sampling on [4N^2,8N^2]',
        stronger_pointwise_same_window=True,physical_receiver_modes='64N^2',
        threshold_receiver_modes=2**62,threshold_earliest_time=2**58,
        current_comparison_is_reduced_state_Markov_theorem=False,
        unresolved_moderate_sizes=True,all_local_perturbations_covered=False,
        integrated_original_trial_joint_arrival_lower='1/589824',
        arrival_is_receiver_set_not_far_single_site_phase_derivative=True)
    r=dict(parent='cb361f8f230371515f1806ab50ca62657a3088be',
        freeze='85fc1bd2c4592a57c814bc5b2879aab65fdb85aa',A=A,B=B,
        numerical_branch_closed=all_present,no_fallback=True,six_event='NOT_RUN_NOT_IN_S9',
        selected_dimension=False,original_CP12_matter_universality=False,
        restored_real_energy_gap=False,continuum_gravity=False,
        distinctive_empirical_prediction=False,external_validation=False)
    (BASE/'summary.json').write_text(json.dumps(r,indent=2)+'\n')
    return r


if __name__=='__main__':print(json.dumps(run(),indent=2))
