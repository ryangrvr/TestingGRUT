"""Build the scoped campaign ledger from saved fresh-worker observations."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json
import enclosures as E

BASE=Path(__file__).resolve().parent


def read(name):return json.loads((BASE/name).read_text())


def combined(rows,weights):
    if not all(x['status']=='COMPLETE' and x['physical_probability_intervals'] for x in rows):return None
    return {event:tuple(E.rational_bounds(sum(w*Q.from_float(x['physical_probability_intervals'][event][j])
        for w,x in zip(weights,rows)))[j] for j in range(2)) for event in ['macro','full']}


def resources(rows):
    return dict(total_sequential_worker_wall_seconds=sum(x['whole_worker_wall_seconds'] for x in rows),
        total_sequential_worker_cpu_seconds=sum(x['whole_worker_cpu_seconds'] for x in rows),
        maximum_worker_peak_resident_MiB=max(x['peak_resident_MiB'] for x in rows))


def main():
    controls=read('controls.json');assert controls['status']=='PASS'
    new=[read('results/s7_N5_root0.json'),read('results/s7_N5_root1.json')]
    old=[read('results/s5_N5_root0.json'),read('results/s5_N5_root1.json')]
    interval=combined(new,[Q(3,5),Q(2,5)])
    passed=bool(interval and interval['macro'][0]>=.5 and all(x['finite_instance_gate'] for x in new))
    q=Q(21,1000)
    ex=sum((q**j/factorial(j) for j in range(17)),Q(0))
    exhi=ex+(q**17/factorial(17))/(1-q/18)
    kappa=Q(1,800)-(exhi-1-q)
    joined=(E.rational_bounds(Q.from_float(interval['macro'][0])*kappa)[0] if interval else None)
    recovery=dict(N=5,physical_receiver_modes=1600,target_time=100,root_weights=['3/5','2/5'],
        gate_pass=passed,macro_probability_interval=(interval['macro'] if interval else None),
        full_assembly_probability_interval=(interval['full'] if interval else None),
        resources=resources(new),maximum_component_state_error=max(x['analytic_state_error_account'] for x in new),
        S4_join=dict(near_response_per_original_trial_lower=joined,probe_time=.005,
            source_variance=.25,conditional_near_response_lower=E.rational_bounds(kappa)[0],
            actual_conditional_receiver_retained=True,rejected_trials_zero_signal=True,
            recording_and_probe_injection_supplied=True,full_physical_probe_simulation=False))
    baseline_interval=combined(old,[Q(3,5),Q(2,5)])
    baseline=dict(resources=resources(old),completed_components=sum(x['status']=='COMPLETE' for x in old),
        probability_intervals=baseline_interval,arithmetic='immutable S5 ordinary double; no roundoff certificate',
        equal_error_speed_comparison=False)
    six=None
    if (BASE/'results/s7_N6_root0.json').exists():
        assert passed
        x=read('results/s7_N6_root0.json')
        six=dict(status=x['status'],stop_reason=x['stop_reason'],target_time=144,
            endpoint_time=x['completed_physical_steps_time'],receiver_modes=2304,
            probability_intervals=x['physical_probability_intervals'],
            one_component_and_proved_mirror_reuse=True,root_weights=['1/2','1/2'],
            target_formation_gate=x['finite_instance_gate'],resources=resources([x]))
    result=dict(recovery=recovery,fresh_S5_baseline=baseline,six_event=six,
        controls_pass=True,unchanged_S5_S6_archives=True,
        general_history_conditioned_action_enrichment_implemented=False,
        shared_group_action_enclosure_implemented=True,
        arithmetic='roundoff-aware under declared IEEE/reduction model; not formal software certification',
        all_size_formation_proved=False,full_independent_disorder_gate_pass=False,
        selected_dimension=False,original_matter_universality=False,
        real_energy_gap_restored=False,distinctive_empirical_prediction=False,
        external_independent_validation=False)
    if (BASE/'WEIGHTED_CONTROLS.json').exists():
        weighted=read('WEIGHTED_CONTROLS.json')
        result['post_campaign_weighted_residual']=dict(small_controls=weighted['status'],
            used_in_S7_trajectories=False,retrospective_gate_change=False,
            six_event_completion_demonstrated=False)
    (BASE/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
