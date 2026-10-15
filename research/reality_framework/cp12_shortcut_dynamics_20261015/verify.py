"""Saved-evidence verification; full reproduction commands are in README.md."""
from pathlib import Path
import argparse,hashlib,json

P=Path(__file__).resolve().parent
FROZEN=P.parent/'cp12_ground_connectivity'
EXPECTED={'frozen_model.py':'e896fce8a3e61948ff41891576b4c3ac028ad26785cafabb30bddef93c5adb9a',
    'exact_destination.py':'1a1ccb4f06602179f14c2834dea3566fada022ca6a23f58370abe0d4f1e5b35b'}


def read(name):return json.loads((P/name).read_text())


def overlaps(a,b):return max(a[0],b[0])<=min(a[1],b[1])+1e-12


def main():
    a=argparse.ArgumentParser();a.add_argument('--write-manifest',action='store_true');args=a.parse_args()
    for name,sha in EXPECTED.items():assert hashlib.sha256((FROZEN/name).read_bytes()).hexdigest()==sha
    checks=read('CHECKS.json');algebra=read('TRANSFER_ALGEBRA.json')
    assert checks['status']==algebra['status']=='PASS'
    assert checks['full_basis_hops_and_Gauss_checks']==24636
    assert algebra['integer_CAR_identity_all_eight_states']
    assert algebra['unregulated_destinations_outside_check_basis']>0
    parent=read('PARENT_TOY_REPRODUCTION.json');assert parent['status']=='PASS'
    assert abs(parent['rows'][2]['full_opposite_prob']-.162306681581)<1e-11
    n4a=read('TRANSPORT_N4_Q5_B128.json');n4b=read('TRANSPORT_N4_Q5_B256.json')
    assert n4b['wait_state_error']<1e-8
    assert n4b['event_probability_interval'][0]>.05
    assert overlaps(n4a['event_probability_interval'],n4b['event_probability_interval'])
    assert abs(sum(n4b['record_outcome_probabilities'])-1)<1e-10
    cutoffs=[];resolved=0
    for small,large in zip(n4a['results'],n4b['results']):
        assert (small['lag'],small['shortcut_present_at_record'])==(large['lag'],large['shortcut_present_at_record'])
        assert overlaps(small['initial_joint_probability_interval'],large['initial_joint_probability_interval'])
        for model in ['full_frozen','recorded_ring_hopping_only']:
            x=small['models'][model];y=large['models'][model]
            assert overlaps(x['receiver_response']['joint_interval'],y['receiver_response']['joint_interval'])
            assert overlaps(x['new_chord_and_ring_joint_interval'],y['new_chord_and_ring_joint_interval'])
            assert y['receiver_response']['phase_shift_identity_error']<1e-12
            if model=='recorded_ring_hopping_only':assert y['new_chord_and_ring_joint_probability']<1e-24
        assert large['state_distance_interval'][0]<=large['full_law_Duhamel_upper']+1e-12
        assert large['new_chord_and_ring_probability_T4_coefficient_conditional_interval'][0]>.13
        if large['shortcut_present_at_record'] and large['lag'] in [.125,.25,.5]:
            assert large['response_change_conditional_interval'][1]<0
            assert large['models']['full_frozen']['receiver_response']['conditional_interval'][1]<0
            assert large['models']['recorded_ring_hopping_only']['receiver_response']['conditional_interval'][0]>0
            resolved+=1
        cutoffs.append(dict(lag=large['lag'],direct_shortcut=large['shortcut_present_at_record'],
            absolute_full_response_cutoff_difference=abs(small['models']['full_frozen']['receiver_response']['joint_value']-
                large['models']['full_frozen']['receiver_response']['joint_value'])))
    n5=read('SPLIT_N5_B50.json');n6=read('SPLIT_N6_B12.json');resp=read('RESPONSE_N5_B50.json')
    for d in [n5,n6,resp]:
        assert d['status']=='PASS'
        c=d['model_checks'];assert c['neutral_core_generator_residual']<1e-10
        assert abs(c['neutral_variance']-2*d['Q'])<1e-9
        assert abs(c['neutral_outside_hopping_square']-2*d['Q'])<1e-9
    for d in [n5,n6]:
        assert d['model_checks']['H_hermiticity']<1e-10
        assert d['model_checks']['V_hermiticity']<1e-10
        for row in d['post_formation_comparisons']:
            assert row['state_distance_interval'][0]<=min(row['exact_full_vs_exact_core_Duhamel_upper'],
                row['exact_full_vs_exact_core_combined_driver_upper'])+1e-10
    assert n5['dimension']==162191 and n6['dimension']==715985 and resp['dimension']==810377
    previous=json.loads((P.parent/'charged_density_dynamics_20261015/RESPONSE_N5_B50.json').read_text())
    reproduction=[]
    for row in resp['results']:
        if row['preparation']=='unfiltered_waiting_state':assert row['response_change_joint_interval'][0]>0
        if row['lag']==.25:
            key='charged_initial_event_response' if row['preparation']=='initial_event_branch' else 'unfiltered_response'
            diff=abs(row['models']['full_frozen']['receiver_response']['joint_value']-previous[key]['joint_susceptibility'])
            assert diff<1e-10;reproduction.append(dict(preparation=row['preparation'],absolute_difference=diff))
            if row['preparation']=='initial_event_branch':
                assert row['response_change_joint_interval'][0]<0<row['response_change_joint_interval'][1]
    summary=dict(status='PASS',frozen_sources=EXPECTED,independent_controls_status=checks['status'],
        exact_transfer_algebra_status=algebra['status'],resolved_shortcut_response_sign_reversals=resolved,
        N4_electric_budget_comparisons=cutoffs,previous_full_N5_response_reproduction=reproduction,
        full_reproduction='Commands in README.md rebuild the matrices and repeat all propagations. This verifier checks saved evidence and integrity only.',
        limitation='Internal analytic/numerical verification using ordinary floating point; no outward rounding, external peer review, continuum limit or autonomous geometric record.')
    if args.write_manifest:
        (P/'VERIFICATION.json').write_text(json.dumps(summary,indent=2)+'\n')
        manifest={f.name:dict(sha256=hashlib.sha256(f.read_bytes()).hexdigest(),bytes=f.stat().st_size)
            for f in sorted(P.iterdir()) if f.is_file() and f.name!='manifest.json'}
        (P/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    else:
        manifest=read('manifest.json')
        for name,info in manifest.items():
            f=P/name;assert f.stat().st_size==info['bytes']
            assert hashlib.sha256(f.read_bytes()).hexdigest()==info['sha256'],name
        actual={f.name for f in P.iterdir() if f.is_file() and f.name!='manifest.json'}
        assert actual==set(manifest)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
