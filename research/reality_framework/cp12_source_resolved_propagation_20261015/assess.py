"""Derive the operational finite-test decisions from saved interval data."""
import hashlib,json
from pathlib import Path
P=Path(__file__).resolve().parent


def read(name):return json.loads((P/name).read_text())


def absolute_ratio(a,b):
    return None if b[0]<=0 else [a[0]/b[1],a[1]/b[0]]


def main():
    files={4:'RESPONSE_N4_Q5_B256_W5_G200000_F2_L12.json',
        5:'RESPONSE_N5_Q6_B256_W5_G200000_F2_L10.json',
        6:'RESPONSE_N6_Q7_B192_W5_G200000_F3_L8.json'}
    formation=[];comparisons=[]
    for n,name in files.items():
        d=read(name);initial=d['initial']
        formation.append(dict(N=n,Q=d['Q'],all_spanning_ring_probability=initial['all_spanning_cycle_probability'],
            longest_chord_class_probability=initial['original_total_formation_probability'],
            longest_chord_class_probability_interval=initial['original_total_formation_interval'],
            fixed_source_endpoint_probability=initial['fixed_source_endpoint_class_probability'],
            source_variance_interval=d['source_variance_interval']))
        if n==4:continue
        for lag in [.125,.25]:
            r=next(r for r in d['rows'] if r['lag']==lag)
            ja,jb=(2,3) if n==5 else (3,2)
            aa=r['fixed_receiver_pairs'][str(ja)];bb=r['fixed_receiver_pairs'][str(jb)]
            a,b=aa['density_response'],bb['density_response']
            margin=a['absolute_interval'][0]-b['absolute_interval'][1]
            comparisons.append(dict(N=n,lag=lag,receiver_a=ja,receiver_b=jb,
                core_distance_a=aa['core_distance'],core_distance_b=bb['core_distance'],
                response_a_interval=a['interval'],response_b_interval=b['interval'],
                ratio_absolute_intervals=absolute_ratio(a['absolute_interval'],b['absolute_interval']),
                absolute_dominance_margin=margin,
                same_core_probability_interval=r['same_core_probability_interval'],
                bypass_ordering_resolved=margin>0,
                survival_benchmark_passed=r['same_core_probability_interval'][0]>.75))
    primary=read('RESPONSE_N5_Q10_B64_W2_G200000_L6.json')
    data=dict(computation_status='PASS',parent_commit='4e7d6c5a91fe96737ab9a8a140ec2b50eb140ebc',
        protocol_sha256=hashlib.sha256((P/'TEST_PROTOCOL.md').read_bytes()).hexdigest(),
        prospective_scope='Six-site comparison fixed locally before that source result; N4/N5 pilots already inspected. Not independent public preregistration.',
        auxiliary_sequence='Q=N+1, original equal coherent neutral preparation in each sector; frozen coefficients unchanged.',
        formation_and_source=formation,fixed_receiver_comparisons=comparisons,
        primary_Q2N_supplement=dict(N=5,Q=10,wait=2,
            exact_spanning_ring_probability=primary['initial']['original_total_formation_probability'],
            source_variance_interval=primary['source_variance_interval'],
            rows=[dict(lag=r['lag'],same_core_probability_interval=r['same_core_probability_interval'],
                near_response_interval=r['shells']['1']['density_response']['interval'],
                far_response_interval=r['shells']['2']['density_response']['interval']) for r in primary['rows']]),
        finite_candidate_metric_verdict='FAIL: the specified six-site recorded-ring branch fails the prospective distance-ordering benchmark with a resolved active-host bypass.',
        formation_verdict='The six-site ring and longest-chord class both fall below the 1% finite-test benchmark. No appreciable spatial-phase scaling is established.',
        emergent_spatial_phase_established=False,growing_distance_front_established=False,
        no_go_for_all_CP12_phases=False,primary_Q2N_large_size_scaling_established=False,
        gravitational_continuum_established=False,distinctive_empirical_prediction_established=False,
        acceptance='Milestone A not achieved. No claim for Milestones B or C. Preserve the mathematics; reconsider the metric/transport architecture under a separate prospective test.',
        limitation='A signed susceptibility is not an arrival probability. Interference and initial correlations can break monotonic ordering in local systems; the conclusion concerns this specified ring metric and active-chord experiment, not all possible collective geometries. Geometry/host recording is a measured intervention with nonzero work and response disturbance.')
    (P/'ASSESSMENT.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))


if __name__=='__main__':main()
