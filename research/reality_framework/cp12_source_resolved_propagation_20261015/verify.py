"""Verify saved evidence, frozen provenance and physical acceptance separately."""
import argparse,hashlib,json,math,platform,subprocess
from pathlib import Path
import numpy as np
import scipy
import matplotlib

P=Path(__file__).resolve().parent;RF=P.parent;ROOT=RF.parent.parent
BASE='4e7d6c5a91fe96737ab9a8a140ec2b50eb140ebc'
PREFIX=P.relative_to(ROOT).as_posix()+'/'


def read(name):return json.loads((P/name).read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def contains(interval,value):return interval[0]<=value<=interval[1]


def finite(value):
    if isinstance(value,dict):
        for v in value.values():finite(v)
    elif isinstance(value,list):
        for v in value:finite(v)
    elif isinstance(value,float):assert math.isfinite(value)


def response_check(r):
    assert contains(r['interval'],r['value'])
    assert r['error_bound']>=0 and r['exact_finite_phase_identity_error']<1e-11
    assert r['zero_excluded']==(abs(r['value'])>r['error_bound'])
    I=r['interval'];A=r['absolute_interval']
    expect=[0. if I[0]<=0<=I[1] else min(abs(v) for v in I),max(abs(v) for v in I)]
    assert np.max(np.abs(np.array(A)-expect))<1e-15


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--refresh',action='store_true',
        help='Refresh the content manifest after intentionally regenerating evidence.')
    args=ap.parse_args()
    frozen={'frozen_model.py':'e896fce8a3e61948ff41891576b4c3ac028ad26785cafabb30bddef93c5adb9a',
        'exact_destination.py':'1a1ccb4f06602179f14c2834dea3566fada022ca6a23f58370abe0d4f1e5b35b'}
    for name,digest in frozen.items():assert sha(RF/'cp12_ground_connectivity'/name)==digest
    dependencies=['cp12_ground_connectivity/frozen_model.py','cp12_ground_connectivity/exact_destination.py',
        'cp12_ground_connectivity/invariant_evolution.py',
        'charged_density_dynamics_20261015/fast_symmetry.py','charged_density_dynamics_20261015/symmetry.cpp',
        'charged_density_dynamics_20261015/stream.py',
        'cp12_shortcut_dynamics_20261015/operators.py','cp12_shortcut_dynamics_20261015/integrals.py',
        'cp12_dynamical_suppression_20261015/boundary.py',
        'cp12_dynamical_suppression_20261015/linear_propagator.py',
        'cp12_dynamical_suppression_20261015/resources.py']
    deps={}
    for name in dependencies:
        p=RF/name;rel=p.relative_to(ROOT).as_posix()
        assert p.read_bytes()==subprocess.check_output(['git','show',f'{BASE}:{rel}'],cwd=ROOT)
        deps[rel]=sha(p)
    changed=subprocess.check_output(['git','diff','--name-only',BASE],cwd=ROOT,text=True).splitlines()
    untracked=subprocess.check_output(['git','ls-files','--others','--exclude-standard'],cwd=ROOT,text=True).splitlines()
    assert all(p.startswith(PREFIX) for p in changed+untracked)
    assert sha(P/'TEST_PROTOCOL.md')=='e33d941f26716e2859f3b20460e467765aa69a33d4d2d3a2fff8642bfda3d10a'
    for p in P.glob('*.py'):compile(p.read_text(),str(p),'exec')
    for p in P.glob('*.json'):finite(json.loads(p.read_text()))
    response_count=0
    for p in P.glob('RESPONSE_*.json'):
        d=json.loads(p.read_text());assert d['status']=='PASS' and d['source_variance_interval'][0]>0
        assert d['build_checks']['H_hermiticity']<1e-8
        assert contains(d['source_variance_interval'],d['source_variance'])
        for row in d['rows']:
            assert contains(row['same_core_probability_interval'],row['same_core_probability'])
            response_check(row['same_core_response']);response_count+=1
            for r in row['shells'].values():
                for key in ['density_response','density_and_same_core_response']:
                    response_check(r[key]);response_count+=1
                if r['density_given_same_core_response_interval'] is not None:
                    assert contains(r['density_given_same_core_response_interval'],r['density_given_same_core_response'])
            for receiver,r in row.get('fixed_receiver_pairs',{}).items():
                assert d['fixed_core_distances'][int(receiver)]==r['core_distance']
                for key in ['density_response','density_and_same_core_response']:
                    response_check(r[key]);response_count+=1
    for p in P.glob('NEUTRAL_*.json'):
        d=json.loads(p.read_text());assert d['status']=='PASS'
        for row in d['rows']:
            assert contains(row['formation_interval'],row['formation_probability'])
            assert contains(row['record_energy_interval'],row['record_energy_change'])
            assert contains(row['record_entropy_interval_bits'],row['record_entropy_bits'])
    for p in P.glob('REFINED_RECORD_*.json'):
        d=json.loads(p.read_text());assert d['status']=='PASS'
        assert abs(sum(r['class_probability'] for r in d['classes'])-d['all_cycle_probability'])<1e-13
        assert sum(r['equivalent_records'] for r in d['classes'])+1==d['complete_outcomes']
        assert contains(d['complete_record_work_interval'],d['complete_record_work'])
        assert contains(d['complete_record_entropy_interval_bits'],d['complete_record_entropy_bits'])
        assert d['raw_rewiring_transition_count']>0
    for p in P.glob('GEOMETRY_*.json'):
        d=json.loads(p.read_text());assert d['status']=='PASS'
        assert abs(sum(v['probability'] for v in d['circulation_edge_count_distribution'].values())-1)<1e-10
        assert d['exact_spanning_cycle_probability']<=d['connected_spanning_circulation_probability']+1e-12
        assert d['connected_spanning_circulation_probability']<=d['all_vertices_covered_probability']+1e-12
        for row in d['all_root_mean_balls_given_spanning_connected_circulation']:
            assert contains(row['conditional_mean_interval'],row['conditional_mean'])
    checks=read('CHECKS.json');refined=read('REFINED_CHECKS.json');krylov=read('KRYLOV_CHECK.json')
    assert all(d['status']=='PASS' for d in [checks,refined,krylov])
    assert checks['raw_record_transition_checks']==68376 and checks['cross_record_H_transitions']==0
    assert checks['rooted_initial_lift_difference']<3e-11 and refined['endpoint_lift_difference']<3e-11
    assert refined['work_identity_difference']<1e-11 and refined['entropy_identity_difference']<1e-11
    for r in refined['labelled_response_checks']:assert r['central_difference']<1e-11
    assert refined['complete_record_effect_on_fixed_receiver'][0]['response_change_interval'][1]<0
    assert len(krylov['rows'])==3
    for r in krylov['rows']:
        assert r['full_basis_difference']<=r['krylov_full_state_bound']+r['independent_full_basis_bound']+5e-10
        assert r['orthogonality_error']<1e-10 and r['boundary_QR_relative_norm_check']<1e-11
    assert krylov['rows'][-1]['N']==6 and krylov['rows'][-1]['full_basis_difference']<1e-10
    neutral=read('NEUTRAL_N6_Q7_B192.json')
    assert neutral['dimension']==1980546 and neutral['outside_dimension']==133348
    assert neutral['rows'][0]['state_error']<3e-9
    verdict=read('ASSESSMENT.json');assert verdict['computation_status']=='PASS'
    six=[r for r in verdict['fixed_receiver_comparisons'] if r['N']==6]
    assert len(six)==2 and all(r['bypass_ordering_resolved'] and r['survival_benchmark_passed'] for r in six)
    quarter=next(r for r in six if r['lag']==.25)
    assert quarter['ratio_absolute_intervals'][0]>550
    assert verdict['formation_and_source'][-1]['all_spanning_ring_probability']<.01
    assert not verdict['emergent_spatial_phase_established'] and not verdict['no_go_for_all_CP12_phases']
    assert not verdict['distinctive_empirical_prediction_established']
    assert (P/'SOURCE_RESPONSE.svg').stat().st_size>1000
    for name in ['REPORT.md','HANDOFF.md','README.md','DERIVATION.md','PROGRAM_TESTS.md','ARCHITECTURE_CHARTER.md']:
        assert len((P/name).read_text())>1000 and 'TO_BE_INSERTED' not in (P/name).read_text()
    files=[p for p in sorted(P.iterdir()) if p.is_file() and p.name not in ['manifest.json','VERIFICATION.json']]
    assert all(p.suffix in ['.md','.py','.json','.svg'] for p in files)
    manifest=dict(parent_commit=BASE,parent_tree=subprocess.check_output(['git','rev-parse',f'{BASE}^{{tree}}'],cwd=ROOT,text=True).strip(),
        packet=PREFIX.rstrip('/'),files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},
        inherited_dependency_sha256=deps,
        exclusions='manifest.json and VERIFICATION.json excluded to avoid self-reference; disposable reproducible numerical caches excluded.',
        law='All inherited sources are byte-identical to parent; no physical hopping, rewiring, preparation or escape diagonal is deleted.')
    if args.refresh:
        (P/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    else:
        assert read('manifest.json')==manifest,'Content differs from saved manifest; inspect before intentionally refreshing.'
    result=dict(numerical_consistency_status='PASS',physical_acceptance_status='FAIL for the specified candidate metric; Milestone A not established',
        parent_commit=BASE,packet_files_hashed=len(files),manifest_sha256=sha(P/'manifest.json'),
        frozen_source_sha256=frozen,inherited_dependencies_unchanged=len(deps),
        response_interval_checks=response_count,independent_labelled_dimension=106914,
        exact_spanning_record_transition_checks=checks['raw_record_transition_checks'],
        six_site_neutral_dimension=neutral['dimension'],six_site_neutral_full_state_bound=neutral['rows'][0]['state_error'],
        six_site_krylov_convergence_difference=krylov['rows'][-1]['full_basis_difference'],
        six_site_distance3_to_distance2_ratio_at_quarter=quarter['ratio_absolute_intervals'],
        scope='Analytic full-H residual and bilinear response bounds evaluated in ordinary floating point with reported allowances, not an outward-rounded formal certificate. No no-go theorem for all phases, no derived continuum and no empirical prediction.',
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__))
    (P/'VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
