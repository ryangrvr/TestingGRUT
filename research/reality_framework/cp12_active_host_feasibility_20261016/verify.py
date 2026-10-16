"""Saved-evidence checks and provenance, distinct from physical acceptance."""
import argparse,hashlib,json,math,platform,subprocess,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import scipy,matplotlib
P=Path(__file__).resolve().parent;RF=P.parent;ROOT=RF.parent.parent
BASE='609db3f19576f1fe6d7713748a7564c91c3a49b7';PREFIX=P.relative_to(ROOT).as_posix()+'/'


def read(name):return json.loads((P/name).read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def inside(I,x):return I[0]-1e-10<=x<=I[1]+1e-10


def finite(x):
    if isinstance(x,dict):
        for v in x.values():finite(v)
    elif isinstance(x,list):
        for v in x:finite(v)
    elif isinstance(x,float):assert math.isfinite(x)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--refresh',action='store_true');a=ap.parse_args()
    deps=['cp12_ground_connectivity/frozen_model.py','cp12_ground_connectivity/exact_destination.py',
        'cp12_ground_connectivity/invariant_evolution.py','charged_density_dynamics_20261015/fast_symmetry.py',
        'charged_density_dynamics_20261015/symmetry.cpp','charged_density_dynamics_20261015/stream.py',
        'cp12_shortcut_dynamics_20261015/operators.py','cp12_shortcut_dynamics_20261015/integrals.py',
        'cp12_dynamical_suppression_20261015/boundary.py','cp12_dynamical_suppression_20261015/linear_propagator.py',
        'cp12_dynamical_suppression_20261015/resources.py']
    deps+=['cp12_source_resolved_propagation_20261015/'+s for s in ['shared.py','source_basis.py','response_run.py','neutral_krylov.py','neutral_prepare.py']]
    hashes={}
    for name in deps:
        path=RF/name;rel=path.relative_to(ROOT).as_posix()
        assert path.read_bytes()==subprocess.check_output(['git','show',f'{BASE}:{rel}'],cwd=ROOT)
        hashes[rel]=sha(path)
    changed=subprocess.check_output(['git','diff','--name-only',BASE],cwd=ROOT,text=True).splitlines()
    untracked=subprocess.check_output(['git','ls-files','--others','--exclude-standard'],cwd=ROOT,text=True).splitlines()
    assert all(p.startswith(PREFIX) for p in changed+untracked)
    for path in P.glob('*.py'):compile(path.read_text(),str(path),'exec')
    for path in P.glob('*.json'):finite(json.loads(path.read_text()))
    response_checks=0
    for n,q,b in [(4,5,256),(5,6,256),(6,7,192)]:
        cat=read(f'CATALOG_N{n}_Q{q}.json');dynamic=read(f'HOST_N{n}_Q{q}_B{b}.json')
        assert sum(r['labelled_count'] for r in cat['classes'])==math.comb(n*(n-1)//2,q)
        assert dynamic['status']=='PASS' and [r['time'] for r in dynamic['rows']]==[0.,.5,1.,2.,5.]
        for r in dynamic['rows']:
            assert inside(r['connectivity_interval'],r['geometry']['connectivity_probability'])
            assert abs(sum(s['probability'] for s in r['host_classes'])-1)<1e-10
            assert inside(r['total_variation_interval'],r['total_variation_from_uniform'])
            assert inside(r['complete_record_work_interval'],r['complete_record_work'])
            assert inside(r['complete_record_entropy_interval_bits'],r['complete_record_entropy_bits'])
        assert dynamic['rows'][0]['total_variation_from_uniform']<1e-12
        assert abs(dynamic['rows'][0]['complete_record_work']-q*(n*(n-1)//2-q)/(10*(n*(n-1)//2)))<1e-12
        for path in P.glob(f'HOST_RESPONSE_N{n}_Q{q}_B{b}_L*.json'):
            d=json.loads(path.read_text());assert d['status']=='PASS' and d['dimension']<=1200000
            assert d['source_variance_interval'][0]>0
            assert d['build_checks']['H_hermiticity']<1e-8
            initial=d['initial'];assert initial['discarded_initial_norm']<=1.1e-10*math.sqrt(initial['one_labelled_host_probability'])
            assert abs(initial['class_probability']-dynamic['rows'][-1]['source_class_probability'])<1e-12
            assert initial['lift_probability_error']<1e-12
            for r in d['rows']:
                assert inside(r['host_retention_interval'],r['host_retention'])
                assert r['host_retention_interval'][1]+1e-10>=r['universal_host_retention_lower_bound']
                assert r['retained_edge_fraction']+1e-10>=r['host_retention']
                assert r['source_neighborhood_retention']+1e-10>=r['host_retention']
                for z in list(r['receivers'].values())+list(r['shells'].values()):
                    rr=z['response'];assert inside(rr['interval'],rr['value'])
                    assert rr['exact_finite_phase_identity_error']<1e-11
                    assert rr['zero_excluded']==(abs(rr['value'])>rr['error_bound']);response_checks+=1
                for dist,shell in r['shells'].items():
                    mean=sum(r['receivers'][str(j)]['response']['value'] for j in shell['vertices'])/len(shell['vertices'])
                    assert abs(mean-shell['response']['value'])<1e-11
    checks=read('CHECKS.json');source=read('RESPONSE_CHECK.json');gate=read('PATH_B_ALGEBRA.json')
    assert all(d['status']=='PASS' for d in [checks,source,gate])
    assert checks['exact_cross_host_norm_square']=='1/720' and checks['raw_rewiring_checks']>1000
    assert checks['host_record_work_identity_difference']<1e-11
    assert checks['primary_fixed_edge_constraints'][1]['hosts']==455
    assert checks['primary_fixed_edge_constraints'][1]['diameters']==[2]
    assert source['host_only_initial_lift_difference']<2e-10
    assert source['record_work_identity_difference']<1e-11 and source['record_entropy_identity_difference']<1e-11
    assert all(r['central_difference']<1e-11 for r in source['fixed_receiver_checks'])
    assert source['complete_record_effects'][0]['response_change_interval'][0]>0
    assert len(gate['rows'])==12
    for r in gate['rows']:
        assert r['target_transfer_probability_T4_coefficient']=='J^2/4'
        assert r['activation_diagonal_detuning']=='-15/8'
        assert r['probe_vacuum_preserves_original_geometry_law']
        assert r['one_geometry_hop_activates_gate']
    base=read('LARGE_N_BASELINE.json');assert base['seed']==20261016 and len(base['rows'])==14
    assert all(inside(r['connectivity_interval'],r['connected']) for r in base['rows'])
    assert base['exact_isolate_moments'][-1]['connected_probability_upper_bound']<.007
    assessment=read('ASSESSMENT.json')
    assert assessment['path_A_spatial_phase']=='NOT ESTABLISHED'
    assert not assessment['distinctive_empirical_prediction'] and not assessment['SCOUT1_reopened']
    assert not assessment['no_go_for_every_CP12_phase']
    assert next(r for r in assessment['conditional_source_tests'] if r['N']==5)['all_adjacent_distance_orders_resolved']
    six_files=[json.loads(p.read_text()) for p in P.glob('HOST_RESPONSE_N6_Q7_B192_L*.json')]
    latest=max(six_files,key=lambda r:r['level']);closure=latest.get('initial_boundary_closure')
    assert closure and closure['new_dimension']<=1200000
    assert closure['initial_residual_norms_after']==[0.,0.] and not closure['physical_law_changed']
    if len(six_files)>1:
        coarse=min(six_files,key=lambda r:r['level'])
        for fine_row in latest['rows']:
            old=next(r for r in coarse['rows'] if r['lag']==fine_row['lag'])
            for j,r in fine_row['receivers'].items():
                I=r['response']['interval'];J=old['receivers'][j]['response']['interval']
                assert max(I[0],J[0])<=min(I[1],J[1])
    ET.parse(P/'ACTIVE_HOST.svg');assert (P/'ACTIVE_HOST.svg').stat().st_size>1000
    for name in ['REPORT.md','HANDOFF.md','README.md','MATHEMATICS.md','PATH_B_DESIGN.md','REVIEW_CONTINUITY.md']:
        assert len((P/name).read_text())>1000 and 'TO_BE_INSERTED' not in (P/name).read_text()
    files=[p for p in sorted(P.iterdir()) if p.is_file() and p.name not in ['manifest.json','VERIFICATION.json']]
    assert all(p.suffix in ['.py','.json','.md','.svg'] for p in files)
    manifest=dict(research_base_commit=BASE,checkpoint_commit='c2751a047c3ad57eb415b39228297b3f796b9287',
        packet=PREFIX.rstrip('/'),files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},
        unchanged_dependency_sha256=hashes,SCOUT1_reference_commit='a2987fed5880da62c874857faa0a93099d23f64d',
        scope='Path A frozen law unchanged. Path B is an explicitly new proposed law and exact supplied-state algebra. Caches and manifest/verification self-references excluded.')
    if a.refresh:(P/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    else:assert read('manifest.json')==manifest,'Evidence differs from saved manifest; inspect before refreshing.'
    result=dict(numerical_consistency='PASS',physical_acceptance='Spatial phase and distinctive empirical prediction NOT ESTABLISHED',
        research_base=BASE,packet_files_hashed=len(files),manifest_sha256=sha(P/'manifest.json'),
        unchanged_dependencies=len(deps),response_interval_checks=response_checks,
        independent_labelled_source_dimension=source['labelled_dimension'],
        exact_rewiring_norm_checks=checks['raw_rewiring_checks'],
        host_source_lift_difference=source['host_only_initial_lift_difference'],
        Path_B_exact_controls=len(gate['rows']),no_large_N_quantum_dynamics_claim=True,
        scope='Analytic full-law residual bounds evaluated in ordinary floating point with stated allowances; not outward-rounded formal certificates or independent external validation.',
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__))
    (P/'VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
