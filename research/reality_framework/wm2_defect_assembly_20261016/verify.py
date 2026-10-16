"""Saved-evidence and provenance check, distinct from WM2 physical acceptance."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET
import model

P=Path(__file__).resolve().parent;ROOT=P.parents[2]
BASE='8741b6e74d571e643e9208187f5d93b5a8e2471c'
PREFIX=P.relative_to(ROOT).as_posix()+'/'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(name):return json.loads((P/name).read_text())
def finite(x):
    if isinstance(x,dict):
        for v in x.values():finite(v)
    elif isinstance(x,list):
        for v in x:finite(v)
    elif isinstance(x,float):assert math.isfinite(x)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--refresh',action='store_true');args=ap.parse_args()
    deps={}
    for name in ['cp12_ground_connectivity/frozen_model.py','cp12_ground_connectivity/exact_destination.py',
        'path_b_protected_transport_20261016/MODEL.md','path_b_protected_transport_20261016/MATHEMATICS.md',
        'path_b_protected_transport_20261016/HANDOFF.md','path_b_protected_transport_20261016/manifest.json']:
        path=P.parent/name;rel=path.relative_to(ROOT).as_posix()
        assert path.read_bytes()==subprocess.check_output(['git','show',f'{BASE}:{rel}'],cwd=ROOT)
        deps[rel]=sha(path)
    old=P.parent/'path_b_protected_transport_20261016';om=json.loads((old/'manifest.json').read_text())
    for name,row in om['files'].items():assert sha(old/name)==row['sha256']
    for name,digest in om['unchanged_dependency_sha256'].items():assert sha(ROOT/name)==digest
    changed=subprocess.check_output(['git','diff','--name-only',BASE],cwd=ROOT,text=True).splitlines()
    untracked=subprocess.check_output(['git','ls-files','--others','--exclude-standard'],cwd=ROOT,text=True).splitlines()
    assert all(x.startswith(PREFIX) for x in changed+untracked)
    for p in P.glob('*.py'):compile(p.read_text(),str(p),'exec')
    for p in P.glob('*.json'):finite(json.loads(p.read_text()))
    e=read('EXACT_CONFIGURATION.json');r=read('RESULTS.json');a=read('ASSESSMENT.json')
    assert r['status']=='PASS' and r['constants']==model.constants(2)
    assert [x['N'] for x in e['word_controls']]==list(range(2,8))
    pairs=0;registers=0;word_cells=0;barriers=0;blocks=0;propagation=0;near=0;pulses=0
    for row in e['register_controls']:
        assert row['finite_register_and_tuple_neighbors_agree'] and row['label_covariance_checked']
        assert row['reverse_primitive_checks'] and row['unseeded_empty_configuration_dark']
        registers+=row['states']
    for row in e['word_controls']:
        assert row['complete_Q_graph_BFS'] and row['all_final_geometry_orders_retained']
        assert row['minimum_virtual_length_minus_d_op']>=2
        assert row['max_enumerated_word_length']==10
        pairs+=row['all_final_protected_pairs_checked']
        for source in row['by_source']:
            for x in source['word_counts']:
                for d,count in x['words_by_d_op'].items():
                    assert int(d)<=x['length']-2 and count>=1
                    word_cells+=1
    for row in e['activation_barriers']:
        for d,k in row['minimum_maximum_defect_count_by_initial_distance'].items():assert k==int(d)
        barriers+=row['long_bond_checks']
    for row in r['Schur']:
        n=row['N'];assert row['constants']==model.constants(n)
        assert row['full_matrix_hermiticity_error']<1e-12
        for x in row['rows']:
            assert x['max_bound_excess']<2e-13 and x['linear_solve_residual']<1e-11
            assert x['near_min']>=.25-1e-13
        if row['spectral']:
            assert row['spectral']['actual_gap']>=13-1e-11
            assert row['spectral']['projector_distance']<=model.constants(n)['spectral_projector_distance_upper']+1e-10
        for x in row['all_geometry_receiver_blocks']:
            assert x['operator_norm']<=x['bound']+2e-13;blocks+=1
        for x in row['arbitrary_near_geometry_packets']:
            assert x['minimum_receiver_action']>=x['uniform_lower']-1e-12;near+=1
        for x in row['full_propagation_blocks']:
            assert x['full_receiver_operator_norm']<=x['Dyson_tail_bound']+1e-11;propagation+=1
    for row in r['source']:
        assert row['source_variance']==.25
        for x in row['rows']:
            assert x['finite_pulse_error']<1e-10 and not x['formation_preparation'];pulses+=1
            assert x['defect_probability']<=model.constants(row['N'])['bare_defect_probability_all_times_upper']+1e-10
    for row in r['formation']+r['formation_robustness']:
        c=model.formation_bound(row['N'],row['U'],row['epsilon'])
        for k,v in c.items():assert row[k]==v
        assert row['sampled_max_probability']<=c['full_assembly_probability_upper_all_times']+1e-12
        assert row['sampled_max_extensive_probability']<=c['extensive_assembly_probability_upper_all_times']+1e-12
        assert abs(row['complete_geometry_record_work'])<=2*row['epsilon']+1e-10
        assert row['complete_geometry_record_entropy_bits_fixed_root']>=-1e-10
    assert r['formation'][-1]['N']==128
    assert r['formation'][-1]['extensive_assembly_probability_upper_all_times']<3e-6
    for row in r['bad_reconnection']:
        assert row['new_word_length']==3 and row['reference_S1_distance']>=4
        assert row['term_qubit_support']==6
    assert a['neutral_seed_formation'].startswith('FAIL:')
    assert a['dimensional_class'].startswith('SUPPLIED')
    for key in ['full_WM2_success','all_WM2_architectures_ruled_out','frozen_CP12_modified',
        'WM1_modified','SCOUT1_reopened','empirical_prediction','earned_quotient_value','large_CP12_dynamics_run']:
        assert not a[key]
    assert a['all_geometry_receiver_bound']
    ET.parse(P/'WM2_CONSTRUCTION.svg')
    for name in ['README.md','REPORT.md','MODEL.md','MATHEMATICS.md','HANDOFF.md','CHARTER.md','SOURCES.md']:
        s=(P/name).read_text();assert len(s)>1000 and 'TO_BE_INSERTED' not in s
    files=[p for p in sorted(P.iterdir()) if p.is_file() and p.name not in ['manifest.json','VERIFICATION.json']]
    assert all(p.suffix in ['.py','.md','.json','.svg'] for p in files)
    manifest=dict(research_base_commit=BASE,packet=PREFIX.rstrip('/'),
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},
        unchanged_dependency_sha256=deps,WM1_packet_files_unchanged=len(om['files']),
        frozen_SCOUT1_commit='a2987fed5880da62c874857faa0a93099d23f64d',
        scope='New WM2-S1 exact construction. Generated event adjacency; supplied path class. Conditional protection proved; neutral formation fails. No full WM2 phase or empirical prediction.')
    if args.refresh:(P/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    else:assert read('manifest.json')==manifest,'Inspect differences before refreshing the manifest.'
    out=dict(numerical_consistency='PASS',physical_acceptance='Full WM2 NOT ESTABLISHED; closed S1 formation FAIL',
        protected_transport='Analytic all-N theorem in full reachable S1 component',
        dimensional_class='Supplied by single-end grammar',
        register_states_checked=registers,source_final_protected_pairs_checked=pairs,
        virtual_word_distance_count_cells=word_cells,activation_barrier_checks=barriers,
        Schur_receiver_operator_blocks=blocks,near_packet_lower_bound_checks=near,
        full_propagation_operator_blocks=propagation,source_pulse_checks=pulses,
        largest_configuration_graph=max(x['full_configuration_states'] for x in e['word_controls']),
        largest_full_Schur_dimension=max(x['dimension'] for x in r['Schur']),
        formation_controls=len(r['formation']),formation_robustness_controls=len(r['formation_robustness']),
        original_neutral_embedding='Explicit CP12 tensor seeded auxiliary module; no hidden reservoir and no original charged-matter locality claim',
        WM1_packet_files_unchanged=len(om['files']),unchanged_dependencies=len(deps),
        packet_files_hashed=len(files),manifest_sha256=sha(P/'manifest.json'),
        arithmetic=r['arithmetic'],environment=r['environment'],external_independent_validation=False)
    (P/'VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))


if __name__=='__main__':main()
