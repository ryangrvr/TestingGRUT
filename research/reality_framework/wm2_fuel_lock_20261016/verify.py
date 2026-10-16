"""Saved-evidence, scope and provenance checks; physical acceptance stays separate."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET
import model

P=Path(__file__).resolve().parent;ROOT=P.parents[2]
BASE='7b1d9234deaeaddea3f34c419bdb35b99efd944b'
PREFIX=P.relative_to(ROOT).as_posix()+'/'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(name):return json.loads((P/name).read_text())
def finite(x):
    if isinstance(x,dict):
        for v in x.values():finite(v)
    elif isinstance(x,list):
        for v in x:finite(v)
    elif isinstance(x,float):assert math.isfinite(x)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--refresh',action='store_true');args=parser.parse_args()
    changed=subprocess.check_output(['git','diff','--name-only',BASE],cwd=ROOT,text=True).splitlines()
    untracked=subprocess.check_output(['git','ls-files','--others','--exclude-standard'],cwd=ROOT,text=True).splitlines()
    assert all(x.startswith(PREFIX) for x in changed+untracked),'Changes outside the separate S2 packet'
    deps={}
    old=P.parent/'wm2_defect_assembly_20261016';om=json.loads((old/'manifest.json').read_text())
    for name,row in om['files'].items():assert sha(old/name)==row['sha256']
    for name,digest in om['unchanged_dependency_sha256'].items():assert sha(ROOT/name)==digest
    for p in sorted(old.iterdir()):
        if not p.is_file():continue
        rel=p.relative_to(ROOT).as_posix()
        assert p.read_bytes()==subprocess.check_output(['git','show',f'{BASE}:{rel}'],cwd=ROOT)
        deps[rel]=sha(p)
    deps.update(om['unchanged_dependency_sha256'])
    wm1=P.parent/'path_b_protected_transport_20261016';wm=json.loads((wm1/'manifest.json').read_text())
    for name,row in wm['files'].items():assert sha(wm1/name)==row['sha256']
    for name,digest in wm['unchanged_dependency_sha256'].items():assert sha(ROOT/name)==digest
    for p in P.glob('*.py'):compile(p.read_text(),str(p),'exec')
    for p in P.glob('*.json'):finite(json.loads(p.read_text()))
    e=read('CONFIGURATION.json');r=read('RESULTS.json');a=read('ASSESSMENT.json')
    assert e['status']==r['status']=='PASS' and r['constants']==model.constants()
    assert [x['N'] for x in e['rows']]==[2,3,4,5]
    for row in e['rows']:
        for flag in ['all_root_spent_masks_retained','independent_register_support_checked',
            'root_moving_label_covariance_checked','energy_changes_checked','reverse_moves_checked',
            'all_P_operational_metric_matches']:assert row[flag]
        assert row['hermiticity_error']==0 and row['minimum_virtual_length_minus_d']==2
    blocks=near=prop=pulses=0;residual=0.
    for row in r['Schur']:
        assert row['dimension']==row['P_dimension']+row['Q_dimension']
        assert row['PQ_norm']<=.25+1e-12 and row['PQ_squared_offdiagonal_max']<1e-14
        if row['spectrum']:
            s=row['spectrum'];assert s['actual_gap']>=42.5-1e-10
            assert s['projector_distance']<=model.constants()['spectral_projector_distance_upper']+1e-9
        for x in row['rows']:
            expected=(1/680)*(6/91)**x['min_distance']
            assert abs(x['bound']-expected)<1e-17
            assert x['Sigma_block_norm']<=x['bound']+2e-13 and x['solve_residual']<1e-11
            residual=max(residual,x['solve_residual']);blocks+=1
        for x in row['near_packet_checks']:
            assert x['minimum_receiver_action']>=model.constants()['near_packet_lower']-1e-12;near+=1
    assert [x['N'] for x in r['formation']]==[2,3,4,8,16,32,64,128]
    ledger=0.
    for x in r['formation']+r['detuned']:
        fc=model.formation_constants(x['N']);assert fc==x['formation_constants']
        assert fc['finite_window_macroscopic_lower']>=.375
        assert abs(x['complete_energy_ledger_error'])<1e-8
        ledger=max(ledger,abs(x['complete_energy_ledger_error']))
        assert x['complete_record_entropy_bits']>=0
        if x['delta']==0:
            assert x['initial_eigenweight_max_error_uniform']<1e-11
            assert x['actual_min_spacing']>=fc['spacing_lower']-1e-12
            assert x['finite_window_macroscopic_probability']>=fc['finite_window_macroscopic_lower']-1e-11
            assert abs(x['infinite_time_full_probability']-1/x['N'])<1e-11
            assert abs(x['infinite_time_macroscopic_probability']-fc['diagonal_ensemble_macroscopic'])<1e-11
            if x['exact_radial_intertwining_error'] is not None:assert x['exact_radial_intertwining_error']<1e-11
        else:
            db=model.detuning_bound(x['N'],x['delta']);assert db==x['detuning_all_time_bound']
            assert x['finite_window_macroscopic_probability']<=db['all_time_macroscopic_upper']+1e-10
    for x in r['detuning_large_N_bounds']:assert x==model.detuning_bound(x['N'],x['delta'])
    for x in r['full_propagation']:
        for row in x['rows']:
            assert row['full_H_receiver_operator_norm']<=row['Dyson_tail_bound']+1e-11;prop+=1
    for x in r['sources']:
        assert x['variance']==.25 and len(x['recorded_word'])==x['assembled_events']
        for row in x['rows']:
            assert row['finite_pulse_error']<1e-10
            assert row['mismatch_probability']<=4*model.constants()['spectral_projector_distance_upper']**2+1e-10
            assert row['bulk_prefix_failure_probability']<=row['full_H_bulk_prefix_bound']+1e-10
            for k,v in row['receiver_response_distance_bounds'].items():assert abs(row['chi_all_geometry'][int(k)])<=v+1e-11
            if row['nearby_full_H_response_lower']>0:
                assert row['chi_all_geometry'][1]>=row['nearby_full_H_response_lower']-1e-11
            pulses+=1
    b=r['branching'];assert b['status']=='PASS'
    for x in b['rows']:
        assert x['protected_primitive_word_length']==5 and x['initial_tree_distance']>5
        assert x['full_fifth_power_unit_coefficient']==-2 and x['lower_powers_zero_through']==4
        assert x['original_tree_and_storage_fully_restored'] and x['all_five_transitions_stay_fuel_consistent']
    for x in b['binary_tree_ball_obstruction']:assert x['leaf_ball_volume_lower']>x['N']/2
    assert a['macroscopic_formation'].startswith('CONDITIONAL PASS')
    assert a['robust_untuned_formation'].startswith('FAIL') and a['naive_branching'].startswith('FAIL')
    assert a['complete_virtual_protection'] and a['source_response_front'] and a['bulk_prefix_persistence']
    for k in ['full_WM2_success','all_WM2_architectures_ruled_out','frozen_CP12_modified','S1_modified',
        'WM1_modified','SCOUT1_reopened','empirical_prediction','earned_quotient_value',
        'common_original_CP12_matter_geometry','large_CP12_dynamics_run','external_independent_validation']:assert not a[k]
    for name in ['README.md','REPORT.md','MODEL.md','MATHEMATICS.md','HANDOFF.md','CHARTER.md','SOURCES.md']:
        s=(P/name).read_text();assert len(s)>1000 and 'TO_BE_INSERTED' not in s
    ET.parse(P/'S2_RESULTS.svg')
    files=[p for p in sorted(P.iterdir()) if p.is_file() and p.name not in ['manifest.json','VERIFICATION.json']]
    assert all(p.suffix in ['.py','.md','.json','.svg'] for p in files)
    manifest=dict(research_base_commit=BASE,packet=PREFIX.rstrip('/'),
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},
        unchanged_dependency_sha256=deps,S1_packet_files_unchanged=len(om['files'])+2,
        WM1_manifest_entries_unchanged=len(wm['files']),
        frozen_SCOUT1_commit='a2987fed5880da62c874857faa0a93099d23f64d',
        scope='New S2 law: resonance-conditional macroscopic formation and complete transport bound. Robust formation and naive branching gates fail. No full spatial phase or empirical prediction.')
    if args.refresh:(P/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    else:assert read('manifest.json')==manifest,'Inspect changes before refreshing content hashes.'
    out=dict(numerical_consistency='PASS',
        physical_acceptance='Conditional resonant module only; robust untuned phase FAIL; dimension NOT EARNED',
        complete_configuration_states_checked=sum(x['complete_reachable_states'] for x in e['rows']),
        source_exit_tests=sum(x['source_exit_pairs_checked'] for x in e['rows']),
        largest_full_one_probe_dimension=max(x['dimension'] for x in r['Schur']),
        Schur_operator_blocks=blocks,near_packet_checks=near,full_H_propagation_blocks=prop,
        finite_pulse_checks=pulses,formation_controls=len(r['formation']),
        detuned_radial_controls=len(r['detuned']),analytic_detuning_controls=len(r['detuning_large_N_bounds']),
        branching_exact_fifth_power_controls=len(b['rows']),
        maximum_Schur_solve_residual=residual,maximum_complete_energy_ledger_error=ledger,
        S1_packet_files_unchanged=len(om['files'])+2,WM1_manifest_entries_unchanged=len(wm['files']),
        unchanged_dependencies=len(deps),packet_files_hashed=len(files),manifest_sha256=sha(P/'manifest.json'),
        environment=r['environment'],arithmetic=r['arithmetic'],external_independent_validation=False)
    (P/'VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))


if __name__=='__main__':main()
