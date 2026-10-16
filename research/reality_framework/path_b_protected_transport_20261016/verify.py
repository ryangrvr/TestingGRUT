"""Evidence/provenance verification, kept distinct from physical acceptance."""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET
import module as wm

P=Path(__file__).resolve().parent; ROOT=P.parents[2]
BASE='678e934367ead8e03cb51e39dd078ca3dc265cae'
PREFIX=P.relative_to(ROOT).as_posix()+'/'


def read(name):return json.loads((P/name).read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def finite(x):
    if isinstance(x,dict):
        for v in x.values():finite(v)
    elif isinstance(x,list):
        for v in x:finite(v)
    elif isinstance(x,float):assert math.isfinite(x)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--refresh',action='store_true'); args=ap.parse_args()
    depnames=['cp12_ground_connectivity/frozen_model.py',
        'cp12_ground_connectivity/exact_destination.py',
        'cp12_active_host_feasibility_20261016/PATH_B_DESIGN.md',
        'cp12_active_host_feasibility_20261016/MATHEMATICS.md',
        'cp12_active_host_feasibility_20261016/HANDOFF.md',
        'cp12_active_host_feasibility_20261016/ASSESSMENT.json',
        'cp12_active_host_feasibility_20261016/manifest.json']
    deps={}
    for name in depnames:
        path=P.parent/name; rel=path.relative_to(ROOT).as_posix()
        assert path.read_bytes()==subprocess.check_output(['git','show',f'{BASE}:{rel}'],cwd=ROOT)
        deps[rel]=sha(path)
    previous=P.parent/'cp12_active_host_feasibility_20261016'
    previous_manifest=json.loads((previous/'manifest.json').read_text())
    for name,row in previous_manifest['files'].items():
        assert sha(previous/name)==row['sha256']
    for name,digest in previous_manifest['unchanged_dependency_sha256'].items():
        assert sha(ROOT/name)==digest
    changed=subprocess.check_output(['git','diff','--name-only',BASE],cwd=ROOT,text=True).splitlines()
    untracked=subprocess.check_output(['git','ls-files','--others','--exclude-standard'],cwd=ROOT,text=True).splitlines()
    assert all(x.startswith(PREFIX) for x in changed+untracked)
    for path in P.glob('*.py'):compile(path.read_text(),str(path),'exec')
    for path in P.glob('*.json'):finite(json.loads(path.read_text()))
    r=read('RESULTS.json'); a=read('ASSESSMENT.json')
    assert r['status']=='PASS' and r['constants']==wm.constants()
    assert len(r['modules'])==11 and len(r['latch'])==7 and len(r['robustness'])==12
    for row in r['modules']:
        assert row['max_degree']<=4
        assert row['near_dressed_coupling']>=r['constants']['dressed_near_lower']
        assert row['intertwining_Frobenius_error']<2e-10
        assert row['orthogonality_Frobenius_error']<2e-10
        assert row['max_Feshbach_bound_excess']<2e-13
        assert row['max_dressed_bound_excess']<2e-13
        if row['kind']=='cycle':assert row['cycle_Green_function_max_error']<1e-13
    for row in r['latch']+r['robustness']:
        c=wm.latch_constants(row['N'],gap=row['Delta'],delta=row['delta'],lam=row['lam'])
        assert row['constants']==c
        for x in row['Schur']:
            assert x['far']+1e-16>=c['shortcut_floor_all_low_energies']
            assert x['near']>=c['nearest_lower_all_low_energies']-1e-13
            assert x['linear_solve_residual']<1e-12
        if row['spectrum']:
            assert row['spectrum']['actual_gap']>=c['full_spectral_gap_lower']
            assert row['spectrum']['canonical_isometry_error']<1e-9
            for x in row['spectrum']['full_dynamics_checks']:
                assert x['geometry_one_probability']<=x['Q_probability']+1e-12
                assert x['Q_probability']<=c['bare_P0_total_Q_probability_upper']+1e-10
    last=r['latch'][-1]; assert last['N']==512 and last['distance']==256
    assert last['Schur'][0]['far']>1e-6
    assert last['constants']['shortcut_floor_all_low_energies']>9.43e-7
    assert all(x['variance']>.24 for x in r['source'])
    assert len(r['source'])==3
    pulses=0
    for row in r['source']:
        for x in row['rows']:
            assert x['finite_pulse_error']<2e-10 and x['norm_error']<2e-10
            assert x['bare_source_high_probability']<=r['constants']['bare_high_probability_upper']+1e-10
            pulses+=1
    words=r['CP12_activation']; assert words['status']=='PASS' and len(words['rows'])==12
    gauss_checks=0
    for row in words['rows']:
        assert row['all_30_word_sum']=='-4' and row['bare_Feshbach_CGTGC_word']=='-1'
        assert sum(Fraction(x) for x in row['nonzero_words'].values())==-4
        assert row['independent_full_H5']['value']=='-4'
        assert row['independent_full_H5']['lower_powers_zero_through']==4
        assert row['independent_full_H5']['vacuum_full_operator_identity_checked']
        gauss_checks+=row['independent_full_H5']['Gauss_checks']
    assert a['WM1_coupled_protection']=='NOT ESTABLISHED'
    assert a['GC1_metric_protection'].startswith('FAIL:')
    for key in ['empirical_prediction','earned_quotient_value','SCOUT1_reopened',
        'frozen_CP12_modified','all_CP12_phases_ruled_out','original_neutral_new_receiver_experiment']:
        assert not a[key]
    ET.parse(P/'PROTECTED_TRANSPORT.svg')
    assert (P/'PROTECTED_TRANSPORT.svg').stat().st_size>1000
    for name in ['README.md','REPORT.md','MATHEMATICS.md','MODEL.md','HANDOFF.md','PROTOCOL.md','SOURCES.md']:
        text=(P/name).read_text(); assert len(text)>1000 and 'TO_BE_INSERTED' not in text
    files=[p for p in sorted(P.iterdir()) if p.is_file() and p.name not in ['manifest.json','VERIFICATION.json']]
    assert all(p.suffix in ['.py','.md','.json','.svg'] for p in files)
    manifest=dict(research_base_commit=BASE,packet=PREFIX.rstrip('/'),
        files={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in files},
        unchanged_dependency_sha256=deps,previous_packet_files_unchanged=len(previous_manifest['files']),
        frozen_SCOUT1_commit='a2987fed5880da62c874857faa0a93099d23f64d',
        scope='New WM1 candidate and separate GC1 control. Static theorem conditional. No full coupled spatial phase, universal matter locality or empirical prediction. Manifest and verification self-references excluded.')
    if args.refresh:(P/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    else:assert read('manifest.json')==manifest,'Inspect differences before refreshing the manifest.'
    output=dict(numerical_consistency='PASS',physical_acceptance='NOT ESTABLISHED',
        WM1_complete_protected_phase='NOT ESTABLISHED',GC1_reference_metric='FAIL',
        static_module='Analytic all-size conditional theorem',
        module_controls=len(r['modules']),geometry_controls=len(r['latch']),
        robustness_controls=len(r['robustness']),source_pulse_controls=pulses,
        exact_activation_controls=len(words['rows']),
        independent_full_H5_controls=len(words['rows']),exact_Gauss_checks=gauss_checks,
        original_neutral_vacuum_identity='Exact operator identity; no new formation simulation claimed.',
        original_neutral_new_receiver_experiment=False,
        unchanged_dependencies=len(deps),previous_packet_files_unchanged=len(previous_manifest['files']),
        inherited_frozen_dependencies_unchanged=len(previous_manifest['unchanged_dependency_sha256']),
        packet_files_hashed=len(files),manifest_sha256=sha(P/'manifest.json'),
        arithmetic=r['arithmetic'],environment=r['environment'],external_independent_validation=False)
    (P/'VERIFICATION.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
