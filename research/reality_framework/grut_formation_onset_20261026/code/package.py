"""Pin actual inputs and audit the limited onset packet."""
import ast,hashlib,re,subprocess,sys
import numpy as np
from common import *

PARENT='34f23f85df3160a86981ae40a159d0c6fffdf95a'
REVIEW='99aedde454c77f2c7dbabaab320a81e95a932fc7'
BUILD='e71730c24fe6a0b28ddd72cd30982c5b39c9aedf'
BASE='research/reality_framework/grut_formation_review_20261025/'

def source_rows():
    pairs=[(REVIEW,BASE+'REVIEW_KINETIC.md'),(REVIEW,BASE+'tracks/kin_review/results/repro.json')]
    pairs += [(BUILD,BASE+x) for x in ['tracks/kin_thermo/code/formation.py','tracks/kin_thermo/code/thermo.py',
        'tracks/kin_source/results/source_interface.json','tracks/kin_rates/results/rates_interface.json',
        'tracks/kin_thermo/results/thermo_interface.json']]
    pairs += [(PARENT,'research/reality_framework/grut_kinetic_continuation_20261026/'+x) for x in ['REPORT.md','RESULTS.json','code/ks_common.py','code/bootstrap.py']]
    pairs += [('2cbc7f1c3a925fea1f3ac38d0b254ca037afd95e','research/reality_framework/grut_reheating_modes_20261026/code/phase_and_orbit.py')]
    rows=[]
    for commit,path in pairs:
        raw=subprocess.check_output(['git','show',commit+':'+path],cwd=REPO)
        rows.append(dict(commit=commit,path=path,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
    return rows

def finite(v):
    if isinstance(v,float):assert np.isfinite(v)
    elif isinstance(v,dict):
        for z in v.values():finite(z)
    elif isinstance(v,list):
        for z in v:finite(z)

def main():
    if '--seal' in sys.argv:
        save('SOURCES.json',dict(internal=source_rows(),review_url='https://github.com/ryangrvr/TestingGRUT/issues/13#issuecomment-6429427680',
            attribution='Claude completed the original builds and focused kinetic review. Codex wrote this finite-response and matched cold-onset calculation. Relative-Floquet construction follows the earlier credited2cbc7f1 implementation. Same-author numerical checks, no external evidence.'))
        c=read('CIRCULAR_RESULTS.json');e=read('ELLIPTIC_RESULTS.json');r=read('ROLL_RESULTS.json');iv=read('INVENTORY_RESULTS.json')
        save('RESULTS.json',dict(accepted_review=REVIEW,accepted_gas_controls_unchanged=True,
            fixed_initial_fragment_charge_claim_withdrawn=True,physical_fragment_distribution=None,
            physical_initial_spectrum=None,physical_f_sym=None,physical_kinetic_threshold=None,
            circular_endpoint=c['examples'][-1],endpoint_circle_charge_over_release=c['endpoint_charge_over_release_charge'],
            endpoint_circle_energy_over_release=c['endpoint_energy_over_release_energy'],
            charge_matched_circle_energy_fraction=c['charge_matched_circular']['fraction_of_release_energy_in_this_circular_state'],
            cold_elliptic_peaks=e['peaks'],cold_seed_roll=[{k:v for k,v in rr.items() if k not in ['times','phi_over_v','X_amplitude_for_supplied_seed']} for rr in r['rows']],
            metastable_energy_cap_eta=iv['eta'],scope='Vacuum onset controls and a necessary thermal inequality, not hot-plasma formation or a final abundance.'))
    pins=read('SOURCES.json')['internal']
    for row in pins:
        raw=subprocess.check_output(['git','show',row['commit']+':'+row['path']],cwd=REPO)
        assert len(raw)==row['bytes'] and hashlib.sha256(raw).hexdigest()==row['sha256']
        if row['commit']==BUILD:
            actual=K.INPUT_ROOT/row['path'].removeprefix(BASE)
            assert hashlib.sha256(actual.read_bytes()).hexdigest()==row['sha256'],actual
        if row['commit']==PARENT and row['path'].endswith(('code/ks_common.py','code/bootstrap.py')):
            assert hashlib.sha256((REPO/row['path']).read_bytes()).hexdigest()==row['sha256']
    v=read('VERIFICATION_RESULTS.json')
    assert len(v['groups'])==7 and all(g['passed'] for g in v['groups'].values())
    assert v['physical_f_sym'] is None
    assert len(read('ELLIPTIC_RESULTS.json')['rows'])==121 and len(read('ROLL_RESULTS.json')['rows'])==2
    for p in ROOT.glob('*.json'):finite(json.loads(p.read_text()))
    for p in (ROOT/'code').glob('*.py'):ast.parse(p.read_text(),filename=str(p))
    for p in ROOT.glob('*.npz'):
        with np.load(p) as z:
            for key in z.files:assert np.all(np.isfinite(z[key]))
    for p in ROOT.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if '://' in link or link.startswith('#'):continue
            assert (p.parent/link.split('#')[0]).exists(),(p,link)
    assert (ROOT/'ONSET.png').stat().st_size>10000
    save('AUDIT_RESULTS.json',dict(passed_numerical_groups=7,source_pins=len(pins),actual_imported_inputs_verified=True,
        JSON_NPZ_finite=True,Python_AST=True,links=True,figure_rendered_and_viewed=True,
        physical_initial_spectrum_derived=False,nonlinear_formation_derived=False,physical_abundance_certified=False,
        independent_reviewer_ran=False))
    if '--seal' in sys.argv:
        rows=[]
        for p in sorted(ROOT.rglob('*')):
            if not p.is_file() or p.name=='MANIFEST.json' or '__pycache__' in p.parts:continue
            rows.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
        save('MANIFEST.json',dict(parent_commit=PARENT,files=rows))
    else:
        for row in read('MANIFEST.json')['files']:
            p=ROOT/row['path'];assert p.stat().st_size==row['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256']
    print('AUDIT:7/7 numerical groups;12 source pins; actual inputs, finite data, syntax, links and manifest pass. Physical outputs remain null.')

if __name__=='__main__':main()
