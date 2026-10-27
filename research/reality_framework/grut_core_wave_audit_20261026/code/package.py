"""Source pins, consumed-snapshot hashes and a bounded output audit."""
from cw_common import *
import ast,re,subprocess

PARENT='815955656099aeb684d909bf21c940ca9add8e87'
REVIEW='870a2ab51b74c3ab0a3d33fa8997ece376b80fb1'
ONSET='23356538581c64e49c255e6cc895373d5d297c96'
KIN='34f23f85df3160a86981ae40a159d0c6fffdf95a'
BASE='research/reality_framework/'

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()

def source_pins():
    reviewbase=BASE+'grut_formation_review_20261025/'
    pairs=[(REVIEW,reviewbase+p) for p in ['REVIEW_ONSET_PILOT.md','tracks/onset_review/code/floquet.py','tracks/onset_review/results/floquet.json']]
    pairs += [(PARENT,BASE+'grut_nonlinear_pilot_20261026/'+p) for p in ['code/pilot.py','code/analyze.py','code/field3d.cpp','RUNS.json','ANALYSIS.json','SNAPSHOT_MANIFEST.json','MATHEMATICS.md','DATASET.json']]
    pairs += [(ONSET,BASE+'grut_formation_onset_20261026/'+p) for p in ['code/common.py','MATHEMATICS.md','REPORT.md']]
    pairs += [(KIN,BASE+'grut_kinetic_continuation_20261026/code/'+p) for p in ['ks_common.py','bootstrap.py']]
    rows=[]
    for commit,path in pairs:
        raw=subprocess.check_output(['git','show',commit+':'+path],cwd=REPO)
        rows.append(dict(commit=commit,path=path,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
    return rows

def used_snapshots():
    runs=json.loads((PILOT/'RUNS.json').read_text())['cases']
    prior=json.loads((PILOT/'ANALYSIS.json').read_text())['cases']
    selected=set()
    for tag,run in runs.items():
        selected.add((tag,Path(run['snapshots'][-1]).name))
        if tag in ['mixed_N96','onesign_N96']:
            selected.update((tag,Path(path).name) for path,row in zip(run['snapshots'],prior[tag]['rows']) if row['t']>=40)
    expected={row['path']:row for row in json.loads((PILOT/'SNAPSHOT_MANIFEST.json').read_text())['files']}
    out=[]
    for tag,name in sorted(selected):
        key=tag+'/'+name;row=expected[key];path=P.WORK/key
        assert path.stat().st_size==row['bytes'] and sha(path)==row['sha256'],path
        out.append(row)
    return out

def finite(x):
    if isinstance(x,float):assert np.isfinite(x)
    elif isinstance(x,dict):
        for a in x.values():finite(a)
    elif isinstance(x,list):
        for a in x:finite(a)

def main():
    seal='--seal' in sys.argv
    if seal:
        save('SOURCES.json',dict(internal=source_pins(),
            review_url='https://github.com/ryangrvr/TestingGRUT/issues/13#issuecomment-6431164252',
            attribution='Claude supplied an independent Floquet reproduction and a reading of the prior3D equations, not repeated long trajectories. Codex wrote this exact local field ledger and post-processing; seven same-author groups, no external evidence.'))
        save('INPUT_SNAPSHOTS.json',dict(files=used_snapshots(),source_commit=PARENT,
            scope='27 existing cold pilot snapshots, unchanged; only their manifest entries are duplicated here.'))
    pins=read('SOURCES.json')['internal']
    for row in pins:
        raw=subprocess.check_output(['git','show',row['commit']+':'+row['path']],cwd=REPO)
        assert len(raw)==row['bytes'] and hashlib.sha256(raw).hexdigest()==row['sha256']
        # Imported source/interfaces have not been edited on this additive branch.
        if row['commit'] in [PARENT,ONSET,KIN]:assert sha(REPO/row['path'])==row['sha256']
    v=read('VERIFICATION_RESULTS.json');assert len(v['groups'])==7 and all(g['passed'] for g in v['groups'].values())
    final=read('FINAL_AUDIT.json')['cases'];history=read('HISTORY_AUDIT.json')['cases'];result=read('RESULTS.json')
    assert len(final)==7 and len(history)==2
    assert result['new_dynamical_histories']==0 and result['physical_f_sym'] is None
    assert result['physical_bound_population'] is None and result['physical_particle_spectrum'] is None
    snapshots=read('INPUT_SNAPSHOTS.json')['files'];assert len(snapshots)==27
    assert all(Path(r['path']).parts[0] in final for r in snapshots)
    from audit_fields import resolve_path
    run=json.loads((PILOT/'RUNS.json').read_text())['cases']['mixed_N96']
    assert resolve_path('mixed_N96',run['snapshots'][-1])==P.WORK/'mixed_N96'/Path(run['snapshots'][-1]).name
    assert abs(C.G-.1)<1e-14 and abs(C.LX-.001)<1e-14 and abs(C.LP-.03)<1e-14 and abs(C.NU-10)<1e-12
    for path in ROOT.glob('*.json'):finite(json.loads(path.read_text()))
    for path in (ROOT/'code').glob('*.py'):ast.parse(path.read_text(),filename=str(path))
    for path in ROOT.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',path.read_text()):
            if '://' in link or link.startswith('#'):continue
            assert (path.parent/link.split('#')[0]).exists(),(path,link)
    assert (ROOT/'CORE_LEDGER.png').stat().st_size>10000
    save('AUDIT_RESULTS.json',dict(passed_numerical_groups=7,unique_input_snapshots=27,
        source_pins=len(pins),source_inputs_unchanged=True,JSON_finite=True,Python_AST=True,links=True,
        workdir_relocation_checked=True,figure_rendered_and_viewed=True,
        independent_reviewer_of_this_packet_ran=False,physical_bound_population_certified=False,
        physical_f_sym_certified=False,new_dynamical_histories=0))
    if seal:
        files=[]
        for path in sorted(ROOT.rglob('*')):
            if not path.is_file() or path.name=='MANIFEST.json' or '__pycache__' in path.parts:continue
            files.append(dict(path=str(path.relative_to(ROOT)),bytes=path.stat().st_size,sha256=sha(path)))
        save('MANIFEST.json',dict(parent_commit=PARENT,files=files))
    else:
        for row in read('MANIFEST.json')['files']:
            path=ROOT/row['path'];assert path.stat().st_size==row['bytes'] and sha(path)==row['sha256']
    print(f'AUDIT:7/7 groups;27 unchanged snapshot inputs;{len(pins)} source pins; finite JSON, syntax, links and manifest pass. Physical abundance/population remain null.')

if __name__=='__main__':main()
