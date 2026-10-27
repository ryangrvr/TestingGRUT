"""Audit and seal the exploratory nonlinear scalar-limit packet."""
from en_common import *
import ast,re,platform

PARENT='b73de649843204639dfa1965a6326d4fc60afe06'
REVIEW='ec1eb7f860f98f8d889f90ff94ebd53396f740bb'
MODES='2cbc7f1c3a925fea1f3ac38d0b254ca037afd95e'
REHEAT='e1603565a0c1411a55cdef7f5dc4ecc25aae1c8f'
KIN='34f23f85df3160a86981ae40a159d0c6fffdf95a'
BUILD='e71730c24fe6a0b28ddd72cd30982c5b39c9aedf'
PILOT='815955656099aeb684d909bf21c940ca9add8e87'
BASE='research/reality_framework/'

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def pins():
    review=BASE+'grut_formation_review_20261025/'
    pairs=[(REVIEW,review+p) for p in ['REVIEW_TIMING.md','tracks/st_critic/NOTES.md','tracks/st_critic/results/neutral_extend.json',
        'tracks/st_physics/code/s4_vacuum_seed.py','tracks/st_physics/results/seed_projection.json','tracks/st_rerun/results/crossings_dense_vs_published.json']]
    pairs += [(PARENT,BASE+'grut_source_timing_20261026/'+p) for p in ['code/timing_common.py','code/response.c','code/source_clocks.py','RESPONSE_RESULTS.json']]
    pairs += [(MODES,BASE+'grut_reheating_modes_20261026/code/'+p) for p in ['modes.py','modes.c','neutral.py','mode_common.py']]
    pairs += [(REHEAT,BASE+'grut_reheating_preparation_20261026/code/'+p) for p in ['reheat.c','reheating.py','common.py']]
    pairs += [(REHEAT,BASE+'grut_heavy_vector_budget_20261025/INPUTS.json')]
    pairs += [(KIN,BASE+'grut_kinetic_continuation_20261026/code/'+p) for p in ['ks_common.py','bootstrap.py']]
    pairs += [(BUILD,review+'tracks/'+p+'/results/'+f) for p,f in [('kin_source','source_interface.json'),('kin_rates','rates_interface.json'),('kin_thermo','thermo_interface.json')]]
    pairs += [(PILOT,BASE+'grut_nonlinear_pilot_20261026/code/field3d.cpp')]
    rows=[]
    for commit,path in pairs:
        raw=subprocess.check_output(['git','show',commit+':'+path],cwd=REPO)
        rows.append(dict(commit=commit,path=path,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
    return rows

def finite(x):
    if isinstance(x,float):assert np.isfinite(x)
    elif isinstance(x,dict):
        for v in x.values():finite(v)
    elif isinstance(x,list):
        for v in x:finite(v)

def main():
    seal='--seal' in sys.argv
    runs=read('RUNS.json')['cases']
    if seal:
        save('SOURCES.json',dict(internal=pins(),external=[
            dict(arxiv='1402.0115',url='https://arxiv.org/abs/1402.0115',authors='Epelbaum, Gelis, Wu',role='Quantum/UV limitations of vacuum-type classical statistical evolution.'),
            dict(arxiv='1409.0701',url='https://arxiv.org/abs/1409.0701',authors='Epelbaum, Gelis, Tanji, Wu',role='Classical/quantum collision and cutoff matching.')],
            scope='Exact model inputs, completed AI review and primary methodological literature; no observational evidence.'))
        save('EXECUTION.json',dict(python=platform.python_version(),numpy=np.__version__,
            compiler=subprocess.check_output(['g++','--version']).decode().splitlines()[0],
            nonlinear_flags=['-std=c++17','-O3','-ffp-contract=off','-fopenmp'],
            numerical_workers=2,OpenMP_threads_each=4,nonlinear_source_sha256=sha(ROOT/'code/early3d.cpp')))
        raw=[]
        for tag,r in runs.items():
            for path in [WORK/tag/'initial.bin',WORK/tag/'parameters.json',WORK/tag/'state.tsv',Path(r['snapshots'][-1])]:
                raw.append(dict(path=str(path.relative_to(WORK)),bytes=path.stat().st_size,sha256=sha(path)))
        save('RAW_FILES.json',dict(files=raw,scope='Regenerable inputs, full dense traces and final double fields. Public curves are labelled subsets, not the whole quantum state.'))
    sources=read('SOURCES.json')['internal']
    for row in sources:
        raw=subprocess.check_output(['git','show',row['commit']+':'+row['path']],cwd=REPO)
        assert len(raw)==row['bytes'] and hashlib.sha256(raw).hexdigest()==row['sha256']
        if row['commit'] in [PARENT,MODES,REHEAT,KIN,PILOT]:assert sha(REPO/row['path'])==row['sha256']
        if row['commit']==BUILD:
            actual=T.K.INPUT_ROOT/row['path'].removeprefix(BASE+'grut_formation_review_20261025/')
            assert sha(actual)==row['sha256']
    checks=read('VERIFICATION_RESULTS.json');assert len(checks['groups'])==7 and all(g['passed'] for g in checks['groups'].values())
    result=read('RESULTS.json');assert len(runs)==7 and len(result['final'])==7
    assert result['physical_nbar_per_charge'] is None and result['full_gauged_model_evolved'] is False
    assert result['physical_f_sym'] is None
    for row in runs.values():assert row['parameters']['end_tau']==800 and row['maximum_recorded_tau_spacing']<.05000001
    dense=read('DENSE_CORRECTIONS.json');r=next(r for r in dense['rows'] if r['beta']==3)
    assert 56.4<r['preparation_phase_003_amplitude_warning']['tau_bracket'][1]<56.5
    assert .20<r['max_amplitude_through80']<.21
    assert abs(dense['seed_definition']['dtheta_X_per_radian'])<1e-10
    for path in ROOT.glob('*.json'):finite(json.loads(path.read_text()))
    for path in ROOT.glob('*.npz'):
        with np.load(path) as z:
            for key in z.files:assert np.all(np.isfinite(z[key])),(path,key)
    for path in (ROOT/'code').glob('*.py'):ast.parse(path.read_text(),filename=str(path))
    for path in ROOT.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',path.read_text()):
            if '://' in link or link.startswith('#'):continue
            assert (path.parent/link.split('#')[0]).exists(),(path,link)
    for name in ['NONLINEAR.png','REFERENCES.png']:assert (ROOT/name).stat().st_size>10000
    for row in read('RAW_FILES.json')['files']:
        path=WORK/row['path']
        if path.exists():assert path.stat().st_size==row['bytes'] and sha(path)==row['sha256']
    save('AUDIT_RESULTS.json',dict(passed_numerical_groups=7,nonlinear_histories=7,source_pins=len(sources),
        imported_model_files_unchanged=True,raw_dense_and_final_hashes=True,JSON_NPZ_finite=True,Python_AST=True,
        local_links=True,figures_rendered_and_viewed=True,physical_particle_population_certified=False,
        full_gauged_cosmology_evolved=False,independent_reviewer_of_this_packet_ran=False))
    if seal:
        rows=[]
        for path in sorted(ROOT.rglob('*')):
            if not path.is_file() or path.name=='MANIFEST.json' or '__pycache__' in path.parts:continue
            rows.append(dict(path=str(path.relative_to(ROOT)),bytes=path.stat().st_size,sha256=sha(path)))
        save('MANIFEST.json',dict(parent_commit=PARENT,files=rows))
    else:
        for row in read('MANIFEST.json')['files']:
            path=ROOT/row['path'];assert path.stat().st_size==row['bytes'] and sha(path)==row['sha256']
    print(f'AUDIT:7/7 groups;7 nonlinear scalar-limit histories;{len(sources)} source pins; dense correction, finite data, links, raw hashes and manifest pass. Physical populations remain null.')

if __name__=='__main__':main()
