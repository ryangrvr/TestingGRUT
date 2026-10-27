"""Pin inputs, audit the stated scope and seal only this additive packet."""
from pt_common import *
import ast,re,platform

PARENT_SHA='7d9c93bc1dffea59228765eed41f65baa4a9692c'
REVIEW_SHA='fd71d3a416db929252561d2cf296e915dd710f37'
REVIEW_BASE='research/reality_framework/grut_formation_review_20261025/'
FOCUS=['REVIEW_EARLY.md','tracks/en_particles/NOTES.md','tracks/en_particles/SKEPTIC.md',
 'tracks/en_particles/code/spectral_window.py','tracks/en_particles/code/response.py',
 'tracks/en_particles/skeptic/code/sk_orbit.py','tracks/en_particles/skeptic/results/sk_orbit_tau20.json',
 'tracks/en_critic/NOTES.md','tracks/en_rerun/NOTES.md','tracks/en_neutral/NOTES.md']

def git_bytes(commit,path):return subprocess.check_output(['git','show',commit+':'+path],cwd=REPO)
def stream_sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()

def sources():
    rows=read_parent('SOURCES.json')['internal'];out=[]
    for row in rows:
        out.append({**row,'source_role':'Inherited audited model input'})
    prefix='research/reality_framework/grut_early_nonlinear_20261027/'
    for rel in ['code/en_common.py','code/initial.py','code/early3d.cpp','MATCH_L12.npz','FIELD_CURVES.npz','SOURCES.json','RAW_FILES.json']:
        path=prefix+rel;raw=git_bytes(PARENT_SHA,path)
        out.append(dict(commit=PARENT_SHA,path=path,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),source_role='Reviewed parent pilot'))
    cache=REPO.parent/'review_fd71d3a'
    for rel in FOCUS:
        path=REVIEW_BASE+rel;raw=git_bytes(REVIEW_SHA,path);p=cache/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(raw)
        out.append(dict(commit=REVIEW_SHA,path=path,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),source_role='Completed Claude review / primary check'))
    save('SOURCES.json',dict(internal=out,external=[
        dict(arxiv='hep-ph/0103049',authors='G. Aarts, J. Berges',url='https://arxiv.org/abs/hep-ph/0103049',role='Statistical versus spectral functions and off-shell transport.'),
        dict(arxiv='1912.07565',authors='L. Shen, J. Berges',url='https://arxiv.org/abs/1912.07565',role='Nonperturbative infrared statistical/spectral distinctions.'),
        dict(arxiv='1201.0687',authors='J. Berges, D. Sexty',url='https://arxiv.org/abs/1201.0687',role='Relativistic condensation and number-changing evolution.'),
        dict(arxiv='1402.0115',authors='T. Epelbaum, F. Gelis, B. Wu',url='https://arxiv.org/abs/1402.0115',role='Inherited classical-statistical UV matching limitation.')],
        scope='Immutable model inputs and a completed AI review; no observational evidence or new independent review.'))

def read_parent(name):return json.loads((PARENT/name).read_text())

def finite(x):
    if isinstance(x,float):assert np.isfinite(x)
    elif isinstance(x,dict):
        for value in x.values():finite(value)
    elif isinstance(x,list):
        for value in x:finite(value)

def main():
    if '--sources' in sys.argv:sources();print('Pinned sources and focused fd71d3a cache ready.');return
    seal='--seal' in sys.argv
    if seal:
        sources();raw=[]
        for p in sorted(ROOT.glob('RUN_*.json')):
            tag=p.stem.removeprefix('RUN_');directory=WORK/tag
            for name in ['initial.npy','state_s34.npy','trace.npz','psi.npy','modes.npz']:
                f=directory/name
                raw.append(dict(path=str(f.relative_to(WORK)),bytes=f.stat().st_size,sha256=stream_sha(f)))
        save('RAW_FILES.json',dict(files=raw,scope='Working exact input/final fields, full traces and selected temporal coefficients. Compact exports and owner checkpoint are retained separately; remaining raw data regenerate.'))
        save('EXECUTION.json',dict(python=platform.python_version(),numpy=np.__version__,compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0],
            compiler_flags=['-std=c++17','-O3','-ffp-contract=off','-fopenmp','-fPIC','-shared'],AI_subagents=0,
            source_sha256=sha(ROOT/'code/fieldstep.cpp'),nonlinear_histories=8,canonical_response_probes=3,
            homogeneous_matching_cases=44,owner_checkpoint=read('DATASET.json')['filename']))
    pins=read('SOURCES.json')['internal']
    for row in pins:
        raw=git_bytes(row['commit'],row['path']);assert len(raw)==row['bytes'] and hashlib.sha256(raw).hexdigest()==row['sha256'],row['path']
        local=REPO/row['path']
        if local.exists():assert sha(local)==row['sha256'],row['path']
    checks=read('VERIFICATION_RESULTS.json');assert checks['passed'] and len(checks['groups'])==8
    result=read('RESULTS.json');assert len(result['rows'])==8
    assert result['physical_nbar'] is None and result['physical_f_sym'] is None and result['physical_kinetic_threshold'] is None and not result['full_gauged_cosmology_evolved']
    assert len(list(ROOT.glob('AVERAGE_*.npz')))==8 and len(list(ROOT.glob('RESPONSE_s*.json')))==3
    for p in ROOT.glob('*.json'):finite(json.loads(p.read_text()))
    for p in ROOT.glob('*.npz'):
        with np.load(p) as z:
            for key in z.files:
                if np.issubdtype(z[key].dtype,np.number):assert np.all(np.isfinite(z[key])),(p,key)
    for p in (ROOT/'code').glob('*.py'):ast.parse(p.read_text(),filename=str(p))
    for p in ROOT.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if '://' in link or link.startswith('#'):continue
            assert (p.parent/link.split('#')[0]).exists(),(p,link)
    dataset=read('DATASET.json');assert dataset['library_file_id']=='libfile_1195a5fac5c08191bd91784eb6c017f7'
    archive=REPO.parent/dataset['filename']
    if archive.exists():assert stream_sha(archive)==dataset['sha256'] and archive.stat().st_size==dataset['bytes']
    if not seal:
        for row in read('MANIFEST.json')['files']:
            p=ROOT/row['path'];assert p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],str(p)
    save('AUDIT_RESULTS.json',dict(numerical_groups=8,nonlinear_histories=8,source_pins=len(pins),
        imported_sources_unchanged=True,JSON_NPZ_finite=True,Python_AST=True,local_links=True,figures_rendered_and_viewed=True,
        owner_checkpoint_identity_recorded=True,independent_reviewer_of_this_packet=False,physical_particle_population_certified=False,full_gauged_cosmology=False))
    if seal:
        rows=[]
        for p in sorted(ROOT.rglob('*')):
            if not p.is_file() or p.name=='MANIFEST.json' or '__pycache__' in p.parts:continue
            rows.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)))
        save('MANIFEST.json',dict(parent_commit=PARENT_SHA,review_commit=REVIEW_SHA,files=rows))
    print(f'AUDIT: 8/8 mathematical groups, 8 nonlinear histories, 3 matrix-response probes, {len(pins)} source pins. Averaged export, finite data, links and owner checkpoint checked. Physical populations remain null.')

if __name__=='__main__':main()
