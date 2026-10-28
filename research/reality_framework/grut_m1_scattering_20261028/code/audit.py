"""Pin inherited inputs and audit the bounded numerical deliverables."""
from common import *
import ast,re,platform

def main():
    repo=RF.parents[1]
    paths=[PARENT/'CANDIDATE_M1.json',PARENT/'PROFILE_M1.csv',PARENT/'DEFORMATION_RESULTS.json',
           PARENT/'code/deformation.py',DYN/'code/radial.cpp',DYN/'code/radial.py',
           DYN/'RADIAL_RESULTS.json',DYN/'EXCITATION_RESULTS.json',DYN/'RADIAL_DATA.npz',DYN/'EXCITATION_DATA.npz',
           RF/'grut_exact_architecture_20261026/ACTION.md',
           RF/'grut_kinetic_continuation_20261026/code/ks_common.py']
    pins=[]
    for p in paths:
        raw=p.read_bytes();pins.append({'path':str(p.relative_to(repo)),'bytes':len(raw),
                'sha256':hashlib.sha256(raw).hexdigest(),
                'git_blob_sha1':hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()})
    prior=ROOT/'SOURCES.json'
    if prior.exists() and json.loads(prior.read_text())['local_pins']!=pins:raise RuntimeError('Inherited input changed')
    kinetic=json.loads((ROOT/'VERTICES_RESULTS.json').read_text())['old_kinetic_input']
    raw=subprocess.check_output(['git','show',kinetic['commit']+':'+kinetic['path']],cwd=repo)
    if hashlib.sha256(raw).hexdigest()!=kinetic['sha256']:raise RuntimeError('Old decay input pin mismatch')
    web=[{'title':'PDG Kinematics, 2025 review','url':'https://pdg.lbl.gov/2025/reviews/rpp2025-rev-kinematics.pdf',
          'read':'primary review retrieved; Eqs49.11,49.18 used for two-body normalization',
          'role':'phase-space convention, not M1 numeric outputs'},
         {'title':'Q-ball Superradiance','authors':'Saffin, Xie, Zhou','arxiv':'2212.03269v2',
          'doi':'10.1103/PhysRevLett.131.111601','url':'https://arxiv.org/html/2212.03269v2',
          'read':'primary HTML Eqs4-17, mode/flux discussion and 3+1D paragraph; abstract and publisher PDF retrieved',
          'role':'prior rotating-background coupled-mode scattering; no paper-code rerun or claimed M1 replication'}]
    save('SOURCES.json',{'parent_commit':'ed73791ccbcb58dc863b153f569e21ae1f29ff2e',
         'local_pins':pins,'git_decay_reference':kinetic,'primary_literature':web,
         'supplied_Argon_dispatch':'user-provided text in this conversation; response and limitations in ARGON_HANDOFF.md',
         'independent_author_review':False})
    for p in (ROOT/'code').glob('*.py'):ast.parse(p.read_text())
    def finite(v):
        if isinstance(v,float) and not np.isfinite(v):raise RuntimeError('Nonfinite JSON')
        if isinstance(v,dict):
            for x in v.values():finite(x)
        elif isinstance(v,list):
            for x in v:finite(x)
    for p in ROOT.glob('*.json'):finite(json.loads(p.read_text()))
    for p in ROOT.glob('*.npz'):
        with np.load(p,allow_pickle=False) as a:
            for key in a.files:
                v=a[key]
                if v.dtype.names:
                    for n in v.dtype.names:
                        if not np.isfinite(v[n]).all():raise RuntimeError('Nonfinite archive')
                elif v.dtype.kind in 'fci' and not np.isfinite(v).all():raise RuntimeError('Nonfinite archive')
    for p in ROOT.glob('*.md'):
        s=re.sub(r'\\\[[\s\S]*?\\\]','',p.read_text())
        for link in re.findall(r'!?\[[^\]\n]+\]\(([^)]+)\)',s):
            if '://' in link or link.startswith('#'):continue
            if not (p.parent/link.split('#')[0]).exists():raise RuntimeError('Missing link '+link)
    save('ENVIRONMENT.json',{'python':platform.python_version(),'numpy':np.__version__,
         'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0],
         'AI_subagents':0,'external_model_runs':0,'nonlinear_incoming_packet_histories':9,
         'regenerated_pulse_histories':5,'conditional_new_stationary_profiles':4})
    files=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts or p.name in ('MANIFEST.json','AUDIT.json'):continue
        raw=p.read_bytes();files.append({'path':str(p.relative_to(ROOT)),'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()})
    save('MANIFEST.json',{'files':files});save('AUDIT.json',{'pins':len(pins),'files':len(files),
        'finite_data':True,'Python_AST':True,'local_links':True,'pass':True,
        'actual_capture_history':None,'actual_final_abundance':None})
    print(json.dumps({'pins':len(pins),'files':len(files),'audit_pass':True}))

if __name__=='__main__':main()
