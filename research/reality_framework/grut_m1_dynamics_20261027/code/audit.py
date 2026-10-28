"""Local provenance, finite data, syntax and document-link audit."""
from common import *
import ast,re
from analysis_io import read_analysis

def main():
    repo=RF.parents[1]
    pins=[]
    sources=[
      'grut_matched_candidate_20261027/CANDIDATE_M1.json',
      'grut_matched_candidate_20261027/PROFILE_M1.csv',
      'grut_matched_candidate_20261027/DEFORMATION_RESULTS.json',
      'grut_matched_candidate_20261027/COLD_MODE_RESULTS.json',
      'grut_matched_candidate_20261027/CONSEQUENCES.json',
      'grut_matched_candidate_20261027/DERIVATION.md',
      'grut_matched_candidate_20261027/VERIFICATION.json',
      'grut_nonlinear_pilot_20261026/code/field3d.cpp',
      'grut_nonlinear_pilot_20261026/code/pilot.py',
      'grut_nonlinear_pilot_20261026/code/analyze.py',
      'grut_daughter_core_20261023/code/radial.cpp',
    ]
    for name in sources:
        p=RF/name;raw=p.read_bytes();pins.append({'path':str(p.relative_to(repo)),
                'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw),
                'git_blob_sha1':hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()})
    sourcefile=ROOT/'SOURCES.json'
    new={'parent_commit':'d4747cbf062325cd54647391d9992802eb5f5c2f','local_pins':pins,
         'new_external_literature_claims':False,'independent_author_review':False,
         'history':'inherited reviewed engines are locally pinned; this changed-potential packet has same-author numerical checks only'}
    if sourcefile.exists():
        prior=json.loads(sourcefile.read_text())
        if prior['local_pins']!=pins:raise RuntimeError('Inherited source changed')
    save('SOURCES.json',new)
    for p in (ROOT/'code').glob('*.py'):ast.parse(p.read_text())
    def finite(v):
        if isinstance(v,float) and not np.isfinite(v):raise RuntimeError('Nonfinite JSON')
        if isinstance(v,dict):
            for x in v.values():finite(x)
        elif isinstance(v,list):
            for x in v:finite(x)
    for p in ROOT.glob('*.json'):finite(json.loads(p.read_text()))
    finite(read_analysis())
    for p in ROOT.glob('*.npz'):
        with np.load(p,allow_pickle=False) as archive:
            for key in archive.files:
                a=archive[key]
                if a.dtype.names:
                    for n in a.dtype.names:
                        if not np.isfinite(a[n]).all():raise RuntimeError('Nonfinite '+str(p))
                elif a.dtype.kind in 'fci' and not np.isfinite(a).all():raise RuntimeError('Nonfinite '+str(p))
    for p in ROOT.glob('*.md'):
        s=re.sub(r'\$\$[\s\S]*?\$\$','',p.read_text());s=re.sub(r'\\\[[\s\S]*?\\\]','',s)
        for link in re.findall(r'!?\[[^\]\n]+\]\(([^)]+)\)',s):
            if '://' in link or link.startswith('#'):continue
            if not (p.parent/link.split('#')[0]).exists():raise RuntimeError('Missing link '+link)
    excluded={'MANIFEST.json','AUDIT.json','ANALYSIS.json'}
    files=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts or p.name in excluded:continue
        raw=p.read_bytes();files.append({'path':str(p.relative_to(ROOT)),'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()})
    save('MANIFEST.json',{'files':files})
    save('AUDIT.json',{'source_pins':len(pins),'manifest_files':len(files),'finite_data':True,
                       'python_AST':True,'local_links':True,'AI_subagents':0,
                       'hot_population_certified':False,'actual_abundance':None,'pass':True})
    print(json.dumps({'pins':len(pins),'files':len(files),'audit_pass':True}))

if __name__=='__main__':main()
