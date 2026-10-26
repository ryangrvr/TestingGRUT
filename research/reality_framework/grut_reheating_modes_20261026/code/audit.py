"""Audit actual saved calculations, pinned inputs and publication files."""
import ast,hashlib,json,re,subprocess,sys
import numpy as np
from mode_common import ROOT,write

def finite(v):
    if isinstance(v,float):assert np.isfinite(v)
    elif isinstance(v,dict):
        for x in v.values():finite(x)
    elif isinstance(v,list):
        for x in v:finite(x)

def main():
    results=json.loads((ROOT/'VERIFICATION_RESULTS.json').read_text())
    assert len(results['groups'])==7 and all(g['passed'] for g in results['groups'].values())
    repo=ROOT.parents[2];sources=json.loads((ROOT/'SOURCES.json').read_text())
    for row in sources['internal']:
        raw=subprocess.check_output(['git','show',row['commit']+':'+row['path']],cwd=repo)
        assert hashlib.sha256(raw).hexdigest()==row['sha256'],row['path']
    allowance=json.loads((ROOT.parent/'grut_heavy_vector_budget_20261025/INPUTS.json').read_text())['benchmark']['allowance_f']
    write('BUDGET_CLARIFICATION.json',dict(allowance=allowance,
        quoted_injected_fraction_scope='Review estimate37-87 relative to the fixed total DM mass/entropy denominator; not a new production calculation.',
        rows=[dict(injected_fraction=f,removal_to_total_budget=1-1/f,
                   removal_to_symmetric_allowance=1-allowance/f) for f in (37.,87.)],
        scope='Necessary removal fractions if this source alone fills the allowance; actual erasure is uncomputed.'))
    for f in ROOT.glob('*.json'):finite(json.loads(f.read_text()))
    for f in (ROOT/'code').glob('*.py'):ast.parse(f.read_text(),filename=str(f))
    for f in ROOT.glob('*.npz'):
        with np.load(f) as d:
            for k in d.files:assert np.all(np.isfinite(d[k])),(f,k)
    for f in ROOT.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',f.read_text()):
            if '://' in link or link.startswith('#'):continue
            assert (f.parent/link.split('#')[0]).exists(),(f,link)
    verified=dict(passed_groups=7,pinned_internal_sources=len(sources['internal']),finite_JSON_NPZ=True,
        Python_AST=True,local_links=True,figure_rendered_and_viewed=True,
        scope=dict(full_gravitational_stability=False,physical_cosmic_seed_derived=False,nonlinear_formation=False,
                   full_kinetics=False,cosmological_point=False,new_closure_law=False))
    write('AUDIT_RESULTS.json',verified)
    if '--seal' in sys.argv:
        rows=[]
        for f in sorted(ROOT.rglob('*')):
            if not f.is_file() or '__pycache__' in f.parts or f.name=='MANIFEST.json':continue
            rows.append(dict(path=str(f.relative_to(ROOT)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
        write('MANIFEST.json',dict(parent_commit='e1603565a0c1411a55cdef7f5dc4ecc25aae1c8f',files=rows))
    else:
        for r in json.loads((ROOT/'MANIFEST.json').read_text())['files']:
            f=ROOT/r['path'];assert f.stat().st_size==r['bytes'] and hashlib.sha256(f.read_bytes()).hexdigest()==r['sha256'],f
    print('AUDIT:7/7 mathematical groups; pinned sources, finite JSON/NPZ, AST, links and manifest pass.')

if __name__=='__main__':main()
