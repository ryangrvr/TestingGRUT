"""Aggregate actual checks and seal the packet after successful calculation.

--seal refreshes MANIFEST.json after scientific groups, source hashes,
finite JSON, Python syntax and local Markdown links pass. It does not
change any physical input or make a failed scientific group pass.
"""
import ast,hashlib,json,re,sys
from pathlib import Path
from common import ROOT,write


def finite(value):
    if isinstance(value,float):
        import math
        assert math.isfinite(value)
    elif isinstance(value,dict):
        for v in value.values():finite(v)
    elif isinstance(value,list):
        for v in value:finite(v)


def main():
    groups={}
    for name in ('REHEATING_VERIFICATION.json','BOUNDARY_RESULTS.json','PORTAL_RESULTS.json'):
        data=json.loads((ROOT/name).read_text());finite(data)
        for key,val in data['groups'].items():
            assert key not in groups
            groups[key]=dict(source=name,**val)
    assert len(groups)==11 and all(g['passed'] for g in groups.values())
    boundary=json.loads((ROOT/'BOUNDARY_RESULTS.json').read_text())
    for m in boundary['models']:
        for r in m['rows']:
            assert abs(r['state_norm']-1)<1e-11 and abs(r['state_purity']-1)<1e-11
            assert r['constraint_residual']<1e-11 and r['physical_effect_commutator']<1e-11
            assert -1e-12<=r['clock_event_probability']<=1+1e-12
    sources=json.loads((ROOT/'SOURCES.json').read_text());repo=ROOT.parents[2]
    for r in sources['internal']:
        p=repo/r['path'];assert p.is_file(),p
        assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],p
    for f in (ROOT/'code').glob('*.py'):ast.parse(f.read_text(),filename=str(f))
    for f in ROOT.glob('*.json'):finite(json.loads(f.read_text()))
    write('VERIFICATION_RESULTS.json',dict(mathematical_groups=groups,passed_groups=len(groups),
        other_checks=dict(internal_source_hashes=len(sources['internal']),finite_JSON=True,
                          Python_AST=True,local_Markdown_links=True,boundary_norm_purity_effect_controls=True),
        physical_claims=dict(complete_cosmological_point=False,derived_initial_phase=False,
            full_fluctuation_stability=False,full_hidden_kinetics=False,complete_core_lifetime=False,
            new_observational_evidence=False,new_closure_law=False),
        authorship='Same-author calculations with distinct numerical methods; Claude review pending.'))
    for f in ROOT.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',f.read_text()):
            if '://' in link or link.startswith('#'):continue
            assert (f.parent/link.split('#')[0]).exists(),(f,link)
    if '--seal' in sys.argv:
        files=[]
        for f in sorted(ROOT.rglob('*')):
            if not f.is_file() or '__pycache__' in f.parts or f.name=='MANIFEST.json':continue
            files.append(dict(path=str(f.relative_to(ROOT)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
        write('MANIFEST.json',dict(parent_commit='04c53f522fe3d275a5f4f0e3caa61eead5c6f940',files=files))
    else:
        manifest=json.loads((ROOT/'MANIFEST.json').read_text())
        for r in manifest['files']:
            f=ROOT/r['path'];assert f.stat().st_size==r['bytes']
            assert hashlib.sha256(f.read_bytes()).hexdigest()==r['sha256'],f
    print('AUDIT: 11/11 mathematical groups; source hashes, finite JSON, AST and local links pass.')


if __name__=='__main__':main()
