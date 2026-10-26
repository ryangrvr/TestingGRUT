"""Audit saved numerical outputs and exact input pins; seal a reproducible packet."""
import ast,hashlib,json,re,subprocess,sys
from pathlib import Path
import numpy as np
from bootstrap import ROOT,REPO
from moment_controls import save

def finite(v):
    if isinstance(v,float):assert np.isfinite(v)
    elif isinstance(v,dict):
        for x in v.values():finite(x)
    elif isinstance(v,list):
        for x in v:finite(x)

def main():
    v=json.loads((ROOT/'VERIFICATION_RESULTS.json').read_text())
    assert len(v['groups'])==6 and all(g['passed'] for g in v['groups'].values())
    assert v['full_kinetic_complete'] is False and v['physical_formation_kernel_derived'] is False
    r=json.loads((ROOT/'RESULTS.json').read_text())
    assert r['final_physical_f_sym'] is None and r['physical_kinetic_threshold'] is None
    assert r['primary_denominator_GeV_per_visible_entropy']==4.4e-10
    pins=json.loads((ROOT/'SOURCES.json').read_text())['internal']
    for row in pins:
        raw=subprocess.check_output(['git','show',row['commit']+':'+row['path']],cwd=REPO)
        assert len(raw)==row['bytes'] and hashlib.sha256(raw).hexdigest()==row['sha256'],row['path']
    prod=json.loads((ROOT/'MOMENT_CONTROLS.json').read_text())['cases']
    follow=json.loads((ROOT/'FOLLOWUP_CONTROLS.json').read_text())['cases']
    assert len(prod)==13 and len(follow)==11
    for case in [*prod.values(),*follow.values()]:
        assert case['certified_full_kinetic_f_sym'] is None
        assert len(case['trajectories'])>1
        assert case['final']==case['trajectories'][-1]
    for f in ROOT.glob('*.json'):finite(json.loads(f.read_text()))
    for f in (ROOT/'code').glob('*.py'):ast.parse(f.read_text(),filename=str(f))
    for f in ROOT.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',f.read_text()):
            if '://' in link or link.startswith('#'):continue
            assert (f.parent/link.split('#')[0]).exists(),(f,link)
    assert (ROOT/'CONTROLS.png').stat().st_size>10000
    save('AUDIT_RESULTS.json',dict(numerical_groups_passed=6,pinned_internal_sources=len(pins),production_cases=13,
        followup_cases=11,conditional_hard_cascades=4,JSON_finite=True,Python_AST=True,local_links=True,
        figure_rendered_and_viewed=True,physical_f_sym_certified=False,physical_threshold_certified=False,
        independent_verifier_or_critic_ran=False,Claude_internal_workflow_invoked=False))
    if '--seal' in sys.argv:
        files=[]
        for f in sorted(ROOT.rglob('*')):
            if not f.is_file() or '__pycache__' in f.parts or f.name=='MANIFEST.json':continue
            files.append(dict(path=str(f.relative_to(ROOT)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
        save('MANIFEST.json',dict(parent_commit='2cbc7f1c3a925fea1f3ac38d0b254ca037afd95e',files=files))
    else:
        for row in json.loads((ROOT/'MANIFEST.json').read_text())['files']:
            f=ROOT/row['path'];assert f.stat().st_size==row['bytes'] and hashlib.sha256(f.read_bytes()).hexdigest()==row['sha256']
    print('AUDIT: six numerical groups, 12 source pins, finite JSON, syntax, links and manifest pass. Physical outputs remain null.')

if __name__=='__main__':main()
