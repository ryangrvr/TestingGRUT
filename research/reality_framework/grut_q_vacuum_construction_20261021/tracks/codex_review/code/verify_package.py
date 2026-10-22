"""Read-only verification of the reconciled construction snapshot.

Checks cryptographic integrity, results, and identities of imported prior
packets. Does not regenerate scientific outputs or call the network.
"""
from pathlib import Path
import hashlib
import json
import subprocess

ROOT=Path(__file__).resolve().parents[3]
REPO=ROOT.parents[2]


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def git(*args):
    return subprocess.check_output(['git',*args],cwd=REPO,text=True).strip()


def main():
    m=json.loads((ROOT/'MANIFEST.json').read_text())
    errors=[]
    for p,want in m['files_sha256'].items():
        f=ROOT/p
        if not f.is_file() or digest(f)!=want:errors.append('hash: '+p)
    actual={str(x.relative_to(ROOT)) for x in ROOT.rglob('*') if x.is_file()
            and '__pycache__' not in x.parts and x.name!='MANIFEST.json'}
    if actual!=set(m['files_sha256']):errors.append('manifest inventory')
    d=json.loads((ROOT/'tracks/codex_review/results/review_results.json').read_text())
    if not d['passed'] or not all(d['checks'].values()) or len(d['checks'])!=17:errors.append('review dispositions')
    f=json.loads((ROOT/'tracks/relate/results/frozen_predictions.json').read_text())
    if not f['checks_passed'] or not all(f['checks'].values()) or len(f['checks'])!=5:errors.append('relation dispositions')
    if f['provenance']['source_results_sha256']!=digest(ROOT/'tracks/v_vacua/results/v3b_sm_check.json'):errors.append('relation source identity')
    if f['provenance']['freeze_code_sha256']!=digest(ROOT/'tracks/relate/code/freeze.py'):errors.append('relation program identity')
    e=json.loads((ROOT/'tracks/codex_review/results/eotwash_replay.json').read_text())
    if abs(e['mass_meV']-3.652606780424738)>1e-10:errors.append('laboratory replay')
    prefix=str(ROOT.relative_to(REPO))+'/'
    archived=git('ls-tree','-r',m['claude_checkpoint'], '--',str(ROOT.relative_to(REPO))).splitlines()
    for line in archived:
        meta,path=line.split('\t',1);expected=meta.split()[2]
        # Git object identity proves archived track files were not rewritten.
        if git('hash-object',str(REPO/path))!=expected:errors.append('changed archived file: '+path)
    for rel in m['preserved_codex_packets']:
        lines=git('ls-tree','-r',m['codex_checkpoint'],'--',rel).splitlines()
        for line in lines:
            meta,path=line.split('\t',1)
            if git('hash-object',str(REPO/path))!=meta.split()[2]:errors.append('changed imported file: '+path)
    # Local report links must exist, or be tracked source files in this sparse checkout.
    import re
    for p in [ROOT/'CONSTRUCTION.md',*list((ROOT/'tracks/codex_review').glob('*.md'))]:
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if '://' in link or link.startswith('#'):continue
            q=(p.parent/link.split('#')[0]).resolve()
            if q.exists():continue
            errors.append('broken report link: '+str(q))
    print(json.dumps({'construction_files':len(m['files_sha256']),
          'archived_Claude_files_checked':len(archived),'review_checks':len(d['checks']),
          'relation_checks':len(f['checks']),'errors':errors,'passed':not errors},indent=2))
    if errors:raise SystemExit(1)


if __name__=='__main__':main()
