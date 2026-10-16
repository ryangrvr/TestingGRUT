"""Verify archived S6 or rerun both frozen protocols in a temporary copy."""
import argparse,hashlib,json,math,os,shutil,subprocess,sys,tempfile
from pathlib import Path
import campaign
import global_residual
import additional_controls
import action_enrichment
import verify

BASE=Path(__file__).resolve().parent
CHARTER_SHA='b7e80adf1abc1494db1e2fa0604a479ba188f8c18f86f95c3c14240142bae0e2'
GLOBAL_CHARTER_SHA='bf9f3c207aa7b4fd5592f5b157f243d97de4e193f97c93e9224c71db8500bf0a'


def charters():
    assert verify.sha(BASE/'CHARTER.md')==CHARTER_SHA
    assert verify.sha(BASE/'GLOBAL_RESIDUAL_CHARTER.md')==GLOBAL_CHARTER_SHA


def manifest():
    files=sorted(p for p in BASE.iterdir() if p.is_file() and p.suffix in ('.md','.py','.json')
                 and p.name not in ('manifest.json','VERIFICATION.json','RERUN.json'))
    return dict(parent_commit='98f1e861c3509c67a06032b7742c295069a9a273',
       original_charter_commit='e627bb5f8287efb47d2ab4c79efbad6fd2984132',
       global_charter_commit='23de49378a64c16c405db87c7b1a703952db9813',
       files={p.name:dict(sha256=verify.sha(p),bytes=p.stat().st_size) for p in files})


def equivalent(a,b):
    if isinstance(a,dict):return a.keys()==b.keys() and all(equivalent(a[k],b[k]) for k in a)
    if isinstance(a,list):return len(a)==len(b) and all(equivalent(x,y) for x,y in zip(a,b))
    if isinstance(a,(int,float)) and not isinstance(a,bool):return math.isclose(a,b,rel_tol=2e-12,abs_tol=2e-12)
    return a==b


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--recompute',action='store_true')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args();charters()
    if args.worker:
        original=campaign.run()
        additional_controls.run()
        extension=global_residual.run()
        action_enrichment.run()
        evidence=dict(original=original,G1=extension,additional_controls=verify.read('ADDITIONAL_CONTROLS.json'),
                      action_controls=verify.read('ACTION_CONTROLS.json'),
                      verification=verify.audit(),external_validation=False,outward_rounded=False)
        campaign.save(BASE/'RERUN.json',evidence);return
    if args.recompute:
        before=manifest()
        with tempfile.TemporaryDirectory(prefix='grut_s6_') as directory:
            root=Path(directory)
            for name in [BASE.name,'wm2_history_formation_20261016','wm2_prepared_disorder_20261016']:
                shutil.copytree(BASE.parent/name,root/name,ignore=shutil.ignore_patterns('__pycache__','*.partial'))
            environment=os.environ.copy();environment['OPENBLAS_NUM_THREADS']='1'
            subprocess.run([sys.executable,str(root/BASE.name/'reproduce.py'),'--worker'],check=True,env=environment)
            evidence=json.loads((root/BASE.name/'RERUN.json').read_text())
        assert manifest()==before,'Archive changed during recomputation'
        if args.output:
            args.output.parent.mkdir(parents=True,exist_ok=True);campaign.save(args.output,evidence)
        print(json.dumps(evidence['verification'],indent=2));return
    assert manifest()==json.loads((BASE/'manifest.json').read_text())
    audited=verify.audit();audited['manifest_sha256']=verify.sha(BASE/'manifest.json')
    assert equivalent(audited,verify.read('VERIFICATION.json'))
    print(json.dumps(audited,indent=2))


if __name__=='__main__':main()
