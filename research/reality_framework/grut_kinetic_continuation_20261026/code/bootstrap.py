"""Locate/extract exact completed review builds without modifying them."""
from pathlib import Path
import os,subprocess

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[2]
WORK=Path(os.environ.get('GRUT_KINETIC_WORKDIR','/tmp/grut_kinetic_continuation_20261026'))
WORK.mkdir(parents=True,exist_ok=True)
COMMIT='e71730c24fe6a0b28ddd72cd30982c5b39c9aedf'
BASE='research/reality_framework/grut_formation_review_20261025'

def inputs():
    override=os.environ.get('GRUT_KINETIC_INPUT_ROOT')
    if override:return Path(override)
    existing=REPO.parent/'review_e71730c'/BASE
    if existing.exists():return existing
    out=WORK/'pinned';out.mkdir(exist_ok=True)
    paths=[BASE+'/'+p for p in ['HANDOFF_KINETIC.md','BRIEF_KINETIC.md','tracks/kin_source','tracks/kin_rates','tracks/kin_thermo','tracks/kin_solve']]
    a=subprocess.Popen(['git','archive',COMMIT,*paths],cwd=REPO,stdout=subprocess.PIPE)
    try:subprocess.run(['tar','-x','-C',str(out)],stdin=a.stdout,check=True)
    finally:a.stdout.close()
    if a.wait()!=0:raise RuntimeError('Fetch pinned review commit before reproducing.')
    return out/BASE

INPUT_ROOT=inputs()
