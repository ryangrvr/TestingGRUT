"""Reproduction driver (en_rerun): run selected Codex production jobs with threads=2 via Codex's run(); merge into RUNS.json."""
import sys,json
sys.path.insert(0,sys.argv[1])
from en_common import *
from run import run
JOBS={'N32':dict(N=32),'N48':dict(N=48),'N64':dict(N=64),'halfstep48':dict(N=48,hmax=.002,sample=.025),
      'seed2_48':dict(N=48,seed=2),'cut10_48':dict(N=48,cutoff=10.),'box16_64':dict(N=64,L=16.)}
for tag in sys.argv[2:]:
    row=run(tag,threads=2,**JOBS[tag])
    cur=read('RUNS.json') if (ROOT/'RUNS.json').exists() else dict(cases={})
    # start from a clean RUNS.json containing only reproduced cases
    cases={k:v for k,v in cur['cases'].items() if v.get('reproduced_by')=='en_rerun'}
    row['reproduced_by']='en_rerun';cases[tag]=row
    save('RUNS.json',dict(cases=cases,physical_f_sym=None,full_model_evolved=False))
    print('RUN',tag,'seconds',row['seconds'],'energy',row['max_recorded_energy_relative_drift'],'Q',row['max_recorded_charge_relative_drift'],flush=True)
