"""Reproduction driver (en_rerun): Codex's early3d engine and initialize(), run in tau segments with restarts
from 64-bit snapshots (so each Bash call stays under 10 minutes). Restart at an output time: clock
tau0 -> segment start, a0 -> a(s_start), same C. Step grid can differ from an unsegmented run at 1e-13.
usage: chunked_run.py <codedir> <tag> <segment_end_tau>   (runs the next segment up to that tau)
       chunked_run.py <codedir> <tag> finalize             (concatenates traces, writes run_result.json, merges RUNS.json)"""
import sys,json,struct,subprocess,time
sys.path.insert(0,sys.argv[1])
from en_common import *
from initial import initialize
JOBS={'N32':dict(N=32),'N48':dict(N=48),'N64':dict(N=64),'halfstep48':dict(N=48,hmax=.002,sample=.025),
      'seed2_48':dict(N=48,seed=2),'cut10_48':dict(N=48,cutoff=10.),'box16_64':dict(N=64,L=16.)}
tag=sys.argv[2];job=dict(N=48,L=12.,seed=1,cutoff=8.,hmax=.004,sample=.05);job.update(JOBS[tag])
d=WORK/(tag+'_chunked');d.mkdir(exist_ok=True);stfile=d/'chunk_state.json'
state=json.loads(stfile.read_text()) if stfile.exists() else None
if state is None:
    q,p,params=initialize(job['N'],job['L'],job['seed'],job['cutoff'])
    write_initial(d/'seg0_initial.bin',q,p,job['L'])
    state=dict(params=params,job=job,segments=[],tau=20.,s=0.,init='seg0_initial.bin')
    stfile.write_text(json.dumps(state,indent=2))
P=state['params'];a_m=P['a_match'];C=P['C_match']
if sys.argv[3]!='finalize':
    end_tau=float(sys.argv[3]);tau0=state['tau'];s0=state['s'];a0=a_m+C*s0
    dur=2*(end_tau-tau0)/(a0+np.sqrt(a0*a0+2*C*(end_tau-tau0)))
    i=len(state['segments']);prefix=d/f'seg{i+1}'
    args=[str(build()),str(d/state['init']),str(prefix),str(job['hmax']),repr(float(dur)),str(job['sample']),'20','4','2','64',repr(float(tau0)),repr(float(a0)),repr(float(C))]
    t=time.monotonic()
    with open(d/f'seg{i+1}.log','w') as f:subprocess.run(args,check=True,stdout=f,stderr=subprocess.STDOUT)
    snaps=sorted(d.glob(f'seg{i+1}_s*.bin'));last=snaps[-1]
    N,L,s,raw=snapshot(last);q=np.array(raw[:3]);p=np.array(raw[3:])
    write_initial(d/f'seg{i+1}_final_initial.bin',q,p,L)
    for x in snaps[:-1]:x.unlink()
    state['segments'].append(dict(prefix=f'seg{i+1}',tau0=tau0,s0=s0,a0=a0,duration=dur,end_tau=end_tau,seconds=time.monotonic()-t,args=args))
    state['tau']=end_tau;state['s']=s0+s;state['init']=f'seg{i+1}_final_initial.bin'
    stfile.write_text(json.dumps(state,indent=2));print('SEGMENT',tag,i+1,tau0,'->',end_tau,'seconds',time.monotonic()-t,flush=True)
else:
    parts=[]
    for i,seg in enumerate(state['segments']):
        tr=np.genfromtxt(str(d/seg['prefix'])+'.tsv',names=True,delimiter='\t')
        arr=np.array(tr.tolist());arr[:,0]+=seg['s0']
        parts.append(arr if i==0 else arr[1:])
    allt=np.vstack(parts);names=tr.dtype.names
    np.savetxt(d/'state.tsv',allt,delimiter='\t',header='\t'.join(names),comments='',fmt='%.17g')
    trace=np.genfromtxt(d/'state.tsv',names=True,delimiter='\t')
    lastsnap=sorted(d.glob(state['segments'][-1]['prefix']+'_s*.bin'))[-1]
    params=dict(P);params.update(hmax=job['hmax'],end_tau=state['tau'],sample_tau=job['sample'],snapshot_tau=20.,
        clock=dict(tau0=20.,a0=a_m,C=C),order=4,threads=2,chunked_restarts=[s['tau0'] for s in state['segments'][1:]])
    result=dict(tag=tag,parameters=params,seconds=sum(s['seconds'] for s in state['segments']),
        max_recorded_energy_relative_drift=float(np.max(abs(trace['rho']/trace['rho'][0]-1))),
        max_recorded_charge_relative_drift=float(np.max(abs(trace['nQ']/trace['nQ'][0]-1))),
        maximum_recorded_tau_spacing=float(np.diff(trace['tau']).max()),
        final={n:float(trace[n][-1]) for n in trace.dtype.names},snapshots=[str(lastsnap)],reproduced_by='en_rerun',chunked=True)
    (d/'run_result.json').write_text(json.dumps(result,indent=2)+'\n')
    # analyze.py reads WORK/<tag>/state.tsv: link the chunked directory under the plain tag name
    link=WORK/tag
    if not link.exists():link.symlink_to(d)
    cur=read('RUNS.json');cases={k:v for k,v in cur['cases'].items() if v.get('reproduced_by')=='en_rerun'};cases[tag]=result
    save('RUNS.json',dict(cases=cases,physical_f_sym=None,full_model_evolved=False))
    print('FINAL',tag,result['max_recorded_energy_relative_drift'],result['max_recorded_charge_relative_drift'],result['maximum_recorded_tau_spacing'])
