"""Run the matched conformal X/phi Hamiltonian; no hidden heat is inserted."""
from en_common import *
from initial import initialize

def run(tag,N=48,L=12.,seed=1,cutoff=8.,end_tau=800.,sample=.05,hmax=.004,threads=4,arrays=None,clock=None,snap=20.):
    directory=WORK/tag;directory.mkdir(exist_ok=True);prefix=directory/'state'
    for stale in directory.glob('state_s*.bin'):stale.unlink()
    if arrays is None:q,p,params=initialize(N,L,seed,cutoff)
    else:q,p=arrays;params=dict(N=N,L=L,verification_initial_state=True)
    if clock is None:clock=dict(tau0=params['tau_match'],a0=params['a_match'],C=params['C_match'])
    tau0,a0,C=clock['tau0'],clock['a0'],clock['C']
    duration=2*(end_tau-tau0)/(a0+np.sqrt(a0*a0+2*C*(end_tau-tau0)))
    init=directory/'initial.bin';write_initial(init,q,p,L);del q,p
    params.update(hmax=hmax,end_tau=end_tau,sample_tau=sample,snapshot_tau=snap,clock=clock,order=4,threads=threads)
    (directory/'parameters.json').write_text(json.dumps(params,indent=2)+'\n')
    args=[str(build()),str(init),str(prefix),str(hmax),str(duration),str(sample),str(snap),'4',str(threads),'64',str(tau0),str(a0),str(C)]
    start=time.monotonic()
    with open(directory/'run.log','w') as f:subprocess.run(args,check=True,stdout=f,stderr=subprocess.STDOUT)
    trace=np.genfromtxt(str(prefix)+'.tsv',names=True,delimiter='\t')
    assert np.all(np.isfinite(np.array(trace.tolist())))
    result=dict(tag=tag,parameters=params,seconds=time.monotonic()-start,
        max_recorded_energy_relative_drift=float(np.max(abs(trace['rho']/trace['rho'][0]-1))),
        max_recorded_charge_relative_drift=float(np.max(abs(trace['nQ']/trace['nQ'][0]-1))) if abs(trace['nQ'][0])>1e-12 else float(np.max(abs(trace['nQ']))),
        maximum_recorded_tau_spacing=float(np.diff(trace['tau']).max()),
        final={name:float(trace[name][-1]) for name in trace.dtype.names},
        snapshots=[str(x) for x in sorted(directory.glob('state_s*.bin'))])
    (directory/'run_result.json').write_text(json.dumps(result,indent=2)+'\n');return result

if __name__=='__main__':
    x=run('benchmark32',N=32,end_tau=80.);save('BENCHMARK.json',x);print(x['seconds'],x['max_recorded_energy_relative_drift'],x['max_recorded_charge_relative_drift'],flush=True)
