"""R1: finite source/current diagnostics and held-out linear responses."""
import argparse,hashlib,itertools,json,time
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import expm_multiply
from shared import HERE,card,write,encode,reconstruction,schwinger,d1_reference
from r1_reference import reference,matrix_diagnostics

def run(n,refined,inputs):
    start=time.monotonic();cfg=card();label=f'N{n}'+('_refined' if refined else '')
    entry=json.loads((HERE/'INPUTS.json').read_text())[label];path=inputs/entry['filename']
    if hashlib.sha256(path.read_bytes()).hexdigest()!=entry['sha256']:raise ValueError('Checkpoint hash mismatch')
    saved=np.load(path,allow_pickle=False);g=saved['ground'];input_norm=float(np.vdot(g,g).real);g=g/np.sqrt(input_norm)
    a=cfg['length']/(n+1);ham,electric,states=schwinger.hamiltonian(n,a**-2,-a/8);h=ham*a/2;del ham
    if not np.array_equal(states,saved['states']):raise ValueError('Basis mismatch')
    E=float(np.vdot(g,h@g).real);residual=float(np.linalg.norm(h@g-E*g))
    _,_,wf,bf,_=reconstruction(cfg,n);stagger=(-1.)**np.arange(n-1)
    F=electric@wf+(bf.sum()-wf@stagger)/4;phase=np.exp(1j*cfg['coupling']*F)
    probes=d1_reference.probe_weights(cfg['length'],n);names=cfg['probe_names'];v=np.array([probes[name]['weights'] for name in names])
    field=np.column_stack([electric@v[i]+probes[name]['constant']-v[i]@stagger/4 for i,name in enumerate(names)]);del electric
    means=np.array([np.vdot(g,field[:,i]*g).real for i in range(len(names))]);xi=(field-means)*g[:,None]
    C=xi.conj().T@xi;K=xi.conj().T@(h@xi-E*xi);K=(K+K.conj().T)/2
    trace=float(h.diagonal().sum());steps=2 if refined else 1;free_error=0.
    def evolve(B,dt):
        nonlocal free_error
        norms=np.sum(abs(B)**2,axis=0)
        for _ in range(steps):
            step=dt/steps;B=expm_multiply(-1j*step*h,B,traceA=-1j*step*trace)*np.exp(1j*E*step)
        free_error=max(free_error,float(np.max(abs(np.sum(abs(B)**2,axis=0)-norms))))
        return B
    responses={'0':encode(C)};current=xi.copy();previous=0.
    for t in cfg['response_times'][1:]:
        current=evolve(current,t-previous);responses[f'{t:g}']=encode(xi.conj().T@current);previous=t
    del current,xi
    # Hopping matrices in a small unitary-state basis, independent of H application.
    def hopping_matrices(B):
        matrices=[]
        for j in range(n-1):
            idx=np.nonzero((((states>>j)&1)==0)&(((states>>(j+1))&1)==1))[0]
            other=np.searchsorted(states,states[idx]^((1<<j)|(1<<(j+1))))
            cross=B[idx,:].conj().T@B[other,:]
            matrices.append((cross+cross.conj().T)/(2*a))
        return np.array(matrices)
    def state_rules(B,coeffs):
        gram=B.conj().T@B;hops=hopping_matrices(B);out={}
        for label,x in coeffs.items():
            norm=x.conj()@gram@x
            if norm.real<=0:raise ValueError('Invalid branch probability')
            hj=np.array([(x.conj()@J@x).real/norm.real for J in hops])
            S=-(v*hj)@v.T
            out[label]={'probability':float(norm.real),'hopping_expectations':hj.tolist(),'fsum':S.tolist()}
        return out,hops
    coeff_single={'vacuum':np.array([1,0],complex),'kick':np.array([0,1],complex),'plus':np.array([1,1],complex)/2,'minus':np.array([1,-1],complex)/2}
    paths=[(0,0),(1,0),(0,1),(1,1)]
    coeff_double={''.join('+' if s==1 else '-' for s in signs):np.array([signs[0]**b[0]*signs[1]**b[1] for b in paths],complex)/4 for signs in itertools.product((1,-1),repeat=2)}
    def observe(B,coeffs):
        out={'gram':encode(B.conj().T@B),'probes':{}}
        for i,name in enumerate(names):
            G=field[:,i];A=B.conj().T@(G[:,None]*B);Q=B.conj().T@((G*G)[:,None]*B);records={}
            for label,x in coeffs.items():
                psi=B@x;p=float(np.vdot(psi,psi).real)
                first=float(np.vdot(psi,G*psi).real/p);second=float(np.vdot(psi,G*G*psi).real/p)
                records[label]={'probability':p,'mean':first,'second_moment':second,'variance':second-first*first}
            out['probes'][name]={'first_matrix':encode(A),'second_matrix':encode(Q),'states':records}
        return out
    pair=np.column_stack([g,phase*g]);four=np.column_stack([pair,phase[:,None]*pair])
    single_rules,single_hops=state_rules(pair,coeff_single);instant_rules,_=state_rules(four,coeff_double)
    twice_rules,_=state_rules(four,{'double_kick':np.array([0,0,0,1],complex)})
    instant={'single':{'0':observe(pair,coeff_single)},'double':observe(four,coeff_double),
      'probe_covariance':{name:{other:float(C[i,j].real) for j,other in enumerate(names)} for i,name in enumerate(names)}}
    hop0=np.array(single_rules['vacuum']['hopping_expectations']);S0=np.array(single_rules['vacuum']['fsum'])
    sumrule_identity=float(np.max(abs(S0-2*K.real)))
    kick_scan={}
    for coupling in cfg['instant_kick_scan']:
        hj=hop0*np.cos(coupling*wf);S=-(v*hj)@v.T
        kick_scan[f'{coupling:g}']={'fsum':S.tolist(),'maximum_link_phase':float(np.max(abs(coupling*wf)))}
    cosine_check=float(np.max(abs(np.array(single_rules['kick']['hopping_expectations'])-hop0*np.cos(cfg['coupling']*wf))))
    twice_cosine_check=float(np.max(abs(np.array(twice_rules['double_kick']['hopping_expectations'])-hop0*np.cos(2*cfg['coupling']*wf))))
    del four,single_hops
    pair=evolve(pair,cfg['backaction_times'][1]);four=np.column_stack([pair,phase[:,None]*pair]);del pair
    four=evolve(four,cfg['backaction_probe_time']-cfg['backaction_times'][1]);delayed_rules,_=state_rules(four,coeff_double);del four
    ref=reference(cfg,n,cfg['modes']);fine=reference(cfg,n,cfg['refined_modes'])
    instant_cfg=dict(cfg,probe_times=[0.],backaction_times=[0.,0.],two_backaction_probe_time=0.)
    instant_ref=d1_reference.reference(instant_cfg,n,cfg['modes'],wf,bf);instant_fine=d1_reference.reference(instant_cfg,n,cfg['refined_modes'],wf,bf)
    primary_responses={t:np.array(mat)[:3,:3,0]+1j*np.array(mat)[:3,:3,1] for t,mat in responses.items()}
    diagnostics=matrix_diagnostics(C[:3,:3].real,K[:3,:3].real,primary_responses,
      minimum_eigenvalue=cfg['gram_minimum_eigenvalue'],condition_limit=cfg['gram_condition_limit'])
    return {'N':n,'refined':refined,'a':a,'dimension':len(g),'input_sha256':entry['sha256'],'input_norm_error':abs(input_norm-1),
      'ground_residual':residual,'ground_energy':E,'probe_names':names,'probe_weights':v.tolist(),'detector_weights':wf.tolist(),'detector_boundary':bf.tolist(),
      'means':means.tolist(),'covariance':C.real.tolist(),'first_moment':K.real.tolist(),'responses':responses,'diagnostics':diagnostics,
      'max_static_imaginary':float(max(np.max(abs(C.imag)),np.max(abs(K.imag)))),'max_free_norm_error':free_error,
      'sumrule_identity_error':sumrule_identity,'kick_cosine_identity_error':cosine_check,'twice_kick_cosine_identity_error':twice_cosine_check,
      'single_current_rules':single_rules,'coincident_double_current_rules':instant_rules,'twice_kick_current_rule':twice_rules,
      'delayed_double_current_rules':delayed_rules,'instant_kick_scan':kick_scan,'instant_observables':instant,
      'reference':ref,'refined_reference':fine,'instant_reference':instant_ref,'instant_refined_reference':instant_fine,
      'seconds':time.monotonic()-start}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--N',type=int,required=True);p.add_argument('--refined',action='store_true')
    p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args();cfg=card()
    if args.N not in cfg['N'] or (args.refined and args.N!=22):raise ValueError('Unregistered point')
    if args.output.exists():raise FileExistsError('Preserve existing output')
    d=run(args.N,args.refined,args.inputs);write(args.output,d);print(json.dumps({k:d[k] for k in ('N','refined','seconds')}),flush=True)
