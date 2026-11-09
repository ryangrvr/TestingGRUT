"""Source-native lattice moments in ground and one/two-measurement states."""
import argparse,hashlib,itertools,json,time
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import expm_multiply
from shared import HERE,card,write,reconstruction,schwinger
from reference import probe_weights,reference,complex_matrix

def run(n,refined,inputs):
    start=time.monotonic();cfg=card();label=f'N{n}'+('_refined' if refined else '')
    if cfg['phase']!=0:raise ValueError('This frozen protocol uses zero detector phase')
    manifest=json.loads((HERE/'INPUTS.json').read_text());entry=manifest[label];path=inputs/entry['filename']
    if hashlib.sha256(path.read_bytes()).hexdigest()!=entry['sha256']:raise ValueError('Input checkpoint hash mismatch')
    data=np.load(path,allow_pickle=False);ground=data['ground'];energy=float(data['eigenvalues'][0]);a=cfg['length']/(n+1)
    ham,electric,states=schwinger.hamiltonian(n,a**-2,-a/8);h=ham*a/2
    if not np.array_equal(states,data['states']):raise ValueError('Input basis mismatch')
    residual=float(np.linalg.norm(h@ground-energy*ground));groundnorm=float(abs(np.vdot(ground,ground).real-1))
    _,_,weights,boundary,rchecks=reconstruction(cfg,n);stag=(-1.)**np.arange(n-1)
    field=electric@weights+(boundary.sum()-weights@stag)/4;phase=np.exp(1j*cfg['coupling']*field)
    y=electric-stag[None,:]/4;probes=probe_weights(cfg['length'],n)
    fields={name:y@np.array(p['weights'])+p['constant'] for name,p in probes.items()}
    # Independent exact energy-square identity in the supplied basis.
    originaldiag=ham.diagonal()*a/2
    smoothdiag=a/2*np.sum(y*y,axis=1)-a*(n-1)/32
    algebra=float(np.max(abs(originaldiag-smoothdiag)))
    del y,electric,ham
    coefficients_single={'vacuum':np.array([1,0],complex),'kick':np.array([0,1],complex),
       'plus':np.array([1,1],complex)/2,'minus':np.array([1,-1],complex)/2}
    paths=[(0,0),(1,0),(0,1),(1,1)]
    coefficients_double={''.join('+' if s==1 else '-' for s in signs):np.array([signs[0]**b[0]*signs[1]**b[1] for b in paths],complex)/4
        for signs in itertools.product((1,-1),repeat=2)}
    max_imag=0.;max_free=0.;trace=float(h.diagonal().sum());substeps=2 if refined else 1
    def evolve(v,dt):
        nonlocal max_free
        norms=np.sum(abs(v)**2,axis=0)
        for _ in range(substeps):
            step=dt/substeps
            v=expm_multiply(-1j*step*h,v,traceA=-1j*step*trace)*np.exp(1j*energy*step)
        max_free=max(max_free,float(np.max(abs(np.sum(abs(v)**2,axis=0)-norms))))
        return v
    def observe(basis,coeffs):
        nonlocal max_imag
        gram=basis.conj().T@basis;result={'gram':complex_matrix(gram),'probes':{}}
        for name,G in fields.items():
            first=basis.conj().T@(G[:,None]*basis);second=basis.conj().T@((G*G)[:,None]*basis);moments={}
            for state,x in coeffs.items():
                # State moments from actual vectors, independently of the reference's matrix contraction.
                v=basis@x;p=np.vdot(v,v)
                if not np.isfinite(p) or p.real<=0:raise ValueError('Invalid branch normalization')
                m=np.vdot(v,G*v)/p;s=np.vdot(v,G*G*v)/p
                max_imag=max(max_imag,abs(p.imag),abs(m.imag),abs(s.imag))
                moments[state]={'probability':float(p.real),'mean':float(m.real),'second_moment':float(s.real),'variance':float(s.real-m.real*m.real)}
            result['probes'][name]={'first_matrix':complex_matrix(first),'second_matrix':complex_matrix(second),'states':moments}
        return result
    pair=np.column_stack([ground,phase*ground]);single={'0':observe(pair,coefficients_single)}
    norm=float(np.vdot(ground,ground).real);ground_expectation=float(np.vdot(ground,h@ground).real/norm)
    kicked=pair[:,1];direct_energy=float(np.vdot(kicked,h@kicked).real/norm-ground_expectation)
    # Independent hopping-cosine calculation of the unitary kick's energy.
    hop=[]
    for j in range(n-1):
        mask=((states>>j)&1)!=((states>>(j+1))&1);idx=np.nonzero(mask)[0]
        other=np.searchsorted(states,states[idx]^((1<<j)|(1<<(j+1))))
        hop.append(float(np.vdot(ground[idx],ground[other]).real/(2*a)))
    cosine_energy=float(np.dot(np.cos(cfg['coupling']*weights)-1,hop)/norm)
    branch_energy={}
    for name,x in coefficients_single.items():
        v=pair@x;p=float(np.vdot(v,v).real)
        branch_energy[name]={'probability':p,'energy_above_vacuum':float(np.vdot(v,h@v).real/p-ground_expectation)}
    means={name:float(np.vdot(ground,G*ground).real) for name,G in fields.items()}
    covariance={name:{other:float(np.vdot(ground,G*fields[other]*ground).real-means[name]*means[other]) for other in fields} for name,G in fields.items()}
    pair=evolve(pair,cfg['backaction_times'][1]);single[f"{cfg['backaction_times'][1]:g}"]=observe(pair,coefficients_single)
    four=np.column_stack([pair,phase[:,None]*pair]);del pair
    four=evolve(four,cfg['two_backaction_probe_time']-cfg['backaction_times'][1])
    single[f"{cfg['two_backaction_probe_time']:g}"]=observe(four[:,:2],coefficients_single)
    double=observe(four,coefficients_double);double.update(unitary_paths=[list(p) for p in paths],backaction_times=cfg['backaction_times'],probe_time=cfg['two_backaction_probe_time'])
    refs=reference(cfg,n,cfg['modes'],weights,boundary);fine=reference(cfg,n,cfg['refined_modes'],weights,boundary)
    return {'N':n,'refined':refined,'dimension':len(states),'a':a,'ground_residual':residual,'ground_norm_error':groundnorm,
      'source_energy_identity_error':algebra,'input_sha256':entry['sha256'],'detector_weights':weights.tolist(),'detector_boundary':boundary.tolist(),
      'probe_weights':probes,'single':single,'double':double,'probe_covariance':covariance,
      'kick_energy':direct_energy,'kick_energy_cosine':cosine_energy,'kick_energy_identity_error':abs(direct_energy-cosine_energy),
      'single_state_energy':branch_energy,'max_free_norm_error':max_free,'max_imaginary_residual':max_imag,
      'references':refs,'refined_reference':fine,'seconds':time.monotonic()-start}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--N',type=int,required=True);p.add_argument('--refined',action='store_true')
    p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    cfg=card()
    if args.N not in cfg['N'] or (args.refined and args.N!=22):raise ValueError('Unregistered point')
    if args.output.exists():raise FileExistsError('Preserve existing output')
    result=run(args.N,args.refined,args.inputs);write(args.output,result)
    print(json.dumps({k:result[k] for k in ('N','refined','seconds')}),flush=True)
