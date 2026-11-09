"""One frozen B2 lattice point, using the inherited Schwinger solver."""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import eigsh,expm_multiply
from model import HERE,base_weights,reconstruction,references,key,distances
from lattice import schwinger

def write(path,data):
    with Path(path).open('x',encoding='utf-8') as f:json.dump(data,f,indent=2,allow_nan=False);f.write('\n')

def run(card,n,mass=0.,refined=False):
    start=time.monotonic();a,edges,v,b=base_weights(card,n)
    _,_,w,boundary,checks=reconstruction(card,n)
    ham,electric,states=schwinger.hamiltonian(n,a**-2,mass-a/8);h=ham*a/2
    tol=card['eigensolver_refined_tol'] if refined else card['eigensolver_tol']
    initial=np.random.default_rng(card['eigensolver_seed_base']+n).standard_normal(h.shape[0])
    energies,vectors=eigsh(h,k=2,which='SA',tol=tol,v0=initial,maxiter=10000)
    order=np.argsort(energies);energies=energies[order];vectors=vectors[:,order]
    ground=vectors[:,0].astype(complex)
    residual=max(float(np.linalg.norm(h@vectors[:,j]-energies[j]*vectors[:,j])) for j in range(2))
    stagger=(-1.)**np.arange(n-1)
    fields={'B1':electric@v+(b-v@stagger)/4,
            'B2':electric@w+(boundary.sum()-w@stagger)/4}
    z=2*schwinger.bits(states,n)-1;s=np.cumsum(z,axis=1)[:,:-1]/2
    identity=float(np.max(abs(fields['B2']-(.25+s@w))))
    trace=float(h.diagonal().sum());steps=2 if refined else 1;free_error=0.
    def evolve(x,dt):
        for _ in range(steps):
            if dt:x=expm_multiply(-1j*(dt/steps)*h,x,traceA=-1j*(dt/steps)*trace)
        return x
    def histories(field):
        nonlocal free_error
        phase=np.exp(1j*(card['phase']+card['coupling']*field))[:,None];out={}
        for times in card['protocols']:
            branches=ground[:,None];labels=[''];previous=0.
            for t in times:
                before=np.linalg.norm(branches)**2;branches=evolve(branches,t-previous)
                free_error=max(free_error,abs(float(np.linalg.norm(branches)**2-before)))
                shifted=phase*branches
                branches=np.concatenate(((branches+shifted)/2,(branches-shifted)/2),axis=1)
                labels=[s+'+' for s in labels]+[s+'-' for s in labels];previous=t
            out[key(times)]={s:float(np.vdot(branches[:,j],branches[:,j]).real) for j,s in enumerate(labels)}
        return out
    output={name:histories(field) for name,field in fields.items()}
    result={'N':n,'mass_R':mass,'refined':refined,'a':a,'dimension':len(ground),
       'energies':energies.tolist(),'eigen_residual':residual,'free_norm_error':free_error,
       'gauss_identity_error':identity,'reconstruction':checks,'weights':w.tolist(),'boundary':boundary.tolist(),
       'histories':output,'field_moments':{name:{'mean':float(np.vdot(ground,f*ground).real),
           'variance':float(np.vdot(ground,f*f*ground).real-np.vdot(ground,f*ground).real**2)} for name,f in fields.items()},
       'probability_norm_error':max(abs(sum(p.values())-1) for group in output.values() for p in group.values()),
       'minimum_probability':min(x for group in output.values() for p in group.values() for x in p.values())}
    if mass==0:
        refs=references(card,n,card['modes']);fine=references(card,n,card['refined_modes'])
        result['references']=refs
        result['reference_mode_doubling']={name:distances(refs[name]['histories'],fine[name]['histories']) for name in ('target','B1','B2')}
    result['seconds']=time.monotonic()-start
    return result,{'ground':ground,'eigenvalues':energies,'states':states,'weights':w,'boundary':boundary}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--N',type=int,required=True);p.add_argument('--mass',type=float,default=0)
    p.add_argument('--refined',action='store_true');p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    card=json.loads((HERE/'protocol_b2_01.json').read_text())
    if args.N not in card['N'] or args.mass not in (0,.5):raise ValueError('Unregistered point')
    if args.output.exists() or args.output.with_suffix('.npz').exists():raise FileExistsError('Preserve existing output')
    result,checkpoint=run(card,args.N,args.mass,args.refined)
    np.savez_compressed(args.output.with_suffix('.npz'),**checkpoint);write(args.output,result)
    print(json.dumps({k:result[k] for k in ('N','mass_R','refined','seconds')}),flush=True)
