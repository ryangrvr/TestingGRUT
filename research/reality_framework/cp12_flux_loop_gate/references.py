"""The same physical coherent mean-energy references as 9d0f2ef."""
from math import sqrt
import json,time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.sparse import coo_matrix,eye
from loops import P,REL,model,resources,observables


def energy_basis(n=4,q=4,budget=256):
    sub,H,B,states,sizes,psi=model(n,q,budget)
    energies,V=eigh(H.toarray())
    assert min(np.diff(energies))>1e-10
    wd=abs(V.T.conj()@psi)**2
    def gibbs(beta):
        weights=np.exp(-beta*(energies-energies[0]));return weights/weights.sum()
    hi=1.
    while gibbs(hi)@energies>q/2:hi*=2
    beta=brentq(lambda x:gibbs(x)@energies-q/2,0.,hi,xtol=1e-13)
    return (sub,H,B,states,sizes,psi),energies,V,dict(spectral_dephased=wd,maximum_entropy=gibbs(beta)),beta


def inclusion(full,source):
    sub,H,B,states,sizes,psi=full;src,Hs,Bs,ss,sz,ps=source
    index={s:i for i,s in enumerate(states)};rr=[];cc=[];vv=[]
    for i,s in enumerate(ss):
        key,sign,size=sub.canonical(s)
        if sign and key in index:
            rr.append(i);cc.append(index[key]);vv.append(sign*sqrt(sz[i]/size))
    J=coo_matrix((vv,(rr,cc)),shape=(len(ss),len(states))).tocsr()
    iso=float(np.linalg.norm((J.T@J-eye(len(states))).data))
    hdiff=float(np.linalg.norm((Hs@J-J@H).data))
    assert iso<1e-11 and hdiff<1e-10
    return J,dict(isometry_error=iso,H_intertwiner_error=hdiff,full_dimension=len(states),source_dimension=len(ss))


def run():
    started=time.monotonic();full,energies,V,weights,beta=energy_basis()
    sub,H,B,states,sizes,psi=full
    names,O,_=observables(4,4,states)
    rows=[]
    previous=json.loads((REL/'MEAN_ENERGY_CONTROL.json').read_text())
    for label,w in weights.items():
        X=V*np.sqrt(w)[None,:];p=np.sum(abs(X)**2,axis=1)
        data=dict(zip(names,map(float,p@O)))
        mean=float(w@energies);second=float(w@energies**2+np.linalg.norm(B@X)**2)
        entropy=float(-sum(x*np.log2(x) for x in w if x>0))
        ref=previous if label=='spectral_dephased' else previous['maximum_entropy_reference']
        assert abs(data['protected_cycle_any']-ref['means']['has_cycle'])<1e-11
        assert abs(mean-2)<1e-11
        rows.append(dict(reference=label,full_mean_energy=mean,full_energy_variance=second-mean*mean,
            global_entropy_bits=entropy,values=data,projected_stationarity=True,
            cutoff_bound_applicable_as_evolution=False))
    out=dict(N=4,Q=4,budget=256,inverse_temperature=beta,rows=rows,
        scope='Exactly the two positive physical comparator states from 9d0f2ef. Mean energy matches the main pure neutral preparation. These are mixed, configuration-coherent finite-support references, not dynamically selected states, full-unbounded-H equilibrium states or simultaneous resource/variance/purity matches. No propagation error attaches to a reference defined as this finite-support physical state; later full-H evolution is separately bounded.',
        seconds=time.monotonic()-started)
    (P/'ENERGY_REFERENCES.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True);return out


if __name__=='__main__':run()
