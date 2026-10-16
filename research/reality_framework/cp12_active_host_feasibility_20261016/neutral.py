"""Original-neutral host distributions and complete host-record resources."""
from math import comb,sqrt
import argparse,json,time
import numpy as np
from scipy.sparse import csr_matrix
from common import *
from graphs import aggregate


def entropy_term(p,L):return -p*np.log2(p/L) if p>0 else 0.


def run(n,q,b,k):
    start=time.monotonic();cat=json.loads((P/f'CATALOG_N{n}_Q{q}.json').read_text())
    m=quotient_model(n,q,b);H=m['H'];V,T,actions,info=reduced(m,k)
    initial=np.zeros(k,complex);initial[0]=1
    trajectory=Propagator(csr_matrix(T),actions,q/2,dt=1/64,batch=16).run(initial,[.5,1.,2.,5.])
    rows=[];classes=cat['classes'];lookup={r['mask']:i for i,r in enumerate(classes)}
    gi=np.fromiter((lookup[s[0]] for s in m['states']),int,count=len(m['states']))
    nv=np.fromiter(((s[1]&((1<<n)-1)).bit_count() for s in m['states']),int,count=len(m['states']))
    parity=(-1.)**nv;M=m['sub'].M
    diag=np.fromiter((q/2+2*v+sum(x*x for x in s[2])/8+
        sum(bool(s[0]>>j&1 and s[1]>>(n+j)&1 and s[2][2*j]==s[2][2*j+1]==0) for j in range(M))*(M-q)/(10*M)
        for s,v in zip(m['states'],nv)),float,count=len(m['states']))
    G=comb(M,q);prior=np.array([r['labelled_count']/G for r in classes])
    def action(x):return H@np.ascontiguousarray(x.real)+1j*(H@np.ascontiguousarray(x.imag))
    physical=[(0.,m['psi'],0.,{})]
    for t,(c,er) in trajectory.items():
        e=float(sum(v[0] for v in er.values())+t*info['full_finite_H_relation_defect_frobenius']+ALLOWANCE)
        physical.append((t,V@c,e,{key:float(v[0]) for key,v in er.items()}))
    del V
    for t,w,e,errors in physical:
        probs=np.bincount(gi,weights=abs(w)**2,minlength=len(classes));assert abs(probs.sum()-1)<1e-10
        geometry=aggregate(n,[(r['mask'],float(p)) for r,p in zip(classes,probs)])
        entropy=0.;entropyI=[0.,0.];host_rows=[]
        for r,p,p0 in zip(classes,probs,prior):
            I=probability_interval(float(p),e);L=r['labelled_count'];values=[entropy_term(v,L) for v in I]
            upper=max(values+[entropy_term(L/np.e,L)] if I[0]<=L/np.e<=I[1] else values)
            entropy+=entropy_term(p,L);entropyI=[entropyI[0]+min(values),entropyI[1]+upper]
            host_rows.append(dict(mask=r['mask'],labelled_count=L,probability=float(p),interval=I,
                uniform_probability=float(p0),selection_ratio=float(p/p0),probability_per_labelled_host=float(p/L)))
        work=float(np.vdot(w,diag*w).real-.5*(np.vdot(w,action(w)).real+np.vdot(parity*w,action(parity*w)).real))
        dmax=q*(M-q)/(10*M);we=dmax*(2*e+e*e)+ALLOWANCE
        tv=float(.5*np.sum(abs(probs-prior)));tv_error=2*e+e*e+ALLOWANCE
        cp=geometry['connectivity_probability'];cI=probability_interval(cp,e)
        selected=host_rows[lookup[cat['source_choice']['canonical_mask']]]
        row=dict(time=t,state_error=e,error_channels=errors,geometry=geometry,
            connectivity_interval=cI,diameter_intervals={d:probability_interval(p,e) for d,p in geometry['diameter_distribution'].items()},
            mean_ball_intervals=[[n*a for a in probability_interval(v/n,e)] for v in geometry['mean_balls']],
            host_classes=host_rows,total_variation_from_uniform=tv,
            total_variation_interval=[max(0.,tv-tv_error),min(1.,tv+tv_error)],
            complete_record_work=work,complete_record_work_interval=[work-we,work+we],
            complete_record_entropy_bits=float(entropy),complete_record_entropy_interval_bits=[max(0.,entropyI[0]),min(float(np.log2(G)),entropyI[1])],
            source_class_probability=selected['probability'],source_class_probability_interval=selected['interval'],
            fixed_source_peripheral_class_probability=selected['probability']*cat['source_choice']['source_orbit_size']/n)
        rows.append(row);np.savez(neutral_path(n,q,b,t),w=w,state_error=e)
        print(json.dumps(dict(stage='host_marginal',N=n,Q=q,time=t,connected=cI,
            diameters=geometry['diameter_distribution'],TV=row['total_variation_interval'],record_work=work)),flush=True)
    data=dict(status='PASS',N=n,Q=q,budget=b,dimension=H.shape[0],outside=m['B'].shape[0],krylov=info,rows=rows,
        preparation='Original coherent neutral ray over all labelled fixed-Q hosts; no circulation, matter or electric filter.',
        record='Every labelled-host projector P_g, all outcomes retained.',
        work_identity='K flips vertex-number parity P; R preserves it. Complete host record work is <D_phys+D_R>-1/2(<H>+<P H P>), retaining raw rewiring hidden on quotient diagonals.',
        error_scope='Full electric boundary, Krylov last column and explicit relation defect; ordinary floating point with reporting allowances.',seconds=time.monotonic()-start)
    (P/f'HOST_N{n}_Q{q}_B{b}.json').write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True);ap.add_argument('--q',type=int,required=True)
    ap.add_argument('--budget',type=int,required=True);ap.add_argument('--columns',type=int,default=128)
    a=ap.parse_args();run(a.n,a.q,a.budget,a.columns)
