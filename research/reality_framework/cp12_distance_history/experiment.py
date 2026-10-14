"""Recorded pre-distance vs post-distance in the same CP12 experiment."""
from pathlib import Path
import argparse,json,time
import numpy as np
from basis import model,distance
from evolution import Evolution,dual

P=Path(__file__).resolve().parent

def run(n,q,budget,wait,lag,steps,post=True):
    start=time.monotonic();sub,H,B,states,sizes,psi=model(n,q,budget)
    E=Evolution(H,B,q/2,steps)
    D=distance(n,states,sub.es);classes=sorted(set(D.tolist()));R=len(classes)
    F=np.column_stack([D==r for r in classes]).astype(float)
    S=np.array([bits&1 for g,bits,field in states],float)
    T=np.array([(bits>>1)&1 for g,bits,field in states],float)
    assert np.max(abs(S*psi))==0
    v,_,ew=E(psi,wait);W=F*v[:,None]
    pr=np.sum(abs(W)**2,axis=0);assert abs(pr.sum()-1)<1e-10
    energy=float(np.vdot(v,H@v).real)
    recorded_energy=float(np.sum(W.conj()*(H@W)).real)
    Z,ec,ef=E(np.column_stack([W,S[:,None]*W,v,S*v]),lag,groups=[range(R),range(R,2*R)])
    X=Z[:,:R];Y=Z[:,R:2*R];x=Z[:,-2];y=Z[:,-1]
    coherence_identity_error=float(max(np.linalg.norm(x-X.sum(axis=1)),np.linalg.norm(y-Y.sum(axis=1))))
    assert coherence_identity_error<1e-10
    ex=ew+ec[:R];ey=ew+ec[R:2*R]
    exf=ew+ef[0];eyf=ew+ef[1]
    # Joint raw signed responses and positive geometry transition weights.
    C=np.array([[float(2*np.vdot(X[:,i],T*F[:,j]*Y[:,i]).imag) for j in range(R)] for i in range(R)])
    Q=np.array([[float(np.linalg.norm(F[:,j]*X[:,i])**2) for j in range(R)] for i in range(R)])
    Qs=np.array([[float(np.linalg.norm(F[:,j]*Y[:,i])**2) for j in range(R)] for i in range(R)])
    assert np.max(abs(Q.sum(axis=1)-pr))<1e-10
    rows=[]
    for i,r in enumerate(classes):
        f=F[:,i]
        a=dual(E,X[:,i],Y[:,i],T,lambda z,f=f:f*z,lambda z,f=f:f*S*z,wait,lag,float(ex[i]),float(ey[i]))
        assert abs(a['value']-C[i].sum())<1e-12
        perr=2*ew+ew**2;plo=max(0.,float(pr[i])-perr);phi=float(pr[i])+perr
        cond=None
        if plo>0:
            cond=[min(a['interval'][0]/plo,a['interval'][0]/phi),max(a['interval'][1]/plo,a['interval'][1]/phi)]
        a.update(distance=r,disconnected=(r==n),probability=float(pr[i]),probability_error=perr,conditional_interval=cond,
                 source_density_at_wait=float(np.linalg.norm(S*W[:,i])**2/pr[i]),target_density_at_wait=float(np.linalg.norm(T*W[:,i])**2/pr[i]))
        rows.append(a);print(json.dumps(dict(stage='past',N=n,lag=lag,**a)),flush=True)
    recorded=[];coherent=[]
    if post:
        for j,r in list(enumerate(classes))+[(None,'all')]:
            O=T if j is None else T*F[:,j]
            # Recording adjoint sums mutually orthogonal distance branches.
            ax=lambda z:np.sum(F*z,axis=1)
            ay=lambda z:S*np.sum(F*z,axis=1)
            a=dual(E,X,Y,O,ax,ay,wait,lag,float(exf),float(eyf))
            a['distance']=r;recorded.append(a)
            b=dual(E,x,y,O,lambda z:z,lambda z:S*z,wait,lag,float(ew+ec[-2]),float(ew+ec[-1]))
            b['distance']=r;coherent.append(b)
            if j is not None:assert abs(a['value']-C[:,j].sum())<1e-12
            else:assert abs(a['value']-C.sum())<1e-12
            print(json.dumps(dict(stage='post',N=n,lag=lag,recorded=a,unrecorded=b)),flush=True)
    change_bound=[]
    for j in range(R):
        bound=0.
        for i in range(R):
            if i==j:continue
            bound+=2*((Q[j,i]**.5+ex[j])*(Qs[j,i]**.5+ey[j])+(Q[i,j]**.5+ex[i])*(Qs[i,j]**.5+ey[i]))
        change_bound.append(float(bound))
    entropy=float(-np.sum(pr[pr>0]*np.log2(pr[pr>0])))
    out=dict(N=n,Q=q,budget=budget,wait=wait,lag=lag,steps=steps,dimension=H.shape[0],outgoing_orbits=B.shape[0],classes=classes,
             initial_energy=q/2,initial_variance=2*q,initial_entropy_bits=0,recording_entropy_bits=entropy,recording_min_entropy_bits=float(-np.log2(max(pr))),
             energy_before_record=energy,energy_after_record=recorded_energy,recording_energy_change=recorded_energy-energy,
             wait_error=ew,recorded_forward_x=exf,recorded_forward_y=eyf,past=rows,post_recorded=recorded,post_unrecorded=coherent,
             joint_responses=C.tolist(),geometry_transitions=Q.tolist(),source_filtered_geometry_transitions=Qs.tolist(),
             past_post_difference_upper_from_transition_weights=change_bound,
             offdiagonal_geometry_probability=float(Q.sum()-np.trace(Q)),
             norm_error=float(max(abs(np.linalg.norm(v)**2-1),np.max(abs(np.sum(abs(X)**2,axis=0)-pr)))),
             zero_wait_kick_is_identity=True,seconds=time.monotonic()-start,
             coherent_branch_sum_error=coherence_identity_error,
             scope='Q=N-1 auxiliary propagation family, not a 3D spatial phase. Gauge-invariant projective distance recording changes the state. Ordinary floating-point evaluation of analytic leakage bounds.')
    (P/f'N{n}_Q{q}_B{budget}_W{wait:g}_T{lag:g}.json').write_text(json.dumps(out,indent=2)+'\n')
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True);ap.add_argument('--q',type=int);ap.add_argument('--budget',type=int,required=True);ap.add_argument('--wait',type=float,default=.1);ap.add_argument('--lag',type=float,default=.1);ap.add_argument('--steps',type=int,default=16);ap.add_argument('--past-only',action='store_true');a=ap.parse_args()
    run(a.n,a.q or a.n-1,a.budget,a.wait,a.lag,a.steps,not a.past_only)
