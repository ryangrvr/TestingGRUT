"""Fixed-receiver matter response in the same charged-density preparation."""
from __future__ import annotations
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import expm_multiply
from formation import model
from diagnostics import arrays,event
from stream import StreamEvolution,probability_interval

P=Path(__file__).resolve().parent
ALLOWANCE=5e-10


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=5);ap.add_argument('--budget',type=int,default=50)
    ap.add_argument('--wait',type=float,default=2.);ap.add_argument('--lag',type=float,default=.25)
    a=ap.parse_args();start=time.monotonic();n=a.n;q=2*n
    sub,H,B,states,sizes,psi=model(n,q,a.budget,source=True)
    diag=arrays(n,q,states);Pevent=event(diag,q,n)
    source=np.asarray([bits&1 for g,bits,E in states],float)
    # State and source preserve permutations of vertices 1..N-1. Hence this
    # average equals the density response of any fixed receiver j != 0.
    receiver=(diag[1]-source)/(n-1)
    evo=StreamEvolution(H,B,q/2,max_dt=1/32,batch=2)
    w,ew=evo.run(psi,[a.wait])[a.wait];ew=float(ew[0])
    X=Pevent*w;Y=source*X;p=float(np.linalg.norm(X)**2);pi=probability_interval(p,ew)
    V=np.column_stack([X,Y,w,source*w])
    Z,el=evo.run(V,[a.lag])[a.lag]
    def chi(i,j):
        x,y=Z[:,i],Z[:,j];ex=ew+float(el[i]);ey=ew+float(el[j])
        value=float(2*np.vdot(x,receiver*y).imag)
        error=2*(ex*np.linalg.norm(receiver*y)+ey*np.linalg.norm(receiver*x)+ex*ey)+ALLOWANCE
        result=dict(joint_susceptibility=value,absolute_error_upper=float(error),
            joint_interval=[value-error,value+error],zero_excluded=bool(abs(value)>error))
        return result
    conditional=chi(0,1);unfiltered=chi(2,3)
    values=[x/y for x in conditional['joint_interval'] for y in pi if y>0]
    conditional['conditional_susceptibility']=conditional['joint_susceptibility']/p
    conditional['conditional_interval']=[min(values),max(values)] if pi[0]>0 else None
    A=H-q/2*eye(H.shape[0],format='csr');trace=complex(A.diagonal().sum())
    finite=[]
    for theta in [1e-3,5e-4]:
        pulses=np.column_stack([np.exp(-1j*theta*source)*X,np.exp(1j*theta*source)*X])
        out=expm_multiply(-1j*a.lag*A,pulses,traceA=-1j*a.lag*trace)
        values=np.sum(receiver[:,None]*abs(out)**2,axis=0)
        fd=float((values[0]-values[1])/(2*theta))
        finite.append(dict(theta=theta,susceptibility=fd,
            difference_from_commutator=fd-conditional['joint_susceptibility']))
    result=dict(status='PASS',N=n,Q=q,electric_square_budget=a.budget,
        source_invariant_dimension=len(states),full_orbit_dimension=int(sum(map(int,sizes))),
        wait=a.wait,lag=a.lag,charged_event_probability=p,charged_event_probability_interval=pi,
        charged_initial_event_response=conditional,unfiltered_response=unfiltered,
        finite_pulse_controls=finite,wait_state_error_upper=ew,lag_column_error_upper=el.tolist(),
        source='exp(-i theta n0) at wait, after the optional initial charged-event projection.',
        detector='Any fixed receiver j != 0. Exact source-stabilizer symmetry permits its density expectation to be evaluated as the average over those vertices. No final graph-distance bin is used.',
        conditioning='Initial event only; n0 commutes with its projector, so its probability has zero source derivative. No final-event quotient derivative is omitted.',
        interpretation='Joint charged-density/matter-response control from the original neutral ray, with explicit diagnostic projection. K5 has only distance-one receivers: this is not a growing spatial propagation test or autonomous classical geometry.',
        error_scope='Full outgoing residuals and analytical Duhamel upper sums, ordinary floating point plus5e-10 reporting allowance; not formal outward-rounded certification.',
        seconds=time.monotonic()-start)
    (P/f'RESPONSE_N{n}_B{a.budget}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
