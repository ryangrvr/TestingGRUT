"""Independent dense unitary dilations test phase, ratio and history bounds."""
import json,time
import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix
from audit import P,Evolution,evolve_paths,interference


def run():
    start=time.monotonic();rng=np.random.default_rng(593241)
    maximum_violation=0.;phase_discrepancy=0.;pinch_discrepancy=0.;cases=0
    for trial in range(24):
        d=10;k=5;rank=1 if trial%2==0 else 3
        A=rng.normal(size=(d,d));HF=(A+A.T)/8+np.diag(np.arange(d)/7)
        H=csr_matrix(HF[:k,:k]);B=csr_matrix(HF[k:,:k])
        X=rng.normal(size=(k,rank))+1j*rng.normal(size=(k,rank));X/=np.linalg.norm(X)
        Xfull=np.vstack([X,np.zeros((d-k,rank),complex)])
        S=rng.integers(0,2,size=d).astype(float)
        C=rng.integers(0,2,size=d).astype(float);C[0]=C[1]=1
        D=rng.uniform(.1,.9,size=d)
        F=rng.integers(0,2,size=d).astype(float);F[0]=1;F[1]=0
        O=dict(density=D[:k],circulation=C[:k],circulation_density=(C*D)[:k])
        OF=dict(density=D,circulation=C,circulation_density=C*D)
        wait=.08+.003*trial;lag=.1+.002*trial
        E=Evolution(H,B,0.,8);v,_,ew=E(X,wait)
        true_wait=expm(-1j*wait*HF)@Xfull
        assert np.linalg.norm(true_wait[:k]-v)<=ew+1e-12
        full,w=evolve_paths(v,ew,S[:k],H,B,O,lag)
        a,wa=evolve_paths(F[:k,None]*v,ew,S[:k],H,B,O,lag)
        b,wb=evolve_paths((1-F[:k,None])*v,ew,S[:k],H,B,O,lag)
        tests=interference(full,a,b,wa,wb,O)
        U=expm(-1j*lag*HF);actual={}
        for label,phase in [('zero',0.),('plus',np.pi/2),('minus',-np.pi/2)]:
            z=U@(np.exp(-1j*phase*S)[:,None]*true_wait)
            z1=U@(np.exp(-1j*phase*S)[:,None]*F[:,None]*true_wait)
            z0=U@(np.exp(-1j*phase*S)[:,None]*(1-F[:,None])*true_wait)
            actual[label]={name:float(np.vdot(z,d[:,None]*z).real) for name,d in OF.items()}
            for name,d in OF.items():
                value=actual[label][name];lo,hi=full['paths'][label]['intervals'][name]
                maximum_violation=max(maximum_violation,lo-value,value-hi)
                assert lo-1e-11<=value<=hi+1e-11
                cross=2*np.vdot(z1,d[:,None]*z0).real
                lo,hi=tests['offdiagonal_history_terms'][label][name]['coherent_minus_pinched_interval']
                maximum_violation=max(maximum_violation,lo-cross,cross-hi)
                assert lo-1e-11<=cross<=hi+1e-11
        x=U@true_wait;y=U@(S[:,None]*true_wait)
        for name,d in OF.items():
            derivative=2*np.vdot(x,d[:,None]*y).imag
            parameter_shift=(actual['plus'][name]-actual['minus'][name])/2
            phase_discrepancy=max(phase_discrepancy,abs(derivative-parameter_shift))
            lo,hi=full['susceptibilities'][name]['interval']
            maximum_violation=max(maximum_violation,lo-derivative,derivative-hi)
            assert lo-1e-11<=derivative<=hi+1e-11
        fp=(actual['plus']['circulation_density']-actual['minus']['circulation_density'])/2
        gp=(actual['plus']['circulation']-actual['minus']['circulation'])/2
        f=actual['zero']['circulation_density'];g=actual['zero']['circulation']
        conditional=(fp*g-f*gp)/g**2
        I=full['density_response_inside_final_circulation']['interval']
        if I is not None:
            maximum_violation=max(maximum_violation,I[0]-conditional,conditional-I[1])
            assert I[0]-1e-11<=conditional<=I[1]+1e-11
        if rank==1:
            z=F[:k,None]*X;zz=X-z;p=float(np.linalg.norm(z)**2)
            rho=X@X.conj().T;pinched=z@z.conj().T+zz@zz.conj().T
            td=float(np.sum(abs(np.linalg.eigvalsh(rho-pinched)))/2)
            pinch_discrepancy=max(pinch_discrepancy,abs(td-np.sqrt(p*(1-p))))
        cases+=1
    assert phase_discrepancy<1e-11 and pinch_discrepancy<1e-11
    out=dict(status='CP12_UNFILTERED_DILATIONS_PASS',unitary_dilation_cases=cases,
        maximum_interval_violation=float(maximum_violation),
        maximum_exact_phase_derivative_discrepancy=float(phase_discrepancy),
        maximum_pinching_trace_distance_discrepancy=float(pinch_discrepancy),
        seconds=time.monotonic()-start,
        scope='Seeded independent finite unitary dilations, with pure and mixed inputs. Tests the full-path positive-observable intervals, exact phase derivative, derivative of conditional density, and offdiagonal history interference bounds. It does not replace the unbounded-H analytic proof or provide interval-arithmetic roundoff certification.')
    (P/'INDEPENDENT_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True);return out


if __name__=='__main__':run()
