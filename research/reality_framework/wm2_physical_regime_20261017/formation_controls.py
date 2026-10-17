"""Independent small matrices and exact rational checks of the current proof.

These diagnostics do not prove the all-N estimates or validate floating-point
hardware. The all-N argument is in FORMATION_ANALYSIS.md.
"""
from pathlib import Path
from itertools import combinations_with_replacement
from fractions import Fraction as Q
from math import factorial,sqrt,pi
import json
import time
import numpy as np
from scipy.special import j1
from scipy.integrate import quad


def histories(M):
    words=[]
    def visit(w):
        words.append(w)
        for v in range(M):
            if v not in w:visit(w+(v,))
    visit(());return words


def operators(M,L,delta,g):
    words=histories(M);wi={w:i for i,w in enumerate(words)}
    children=[]
    for v in range(M):
        x=np.zeros((len(words),len(words)),complex)
        for w,j in wi.items():
            if v not in w:x[wi[w+(v,)],j]=1
        children.append(x)
    S=sum(children);K=np.diag([len(w) for w in words])
    D=np.diag([sum(delta[v] for v in w) for w in words])
    occ=[o for k in range(M+2) for o in combinations_with_replacement(range(L),k)]
    oi={o:i for i,o in enumerate(occ)};ann=[]
    for v in range(L):
        x=np.zeros((len(occ),len(occ)),complex)
        for o,j in oi.items():
            if v in o:
                a=list(o);a.remove(v);x[oi[tuple(a)],j]=sqrt(o.count(v))
        ann.append(x)
    field=sum(-(ann[v].conj().T@ann[v+1]+ann[v+1].conj().T@ann[v]) for v in range(L-1))
    A=np.kron(np.eye(len(words)),ann[0]);Sp=np.kron(S,np.eye(len(occ)))
    H=np.kron(D,np.eye(len(occ)))+np.kron(np.eye(len(words)),field)-g/(M+1)*(Sp@A.conj().T+Sp.conj().T@A)
    initial=np.zeros(len(H),complex);initial[0]=1
    return words,children,S,K,D,A,Sp,H,initial


def matrix_checks():
    rows=[]
    for M,L,delta,g in [(2,3,[-.75,.63],.9),(3,2,[-.25,.11,.75],1.1)]:
        words,children,S,K,D,A,Sp,H,x0=operators(M,L,delta,g)
        identity=0.
        for v,a in enumerate(children):
            for w,b in enumerate(children):
                want=np.diag([int(v not in u) for u in words]) if v==w else np.zeros_like(S)
                identity=max(identity,float(np.linalg.norm(a.conj().T@b-want)))
        assert identity==0 and np.linalg.norm(S.conj().T@S-(M*np.eye(len(K))-K))==0
        eig,V=np.linalg.eigh(H)
        def evolve(x,t):return V@(np.exp(-1j*t*eig)*(V.conj().T@x))
        def heis(op,t):return evolve(op@evolve(x0,t),-t)
        h=np.diag(-np.ones(L-1),1)+np.diag(-np.ones(L-1),-1)
        eh,vh=np.linalg.eigh(h)
        def C(t):return np.dot(abs(vh[0])**2,np.exp(-1j*t*eh))
        points,weights=np.polynomial.legendre.leggauss(72)
        for t in [.2,1.3,4.]:
            integral=np.zeros(len(H),complex)
            for point,weight in zip(points,weights):
                u=t*(point+1)/2;integral+=weight*C(u)*heis(Sp,t-u)*t/2
            left=heis(A,t);right=1j*g/(M+1)*integral
            error=float(np.linalg.norm(left-right));assert error<3e-12
            psi=evolve(x0,t);number=np.kron(K,np.eye(len(H)//len(K)))
            direct=float((1j*np.vdot(psi,(H@number-number@H)@psi)).real)
            current=2*g/(M+1)*float(np.vdot(heis(Sp,t),left).imag)
            assert abs(direct-current)<3e-12
            rows.append(dict(N=M+1,L=L,t=t,full_tensor_dimension=len(H),
                exact_word_identity_error=identity,receiver_equation_error=error,
                direct_count_current=direct,Heisenberg_count_current=current))
        u=.7;Aw=sum(np.exp(-1j*d*u)*a for d,a in zip(delta,children))
        creation=np.kron(S@Aw-Aw@S,np.eye(len(H)//len(S)))@A.conj().T@x0
        absorption=np.kron(S.conj().T@Aw-Aw@S.conj().T,np.eye(len(H)//len(S)))@A@x0
        assert np.linalg.norm(creation)>1e-3 and np.linalg.norm(absorption)==0
        rows.append(dict(N=M+1,creation_commutator_on_seed=float(np.linalg.norm(creation)),
            absorption_commutator_on_seed=float(np.linalg.norm(absorption)),
            dropping_creation_commutator_is_invalid=True))
        # Check the detuning cancellation and time arguments of (11).
        t,u,v,hstep=1.3,.7,.2,1e-5
        Dfull=np.kron(D,np.eye(len(H)//len(D)))
        fieldfull=H-Dfull+g/(M+1)*(Sp@A.conj().T+Sp.conj().T@A)
        interaction=H-Dfull-fieldfull
        def phased(z):return np.kron(sum(np.exp(-1j*d*z)*a for d,a in zip(delta,children)),np.eye(len(H)//len(D)))
        def retarded(y):
            op=phased(u-y);return evolve(op@evolve(x0,t-y),-(t-y))
        derivative=(retarded(v+hstep)-retarded(v-hstep))/(2*hstep)
        op=phased(u-v);commutator=-1j*(interaction@op-op@interaction)
        expected=evolve(commutator@evolve(x0,t-v),-(t-v))
        mismatch=float(np.linalg.norm(derivative-expected));assert mismatch<2e-8
        rows.append(dict(N=M+1,retarded_derivative_error=mismatch,
            D_phase_cancels=True,both_interaction_terms_retained=True))
    return rows


def kernels():
    t=np.r_[np.linspace(.001,2,100),np.geomspace(2,1e4,250)]
    c=abs(j1(2*t)/t);bound=np.minimum(1.,2/t**1.5)
    assert np.all(c<=bound+1e-14)
    moments=[]
    for L in [2,3,6]:
        h=np.diag(np.ones(L-1,dtype=object),1)+np.diag(np.ones(L-1,dtype=object),-1)
        full=np.diag(np.ones(2*L,dtype=object),1)+np.diag(np.ones(2*L,dtype=object),-1)
        a=np.eye(L,dtype=object);b=np.eye(2*L+1,dtype=object)
        for n in range(2*L):
            assert a[0,0]==b[0,0];a=a@h;b=b@full
        assert abs(int(a[0,0]-b[0,0]))==1
        moments.append(dict(L=L,matching_orders=list(range(2*L)),first_difference_order=2*L))
    rates=[]
    for d in [-.75,-.23,0.,.61,.75]:
        for eta in [.1,.02]:
            observed=quad(lambda x:sqrt(max(0.,4-x*x))/(2*pi)*eta/(eta*eta+(x+d)**2),-2,2,epsabs=1e-11)[0]
            z=eta+1j*d;analytic=((np.sqrt(z*z+4)-z)/2).real
            assert abs(observed-analytic)<2e-9
            rates.append(dict(delta=d,eta=eta,Abel_real_integral=observed,analytic=analytic,
                limiting_rate=sqrt(4-d*d)))
    return dict(kernel_grid_points=len(t),maximum_ratio=float(np.max(c/bound)),
        closed_walk_moments=moments,Abel_rate_checks=rates,
        finite_exponent_status='exact walk/Taylor argument; tiny physical tail below floating resolution')


def constants():
    exponential=sum((Q(729,125)**k)/factorial(k) for k in range(41))
    assert exponential>256
    N=2**28
    failure=Q(1,128)+Q(35,72)+Q(146,N)+Q(1,2**100)
    assert failure<Q(1,2)
    assert Q(4,9)+Q(2,3)+Q(4,5)<2
    assert Q(187,40)<5
    assert Q(55,16)>Q(81,25)
    assert Q(13,10)**2>1+Q(33,40)**2
    assert 20*Q(121,100)*3<73 and 128*Q(121,100)<155
    return dict(N0=N,formation_lower_display=float(1-failure),
        formation_lower_exact=str(1-failure),claimed_uniform_lower='1/2',
        exponential_Taylor_lower_exceeds_256=True,
        small_sizes_uncovered=True,hardware_arithmetic_certification=False)


def run():
    start=time.monotonic();r=dict(matrix_identities=matrix_checks(),kernels=kernels(),constants=constants())
    r.update(status='PASS',elapsed_seconds=time.monotonic()-start,scope='implementation diagnostics; proof is analytic, no external validation')
    return r


if __name__=='__main__':
    r=run();Path(__file__).with_name('FORMATION_CONTROLS.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))
