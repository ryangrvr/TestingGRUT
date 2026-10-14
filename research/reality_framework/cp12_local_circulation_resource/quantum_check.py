"""Complete small frozen-law matrices: noncommuting Gibbs and quench tests.
Imports only the SHA-checked published CP-12 builder. These finite-cutoff
controls test identities, not the analytic unregulated thermodynamic theorem.
"""
from pathlib import Path
from math import comb, log
from itertools import combinations
import hashlib,json,sys,time,os
import numpy as np
from scipy.linalg import eigh
from scipy.special import logsumexp
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply

P=Path(__file__).resolve().parent
S=P.parent/'cp12_ground_connectivity'
assert hashlib.sha256((S/'frozen_model.py').read_bytes()).hexdigest()=='e896fce8a3e61948ff41891576b4c3ac028ad26785cafabb30bddef93c5adb9a'
sys.path.insert(0,str(S));from frozen_model import build,wilson_vector
sys.path.insert(0,str(P));from check import cycles

def kl(p,q):
    ans=0.
    if p:ans+=p*log(p/q)
    if p<1:ans+=(1-p)*log((1-p)/(1-q))
    return ans

def run(n,q,L):
    start=time.monotonic();Hs,states,lookup,_=build(n,q,L);H=Hs.toarray();dim=len(states);M=comb(n,2)
    A=np.array([q/2+2*(bits&((1<<n)-1)).bit_count()+sum(x*x for x in E)/8 for _,bits,E in states])
    V=H-np.diag(A);ev,U=eigh(H,check_finite=False)
    vmin=-2*q;vmax=2*q+2*.1*q*(M-q)/M;osc=vmax-vmin
    ve=eigh(V,eigvals_only=True,check_finite=False)
    assert ve[0]>=vmin-1e-10 and ve[-1]<=vmax+1e-10
    es=list(combinations(range(n),2));flags={}
    for g,_,_ in states:
        if g not in flags:flags[g]=bool(cycles(n,[e for i,e in enumerate(es) if g>>i&1],3))
    target=np.array([flags[g] for g,_,_ in states],bool)
    # Noncommutation is tested explicitly; diagonal/classical Gibbs substitution would be invalid.
    comm=np.linalg.norm((target[:,None].astype(int)-target[None,:].astype(int))*H)
    assert comm>1e-6
    gibbs=[];quenches=[]
    sea_bits=sum(1<<(n+e) for e in range(M));sea=[i for i,(_,bits,E) in enumerate(states) if bits==sea_bits and not any(E)]
    count=comb(M,q);assert len(sea)==count
    W0=np.zeros((dim,count),complex);W0[sea,np.arange(count)]=1/np.sqrt(count)
    eigW=U.T@W0;initialV=float(np.trace(W0.conj().T@V@W0).real)
    coherent=W0.sum(axis=1);assert abs(np.linalg.norm(coherent)-1)<1e-12
    mixture_energy=float(np.trace(W0.conj().T@H@W0).real)
    coherent_energy=float(np.vdot(coherent,H@coherent).real)
    assert abs(mixture_energy-(q/2+.1*q*(M-q)/M))<1e-11
    assert abs(coherent_energy-q/2)<1e-11
    for beta in [.3,1.,3.]:
        za=logsumexp(-beta*A);zh=logsumexp(-beta*ev)
        sw=np.exp(-beta*A-za);w=np.exp(-beta*ev-zh)
        diagonal=(abs(U)**2)@w;p=float(diagonal@target);qref=float(sw@target)
        meanA=float(diagonal@A);meanH=float(w@ev);meanV=meanH-meanA
        D=-beta*meanV+za-zh
        direct=float(np.sum(w*np.log(w)))+beta*meanA+za
        binary=kl(p,qref)
        assert abs(D-direct)<2e-11 and D>=binary-1e-10 and D<=beta*osc+1e-10
        gibbs.append(dict(beta=beta,probability=p,reference_probability=qref,
                          relative_entropy=D,binary_relative_entropy=binary,
                          entropy_identity_error=abs(D-direct),analytic_D_upper=beta*osc))
        entropy=log(count);d0=beta*q/2+za-entropy
        coherent_d0=beta*q/2+za
        assert abs(coherent_d0-d0-log(count))<1e-12
        rows=[]
        for T in [0.,.5,2.,8.]:
            W=U@(np.exp(-1j*T*ev)[:,None]*eigW)
            diag=np.sum(abs(W)**2,axis=1);pa=float(diag@target);meanA=float(diag@A)
            meanV=float(np.trace(W.conj().T@V@W).real)
            D=-entropy+beta*meanA+za
            identity=d0+beta*(initialV-meanV)
            assert abs(D-identity)<3e-10 and D<=d0+beta*osc+1e-10
            assert D>=kl(pa,qref)-1e-10
            rows.append(dict(T=T,probability=pa,relative_entropy=D,
                             binary_relative_entropy=kl(pa,qref),
                             conservation_identity_error=abs(D-identity)))
        quenches.append(dict(beta=beta,initial_D=d0,coherent_initial_D=coherent_d0,
                            coherence_entropy_difference=log(count),all_time_D_upper=d0+beta*osc,rows=rows))
    wilson=None
    if n==4 and q==4:
        cycle=[0,1,2,3];edges=[tuple(sorted(e)) for e in zip(cycle,cycle[1:]+cycle[:1])]
        mask=sum(1<<es.index(e) for e in edges);s=(mask,sea_bits,(0,)*(2*M));v=np.zeros(dim);v[lookup[s]]=1
        rr=[];cc=[];vv=[]
        for j,state in enumerate(states):
            for c in [cycle,list(reversed(cycle))]:
                z=wilson_vector(n,c,state)
                if z in lookup:rr.append(lookup[z]);cc.append(j);vv.append(1/np.sqrt(2))
        O=coo_matrix((vv,(rr,cc)),shape=Hs.shape).tocsr()
        assert np.linalg.norm((O-O.T).data)<1e-12
        Ov=O@v;O2v=O@Ov
        moment=float(2*np.vdot(Ov,Hs@Ov).real-2*np.vdot(O2v,Hs@v).real)
        expected=2*len(cycle)*(.25-.1*(M-q)/M)
        assert abs(moment-expected)<2e-12 and abs(np.linalg.norm(Ov)-1)<1e-12
        rows=[];eps=1e-5
        plus=expm_multiply(1j*eps*O,v);minus=expm_multiply(-1j*eps*O,v)
        for T in [.002,.001,.0005]:
            vp=expm_multiply(-1j*T*Hs,plus);vm=expm_multiply(-1j*T*Hs,minus)
            response=float((np.vdot(vp,O@vp).real-np.vdot(vm,O@vm).real)/(2*eps))
            rows.append(dict(T=T,response=response,response_over_T=response/T,error=abs(response/T-expected)))
        assert rows[-1]['error']<2e-5
        assert rows[-1]['error']<.3*rows[-2]['error']
        wilson=dict(normalized_Wilson_cosine_double_commutator=moment,exact_expected=expected,
                    impulse_epsilon=eps,impulse_rows=rows,
                    scope='The first moment is unregulated-exact; finite-time impulse rows use the stated cutoff.')
    out=dict(N=n,Q=q,L=L,full_dimension=dim,graph_projector_H_commutator_norm=comm,
             interaction_spectrum=[float(ve[0]),float(ve[-1])],analytic_interaction_interval=[vmin,vmax],
             gibbs=gibbs,neutral_mixture_quenches=quenches,
             neutral_preparation_full_energies=dict(mixture=mixture_energy,coherent=coherent_energy),
             wilson_response=wilson,seconds=time.monotonic()-start)
    print(json.dumps(dict(N=n,Q=q,dimension=dim,seconds=out['seconds'],status='PASS')),flush=True)
    return out

def main():
    out=dict(status='CP12_LOCAL_CIRCULATION_QUANTUM_CONTROLS_PASS',models=[run(4,3,1),run(4,4,1)],
             scope='Fresh full-matrix Gibbs and unitary calculations at finite cutoff; no large-N or zero-temperature phase inference.')
    (Path(os.environ.get('CP12_OUTPUT_DIR',str(P)))/'quantum_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
if __name__=='__main__':main()
