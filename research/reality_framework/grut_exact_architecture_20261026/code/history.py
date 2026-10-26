"""Exact finite internal-clock/record construction, with a supplied circuit.

No spatial geometry, Born-rule derivation, outcome selection or cosmic-state
selection is claimed. There is one fixed propagation constraint on 64 states.
"""
import math
import numpy as np
from common import write

I2=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Z=np.diag([1.,-1.]).astype(complex)
H2=np.array([[1,1],[1,-1]],complex)/math.sqrt(2)
N,D=8,8


def gate(bit,op):
    ans=np.array([[1.]],complex)
    for b in range(3):ans=np.kron(ans,op if b==bit else I2)
    return ans


def cnot(control,target):
    ans=np.zeros((D,D),complex)
    for j in range(D):
        bits=[(j>>(2-b))&1 for b in range(3)]
        bits[target]^=bits[control]
        out=sum(v<<(2-b) for b,v in enumerate(bits))
        ans[out,j]=1
    return ans


def build():
    had=gate(0,H2);a=cnot(0,1);b=cnot(0,2)
    ph=gate(0,np.diag(np.exp(-.37j*np.array([1.,-1.]))))
    steps=[had,a,b,ph,ph.conj().T,b.conj().T,a.conj().T,had.conj().T]
    W=[np.eye(D,dtype=complex)]
    for u in steps:W.append(u@W[-1])
    C=np.zeros((N*D,N*D),complex)
    for t,u in enumerate(steps):
        nxt=(t+1)%N;sl=slice(t*D,(t+1)*D);sn=slice(nxt*D,(nxt+1)*D)
        C[sl,sl]+=.5*np.eye(D);C[sn,sn]+=.5*np.eye(D)
        C[sn,sl]-=.5*u;C[sl,sn]-=.5*u.conj().T
    V=np.concatenate(W[:N],axis=0)/math.sqrt(N)
    dressing=np.zeros_like(C)
    for t in range(N):dressing[t*D:(t+1)*D,t*D:(t+1)*D]=W[t]
    return C,V,W,steps,dressing


def reduced(rho,keep):
    # Independent direct basis trace; eight states keep this transparent.
    n=len(keep);out=np.zeros((2**n,2**n),complex)
    rest=[b for b in range(3) if b not in keep]
    for i in range(D):
        ib=[(i>>(2-b))&1 for b in range(3)]
        ii=sum(ib[b]<<(n-1-j) for j,b in enumerate(keep))
        for j in range(D):
            jb=[(j>>(2-b))&1 for b in range(3)]
            if any(ib[b]!=jb[b] for b in rest):continue
            jj=sum(jb[b]<<(n-1-k) for k,b in enumerate(keep))
            out[ii,jj]+=rho[i,j]
    return out


def entropy(rho):
    w=np.linalg.eigvalsh(rho)
    return float(-sum(x*math.log2(x) for x in w if x>1e-14))


def classical_mi(rho,a,b):
    tab=np.zeros((2,2))
    for j,p in enumerate(np.real(np.diag(rho))):tab[(j>>(2-a))&1,(j>>(2-b))&1]+=p
    pa=tab.sum(axis=1);pb=tab.sum(axis=0)
    return float(sum(tab[i,j]*math.log2(tab[i,j]/(pa[i]*pb[j])) for i in range(2) for j in range(2) if tab[i,j]>1e-14))


def seed(p):
    s=H2@np.array([math.sqrt(1-p),math.sqrt(p)],complex)
    return np.kron(s,np.array([1,0,0,0],complex))


def main():
    C,V,W,steps,dressing=build();P=V@V.conj().T
    eig=np.linalg.eigvalsh(C)
    lap=np.eye(N)
    for t in range(N):lap[(t+1)%N,t]-=.5;lap[t,(t+1)%N]-=.5
    reference=np.kron(lap,np.eye(D))
    spectral=np.sort(np.repeat([1-math.cos(2*math.pi*k/N) for k in range(N)],D))
    projector_res=np.linalg.norm(P@P-P)
    c_v=np.linalg.norm(C@V)
    spectrum_error=float(np.max(abs(eig-spectral)))
    dressed_error=np.linalg.norm(C-dressing@reference@dressing.conj().T)

    rows=[];prob_error=0.;effect_comm=0.;effect_sum=0.;state_res=0.
    primary=None;dephased_witness=None
    for p in (0.,.1,.5,.9,1.):
        v=seed(p);history=V@v;rho=np.outer(history,history.conj())
        state_res=max(state_res,np.linalg.norm(C@history))
        ticks=[]
        for t in range(N):
            sl=slice(t*D,(t+1)*D);clock_prob=float(np.real(np.trace(rho[sl,sl])))
            conditional=rho[sl,sl]/clock_prob
            direct=np.outer(W[t]@v,(W[t]@v).conj())
            prob_error=max(prob_error,np.linalg.norm(conditional-direct),abs(clock_prob-1/N))
            effects=[];outcomes=[]
            for outcome in (0,1):
                proj=gate(1,np.diag([1.,0.]) if outcome==0 else np.diag([0.,1.]))
                # Physical-subspace POVM; raw clock effects need not commute C.
                A=np.zeros_like(C);A[sl,sl]=proj
                E=N*P@A@P
                effects.append(E);effect_comm=max(effect_comm,np.linalg.norm(C@E-E@C))
                outcomes.append(float(np.real(np.trace(rho@E))))
                ev=np.linalg.eigvalsh(V.conj().T@E@V)
                assert min(ev)>-1e-12 and max(ev)<1+1e-12
            effect_sum=max(effect_sum,np.linalg.norm(effects[0]+effects[1]-P))
            ticks.append(dict(tick=t,clock_probability=clock_prob,record1_probability_one=outcomes[1],
                              matter_entropy_bits=entropy(reduced(conditional,[0])),
                              pointer_MI_matter_R1=classical_mi(conditional,0,1),
                              pointer_MI_matter_R2=classical_mi(conditional,0,2),
                              pointer_MI_R1_R2=classical_mi(conditional,1,2),
                              global_XXX_witness=float(np.real(np.trace(conditional@gate(0,X)@gate(1,X)@gate(2,X))))))
            if p==.5 and t==3:
                mixture=np.diag(np.diag(conditional))
                dephased_witness=float(np.real(np.trace(mixture@gate(0,X)@gate(1,X)@gate(2,X))))
        prob_error=max(prob_error,abs(ticks[3]['record1_probability_one']-p))
        rows.append(dict(input_branch_probability=p,constraint_residual=float(np.linalg.norm(C@history)),
                         global_entropy_bits=entropy(rho),ticks=ticks))
        if p==.5:primary=rows[-1]

    # The maximum-entropy state on the physical subspace has conditional I/D
    # for every clock tick, regardless of the particular supplied unitary path.
    rho_max=P/D
    max_state_error=0.;max_MI=0.
    for t in range(N):
        sl=slice(t*D,(t+1)*D);r=N*rho_max[sl,sl]
        max_state_error=max(max_state_error,np.linalg.norm(r-np.eye(D)/D))
        max_MI=max(max_MI,abs(classical_mi(r,0,1)),abs(classical_mi(r,1,2)))

    # A state filter genuinely restricts admissible preparations, but leaves
    # the original physical Hilbert space and its operators unchanged. It
    # is an additional supplied law, not a consequence of closure.
    v0=seed(.5);filter_on_phys=np.eye(D)-np.outer(v0,v0.conj())
    select=V@filter_on_phys@V.conj().T
    sel_eigs=np.linalg.eigvalsh(filter_on_phys)
    selector_comm=float(np.linalg.norm(C@select-select@C))
    wrong=V@seed(.1)
    rejection=float(np.real(wrong.conj()@select@wrong))
    groups={
        'constraint_and_exact_spectrum':dict(passed=c_v<1e-12 and spectrum_error<1e-12 and dressed_error<1e-12 and projector_res<1e-12,
            constraint_isometry_residual=c_v,spectrum_error=spectrum_error,dressing_residual=dressed_error,
            physical_projector_residual=projector_res,kernel_dimension=int(sum(abs(eig)<1e-10)),gap=float(eig[D])),
        'physical_POVM_and_clock_conditioning':dict(passed=prob_error<1e-12 and effect_comm<1e-12 and effect_sum<1e-12,
            probability_and_conditional_error=prob_error,effect_constraint_commutator=effect_comm,effect_completeness_error=effect_sum),
        'records_coherence_and_erasure':dict(passed=abs(primary['ticks'][3]['pointer_MI_R1_R2']-1)<1e-12 and abs(primary['ticks'][3]['global_XXX_witness']-1)<1e-12 and abs(dephased_witness)<1e-12 and abs(primary['ticks'][7]['matter_entropy_bits'])<1e-12,
            record_window_ticks=[3,4,5],R1_R2_information_at_tick3=primary['ticks'][3]['pointer_MI_R1_R2'],
            pure_vs_dephased_XXX_at_tick3=[primary['ticks'][3]['global_XXX_witness'],dephased_witness],
            erasure_entropy_at_tick7=primary['ticks'][7]['matter_entropy_bits']),
        'state_selection_countermodels':dict(passed=max_state_error<1e-12 and max_MI<1e-12 and sum(abs(sel_eigs)<1e-10)==1 and rejection>1e-3 and state_res<1e-12 and selector_comm<1e-12,
            maximum_entropy_conditional_state_error=max_state_error,maximum_entropy_pointer_MI=max_MI,
            state_filter_zero_eigenvalues_on_original_physical_space=int(sum(abs(sel_eigs)<1e-10)),excluded_p01_filter_expectation=rejection,
            state_filter_constraint_commutator=selector_comm,
            supplied_state_probabilities=[r['input_branch_probability'] for r in rows])
    }
    result=dict(clock_readings=N,matter_and_record_dimension=D,total_dimension=N*D,physical_dimension=D,
                gate_order=['H_S','CNOT_S_R1','CNOT_S_R2','exp(-i*.37 Z_S)','inverse_phase','inverse_copy_R2','inverse_copy_R1','inverse_H'],
                cyclic_gate_product_residual=float(np.linalg.norm(W[-1]-np.eye(D))),groups=groups,
                histories=rows,max_entropy_physical_state_entropy_bits=math.log2(D),
                scope='Exact supplied finite history constraint and physical POVMs; no gravity, selected initial state, Born rule or unique outcome derived.')
    write('HISTORY_RESULTS.json',result)
    print('HISTORY',sum(bool(g['passed']) for g in groups.values()),'/',len(groups),'kernel',D)
    if not all(g['passed'] for g in groups.values()):raise SystemExit(1)


if __name__=='__main__':main()
