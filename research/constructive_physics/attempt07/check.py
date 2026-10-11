"""Exact state-vector and independent branch controls for frozen CP-07.
NumPy only. Run: python check.py. No parameter fitting or Monte Carlo.
"""
import itertools, json, math
from pathlib import Path
import numpy as np
I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1,-1]).astype(complex)
P=[np.diag([1,0]),np.diag([0,1])]
theta0=math.pi/4; g=math.pi/12; gamma0=math.pi/4; phi0=math.pi/3

def gate(gamma,phi):
    # Tensor order: relation, matter i, matter j.
    R=np.zeros((8,8),complex)
    for i,j in itertools.product(range(2),repeat=2):
        th=theta0+g*(1-2*i)*(1-2*j)
        rot=math.cos(th)*I-1j*math.sin(th)*Y
        R+=np.kron(rot,np.kron(P[i],P[j]))
    phase=np.kron(np.diag([1,np.exp(1j*phi)]),np.eye(4))
    V=np.kron(P[0],np.eye(4))+np.kron(P[1],math.cos(gamma)*np.eye(4)-1j*math.sin(gamma)*np.kron(Z,X))
    return V@phase@R

def apply(v,U,axes,nq):
    axes=list(axes); rest=[a for a in range(nq) if a not in axes]
    order=axes+rest
    t=v.reshape([2]*nq).transpose(order)
    t=(U@t.reshape(2**len(axes),-1)).reshape([2]*nq)
    return t.transpose(np.argsort(order)).reshape(-1)

def evolve(n,gamma=gamma0,phi=phi0,initial=None):
    edges=[(i,j) for j in range(1,n) for i in range(j)]
    nq=n+len(edges)
    v=np.zeros(2**nq,complex)
    if initial is None:v[0]=1
    else:v= np.kron(initial, np.eye(1,2**len(edges),0).ravel())
    U=gate(gamma,phi)
    for a,(i,j) in enumerate(edges):v=apply(v,U,[n+a,i,j],nq)
    return v,edges

def independent_branches(n,gamma=gamma0,phi=phi0):
    # Kraus branch multiplication in matter space, independently of edge gates.
    edges=[(i,j) for j in range(1,n) for i in range(j)]
    m=np.zeros(2**n,complex);m[0]=1
    branches=[m]; nq=n
    configs=np.arange(2**n)
    for i,j in edges:
        zi=1-2*((configs>>(n-1-i))&1)
        zj=1-2*((configs>>(n-1-j))&1)
        th=theta0+g*zi*zj
        vm=math.cos(gamma)*np.eye(4)-1j*math.sin(gamma)*np.kron(Z,X)
        new=[]
        for old in branches:
            new.append(np.cos(th)*old)
            new.append(np.exp(1j*phi)*apply(np.sin(th)*old,vm,[i,j],nq))
        branches=new
    return np.stack(branches,axis=1).reshape(-1)

def reduced(v,n):
    m=v.reshape(2**n,-1)
    return m@m.conj().T,m.conj().T@m,np.sum(abs(m)**2,axis=0)

def statistics(probs,edges,n):
    hsum=rsum=0.; hist={}
    for k,p in enumerate(probs):
        if p<1e-18:continue
        C=np.zeros((n,n),bool)
        for a,(i,j) in enumerate(edges):C[i,j]=bool((k>>(len(edges)-1-a))&1)
        for t in range(n):C |= C[:,t,None]&C[None,t,:]
        h=np.ones(n,int)
        for j in range(n):
            prev=np.flatnonzero(C[:,j]);h[j]=1+max(h[prev],default=0)
        height=int(h.max());hist[height]=hist.get(height,0)+float(p)
        hsum+=p*height;rsum+=p*C.sum()
    return {'mean_height':float(hsum),'mean_relations':float(rsum),'height_distribution':hist}

def test():
    err_unit=float(np.linalg.norm(gate(gamma0,phi0).conj().T@gate(gamma0,phi0)-np.eye(8)))
    assert err_unit<1e-12
    rows=[];max_branch=max_inverse=max_phase=0.
    for n in [3,4,5]:
        v,e=evolve(n); nq=n+len(e)
        assert abs(np.vdot(v,v)-1)<1e-12
        b=independent_branches(n); be=float(np.linalg.norm(v-b));max_branch=max(max_branch,be)
        assert be<1e-12
        inverse=v.copy();U=gate(gamma0,phi0)
        for a in reversed(range(len(e))):inverse=apply(inverse,U.conj().T,[n+a,*e[a]],nq)
        target=np.zeros_like(v);target[0]=1
        ie=float(np.linalg.norm(inverse-target));max_inverse=max(max_inverse,ie);assert ie<1e-12
        rm,rg,p=reduced(v,n)
        vphase,_=evolve(n,phi=0);rm0,rg0,p0=reduced(vphase,n)
        pe=float(np.max(abs(p-p0)));max_phase=max(max_phase,pe)
        assert pe<1e-12 and np.linalg.norm(rm-rm0)<1e-12
        vcontrol,_=evolve(n,gamma=0);_,_,pc=reduced(vcontrol,n)
        tv=float(np.sum(abs(p-pc))/2);assert tv>0.01
        rows.append({'N':n,'full_qubits':nq,'geometry_total_variation_vs_gamma0':tv,
                     'joint_geometry_coherence_change_with_phi':float(np.linalg.norm(rg-rg0)),
                     'growth':statistics(p,e,n),'gamma0_control':statistics(pc,e,n)})
    # Independently predicted one-edge Ramsey response: not a history weight.
    v2,_=evolve(2); v20,_=evolve(2,phi=0)
    _,edge2,_=reduced(v2,2); _,edge20,_=reduced(v20,2)
    x_phi=float(np.trace(edge2@X).real); x_zero=float(np.trace(edge20@X).real)
    expected=math.sin(2*(theta0+g))*math.cos(gamma0)
    assert abs(x_phi-expected*math.cos(phi0))<1e-12
    assert abs(x_zero-expected)<1e-12
    # Recombine the same relation coherently with an inverse preparation rotation.
    th=theta0+g; inverse_rot=math.cos(th)*I+1j*math.sin(th)*Y
    recombined=apply(v2,inverse_rot,[2],3)
    p_return=float(np.sum(abs(recombined.reshape(4,2)[:,0])**2))
    expected_return=math.cos(th)**4+math.sin(th)**4+2*(math.cos(th)*math.sin(th))**2*math.cos(gamma0)*math.cos(phi0)
    assert abs(p_return-expected_return)<1e-12
    # Numerical affine-mixture test on the 2D input subspace of first matter qubit.
    n=3; init0=np.zeros(8,complex);init0[0]=1;init1=np.zeros(8,complex);init1[4]=1
    F=np.stack([evolve(n,initial=i)[0] for i in [init0,init1]],axis=1)
    assert np.linalg.norm(F.conj().T@F-np.eye(2))<1e-12
    r1=np.array([[.5,.5],[.5,.5]],complex);r2=np.array([[.2,0],[0,.8]],complex);p=.37
    lhs=F@(p*r1+(1-p)*r2)@F.conj().T
    rhs=p*(F@r1@F.conj().T)+(1-p)*(F@r2@F.conj().T)
    mixture=float(np.linalg.norm(lhs-rhs));assert mixture<1e-12
    # An explicit observer copies the first relation. Tracing its record removes
    # cross-sector geometry coherence; inverse copy restores the entire state.
    v,e=evolve(3);nq=3+len(e); state=np.kron(v,np.array([1,0]))
    cnot=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]],complex)
    copied=apply(state,cnot,[3,nq],nq+1)
    branches=copied.reshape(-1,2)
    rho_copied=branches@branches.conj().T
    rho_before=np.outer(v,v.conj())
    record_change=float(np.linalg.norm(rho_copied-rho_before));assert record_change>.01
    restored=apply(copied,cnot,[3,nq],nq+1)
    restore=float(np.linalg.norm(restored-state));assert restore<1e-12
    result={'verdict':'CONDITIONAL_QUANTUM_MECHANISM; NOT_A_CONTINUUM_GRAVITY_THEORY',
            'unitarity_error':err_unit,'maximum_independent_branch_error':max_branch,
            'maximum_inverse_error':max_inverse,'maximum_phase_probability_error':max_phase,
            'mixture_error':mixture,'record_reduced_state_change':record_change,
            'record_inverse_error':restore,'cases':rows,
            'one_edge_X_at_phi':x_phi,'one_edge_X_at_zero_phi':x_zero,
            'coherent_recombination_return_probability':p_return,
            'fixed_edge_probability_lower_bound':.25,
            'exact_mean_height_lower_bound':'1+(N-1)/4'}
    print(json.dumps(result,indent=2))
    Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('CP07_EXACT_QUANTUM_GROWTH_CHECKS_PASS')
if __name__=='__main__':test()
