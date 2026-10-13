"""Connection-tensor rank and existing composite gauge-tensor probes.

The weighted-Laplacian part audits a supplied operational metric, not
the CP-12 matter dispersion. The 126-state star block is exact with
unregulated electric fields; all retain the actual rewiring diagonal.
"""
from pathlib import Path
import itertools
import json
import math
import sys

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import diags
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parents[1]/"constructive_physics/attempt12"))
from model import build, gauss, GAMMA


def symmetric_components(A):
    return np.array([A[0,0],A[1,1],A[2,2],math.sqrt(2)*A[1,2],
                     math.sqrt(2)*A[0,2],math.sqrt(2)*A[0,1]])


def tt_basis(direction):
    k=np.array(direction,dtype=float);k/=np.linalg.norm(k)
    candidate=np.eye(3)[np.argmin(abs(k))]
    p=np.cross(k,candidate);p/=np.linalg.norm(p);q=np.cross(k,p)
    return [(np.outer(p,p)-np.outer(q,q))/math.sqrt(2),
            (np.outer(p,q)+np.outer(q,p))/math.sqrt(2)]


def metric_controls():
    # Unit tetrahedral nearest-neighbour directions in an explicitly supplied embedding.
    d=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]],float)/math.sqrt(3)
    directions=[np.outer(v,v)/2 for v in d]
    linear=np.column_stack([symmetric_components(A) for A in directions])
    pairs=list(itertools.combinations(range(4),2))
    composite=np.column_stack([symmetric_components(np.outer(d[a]-d[b],d[a]-d[b]))
                               for a,b in pairs])
    assert np.linalg.matrix_rank(linear)==4 and np.linalg.matrix_rank(composite)==6
    missing=[np.diag([1.,-1.,0.])/math.sqrt(2),np.diag([1.,1.,-2.])/math.sqrt(6)]
    assert max(abs(np.trace(A @ T)) for A in directions for T in missing)<1e-14
    rows=[]
    for vector in [[0,0,1],[1,0,0],[0,1,0],[1,1,1],[1,2,3]]:
        polarizations=tt_basis(vector)
        response=np.array([[np.trace(T @ A) for A in directions] for T in polarizations])
        rows.append(dict(direction=vector,linear_tt_rank=int(np.linalg.matrix_rank(response,tol=1e-12)),
                         singular_values=np.linalg.svd(response,compute_uv=False).tolist()))
    def stiffness(weights):
        W=weights.sum();mean=weights @ d
        return .5*(sum(w*np.outer(v,v) for w,v in zip(weights,d))-np.outer(mean,mean)/W)
    def band(k,weights):
        return float(weights.sum()-abs(np.sum(weights*np.exp(1j*(d @ k)))))
    finite=[]
    for weights in [np.ones(4),np.array([1.1,.9,1.05,.95])]:
        for magnitude in [.1,.05,.025]:
            k=magnitude*np.array([1.,2.,3.])/math.sqrt(14)
            exact=band(k,weights); predicted=float(k @ stiffness(weights) @ k)
            finite.append(dict(weights=weights.tolist(),k_magnitude=magnitude,exact_bloch_band=exact,
                               quadratic_prediction=predicted,scaled_fourth_order_remainder=(exact-predicted)/magnitude**4))
    # The Jacobian at the symmetric point is the four dyads even after the internal
    # A/B sublattice displacement has been eliminated.
    derivative_errors=[]
    eps=1e-5
    for a in range(4):
        w=np.ones(4);w[a]+=eps;up=stiffness(w);w[a]-=2*eps;down=stiffness(w)
        derivative_errors.append(float(np.linalg.norm((up-down)/(2*eps)-directions[a])))
    assert max(derivative_errors)<1e-9
    return dict(linear_connection_metric_rank=4,pair_tensor_rank=6,
                missing_diagonal_shear_generators=[T.tolist() for T in missing],tt_controls=rows,
                stiffness_jacobian_errors=derivative_errors,bloch_controls=finite,
                scope="Supplied weighted diamond Laplacian. Not a derived wave equation or emergent Einstein metric.")


def spectral_controls():
    n=5; pairs=list(itertools.combinations(range(n),2));M=len(pairs)
    occupied=[pairs.index((0,i)) for i in range(1,n)];mask=sum(1 << e for e in occupied)
    # On this star every leaf has E_leaf=n_leaf in {0,1}; E_center=hole-n_leaf
    # lies in {-1,0,1}. Thus L=1 includes the complete UNREGULATED star block.
    H,states,lookup,_=build(n,4,1,masks=[mask],rewire=False)
    a=GAMMA*(M-4)/M
    removable=np.array([sum(bits >> (n+e) & 1 and E[2*e]==E[2*e+1]==0
                           for e in occupied) for _,bits,E in states],int)
    H=(H+diags(a*removable)).tocsr()
    H2,states2,_,_=build(n,4,2,masks=[mask],rewire=False)
    assert states2==states and H.shape==(126,126)
    assert all(not np.any(gauss(n,s)) for s in states)
    energies,vectors=eigh(H.toarray()); ground=vectors[:,0];gap=energies-energies[0]
    u=np.array([[float(E[2*e]**2) for e in occupied] for _,_,E in states])
    pair_ids=list(itertools.combinations(range(4),2))
    pair_values=np.column_stack([u[:,i]*u[:,j] for i,j in pair_ids])
    E_basis=np.array([[1,-1,0,0,-1,1],[1,1,-2,-2,1,1]],float)
    E_basis[0]/=2;E_basis[1]/=math.sqrt(12)
    T_basis=np.array([[1,0,0,-1],[0,1,-1,0],[1,-1,-1,1]],float)
    # Explicitly orthonormal basis of the sum-zero four-leg representation.
    T_basis=np.linalg.qr(T_basis.T)[0].T
    probes_E=pair_values @ E_basis.T;probes_T=u @ T_basis.T
    matrix_E=vectors.T @ (probes_E*ground[:,None]);matrix_T=vectors.T @ (probes_T*ground[:,None])
    matrix_E[0]=0;matrix_T[0]=0
    strength_E=np.sum(matrix_E**2,axis=1);strength_T=np.sum(matrix_T**2,axis=1)
    shells=[]
    for index in range(1,len(gap)):
        if shells and abs(gap[index]-shells[-1]['gap'])<1e-9:
            shells[-1]['degeneracy']+=1;shells[-1]['E_weight']+=float(strength_E[index]);shells[-1]['T_weight']+=float(strength_T[index])
        else:shells.append(dict(gap=float(gap[index]),degeneracy=1,E_weight=float(strength_E[index]),T_weight=float(strength_T[index])))
    first_E=next(x for x in shells if x['E_weight']>1e-12)
    first_T=next(x for x in shells if x['T_weight']>1e-12)
    assert strength_E.sum()>1e-8 and strength_T.sum()>1e-8
    # Independently construct every simultaneous leaf/mode/rotor permutation.
    # Include the fermionic reordering sign, including the filled absent edges.
    group=[]
    for perm in itertools.permutations(range(1,5)):
        vertex=[0]+list(perm); mode=vertex+[n+pairs.index(tuple(sorted((vertex[i],vertex[j])))) for i,j in pairs]
        target=[]; signs=[]
        for _,bits,E in states:
            occupied_modes=[mode[j] for j in range(n+M) if bits >> j & 1]
            parity=sum(occupied_modes[i]>occupied_modes[j] for i in range(len(occupied_modes)) for j in range(i+1,len(occupied_modes)))
            nb=sum(1 << j for j in occupied_modes); ne=[0]*(2*M)
            for e,(i,j) in enumerate(pairs):
                ni,nj=vertex[i],vertex[j];f=pairs.index(tuple(sorted((ni,nj))))
                ne[2*f+(ni>nj)]=E[2*e];ne[2*f+(nj>ni)]=E[2*e+1]
            target.append(lookup[(mask,nb,tuple(ne))]);signs.append((-1)**parity)
        U=coo_matrix((signs,(target,np.arange(len(states)))),shape=H.shape).toarray().astype(float)
        assert np.linalg.norm(U.T @ U-np.eye(len(states)))<1e-12
        assert np.linalg.norm(U @ H.toarray()-H.toarray() @ U)<1e-11
        ground_phase=float(ground @ U @ ground)
        assert abs(abs(ground_phase)-1)<1e-10
        # A common phase of a symmetry operation does not affect observables;
        # choose the convention where the unique ground state is invariant.
        U/=ground_phase
        remaining=set(range(1,5));lengths=[]
        while remaining:
            start=min(remaining);j=start;length=0
            while j in remaining:remaining.remove(j);length+=1;j=vertex[j]
            lengths.append(length)
        kind=tuple(sorted(lengths,reverse=True))
        shell_E=np.flatnonzero(abs(gap-first_E['gap'])<1e-9)
        shell_T=np.flatnonzero(abs(gap-first_T['gap'])<1e-9)
        char_E=float(np.trace(vectors[:,shell_E].T @ U @ vectors[:,shell_E]))
        char_T=float(np.trace(vectors[:,shell_T].T @ U @ vectors[:,shell_T]))
        expected_E={(1,1,1,1):2,(2,1,1):0,(2,2):2,(3,1):-1,(4,):0}[kind]
        expected_T={(1,1,1,1):3,(2,1,1):1,(2,2):-1,(3,1):0,(4,):-1}[kind]
        assert abs(char_E-expected_E)<1e-9 and abs(char_T-expected_T)<1e-9
        group.append(dict(permutation=vertex,cycle_type=list(kind),E_character=char_E,T_character=char_T,
                          expected_E_character=expected_E,expected_T_character=expected_T))
    d=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]],float)/math.sqrt(3)
    polarizations=[np.diag([1.,-1.,0.])/math.sqrt(2),
                   np.array([[0,1,0],[1,0,0],[0,0,0]],float)/math.sqrt(2)]
    tensor_probes=[]
    for kind in ['linear','composite']:
        for name,T in zip(['plus','cross'],polarizations):
            weights=np.array([v @ T @ v for v in d]) if kind=='linear' else np.array([(d[i]-d[j]) @ T @ (d[i]-d[j]) for i,j in pair_ids])
            values=(u if kind=='linear' else pair_values) @ weights
            amplitude=vectors.T @ (values*ground);amplitude[0]=0
            active=np.flatnonzero(amplitude**2>1e-12)
            tensor_probes.append(dict(kind=kind,polarization=name,total_weight=float(amplitude @ amplitude),
                first_active_gap=float(gap[active[0]]) if len(active) else None,
                static_susceptibility=float(2*np.sum(amplitude[1:]**2/gap[1:]))))
    assert tensor_probes[0]['total_weight']<1e-24
    assert next(x for x in tensor_probes if x['kind']=='composite' and x['polarization']=='plus')['total_weight']>1e-8
    # Spectral matrices are genuinely independent probes, not proportional operators.
    static_E=2*(matrix_E[1:].T/gap[1:]) @ matrix_E[1:]
    static_T=2*(matrix_T[1:].T/gap[1:]) @ matrix_T[1:]
    assert np.linalg.norm(static_E-np.eye(2)*np.trace(static_E)/2)<1e-9
    assert np.linalg.norm(static_T-np.eye(3)*np.trace(static_T)/3)<1e-9
    responses=[]
    for channel,values in [('E',probes_E[:,0]),('T',probes_T[:,0])]:
        for T in [.1,.4,1.,2.]:
            evolved=expm_multiply(-1j*T*H,np.column_stack([ground,values*ground]))
            response=-2*float(np.imag(np.vdot(evolved[:,0],values*evolved[:,1])))
            spectral=2*float(np.sum((vectors.T @ (values*ground))[1:]**2*np.sin(gap[1:]*T)))
            assert abs(response-spectral)<1e-10
            eps=1e-5
            plus=expm_multiply(-1j*T*H,np.exp(1j*eps*values)*ground)
            minus=expm_multiply(-1j*T*H,np.exp(-1j*eps*values)*ground)
            derivative=float((values @ abs(plus)**2-values @ abs(minus)**2)/(2*eps))
            assert abs(derivative-response)<2e-8
            responses.append(dict(channel=channel,T=T,response=response,spectral_response=spectral,
                                  impulse_derivative=derivative,impulse_error=abs(derivative-response)))
    # Exact neutral excitation with unchanged exterior boundary variables.
    sea=sum(1 << (n+e) for e in range(M));E0=(0,)*(2*M)
    bits=sea ^ (1 << (n+occupied[0])) ^ 1;E1=list(E0);E1[2*occupied[0]]=1
    neutral=lookup[(mask,sea,E0)];pair=lookup[(mask,bits,tuple(E1))]
    boundary0=[(sea >> i & 1,E0[2*e+1]) for i,e in zip(range(1,n),occupied)]
    boundary1=[(bits >> i & 1,E1[2*e+1]) for i,e in zip(range(1,n),occupied)]
    assert boundary0==boundary1 and abs(H[pair,pair]-H[neutral,neutral]-(2+.125-a))<1e-12
    return dict(status="CP12_TENSOR_CHANNEL_CHECKS_PASS",dimension=len(states),ground_energy=float(energies[0]),
                actual_rewiring_diagonal_per_removable_edge=a,
                complete_star_basis_at_L1=True,first_composite_E_shell=first_E,first_linear_T_shell=first_T,
                total_E_spectral_weight=float(strength_E.sum()),total_T_spectral_weight=float(strength_T.sum()),
                simultaneous_permutation_controls=group,tensor_polarization_controls=tensor_probes,
                static_E_matrix=static_E.tolist(),static_T_matrix=static_T.tolist(),shells=shells,responses=responses,
                identical_boundary_control=dict(neutral_energy_expectation=float(H[neutral,neutral]),
                    neutral_pair_energy_expectation=float(H[pair,pair]),difference=float(H[pair,pair]-H[neutral,neutral]),
                    boundary_variables=boundary0),
                scope="Exact conditional star at frozen couplings; composite tensor RESPONSE, not a graviton or an emergent metric.")


if __name__=="__main__":
    out=dict(metric=metric_controls(),physical_star=spectral_controls())
    (P/"tensor_channel_results.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"metric_ranks":[out['metric']['linear_connection_metric_rank'],out['metric']['pair_tensor_rank']],
                      **{k:out['physical_star'][k] for k in ['first_composite_E_shell','first_linear_T_shell','total_E_spectral_weight','identical_boundary_control']}},indent=2))
