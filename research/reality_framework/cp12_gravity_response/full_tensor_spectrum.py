"""Full fluctuating-graph CP-12 spectrum for independent tensor probes.

All graphs at N=5,Q=4, all physical fermions and regulated electric fields.
The event reference/tetrahedral directions are supplied, not generated.
Probe values are bounded independently of the electric cutoff.
"""
from pathlib import Path
import argparse
import itertools
import json
import math
import sys

import numpy as np
from scipy.sparse.linalg import eigsh, cg, LinearOperator, expm_multiply

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parents[1]/"constructive_physics/attempt12"))
from model import build


def permutation_map(n,states,lookup,permutation):
    edges=list(itertools.combinations(range(n),2));M=len(edges)
    target_edges=[edges.index(tuple(sorted((permutation[i],permutation[j])))) for i,j in edges]
    mode=permutation+[n+e for e in target_edges]
    flipped=[permutation[i]>permutation[j] for i,j in edges]
    destination=[];signs=[]
    for mask,bits,E in states:
        occupied=[mode[j] for j in range(n+M) if bits >> j & 1]
        parity=sum(occupied[i]>occupied[j] for i in range(len(occupied)) for j in range(i+1,len(occupied)))
        # CP-12 is an oriented incidence coupling. Relabelling a reversed edge
        # also changes the sign of its edge fermion; omitting this is incorrect.
        parity+=sum(flip and bits >> (n+e) & 1 for e,flip in enumerate(flipped))
        newbits=sum(1 << j for j in occupied);newmask=0;newE=[0]*(2*M)
        for e,f in enumerate(target_edges):
            if mask >> e & 1:newmask |= 1 << f
            flip=flipped[e];newE[2*f+flip]=E[2*e];newE[2*f+(not flip)]=E[2*e+1]
        destination.append(lookup[(newmask,newbits,tuple(newE))]);signs.append(-1 if parity%2 else 1)
    return np.array(destination),np.array(signs)


def run(cutoff, eigenvectors=28):
    n=5;q=4; pairs=list(itertools.combinations(range(n),2));M=len(pairs)
    H,states,lookup,_=build(n,q,cutoff)
    print(json.dumps(dict(stage="built",L=cutoff,dimension=len(states),nonzeros=H.nnz)),flush=True)
    # Deterministic independent starting vector. No known spectrum is encoded.
    rng=np.random.default_rng(19312+cutoff)
    values,vectors=eigsh(H,k=eigenvectors,which='SA',tol=2e-11,maxiter=10000,v0=rng.normal(size=len(states)))
    order=np.argsort(values);values=values[order];vectors=vectors[:,order]
    ground=vectors[:,0];gaps=values-values[0]
    residuals=np.linalg.norm(H @ vectors-vectors*values,axis=0)
    assert max(residuals)<1e-7
    representation=[]
    shell=np.flatnonzero(abs(gaps-gaps[1])<1e-8)
    assert len(shell)==5
    conjugacy=[([0,1,2,3,4],[1,1,1,1,1]),([1,0,2,3,4],[2,1,1,1]),
               ([1,0,3,2,4],[2,2,1]),([1,2,0,3,4],[3,1,1]),
               ([1,2,0,4,3],[3,2]),([1,2,3,0,4],[4,1]),([1,2,3,4,0],[5])]
    for perm,kind in conjugacy:
        destination,signs=permutation_map(n,states,lookup,perm)
        def action(x):
            out=np.empty_like(x)
            out[destination]=signs[:,None]*x if x.ndim==2 else signs*x
            return out
        assert len(np.unique(destination))==len(states)
        commutator=float(np.linalg.norm(action(H @ vectors[:,shell])-H @ action(vectors[:,shell])))
        ground_error=float(np.linalg.norm(action(ground)-ground))
        character=float(np.trace(vectors[:,shell].T @ action(vectors[:,shell])))
        fixed=kind.count(1);twocycles=kind.count(2)
        expected=fixed*(fixed-3)/2+twocycles
        assert commutator<1e-8 and ground_error<1e-7 and abs(character-expected)<1e-7
        representation.append(dict(permutation=perm,cycle_type=kind,lowest_shell_character=character,
             unordered_pair_minus_vertex_character=expected,commutator_residual=commutator,ground_invariance_error=ground_error))
    star_edges=[pairs.index((0,i)) for i in range(1,n)]
    u=np.array([[float(E[2*e]!=0) for e in star_edges] for _,_,E in states])
    pair_ids=list(itertools.combinations(range(4),2))
    p=np.column_stack([u[:,i]*u[:,j] for i,j in pair_ids])
    weights_E=np.array([[1,-1,0,0,-1,1],[1,1,-2,-2,1,1]],float)
    weights_E[0]/=2;weights_E[1]/=math.sqrt(12)
    weights_T=np.linalg.qr(np.array([[1,0,0,-1],[0,1,-1,0],[1,-1,-1,1]],float).T)[0].T
    probes=dict(E=p @ weights_E.T,T=u @ weights_T.T)
    boundary=np.array([any(abs(z)==cutoff for z in E) for _,_,E in states])
    rows=[]
    A=LinearOperator(H.shape,matvec=lambda x:H @ x-values[0]*x+ground*np.vdot(ground,x),dtype=float)
    for name,operators in probes.items():
        mean=ground**2 @ operators; rhs=(operators-mean)*ground[:,None]
        coefficients=vectors.T @ rhs
        low_shells=[]
        for index in range(1,len(values)):
            weight=float(np.sum(coefficients[index]**2))
            if low_shells and abs(gaps[index]-low_shells[-1]['gap'])<1e-8:
                low_shells[-1]['degeneracy']+=1;low_shells[-1]['weight']+=weight
                low_shells[-1]['electric_boundary_weight']+=float(boundary @ (vectors[:,index]**2))
            else:
                low_shells.append(dict(gap=float(gaps[index]),degeneracy=1,weight=weight,
                    electric_boundary_weight=float(boundary @ (vectors[:,index]**2))))
        static=np.zeros((operators.shape[1],operators.shape[1])); solve_residuals=[]
        for index in range(operators.shape[1]):
            solution,info=cg(A,rhs[:,index],rtol=1e-11,atol=1e-13,maxiter=5000)
            assert info==0
            solve_residuals.append(float(np.linalg.norm(A @ solution-rhs[:,index])))
            static[:,index]=2*rhs.T @ solution
        assert max(solve_residuals)<1e-8
        assert np.linalg.norm(static-static.T)<1e-8
        assert np.linalg.norm(static-np.eye(len(mean))*np.trace(static)/len(mean))<1e-7
        responses=[]
        for T in [.1,.4,1.,2.,5.]:
            evolved=expm_multiply(-1j*T*H,np.column_stack([ground,operators[:,0]*ground]))
            R=-2*float(np.imag(np.vdot(evolved[:,0],operators[:,0]*evolved[:,1])))
            item=dict(T=T,response=R)
            if T==.4:
                eps=1e-5
                plus=expm_multiply(-1j*T*H,np.exp(1j*eps*operators[:,0])*ground)
                minus=expm_multiply(-1j*T*H,np.exp(-1j*eps*operators[:,0])*ground)
                impulse=float((operators[:,0] @ abs(plus)**2-operators[:,0] @ abs(minus)**2)/(2*eps))
                assert abs(impulse-R)<2e-8
                item.update(independent_impulse_derivative=impulse,impulse_error=abs(impulse-R))
            responses.append(item)
        total_weight=float(np.sum(rhs**2));captured=float(np.sum(coefficients[1:]**2))
        first=next((shell for shell in low_shells if shell['weight']>1e-11),None)
        rows.append(dict(channel=name,probe_components=len(mean),probe_mean=mean.tolist(),
                variance_per_component=np.sum(rhs**2,axis=0).tolist(),total_variance=total_weight,
                captured_low_shell_weight=captured,captured_variance_fraction=captured/total_weight,
                first_active_shell=first,low_shells=low_shells,static_susceptibility=static.tolist(),
                linear_solve_residuals=solve_residuals,responses=responses))
        print(json.dumps(dict(stage="channel",L=cutoff,channel=name,first_active_shell=first,
                              static_per_component=float(np.trace(static)/len(mean)))),flush=True)
    out=dict(status="CP12_FULL_TENSOR_SPECTRUM_PASS",N=n,Q=q,L=cutoff,dimension=len(states),
             ground_energy=float(values[0]),full_gap=float(gaps[1]),low_eigenvalues=values.tolist(),
             eigen_residuals=residuals.tolist(),ground_boundary_weight=float(boundary @ ground**2),
             lowest_shell_permutation_representation=representation,
             channels=rows,scope="Complete finite regulated CP-12, independently calibrated vertex-0 tensor probes. No continuum pole, label-gauge quotient, or physical metric identification.")
    (P/f"full_tensor_L{cutoff}.json").write_text(json.dumps(out,indent=2)+"\n")
    return out


if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument('--cutoff',type=int,default=2)
    ap.add_argument('--eigenvectors',type=int,default=28);args=ap.parse_args()
    run(args.cutoff,args.eigenvectors)
