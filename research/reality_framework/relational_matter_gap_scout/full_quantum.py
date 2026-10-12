#!/usr/bin/env python3
"""Full finite quantum hard-core matter on a cyclic phase lattice.

N identical hard-core particles on q=2N bins, no independent graph qubits.
Fixed circumference L=N*a, where a=sqrt(g/P) was separately derived
in LAW.md. NOT a quantum operator for a variable circumference.
Requires numpy/scipy. Run python full_quantum.py.
"""
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

J,G,P,A=0.2,1.,1.,1.

def build(n):
    q=2*n;dx=A/2;L=n*A
    masks=[];gaplist=[];onsite=[]
    for pos in itertools.combinations(range(q),n):
        mask=sum(1<<int(i) for i in pos)
        masks.append(mask)
        gaps=np.diff(np.r_[pos,pos[0]+q]).astype(float)*dx
        gaplist.append(gaps)
        onsite.append(G*np.sum(1/gaps)+P*L+2*J*n)
    lookup={m:i for i,m in enumerate(masks)}
    row=list(range(len(masks)));col=list(range(len(masks)))
    data=list(onsite)
    for i,mask in enumerate(masks):
        for s in range(q):
            if mask>>s&1:
                for dest in ((s-1)%q,(s+1)%q):
                    if not(mask>>dest&1):
                        new=mask^(1<<s)^(1<<dest)
                        row.append(lookup[new]);col.append(i);data.append(-J)
    H=coo_matrix((data,(row,col)),shape=(len(masks),len(masks))).tocsr()
    assert np.linalg.norm((H-H.conj().T).data)<1e-12
    # Simultaneous rotation is an exact global symmetry. Calling it a
    # gauge redundancy would be an ADDITIONAL physical identification;
    # microscopic hopping and gap potential commute with the shift.
    permutation=[]
    for mask in masks:
        shifted=((mask<<1)&((1<<q)-1))|(mask>>(q-1))
        permutation.append(lookup[shifted])
    S=coo_matrix((np.ones(len(masks)),
        (permutation,np.arange(len(masks)))),shape=H.shape).tocsr()
    assert np.linalg.norm((H@S-S@H).data)<1e-10
    return H,masks,np.array(gaplist),S

def analyze(n):
    H,masks,gaps,S=build(n)
    k=min(16,H.shape[0]-2)
    values,vectors=eigsh(H,k=k,which="SA",tol=1e-10,
                         v0=np.ones(H.shape[0]))
    ind=np.argsort(values);E=values[ind];V=vectors[:,ind]
    g=V[:,0];prob=abs(g)**2
    residual=float(np.linalg.norm(H@g-E[0]*g))
    assert residual<1e-8
    ground_shift=float(np.linalg.norm(S@g-g))
    assert ground_shift<1e-7
    q=2*n
    # A PHYSICAL density Fourier mode, not inferred from energy rankings.
    occ=np.array([[int((mask>>s)&1) for s in range(q)]
                  for mask in masks],float)
    first_density=np.array([np.sum(x*np.exp(-2j*np.pi*np.arange(q)/q))
                            /np.sqrt(q) for x in occ])
    psi=first_density*g
    w=np.abs(V[:,1:].conj().T@psi)**2
    first_excited_shift=float(np.vdot(V[:,1],S@V[:,1]).real)
    first_excited_density_weight=float(abs(np.vdot(V[:,1],psi))**2)
    assert abs(first_excited_shift+1)<1e-7
    assert first_excited_density_weight<1e-12
    active=np.flatnonzero(w>1e-6)
    assert len(active)>=2
    density_gap=float(np.min(E[1:][active])-E[0])
    assert abs(E[3]-E[4])<1e-8
    perfect=float(sum(prob[i] for i,ds in enumerate(gaps)
                      if np.max(abs(ds-A))<1e-12))
    mean_g=float((prob[:,None]*gaps).sum()/n)
    var_g=float((prob[:,None]*gaps**2).sum()/n-mean_g**2)
    return dict(N=n,phase_bins=q,hilbert_dimension=H.shape[0],
        ground_energy=float(E[0]),ground_energy_per_particle=float(E[0]/n),
        lowest_gap=float(E[1]-E[0]),first_density_coupled_gap=density_gap,
        first_excited_shift_character=first_excited_shift,
        first_excited_density_weight=first_excited_density_weight,
        first_density_spectral_weight=float(np.sum(prob*abs(first_density)**2)),
        first_15_density_weight=float(w.sum()),
        ground_gap_variance=var_g,mean_neighbor_gap=mean_g,
        perfectly_even_occupancy_probability=perfect,
        eigen_residual=residual,ground_shift_residual=ground_shift,
        lowest_six_energies=E[:6].tolist())

def main():
    rows=[]
    for n in (3,4,5,6,7):
        r=analyze(n);rows.append(r);print(json.dumps(r),flush=True)
    path=Path(__file__).with_name("full_quantum_results.json")
    path.write_text(json.dumps(rows,indent=2)+"\n")
    print("SCOUT_RG01_FULL_QUANTUM_FIXED_SIZE_PASS")

if __name__=="__main__":
    main()
