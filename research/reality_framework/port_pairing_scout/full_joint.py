#!/usr/bin/env python3
"""SCOUT-RG02 full, non-adiabatic graph+fermion quantum benchmark.

N vertices, two persistent fermionic port modes per vertex, exactly N
fermions, every perfect pairing of ports as fluctuating graph. Graph
reconnection acts as identity on the persistent port-fermion Fock space.
N=3,4 signs +/- and N=5 favorable sign tested. Never a lattice input.

OPENBLAS_NUM_THREADS=1 python full_joint.py    (numpy/scipy)
"""
import itertools, json, math
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix, diags, eye
from scipy.sparse.linalg import eigsh

J=0.1
def matchings(vertices):
    if not vertices:
        yield ();return
    a=vertices[0]
    for j in range(1,len(vertices)):
        for tail in matchings(vertices[1:j]+vertices[j+1:]):
            yield ((a,vertices[j]),)+tail

def canonical(p):
    return tuple(sorted(tuple(sorted(edge)) for edge in p))

def cycle_lengths(p,n):
    parent=list(range(n))
    def root(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    for a,b in p:
        x,y=root(a//2),root(b//2)
        if x!=y:parent[x]=y
    sizes={}
    for i in range(n):sizes[root(i)]=sizes.get(root(i),0)+1
    return tuple(sorted(sizes.values()))

def fermion_basis(n):
    bits=[sum(1<<i for i in comb) for comb in
          itertools.combinations(range(2*n),n)]
    return bits,{v:i for i,v in enumerate(bits)}

def fermion_hopping(n,matching,fbits,index,relative_sign=-1.):
    # Persistent fermion modes on physical ports. Normal intra-vertex
    # hopping +1; connection hopping relative_sign (+/- 1).
    terms=[(2*i,2*i+1,1.) for i in range(n)]
    terms +=[(a,b,relative_sign) for a,b in matching]
    rows=[];cols=[];vals=[]
    for j,bits in enumerate(fbits):
        for a,b,hopping in terms:
            if ((bits>>a)&1)==((bits>>b)&1):continue
            nxt=bits^(1<<a)^(1<<b)
            low,high=min(a,b),max(a,b)
            between=((1<<high)-1)^((1<<(low+1))-1)
            fermionic_sign=-1 if (bits&between).bit_count()%2 else 1
            rows.append(index[nxt]);cols.append(j)
            vals.append(hopping*fermionic_sign)
    return coo_matrix((vals,(rows,cols)),
               shape=(len(fbits),len(fbits))).tocsr()

def assemble(n,relative_sign=-1.,gamma=J):
    ms=list(matchings(tuple(range(2*n))))
    lookup={canonical(p):i for i,p in enumerate(ms)}
    fb,fbindex=fermion_basis(n)
    nm=len(ms);nf=len(fb);dimension=nm*nf
    rows=[];cols=[];vals=[]
    for i,p in enumerate(ms):
        # All two-bond alternative pairings; kinetic positive Laplacian.
        for a,b in itertools.combinations(range(n),2):
            u,v=p[a];x,y=p[b]
            for replacement in (((u,x),(v,y)),((u,y),(v,x))):
                new=list(p);new[a],new[b]=replacement
                destination=lookup[canonical(new)]
                for k in range(nf):
                    rows.append(i*nf+k);cols.append(destination*nf+k)
                    vals.append(-gamma/n)
        F=fermion_hopping(n,p,fb,fbindex,relative_sign).tocoo()
        rows.extend((F.row+i*nf).tolist())
        cols.extend((F.col+i*nf).tolist())
        vals.extend(F.data.tolist())
    H=coo_matrix((vals,(rows,cols)),shape=(dimension,dimension)).tocsr()
    # Exactly common reference shift: for sign -1 the ground of each
    # filled-band block >= -4N/pi. Graph kinetic is positive.
    H=H+diags(np.full(dimension,gamma*(n-1)+4*n/math.pi))
    assert np.linalg.norm((H-H.getH()).data)<1e-10
    return H,ms,fb,fbindex

def fermion_permutation_parity(occupied, mapping):
    ordered=[mapping[i] for i in occupied]
    inversions=sum(ordered[i]>ordered[j] for i in range(len(ordered))
                   for j in range(i+1,len(ordered)))
    return -1 if inversions%2 else 1

def signed_vertex_swap_commutator(n,H,ms,fb,index):
    """True CAR wedge representation, not unsigned computational relabelling."""
    assert n<=4
    lookup={canonical(p):i for i,p in enumerate(ms)}
    vertex=list(range(n))
    vertex[0],vertex[1]=vertex[1],vertex[0]
    mapping=[2*vertex[i//2]+i%2 for i in range(2*n)]
    rows=[];cols=[];data=[]
    nf=len(fb)
    for a,p in enumerate(ms):
        j=lookup[canonical(tuple((mapping[u],mapping[v]) for u,v in p))]
        for k,bits in enumerate(fb):
            occ=[x for x in range(2*n) if bits>>x&1]
            targetbits=sum(1<<mapping[x] for x in occ)
            rows.append(j*nf+index[targetbits]);cols.append(a*nf+k)
            data.append(fermion_permutation_parity(occ,mapping))
    T=coo_matrix((data,(rows,cols)),shape=H.shape).tocsr()
    assert np.linalg.norm((T@T.getH()-eye(H.shape[0])).data)<1e-10
    return float(np.linalg.norm((H@T-T@H).data))

def run(n,relative_sign):
    H,ms,fb,index=assemble(n,relative_sign)
    E,V=eigsh(H,k=4,which="SA",v0=np.ones(H.shape[0]),
              tol=1e-9,maxiter=25000)
    order=np.argsort(E);E=E[order];g=V[:,order[0]]
    residual=float(np.linalg.norm(H@g-E[0]*g))
    assert residual<1e-8
    nf=len(fb)
    blocks=g.reshape(-1,nf)
    probabilities=np.sum(abs(blocks)**2,axis=1)
    lengths=[cycle_lengths(p,n) for p in ms]
    connected=float(sum(probabilities[i] for i,p in enumerate(lengths)
                         if p==(n,)))
    singular=np.linalg.svd(blocks,compute_uv=False)
    eigen=np.abs(singular)**2;eigen/=sum(eigen)
    entropy=float(-sum(p*math.log(p) for p in eigen if p>1e-14))
    purity=float(sum(eigen**2))
    comm=(signed_vertex_swap_commutator(n,H,ms,fb,index)
          if n<=4 else None)
    if comm is not None:assert comm<1e-10
    result=dict(N=n,relative_sign=relative_sign,dimension=H.shape[0],
        graph_pairings=len(ms),fermion_states=nf,
        lowest_energies=E.tolist(),
        connected_probability=connected,
        largest_component_fraction=float(sum(
            probabilities[i]*max(p)/n for i,p in enumerate(lengths))),
        geometry_matter_entanglement_nats=entropy,
        reduced_geometry_purity=purity,
        schmidt_rank=int(sum(eigen>1e-10)),
        vertex_swap_commutator=comm,
        eigen_residual=residual)
    print(json.dumps(result),flush=True)
    return result

def main():
    rows=[run(n,sign) for n,sign in (
        (3,-1.),(3,1.),(4,-1.),(4,1.),(5,-1.))]
    by={(r['N'],r['relative_sign']):r for r in rows}
    assert abs(by[(5,-1.)]['connected_probability']
               -0.96894089748856)<1e-9
    assert by[(3,1.)]['connected_probability']<1e-4
    assert by[(4,1.)]['connected_probability']<1e-6
    assert .39<by[(5,-1.)]['reduced_geometry_purity']<.42
    Path(__file__).with_name("full_joint_results.json").write_text(
        json.dumps(rows,indent=2)+"\n")
    print("SCOUT_RG02_FULL_JOINT_FERMIONS_PASS")

if __name__=="__main__":main()
