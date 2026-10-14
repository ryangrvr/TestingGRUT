"""Independent exact combinatorial and full physical-matrix controls."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import sqrt
import hashlib, json, sys
import numpy as np
from scipy.sparse import diags, coo_matrix, eye
from scipy.sparse.linalg import eigsh, expm_multiply

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'cp12_ground_connectivity'))
sys.path.insert(0,str(P.parent/'cp12_unmasked_formation'))
from frozen_model import build, gauss
from exact_destination import transitions
from run import model
from dynamics import resources, rewire_matrix


def bridge_mask(n,g,es):
    def components(h):
        adj=[[] for _ in range(n)]
        for e,(a,b) in enumerate(es):
            if h>>e&1:adj[a].append(b);adj[b].append(a)
        seen=set();count=0
        for a in range(n):
            if a in seen:continue
            count+=1;seen.add(a);todo=[a]
            for v in todo:
                for w in adj[v]:
                    if w not in seen:seen.add(w);todo.append(w)
        return count
    c=components(g)
    return sum(1<<e for e in range(len(es)) if g>>e&1 and components(g^(1<<e))>c)


def destination_control():
    n=6;es=list(combinations(range(n),2));M=len(es);scale=__import__('math').lcm(8,10*M)
    pairs=[(0,1),(1,2),(0,2),(2,3),(3,4),(4,5)]
    g=sum(1<<es.index(e) for e in pairs)
    sea=sum(1<<(n+e) for e in range(M));field=[0]*(2*M)
    for a,b in [(0,1),(1,2),(2,0)]:
        e=es.index(tuple(sorted((a,b))));sign=1 if a<b else -1
        field[2*e]=sign;field[2*e+1]=-sign
    s=(g,sea,tuple(field));assert not np.any(gauss(n,s))
    target_pairs=[(0,3),(0,4),(0,5)]
    counts={pair:F(0) for pair in target_pairs}
    for dest,c in transitions(n,s,scale):
        for pair in target_pairs:
            f=es.index(pair)
            if dest[0]>>f&1:counts[pair]+=F(c*c,scale*scale)
    assert len(set(counts.values()))==1
    assert next(iter(counts.values()))==F(3,100*M*M)
    return dict(N=n,Q=6,prepared_comparator_only=True,flux_protected_triangle=True,
                removable_edges=3,rows=[dict(candidate=list(pair),initial_distance=d,
                edge_creation_coupling_square=str(counts[pair]))
                for pair,d in zip(target_pairs,[2,3,4])],
                scope='Exact instantaneous cross coupling only; finite-time higher-order topology dependence is not excluded. Supplied comparator, not an autonomous spatial state.')


def energy_envelope(n,q,energy_density=None):
    nu=q/n;a=F(15,8);eps=F(1,4)
    c=(sqrt(float(a*a)+8*nu)-float(a))/2
    excess=0. if energy_density is None else energy_density-q/(2*n)
    return dict(N=n,Q=q,density=nu,negative_kinetic_mass_allowance_per_vertex=c,
                neutral_energy_removable_fraction_lower=max(0.,1-c/(float(eps)*nu)),
                general_removable_per_vertex_lower=nu-(excess+c)/float(eps))


def physical_control():
    n=4;q=4;L=3;budget=8
    Hfull,full,_,_=build(n,q,L)
    es=list(combinations(range(n),2));M=len(es)
    rr=[];kk=[];nv=[];e2=[];vacuum=0
    for s in full:
        g,b,E=s
        r=sum(bool(g>>e&1) and bool(b>>(n+e)&1) and E[2*e]==E[2*e+1]==0 for e in range(M))
        k=sum(bool(g>>e&1) and bool(b>>(n+e)&1) and (E[2*e]!=0 or E[2*e+1]!=0) for e in range(M))
        v=(b&((1<<n)-1)).bit_count();f=sum(x*x for x in E)
        assert r+k+v==q and f>=2*k+v
        if v==0:
            bridges=bridge_mask(n,g,es)
            assert all(not(bridges>>e&1) or E[2*e]==E[2*e+1]==0 for e in range(M))
            vacuum+=1
        rr.append(r);kk.append(k);nv.append(v);e2.append(f)
    rr=np.asarray(rr);e2=np.asarray(e2)
    # R commutes with r on the independently enumerated physical space.
    mask=np.array([s[0] for s in full]);off=Hfull.tocoo()
    rewired=(mask[off.row]!=mask[off.col])
    assert np.all(rr[off.row[rewired]]==rr[off.col[rewired]])
    # Quadratic-fermion lower bound, including the full rewiring Laplacian.
    c=energy_envelope(n,q)['negative_kinetic_mass_allowance_per_vertex']
    positive=Hfull-diags(q/2+float(F(1,4))*(q-rr)-n*c)
    lower=float(eigsh(positive,k=1,which='SA',return_eigenvectors=False,tol=1e-11)[0])
    assert lower>=-1e-9
    sub,H,B,states,sizes,psi=model(n,q,budget)
    row=[];col=[];data=[];index={s:a for a,s in enumerate(states)}
    for i,s in enumerate(full):
        if sum(x*x for x in s[2])>budget:continue
        key,sign,size=sub.canonical(s)
        j=index.get(key)
        if sign and j is not None:row.append(i);col.append(j);data.append(sign/sqrt(size))
    T=coo_matrix((data,(row,col)),shape=(len(full),len(states))).tocsr()
    iso=T.T@T-eye(len(states),format='csr')
    kept=(e2<=budget).astype(float)
    inside=diags(kept)@Hfull@T-T@H
    outside=diags(1-kept)@Hfull@T
    gram=outside.T@outside-B.T@B
    assert np.linalg.norm(iso.data)<1e-10
    assert np.linalg.norm(inside.data)<1e-9
    assert np.linalg.norm(gram.data)<1e-9
    r,k,v,f=resources(n,q,sub,states)
    # Compare a derivative against an independent full-space commutator.
    w=expm_multiply(-.31j*H,psi)
    rate=-2*np.vdot(H@w,r*w).imag
    wf=T@w
    fullrate=-2*np.vdot(Hfull@wf,rr*wf).imag
    assert abs(rate-fullrate)<1e-10
    # Independent graph/environment partial trace checks the orbit expansion.
    from records import expansion
    ar,ac,aj,aw,_,G,C=expansion(sub,states,sizes,r)
    A=coo_matrix((aw*w[aj],(ar,ac)),shape=(G,C)).tocsr()
    ev=np.linalg.eigvalsh((A@A.conj().T).toarray())
    gd={};ed={};fr=[];fc=[]
    for g,b,E in full:
        if g not in gd:gd[g]=len(gd)
        if (b,E) not in ed:ed[b,E]=len(ed)
        fr.append(gd[g]);fc.append(ed[b,E])
    Af=coo_matrix((wf,(fr,fc)),shape=(len(gd),len(ed))).tocsr()
    evf=np.linalg.eigvalsh((Af@Af.conj().T).toarray())
    trace_difference=float(np.max(abs(ev-evf)))
    assert trace_difference<1e-11
    from evolution import Evolution
    _,_,err=Evolution(H,B,q/2,16)(psi,.31)
    reference=expm_multiply(-.31j*Hfull,T@psi)
    Fref=float(abs(reference)**2@e2/8)
    Fproj=float(abs(w)**2@f/8)
    K=2*n+sqrt(2*n*q)+2*.1*q*(M-q)/M
    energy_bound=K*(2*err+err**2)
    assert abs(Fref-Fproj)<=energy_bound+1e-12
    O=np.column_stack([r/q,k/q]);OF=np.column_stack([rr/q,np.asarray(kk)/q])
    _,back,_=Evolution(H,B,q/2,16)(O*w[:,None],.31,back=True)
    actual_counts=abs((abs(reference)**2@OF)-(abs(w)**2@O))
    count_bounds=2*err*back+err**2
    assert np.all(actual_counts<=count_bounds+1e-12)
    rate_schur=[]
    for count in [rr,np.asarray(kk)]:
        comm=Hfull@diags(count)-diags(count)@Hfull
        norm=float(np.max(np.asarray(abs(comm).sum(axis=1))))/q
        assert norm<=2+1e-12
        rate_schur.append(norm)
    dr=[];dc=[];dv=[];deleted={}
    for j,(g,b,Efields) in enumerate(full):
        for e in range(M):
            if g>>e&1 and b>>(n+e)&1 and Efields[2*e]==Efields[2*e+1]==0:
                key=(g^(1<<e),b,Efields)
                if key not in deleted:deleted[key]=len(deleted)
                dr.append(deleted[key]);dc.append(j);dv.append(1.)
    D=coo_matrix((dv,(dr,dc)),shape=(len(deleted),len(full))).tocsr()
    off=Hfull.tocoo();which=mask[off.row]!=mask[off.col]
    Rfull=coo_matrix((off.data[which],(off.row[which],off.col[which])),shape=Hfull.shape).tocsr()+diags(.1/M*rr*(M-q))
    incidence=D.T@D-diags((M-q+1)*rr)+M/.1*Rfull
    incidence_error=float(np.linalg.norm(incidence.data))
    assert incidence_error<1e-9
    Rfast=rewire_matrix(n,q,sub,states,sizes,H=H)
    Rslow=rewire_matrix(n,q,sub,states,sizes)
    fast_error=float(np.linalg.norm((Rfast-Rslow).data))
    physical_R_error=float(np.linalg.norm((Rfast-T.T@Rfull@T).data))
    assert fast_error<1e-10 and physical_R_error<1e-10
    return dict(N=n,Q=q,L=L,electric_budget=budget,full_dimension=len(full),
                invariant_dimension=len(states),vacuum_states_checked=vacuum,
                kinetic_envelope_smallest_eigenvalue=lower,
                isometry_residual=float(np.linalg.norm(iso.data)),
                Hamiltonian_intertwiner_residual=float(np.linalg.norm(inside.data)),
                outgoing_Gram_residual=float(np.linalg.norm(gram.data)),
                rate_difference=float(abs(rate-fullrate)),
                independent_graph_trace_spectrum_difference=trace_difference,
                independent_electric_energy_difference=abs(Fref-Fproj),
                electric_energy_inferred_error_upper=energy_bound,
                detector_specific_count_actual_errors=actual_counts.tolist(),
                detector_specific_count_error_bounds=count_bounds.tolist(),
                normalized_count_commutator_Schur_bounds=rate_schur,
                full_physical_deletion_rewiring_identity_residual=incidence_error,
                fast_rewiring_versus_transition_assembly_residual=fast_error,
                fast_rewiring_versus_full_physical_residual=physical_R_error)


def main():
    control=physical_control()
    out=dict(status='CP12_AUTONOMOUS_FLUX_EXACT_CONTROLS_PASS',
             exact_frozen_density_threshold='49/2',
             energy_envelopes=[energy_envelope(200,int(200*nu)) for nu in [1,2,3,10,24.5,30,50]],
             destination_control=destination_control(),physical_matrix_control=control,
             qualification='Full-law algebraic proofs are in DERIVATION.md; matrix evaluations use ordinary floating point.')
    (P/'EXACT_CONTROLS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)


if __name__=='__main__':main()
