"""Exact graph/word analysis precedes small matrix and formation controls."""
from __future__ import annotations
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
import json
import math
import platform
import heapq
import numpy as np
import scipy
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh,splu,expm_multiply
from scipy import sparse
from scipy.special import gammainc
import model

P=Path(__file__).resolve().parent


def register_check(n):
    states,index=model.basis(n);checks=0
    assert not list(model.register_neighbors(n,(frozenset(),None,frozenset(range(n)),None)))
    for w,p in states:
        expected={(model.register_encode(n,*target),kind)
            for target,kind in model.neighbors(n,w,p)}
        actual=set(model.register_neighbors(n,model.register_encode(n,w,p)))
        assert expected==actual
        for target,kind in actual:
            assert (model.register_encode(n,w,p),kind) in set(model.register_neighbors(n,target))
        perm=lambda a:(a+1)%n
        renamed=(tuple(perm(a) for a in w),perm(p))
        renamed_neighbors={(tuple(perm(a) for a in target[0]),perm(target[1]),kind)
            for target,kind in model.neighbors(n,w,p)}
        direct_renamed={(target[0],target[1],kind)
            for target,kind in model.neighbors(n,*renamed)}
        assert renamed_neighbors==direct_renamed
        checks+=1
    # No head and no edges: every assembly transition vanishes. This is a
    # separate input check, not the seeded initial state used in formation.
    return dict(N=n,states=checks,finite_register_and_tuple_neighbors_agree=True,
        reverse_primitive_checks=True,label_covariance_checked=True,
        unseeded_empty_configuration_dark=True)


def words_check(n,max_length=10):
    original=tuple(range(n));full=[w for w in model.words(n) if len(w)==n]
    pairs=0;finite_virtual=0;min_margin=None;max_constant_word_distance=0
    by_source=[]
    for source in range(n):
        d_full=model.config_distance(n,(original,source))
        d_q=model.config_distance(n,(original[:-1],source),virtual=True)
        for w in full:
            for receiver in range(n):
                d_op=abs(source-receiver)  # positions in this formed original word
                # Compare to *every* final fully assembled ordering, not only original.
                full_length=d_full.get((w,receiver),math.inf)
                assert full_length>=d_op
                vlength=2+d_q.get((w[:-1],receiver),math.inf)
                if math.isfinite(vlength):
                    assert vlength>=d_op+2
                    finite_virtual+=1
                    margin=vlength-d_op
                    min_margin=margin if min_margin is None else min(min_margin,margin)
                pairs+=1
        vector={(original[:-1],source):1};rows=[]
        for middle in range(max_length-1):
            counts=defaultdict(int)
            for (w,p),count in vector.items():
                if len(w)==n-1:
                    d=abs(source-p);assert middle>=d
                    counts[d]+=count
                    max_constant_word_distance=max(max_constant_word_distance,d)
            rows.append(dict(length=middle+2,words_by_d_op={str(k):v for k,v in sorted(counts.items())}))
            nxt=defaultdict(int)
            for c,count in vector.items():
                for target,kind in model.neighbors(n,*c,virtual=True):nxt[target]+=count
            vector=nxt
        by_source.append(dict(source=source,word_counts=rows))
    return dict(N=n,full_configuration_states=len(model.basis(n)[0]),
        protected_configurations=len(full)*n,all_final_protected_pairs_checked=pairs,
        finite_virtual_pairs=finite_virtual,minimum_virtual_length_minus_d_op=min_margin,
        max_enumerated_word_length=max_length,by_source=by_source,
        complete_Q_graph_BFS=True,all_final_geometry_orders_retained=True)


def barrier_check(n):
    original=tuple(range(n));height={original:0};queue=[(0,original)]
    while queue:
        cost,w=heapq.heappop(queue)
        if cost!=height[w]:continue
        for (v,p),kind in model.neighbors(n,w,0):
            if kind!='assembly':continue
            candidate=max(cost,n-len(v))
            if candidate<height.get(v,math.inf):
                height[v]=candidate;heapq.heappush(queue,(candidate,v))
    count=0;mins={}
    for w,k in height.items():
        for i,j in zip(w[:-1],w[1:]):
            d=abs(i-j)
            if d>1:
                assert k>=d
                mins[d]=min(mins.get(d,math.inf),k);count+=1
    assert all(d==k for d,k in mins.items())
    return dict(N=n,all_geometry_states=len(height),long_bond_checks=count,
        minimum_maximum_defect_count_by_initial_distance={str(d):k for d,k in sorted(mins.items())},
        method='Complete minimax Dijkstra over every partial and full assembly; not a selected reconnection path.')


def schur_check(n):
    m=model.matrices(n);h=m['H'];pi=m['P'];qi=m['Q'];states=m['states'];idx=m['index']
    hp=h[pi,:][:,pi].tocsr();hq=h[qi,:][:,qi].tocsc();b=h[qi,:][:,pi].tocsc()
    pindex={int(k):i for i,k in enumerate(pi)};c=model.constants(n)
    assert np.max(np.abs((b.T@b-sparse.eye(len(pi))/n).data),initial=0)<1e-14
    original=tuple(range(n));rows=[]
    for E in [-.5,0.,.5]:
        lu=splu(hq-E*sparse.eye(len(qi),format='csc'))
        for source in range(n):
            k=pindex[idx[original,source]];rhs=b[:,k].toarray().ravel();sol=lu.solve(rhs)
            sigma=-np.asarray(b.T@sol).ravel();effective=hp[:,k].toarray().ravel()+sigma
            rc=model.constants(n,Emax=E);violations=[];byd=defaultdict(list)
            for row,physical in enumerate(pi):
                w,j=states[physical];d=abs(source-j)
                bound=rc['Sigma_prefactor']*rc['Sigma_ratio']**d
                violations.append(abs(sigma[row])-bound);byd[d].append(abs(float(sigma[row])))
            excess=max(violations);assert excess<2e-13
            near=[pindex[idx[original,j]] for j in [source-1,source+1] if 0<=j<n]
            assert min(abs(effective[j]) for j in near)>=.25-1e-13
            rows.append(dict(E=E,source=source,max_bound_excess=excess,
                linear_solve_residual=float(np.linalg.norm((hq-E*sparse.eye(len(qi)))@sol-rhs)),
                max_abs_Sigma_by_d={str(d):max(v) for d,v in sorted(byd.items())},
                near_min=float(min(abs(effective[j]) for j in near))))
    spectral=None;propagation=[]
    if n<=5:
        ev,u=eigh(h.toarray());pdim=len(pi)
        low=u[:,:pdim];gram=low[pi,:]@low[pi,:].T
        angle=float(math.sqrt(max(0.,1-np.linalg.eigvalsh(gram).min())))
        assert ev[pdim]-ev[pdim-1]>=c['protected_spectral_gap_lower']-1e-11
        assert angle<=c['spectral_projector_distance_upper']+1e-10
        spectral=dict(low_max=float(ev[pdim-1]),high_min=float(ev[pdim]),
            actual_gap=float(ev[pdim]-ev[pdim-1]),projector_distance=angle)
        receiver=[i for i,(w,p) in enumerate(states) if p==1]
        for time in [.1,.25,1.]:
            evolution=(u*np.exp(-1j*ev*time))@u.T
            for distance in range(1,n):
                source=[i for i,(w,p) in enumerate(states) if len(w)==n and p==0 and w.index(1)>=distance]
                norm=float(np.linalg.norm(evolution[np.ix_(receiver,source)],ord=2))
                x=(2+.5)*time; bound=float(math.exp(x)*gammainc(distance,x))
                assert norm<=bound+1e-11
                propagation.append(dict(t=time,min_distance=distance,
                    full_receiver_operator_norm=norm,Dyson_tail_bound=bound))
    # Receiver projectors sum ALL final full geometries. Source columns
    # include arbitrary coherent geometry packets with separation >=d.
    full_sigma=-b.T@splu(hq).solve(b.toarray())
    p_states=[states[k] for k in pi];blocks=[];near_packets=[]
    pairs=[(i,j) for i in range(n) for j in range(n) if i!=j] if n<=5 else [(0,1)]
    for i,j in pairs:
        receiver=[k for k,(w,p) in enumerate(p_states) if p==j]
        for d in range(1,n):
            source=[k for k,(w,p) in enumerate(p_states)
                if p==i and abs(w.index(i)-w.index(j))>=d]
            if not source:continue
            norm=float(np.linalg.norm(full_sigma[np.ix_(receiver,source)],ord=2))
            bound=model.constants(n,Emax=0)['Sigma_prefactor']*model.constants(n,Emax=0)['Sigma_ratio']**d
            assert norm<=bound+2e-13
            blocks.append(dict(source_event=i,receiver_event=j,min_distance=d,
                source_columns=len(source),all_geometry_receiver_rows=len(receiver),
                operator_norm=norm,bound=bound))
        source=[k for k,(w,p) in enumerate(p_states) if p==i and abs(w.index(i)-w.index(j))==1]
        if source:
            mat=(hp.toarray()+full_sigma)[np.ix_(receiver,source)]
            minimum=float(np.linalg.svd(mat,compute_uv=False).min())
            lower=.25-model.constants(n,Emax=0)['Sigma_prefactor']
            assert minimum>=lower-1e-12
            near_packets.append(dict(source_event=i,receiver_event=j,
                minimum_receiver_action=minimum,uniform_lower=lower))
    return dict(N=n,dimension=h.shape[0],protected_dimension=len(pi),
        Q_dimension=len(qi),full_matrix_hermiticity_error=float(sparse.linalg.norm(h-h.T)),
        constants=c,spectral=spectral,rows=rows,all_geometry_receiver_blocks=blocks,
        arbitrary_near_geometry_packets=near_packets,full_propagation_blocks=propagation)


def source_check(n):
    m=model.matrices(n);h=m['H'];idx=m['index'];states=m['states'];original=tuple(range(n))
    initial=np.zeros(len(states));initial[idx[original,0]]=initial[idx[original,1]]=1/math.sqrt(2)
    number=np.array([p==0 for w,p in states],float);y0=number*initial
    assert abs(np.vdot(initial,y0).real-.5)<1e-14
    rows=[]
    for t in [.01,.1,.25,1.]:
        x=expm_multiply(-1j*t*h,initial);y=expm_multiply(-1j*t*h,y0)
        chi=np.zeros(n)
        for a,(w,p) in enumerate(states):chi[p]+=2*(x[a].conjugate()*y[a]).imag
        xp=expm_multiply(-1j*t*h,initial+(-1j-1)*y0)
        xm=expm_multiply(-1j*t*h,initial+(1j-1)*y0)
        fd=np.zeros(n)
        for a,(w,p) in enumerate(states):fd[p]+=(abs(xp[a])**2-abs(xm[a])**2)/2
        assert max(abs(fd-chi))<1e-10
        defects=float(np.vdot(x[m['Q']],x[m['Q']]).real)
        assert defects<=model.constants(n)['bare_defect_probability_all_times_upper']+1e-10
        original_indices=[i for i,(w,p) in enumerate(states) if w==original]
        retention=float(np.vdot(x[original_indices],x[original_indices]).real)
        assert retention>=max(0.,1-t*t/n)-1e-10
        rows.append(dict(t=t,chi=chi.tolist(),finite_pulse_error=float(max(abs(fd-chi))),
            defect_probability=defects,original_geometry_probability=retention,
            formation_preparation=False))
    return dict(N=n,source_variance=.25,
        preparation='A recorded, fully assembled branch with probe (|root>+|next>)/sqrt(2); NOT claimed produced with appreciable original-neutral probability.',rows=rows)


def formation(n,U=16.,eps=1.):
    # Exact symmetry-reduced invariant radial Hamiltonian, not a regulator.
    k=np.arange(n-1,-1,-1);b=-eps*np.sqrt(np.arange(n-1,0,-1)/n)
    h=np.diag(U*k)+np.diag(b,1)+np.diag(b,-1)
    ev,u=eigh(h);initial=np.zeros(n);initial[0]=1
    coeff=u.T@initial
    times=np.linspace(0,4*np.pi/U,257)
    waves=(u@(coeff[:,None]*np.exp(-1j*ev[:,None]*times))).T
    probs=abs(waves[:,-1])**2;bound=model.formation_bound(n,U,eps)
    assert probs.max()<=bound['full_assembly_probability_upper_all_times']+1e-12
    extended=np.sum(abs(waves[:,bound['extensive_length_threshold']-1:])**2,axis=1)
    assert extended.max()<=bound['extensive_assembly_probability_upper_all_times']+1e-12
    pick=int(np.argmax(probs));v=waves[pick];p=abs(v)**2
    multiplicities=np.array([math.factorial(n-1)//math.factorial(n-m) for m in range(1,n+1)],dtype=object)
    positive=p>0
    entropy=float(-sum(float(p[i])*math.log2(float(p[i])) for i in range(n) if positive[i])
        +sum(float(p[i])*math.log2(int(multiplicities[i])) for i in range(n)))
    kinetic=np.diag(b,1)+np.diag(b,-1)
    work=-float(np.vdot(v,kinetic@v).real)
    energy=float(np.vdot(v,h@v).real)
    assert abs(energy-U*(n-1))<1e-9
    assert abs(work)<=2*eps+1e-10
    exact_reduction_error=None
    if n<=5:
        full=model.matrices(n,U,eps,probe=False);s=full['states']
        w=np.zeros((len(s),n))
        for i,(word,_) in enumerate(s):w[i,len(word)-1]=1/math.sqrt(int(multiplicities[len(word)-1]))
        exact_reduction_error=float(np.linalg.norm(full['H']@w-w@h))
        assert exact_reduction_error<1e-11
    return dict(**bound,sample_time_window=[0.,float(times[-1])],sample_points=len(times),
        sampled_max_probability=float(probs[pick]),time_of_sampled_max=float(times[pick]),
        sampled_max_extensive_probability=float(extended.max()),
        probability_notice='Tiny sampled probabilities can be an unresolved double-precision floor; the all-time analytic bound decides the formation gate.',
        complete_geometry_record_work=work,complete_geometry_record_entropy_bits_fixed_root=entropy,
        extra_uniform_root_record_entropy_bits=math.log2(n),
        auxiliary_energy_check=energy,exact_radial_intertwining_error=exact_reduction_error)


def reconnection_control(n):
    original=tuple(range(n));j=n//2
    target=(0,)+tuple(range(j,n))+tuple(range(1,j))
    old=set(zip(original[:-1],original[1:]));new=set(zip(target[:-1],target[1:]))
    removed=sorted(old-new);added=sorted(new-old)
    assert len(removed)==len(added)==2 and (0,j) in new
    # A finite-support extra term can perform this block rotation, then
    # one current-edge matter hop and the inverse restore the old graph.
    return dict(N=n,reference_S1_distance=j,new_word_length=3,
        removed_bonds=removed,added_bonds=added,
        changed_head=[original[-1],target[-1]],
        term_qubit_support=6,
        status='EXCLUDED from S1 by its end-defect grammar; adding this move collapses the original protected event metric and requires recomputing d_op.')


def main():
    # Exact register and graph tests run first, before Schur or dynamics.
    registers=[register_check(n) for n in range(2,7)]
    words=[words_check(n) for n in range(2,8)]
    partial=dict(register_controls=registers,word_controls=words,
        activation_barriers=[barrier_check(n) for n in range(3,8)])
    (P/'EXACT_CONFIGURATION.json').write_text(json.dumps(partial,indent=2)+'\n')
    print('Exact configuration and word checks passed',flush=True)
    r=dict(status='PASS',environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        constants=model.constants(2),
        Schur=[schur_check(n) for n in range(2,7)],
        source=[source_check(n) for n in [3,4,5]],
        formation=[formation(n) for n in [2,3,4,5,6,8,12,16,32,64,128]],
        formation_robustness=[formation(16,U=u) for u in [4.,8.,16.,32.]],
        bad_reconnection=[reconnection_control(n) for n in [8,16,32,64,128]],
        acceptance='S1 has generated event adjacency and protected conditional transport. Its one-dimensional class is supplied by the grammar and appreciable closed-system formation fails. Full WM2 NOT established.',
        arithmetic='Exact combinatorial graph/word enumeration; analytic all-N proofs; ordinary floating-point finite matrices and radial checks, not outward-rounded certificates.')
    (P/'RESULTS.json').write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=r['status'],Schur_controls=len(r['Schur']),source_controls=len(r['source']),formation_controls=len(r['formation']))),flush=True)


if __name__=='__main__':main()
