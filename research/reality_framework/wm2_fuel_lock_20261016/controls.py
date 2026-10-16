"""Complete small fuel graph, exact radial formation and branching controls."""
from collections import deque,defaultdict
from pathlib import Path
import json
import math
import platform
import numpy as np
import scipy
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import splu,expm_multiply
from scipy.special import gammainc
import model
import branching

P=Path(__file__).resolve().parent


def graph_check(n):
    m=model.matrices(n);states=m['states'];index=m['index'];h=m['H'];source=((0,),1,0)
    visited={source};queue=deque([source])
    while queue:
        c=queue.popleft()
        encoded=model.register_encode(n,*c)
        expected={(model.register_encode(n,*r),kind,amp) for r,kind,amp in model.neighbors(n,*c)}
        assert set(model.register_neighbors(n,encoded))==expected
        # Move the root as well: no event-label order enters the primitive law.
        perm=lambda v:(v+1)%n
        relabel=lambda s:(tuple(perm(v) for v in s[0]),
            sum(1<<perm(v) for v in range(n) if s[1]>>v&1),perm(s[2]))
        assert {(relabel(r),kind,amp) for r,kind,amp in model.neighbors(n,*c)}==set(model.neighbors(n,*relabel(c)))
        for target,kind,amp in model.neighbors(n,*c):
            assert c in [r for r,k,a in model.neighbors(n,*target)]
            change=model.diagonal(n,*target[:2])-model.diagonal(n,*c[:2])
            if kind.startswith('stored') or kind=='matter':assert change==0
            else:assert change in [-80.,-48.,48.,80.]
            if target not in visited:visited.add(target);queue.append(target)
    assert len(visited)==len(states)
    pairs=0;minimum=math.inf;rows=[];popappend_checked=0
    for w in model.words(n):
        q=model.spent(w)
        for source_event in w:
            start=(w,q,source_event);distance={};queue=deque()
            for dest,kind,amp in model.neighbors(n,*start):
                if model.mismatch(n,dest[0],dest[1]):
                    distance[dest]=1;queue.append(dest)
            virtual_min={}
            while queue:
                c=queue.popleft();length=distance[c]
                for dest,kind,amp in model.neighbors(n,*c):
                    if not model.mismatch(n,dest[0],dest[1]):
                        j=dest[2]
                        if j in w:
                            d=abs(w.index(source_event)-w.index(j))
                            assert length+1>=d+2
                            minimum=min(minimum,length+1-d);pairs+=1
                            virtual_min[j]=min(virtual_min.get(j,math.inf),length+1)
                    elif dest not in distance:
                        distance[dest]=length+1;queue.append(dest)
            # Operational metric: ALL consistent primitive operations, not
            # just matter hops on a frozen graph. It still equals the path metric.
            pd={start:0};queue=deque([start]);event_min={source_event:0}
            while queue:
                c=queue.popleft()
                for dest,kind,amp in model.neighbors(n,*c):
                    if model.mismatch(n,dest[0],dest[1]):continue
                    if dest not in pd:
                        pd[dest]=pd[c]+1;queue.append(dest)
                        event_min[dest[2]]=min(event_min.get(dest[2],math.inf),pd[dest])
            for j in w:assert event_min[j]==abs(w.index(source_event)-w.index(j))
            rows.append(dict(initial_word=list(w),source=source_event,
                minimum_virtual_length_by_receiver={str(j):int(k) for j,k in sorted(virtual_min.items())}))
    assert minimum==2
    return dict(N=n,complete_reachable_states=len(states),all_root_spent_masks_retained=True,
        independent_register_support_checked=True,energy_changes_checked=True,
        root_moving_label_covariance_checked=True,
        reverse_moves_checked=True,source_exit_pairs_checked=pairs,
        minimum_virtual_length_minus_d=minimum,all_P_operational_metric_matches=True,
        rows=rows,hermiticity_error=float(sparse.linalg.norm(h-h.T)))


def schur(n):
    m=model.matrices(n);h=m['H'];pi=m['P'];qi=m['Q'];states=m['states']
    hp=h[pi,:][:,pi].toarray();q=h[qi,:][:,qi].tocsc();b=h[qi,:][:,pi].tocsc()
    c=model.constants();bb=(b.T@b).toarray()
    off=bb-np.diag(np.diag(bb));assert abs(off).max(initial=0)<1e-14
    bnorm=math.sqrt(np.linalg.eigvalsh(bb).max());assert bnorm<=.25+1e-13
    p_states=[states[i] for i in pi];rows=[];near=[]
    for E in [-2.5,0.,2.5]:
        lu=splu(q-E*sparse.eye(len(qi),format='csc'));sol=lu.solve(b.toarray())
        sigma=-b.T@sol;effective=hp+sigma
        residual=float(np.linalg.norm((q-E*sparse.eye(len(qi)))@sol-b.toarray()))
        assert residual<1e-11
        pairs=[(i,j) for i in range(n) for j in range(n) if i!=j] if n<=4 else [(0,1),(1,0)]
        for i,j in pairs:
            target=[a for a,(w,f,p) in enumerate(p_states) if p==j]
            for d in range(1,n):
                source=[a for a,(w,f,p) in enumerate(p_states)
                    if p==i and i in w and j in w and abs(w.index(i)-w.index(j))>=d]
                if not source:continue
                norm=float(np.linalg.norm(sigma[np.ix_(target,source)],ord=2))
                bound=c['Sigma_C0']*c['Sigma_ratio']**d
                assert norm<=bound+2e-13
                if d>1:assert abs(hp[np.ix_(target,source)]).max(initial=0)==0
                rows.append(dict(E=E,source=i,receiver=j,min_distance=d,
                    source_columns=len(source),all_consistent_receiver_rows=len(target),
                    Sigma_block_norm=norm,bound=bound,solve_residual=residual))
            source=[a for a,(w,f,p) in enumerate(p_states)
                if p==i and i in w and j in w and abs(w.index(i)-w.index(j))==1]
            if source:
                minimum=float(np.linalg.svd(effective[np.ix_(target,source)],compute_uv=False).min())
                assert minimum>=c['near_packet_lower']-1e-12
                near.append(dict(E=E,source=i,receiver=j,minimum_receiver_action=minimum))
    spectrum=None
    if n<=4:
        ev,u=eigh(h.toarray());pd=len(pi);low=u[:,:pd];gram=low[pi,:]@low[pi,:].T
        angle=math.sqrt(max(0.,1-np.linalg.eigvalsh(gram).min()))
        assert ev[pd]-ev[pd-1]>=c['protected_spectral_gap_lower']-1e-10
        assert angle<=c['spectral_projector_distance_upper']+1e-9
        spectrum=dict(low_max=float(ev[pd-1]),high_min=float(ev[pd]),
            actual_gap=float(ev[pd]-ev[pd-1]),projector_distance=angle)
    return dict(N=n,dimension=h.shape[0],P_dimension=len(pi),Q_dimension=len(qi),
        PQ_norm=bnorm,PQ_squared_offdiagonal_max=float(abs(off).max(initial=0)),
        spectrum=spectrum,rows=rows,near_packet_checks=near)


def formation(n,delta=0.):
    h=model.radial(n,delta=delta);ev,u=eigh(h);c=u[0,:]
    fc=model.formation_constants(n);T=fc['observation_window'];L=fc['threshold']
    differences=ev[:,None]-ev[None,:]
    avgE=(c[:,None]*c[None,:])*np.exp(-.5j*differences*T)*np.sinc(differences*T/(2*np.pi))
    rho=u@avgE@u.T;pop=np.real(np.diag(rho))
    diagonal_pop=(u*u)@(c*c)
    averaged=float(pop[L-1:].sum());infinite=float(diagonal_pop[L-1:].sum())
    if delta==0:
        assert abs(c*c-1/n).max()<1e-11
        assert np.diff(ev).min()>=fc['spacing_lower']-1e-12
        assert averaged>=fc['finite_window_macroscopic_lower']-1e-11
        assert abs(diagonal_pop-1/n).max()<1e-11
    db=model.detuning_bound(n,delta) if delta else None
    if db:assert averaged<=db['all_time_macroscopic_upper']+1e-10
    energy=float(np.trace(h@rho).real);assert abs(energy)<1e-10
    mean_m=float(np.arange(1,n+1)@pop)
    kinetic=h-np.diag(delta*np.arange(n));work=-float(np.trace(kinetic@rho).real)
    kinetic_energy=float(np.trace(kinetic@rho).real)
    seed_free=16*(n-1);binding=16*(n-mean_m);storage_gain=(16+delta)*(mean_m-1)
    ledger_error=seed_free-binding-storage_gain-kinetic_energy
    assert abs(ledger_error)<1e-8
    multiplicities=[math.factorial(n-1)//math.factorial(n-m) for m in range(1,n+1)]
    entropy=-sum(float(p)*math.log2(float(p)) for p in pop if p>0)
    entropy+=sum(float(p)*math.log2(d) for p,d in zip(pop,multiplicities))+math.log2(n)
    reduction=None
    if n<=4:
        full=model.matrices(n,delta=delta,probe=False);w=np.zeros((len(full['states']),n))
        for i,(word,q,p) in enumerate(full['states']):w[i,len(word)-1]=1/math.sqrt(multiplicities[len(word)-1])
        reduction=float(np.linalg.norm(full['H']@w-w@h));assert reduction<1e-11
    return dict(N=n,delta=delta,formation_constants=fc,
        initial_eigenweight_max_error_uniform=float(abs(c*c-1/n).max()),
        actual_min_spacing=float(np.diff(ev).min()),
        finite_window_macroscopic_probability=averaged,
        infinite_time_macroscopic_probability=infinite,
        infinite_time_full_probability=float(diagonal_pop[-1]),
        finite_window_full_probability=float(pop[-1]),
        detuning_all_time_bound=db,mean_assembled_events=mean_m,
        initial_binding_energy=seed_free,remaining_binding_energy=binding,
        initial_storage_energy=16+delta,total_initial_auxiliary_energy=16*n+delta,
        stored_energy_gain=storage_gain,residual_kinetic_energy=kinetic_energy,
        complete_energy_ledger_error=ledger_error,
        mean_complete_record_work=work,complete_record_entropy_bits=entropy,
        exact_radial_intertwining_error=reduction,
        arithmetic_notice='Tiny detuned probabilities can be a numerical floor; the analytic all-time upper bound decides scaling.')


def dyson_tail(d,z):
    return 1. if d==0 else min(1.,math.exp(z)*gammainc(d,z))


def propagation(n):
    m=model.matrices(n);h=m['H'];states=m['states'];pi=m['P']
    initial=np.eye(h.shape[0])[:,pi];rows=[]
    for t in [.05,.25]:
        u=expm_multiply(-1j*t*h,initial)
        target=[a for a,(w,q,p) in enumerate(states) if p==1]
        for d in range(1,n):
            source=[a for a,k in enumerate(pi) if states[k][2]==0 and 1 in states[k][0]
                and abs(states[k][0].index(0)-states[k][0].index(1))>=d]
            norm=float(np.linalg.norm(u[np.ix_(target,source)],ord=2))
            bound=dyson_tail(d,3*t);assert norm<=bound+1e-11
            rows.append(dict(t=t,min_distance=d,source_columns=len(source),
                all_geometry_and_storage_receiver_rows=len(target),
                full_H_receiver_operator_norm=norm,Dyson_tail_bound=bound))
    return dict(N=n,rows=rows)


def source(n,length):
    m=model.matrices(n);h=m['H'];states=m['states'];idx=m['index']
    w=tuple(range(length));q=model.spent(w);initial=np.zeros(len(states))
    initial[idx[w,q,0]]=initial[idx[w,q,1]]=1/math.sqrt(2)
    ni=np.array([p==0 for w,q,p in states]);y0=initial*ni;rows=[]
    for t in [.005,.05,.25,1.]:
        x=expm_multiply(-1j*t*h,initial);y=expm_multiply(-1j*t*h,y0)
        xp=expm_multiply(-1j*t*h,initial+(-1j-1)*y0);xm=expm_multiply(-1j*t*h,initial+(1j-1)*y0)
        chi=np.zeros(n);fd=np.zeros(n)
        for a,(word,f,p) in enumerate(states):
            chi[p]+=2*(x[a].conjugate()*y[a]).imag
            fd[p]+=(abs(xp[a])**2-abs(xm[a])**2)/2
        error=float(abs(chi-fd).max());assert error<1e-10
        near_lower=.25*t-math.expm1(6*t)+6*t
        if near_lower>0:assert chi[1]>=near_lower-1e-11
        distance_bounds={}
        for k in w:
            bound=dyson_tail(abs(k),3*t)*dyson_tail(abs(k-1),3*t)
            assert abs(chi[k])<=bound+1e-11
            distance_bounds[str(k)]=float(bound)
        bad=float(np.vdot(x[m['Q']],x[m['Q']]).real)
        assert bad<=4*model.constants()['spectral_projector_distance_upper']**2+1e-10
        bprefix=math.ceil(length/2);steps=length-bprefix+1
        bad_prefix=[i for i,(word,f,p) in enumerate(states) if word[:bprefix]!=w[:bprefix]]
        prefix_probability=float(np.vdot(x[bad_prefix],x[bad_prefix]).real)
        tau=2*(1+.25)/math.sqrt(n);z=tau*t
        prefix_bound=min(1.,(math.exp(z)*gammainc(steps,z))**2)
        assert prefix_probability<=prefix_bound+1e-10
        rows.append(dict(t=t,chi_all_geometry=chi.tolist(),finite_pulse_error=error,
            nearby_full_H_response_lower=near_lower,
            receiver_response_distance_bounds=distance_bounds,
            mismatch_probability=bad,bulk_prefix_length=bprefix,
            bulk_prefix_failure_probability=prefix_probability,
            full_H_bulk_prefix_bound=prefix_bound))
    return dict(N=n,recorded_word=list(w),assembled_events=length,variance=.25,
        preparation='Recorded formed component and injected neutral probe; all later storage/geometry states retained',rows=rows)


def main():
    graph=[graph_check(n) for n in [2,3,4,5]]
    (P/'CONFIGURATION.json').write_text(json.dumps(dict(status='PASS',rows=graph),indent=2)+'\n')
    print('Complete storage-configuration and virtual-word checks passed',flush=True)
    data=dict(status='PASS',environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        constants=model.constants(),Schur=[schur(n) for n in [2,3,4,5]],
        formation=[formation(n) for n in [2,3,4,8,16,32,64,128]],
        detuned=[formation(n,delta=d) for n in [16,32,64,128] for d in [.01,.1,.25,1.]],
        detuning_large_N_bounds=[model.detuning_bound(n,d) for n in [128,512,2048,8192,32768] for d in [.01,.1,.25,1.]],
        full_propagation=[propagation(n) for n in [2,3,4]],
        sources=[source(n,m) for n,m in [(3,3),(4,4),(4,2),(4,3),(5,3)]],
        branching=branching.run(),
        arithmetic='Analytic all-size bounds, complete finite graphs, and ordinary floating-point matrices; not outward-rounded certificates.',
        decision='Resonant macroscopic formation and conditional protected transport proved. Robust phase, irreversible lock, selected dimension and universal matter NOT established.')
    (P/'RESULTS.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status='PASS',Schur=len(data['Schur']),formation=len(data['formation']),detuned=len(data['detuned']))),flush=True)


if __name__=='__main__':main()
