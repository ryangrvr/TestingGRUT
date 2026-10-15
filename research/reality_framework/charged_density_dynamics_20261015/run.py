"""Execute charged-loop density, coherent currents and two-time histories."""
from __future__ import annotations
from pathlib import Path
from math import comb,log2,sqrt
import argparse,hashlib,json,sys,time
import numpy as np
from scipy.sparse import coo_matrix,eye
from scipy.sparse.linalg import expm_multiply
from formation import model
from diagnostics import arrays,event
from stream import StreamEvolution,probability_interval

P=Path(__file__).resolve().parent
ALLOWANCE=5e-10


def entropy_error(error,outcomes):
    tv=min(1.,2*error+error**2)
    if tv>=1-1/outcomes:return log2(outcomes)
    h=0. if tv==0 else -tv*log2(tv)-(1-tv)*log2(1-tv)
    return h+tv*log2(outcomes-1)


def graph_record(n,q,sub,H,states,diag,w,mask,error):
    """Exact graph/event pinching expectation from invariant amplitudes.

    Hopping changes nv by one and preserves graph. Rewiring preserves nv.
    Thus their invariant matrix elements can be separated without expanding
    every graph record into a larger physical basis.
    """
    es=sub.es;M=sub.M;pr=abs(w)**2;nv=diag[1]
    coo=H.tocoo();hopping=nv[coo.row]!=nv[coo.col]
    sameflag=mask[coo.row]==mask[coo.col]
    Kpinch=coo_matrix((coo.data[hopping&sameflag],
        (coo.row[hopping&sameflag],coo.col[hopping&sameflag])),shape=H.shape).tocsr()
    # The physical rewiring escape diagonal is retained after graph pinching,
    # including transitions folded into a diagonal quotient matrix element.
    diagonal=[]
    for g,bits,E in states:
        r=sum(bool(g>>e&1) and bool(bits>>(n+e)&1) and
              E[2*e]==E[2*e+1]==0 for e in range(M))
        diagonal.append(q/2+2*(bits&((1<<n)-1)).bit_count()+sum(x*x for x in E)/8+.1*r*(M-q)/M)
    after=float(np.dot(diagonal,pr)+np.vdot(w,Kpinch@w).real)
    before=float(np.vdot(w,H@w).real)
    # Use Phi(H)-H for the uncertainty, not an unbounded-energy norm claim.
    operator_bound=2*(sqrt(2*q*min(n,q))+.1*q*(M-q)/M)
    energy_error=operator_bound*(2*error+error**2)+ALLOWANCE
    record=[];entropy=0.
    for g in sorted({s[0] for s in states}):
        multiplicity=len(sub.group)//len(sub.graph_info(g)[1])
        group=np.array([s[0]==g for s in states])
        for flag in [False,True]:
            prob=float(pr[group&(mask==flag)].sum())
            if prob:entropy-=prob*log2(prob/multiplicity)
            record.append(dict(graph_representative=g,flag=flag,labelled_graphs=multiplicity,
                combined_probability=prob,probability_per_label=prob/multiplicity))
    return dict(energy_change=after-before,energy_change_error=energy_error,
        entropy_bits=entropy,entropy_error_bits=entropy_error(error,2*comb(M,q)),
        outcomes=record,
        interpretation='External complete graph/event recording is a diagnostic intervention. Its pure-state branch entropy and system-energy change are accounted; autonomous record formation and apparatus reset costs are not derived.')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True)
    ap.add_argument('--budget',type=int,required=True);ap.add_argument('--q',type=int)
    ap.add_argument('--times',default='0.25,0.5,1,2');ap.add_argument('--wait',type=float,default=.5)
    ap.add_argument('--lag',type=float,default=.25);ap.add_argument('--a',type=float,default=.25)
    ap.add_argument('--b',type=float,default=.25);ap.add_argument('--dt',type=float,default=1/32)
    ap.add_argument('--source-current',action='store_true');ap.add_argument('--history',action='store_true')
    ap.add_argument('--reference-builder',action='store_true')
    args=ap.parse_args();started=time.monotonic();n=args.n;q=args.q or 2*n
    builder=model
    if args.reference_builder:
        import importlib.util
        spec=importlib.util.spec_from_file_location('density_reference_builder',P.parent/'cp12_unmasked_formation/run.py')
        reference=importlib.util.module_from_spec(spec);spec.loader.exec_module(reference)
        builder=reference.model
    sub,H,B,states,sizes,psi=builder(n,q,args.budget)
    diag=arrays(n,q,states);mask=event(diag,q,n,args.a,args.b)
    readouts=sorted({0.,args.wait,args.wait+args.lag}|{float(t) for t in args.times.split(',')})
    evo=StreamEvolution(H,B,q/2,max_dt=args.dt)
    direct=evo.run(psi,readouts)
    print(json.dumps(dict(stage='direct_complete',N=n,budget=args.budget,dimension=len(states))),flush=True)
    ip=np.flatnonzero(mask);iq=np.flatnonzero(~mask)
    source=None
    if args.source_current:
        # All omitted complement destinations are covered by this B. Its
        # extra rows for omitted target states only make the bound looser.
        D=H[iq][:,iq];BD=B[:,iq]
        source=StreamEvolution(D,BD,q/2,max_dt=args.dt).run(psi[iq],readouts)
        C=H[ip][:,iq]
    rows=[]
    normC=sqrt(2*q*min(n,q)) # [R,P]=0 for circulation and nv thresholds.
    for t,(w,errors) in direct.items():
        error=float(errors[0]);prob=abs(w)**2;p=float(prob[mask].sum())
        row=dict(time=t,event_probability=p,probability_interval=probability_interval(p,error),
            state_error_upper=error,norm_error=float(abs(prob.sum()-1)),
            mean_circulation_fraction=float(diag[0]@prob/q),mean_matter_density=float(diag[1]@prob/n),
            mean_balanced_circulation_fraction=float(diag[2]@prob/q),
            mean_circulation_vertex_fraction=float(diag[3]@prob/n),
            mean_largest_circulating_component_fraction=float(diag[4]@prob/n),
            mean_circulating_cycle_rank=float(diag[5]@prob),mean_graph_degree_variance=float(diag[6]@prob),
            threshold_sensitivity=[])
        for a in [.25,.5,.75,1.]:
            for b in [.2,.25,.5]:
                m=event(diag,q,n,a,b);value=float(prob[m].sum())
                row['threshold_sensitivity'].append(dict(a=a,b=b,c_min=int(np.ceil(a*q-1e-12)),
                    nv_min=int(np.ceil(b*n-1e-12)),p=value,interval=probability_interval(value,error)))
        if source is not None:
            x,z=w[ip],w[iq];zd,ed=source[t];ed=float(ed[0])
            Cz=C@z;Czd=C@zd;Cdx=C.T.conj()@x
            flow=float(2*np.vdot(x,Cz).imag)
            source_flow=float(2*np.vdot(x,Czd).imag)
            # Full observable actions include all omitted electric directions.
            # Their target/complement flags are bounded by contraction, so no
            # unavailable outer-state labels are silently assigned.
            BPx=float(np.linalg.norm(B@(mask*w)))
            Bw=float(np.linalg.norm(B@w))
            Jaction=float(np.sqrt(np.linalg.norm(Cz)**2+np.linalg.norm(Cdx)**2+(BPx+Bw)**2))
            faction=float(np.sqrt(np.linalg.norm(Czd)**2+np.linalg.norm(BD@zd)**2))
            dxaction=float(np.sqrt(np.linalg.norm(Cdx)**2+BPx**2))
            total_error=min(normC*(2*error+error**2),2*error*Jaction+normC*error**2)+ALLOWANCE
            source_error=min(2*normC*(error+ed+error*ed),
                2*(error*faction+ed*dxaction+normC*error*ed))+ALLOWANCE
            row['current']=dict(total=flow,initial_source=source_flow,
                coherent_return=flow-source_flow,total_error_upper=total_error,
                source_error_upper=source_error,return_error_upper=total_error+source_error,
                full_current_operator_action_upper=Jaction,
                full_source_action_upper=faction,full_return_leg_action_upper=dxaction,
                complement_source_state_error=ed,
                method='Source propagated with the compressed complement generator; return equals exact total minus source, not an independently integrated quadrature.')
        rows.append(row)
        print(json.dumps(dict(stage='time',N=n,budget=args.budget,time=t,p=p,error=error,
            interval=row['probability_interval'],current=row.get('current'))),flush=True)
    # Independent finite differences test the current in the finite model.
    finite_difference=[]
    A=H-q/2*eye(H.shape[0],format='csr')
    for t in [args.wait]:
        w,_=direct[t];flow=float(2*np.vdot(mask*w,H@((~mask)*w)).imag)
        for h in [1e-3,5e-4]:
            plus=expm_multiply(-1j*h*A,w,traceA=-1j*h*complex(A.diagonal().sum()))
            minus=expm_multiply(1j*h*A,w,traceA=1j*h*complex(A.diagonal().sum()))
            fd=float((np.linalg.norm(plus[mask])**2-np.linalg.norm(minus[mask])**2)/(2*h))
            finite_difference.append(dict(time=t,delta=h,current=flow,finite_difference=fd,error=abs(fd-flow)))
    histories=None
    w,ew=direct[args.wait];ew=float(ew[0])
    if args.history:
        V=np.column_stack([mask*w,(~mask)*w])
        lag=evo.run(V,[args.lag])[args.lag]
        Z,el=lag;eh=ew+float(el[0]);ee=ew+float(el[1])
        hit=mask*Z[:,0];enter=mask*Z[:,1]
        hitprob=float(np.linalg.norm(hit)**2);enterprob=float(np.linalg.norm(enter)**2)
        overlap=complex(np.vdot(hit,enter));pwait=float(np.linalg.norm(mask*w)**2)
        pi=probability_interval(pwait,ew);ji=probability_interval(hitprob,eh)
        conditional=[ji[0]/pi[1],min(1.,ji[1]/pi[0])] if pi[0]>0 else [0.,1.]
        nextw,_=direct[args.wait+args.lag]
        incoherent=hitprob+enterprob;unmeasured=float(np.linalg.norm(mask*nextw)**2)
        assert np.linalg.norm(Z.sum(axis=1)-nextw)<1e-9
        assert abs(unmeasured-incoherent-2*overlap.real)<1e-9
        histories=dict(wait=args.wait,lag=args.lag,initial_event_probability=pwait,
            sequential_joint_probability=hitprob,sequential_joint_interval=ji,
            conditional_event_survival=hitprob/pwait if pwait else None,conditional_interval=conditional,
            initially_outside_then_inside_probability=enterprob,
            unmeasured_later_probability=unmeasured,
            initial_flag_recorded_later_probability=incoherent,
            interference_in_later_probability=2*overlap.real,
            complex_history_overlap=dict(real=overlap.real,imag=overlap.imag),
            normalized_history_coherence=abs(overlap)/sqrt(hitprob*enterprob) if hitprob*enterprob else None,
            overlap_error_upper=eh*np.linalg.norm(enter)+ee*np.linalg.norm(hit)+eh*ee+ALLOWANCE,
            interpretation='The sequential probability explicitly projects at wait. The unmeasured probability retains both history branches and their interference; event survival is not survival of an identified graph or autonomous classical history.')
    recording=graph_record(n,q,sub,H,states,diag,w,mask,ew)
    M=sub.M;CN=.1*sqrt(q*(M-q))/M
    recording['conditional_same_graph_survival_lower_bound']=max(0.,1-CN**2*args.lag**2)
    recording['same_graph_bound_scope']='After the explicit graph record, each supplied graph branch obeys the exact single-host leakage bound. This generic normalization effect also preserves random graphs; it does not select spatial geometry.'
    result=dict(status='PASS',parent='1e68b1a7d29f3f7b36350aaf7d72700175ad6182',N=n,Q=q,
        electric_square_budget=args.budget,dimension=len(states),full_orbit_dimension=int(sum(map(int,sizes))),
        basis_builder='reference Python' if args.reference_builder else 'integer C++ symmetry acceleration',
        outgoing_orbits=B.shape[0],event=dict(a=args.a,b=args.b,c_min=int(np.ceil(args.a*q-1e-12)),
            nv_min=int(np.ceil(args.b*n-1e-12)),definition='Opposite-sign half-leg fields on edges belonging to directed cycles.'),
        source_current_computed=args.source_current,rows=rows,finite_difference_checks=finite_difference,
        two_time_histories=histories,complete_graph_event_record=recording,
        full_law_omitted_states=True,numerical_allowance=ALLOWANCE,max_dt=args.dt,
        error_scope='Analytical Duhamel upper sums evaluated in ordinary floating point; no formal outward-rounding certificate. Unbounded diagonal energy is removed from the recording-energy uncertainty operator.',
        geometry_scope='At Q=2N, K5 and all N6 hosts have saturated diameter at most two. No continuum dimension or macroscopic causal cone can be inferred.',
        seconds=time.monotonic()-started)
    path=P/f'N{n}_Q{q}_B{args.budget}_A{args.a:g}_V{args.b:g}.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(stage='saved',path=path.name,status='PASS',seconds=result['seconds'])),flush=True)


if __name__=='__main__':main()
