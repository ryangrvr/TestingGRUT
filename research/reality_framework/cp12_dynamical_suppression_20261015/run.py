"""Full frozen-H action-ratio and mechanism tests on neutral-origin branches."""
from pathlib import Path
from math import sqrt
import argparse,gc,json,sys,time
import numpy as np
from boundary import quotient_model,cache_root,RF,load_model
from mechanisms import statistics,action_operators
from action_integrals import ActionEvolution,ALLOWANCE
from hopping_bound import hopping_bound
from linear_propagator import MatrixFreePropagator
from resources import include_binary_record_entropy
sys.path.insert(0,str(RF/'cp12_shortcut_dynamics_20261015'))
from integrals import Propagator
sys.path.insert(0,str(RF/'charged_density_dynamics_20261015'))
from stream import probability_interval


def ratio(I,J):
    if J[0]<=0:return [0.,None]
    return [I[0]/J[1],I[1]/J[0]]


def main():
    a=argparse.ArgumentParser();a.add_argument('--n',type=int,required=True);a.add_argument('--q',type=int)
    a.add_argument('--budget',type=int,default=12);a.add_argument('--model')
    a.add_argument('--waits',default='1');a.add_argument('--lags',default='.125,.25,.5')
    a.add_argument('--dt',type=float,default=1/32);a.add_argument('--pilot-only',action='store_true')
    a.add_argument('--out',required=True);args=a.parse_args();start=time.monotonic()
    n=args.n;q=args.q or 2*n;name=args.model or f'N{n}_B{args.budget}'
    hopbound,hop_squares,hop_minors=hopping_bound(n,q)
    m=load_model(n,q,args.budget,args.model);H,B,psi=(m[k] for k in ['H','B','psi'])
    print(json.dumps(dict(stage='loaded_model',model=name,dimension=H.shape[0])),flush=True)
    # Build observables once, then release millions of Python basis labels
    # before sparse propagation. Labels remain reproducibly in the cache.
    from operators import flow
    if not args.pilot_only:
        ds,db=statistics(m,n,q,name);P=ds['charged_event'].astype(bool);actions=action_operators(m,n)
    else:
        P=np.array([flow(n,s)[0].bit_count()>=int(np.ceil(q/4)) and
            (s[1]&((1<<n)-1)).bit_count()>=int(np.ceil(n/4)) for s in m['states']])
    for key in ['states','outer_states','sub','sizes','V','BV']:m.pop(key,None)
    gc.collect()
    print(json.dumps(dict(stage='prepared_observables',model=name)),flush=True)
    waits=sorted(float(x) for x in args.waits.split(','));lags=sorted(float(x) for x in args.lags.split(','))
    propagator=MatrixFreePropagator if n>=6 else Propagator
    neutral=propagator(H,{'electric':B},q/2,dt=args.dt).run(psi,waits)
    rows=[]
    for wait,(w,er) in neutral.items():
        ew=float(er['electric'][0]);X=P*w;p=float(np.linalg.norm(X)**2);pI=probability_interval(p,ew)
        de=-2*float(np.vdot(X,H@((1-P)*w)).real)
        workerr=2*sqrt(2*q*min(n,q))*(2*ew+ew*ew)+ALLOWANCE
        row=dict(wait=wait,formation_probability=p,formation_interval=pI,waiting_state_error=ew,
            record_energy_change=de,record_energy_interval=[de-workerr,de+workerr])
        print(json.dumps(dict(stage='neutral_formation',model=name,**row)),flush=True)
        if not args.pilot_only:
            results=ActionEvolution(H,B,actions,ds,db,q,q/2,dt=args.dt,
                hopping_norm_bound=hopbound,matrix_free=(n>=6)).run(X,lags,initial_error=ew)
            row['probes']=[]
            for lag,data in results.items():
                ints=data['integrals'];Iout=ints['outside']['full_law_interval'];Iin=ints['inside']['full_law_interval']
                outref=ints['dephased_outside_action']['full_law_interval'];inref=ints['dephased_inside_action']['full_law_interval']
                Z=ints['neutral_edges']['full_law_interval'];av=ints['neutral_edge_vacancies']['full_law_interval']
                availability=ratio(av,[2*Z[0],2*Z[1]])
                delta2=ratio(ints['outside_detuning_square_capped8']['full_law_interval'],outref)
                survival=ratio(data['final']['charged_event']['interval'],pI)
                probe=dict(lag=lag,joint_integrals=ints,
                    conditional_integrals={key:ratio(value['full_law_interval'],pI) for key,value in ints.items()},
                    eta=ints['outside']['value']/ints['inside']['value'],eta_interval=ratio(Iout,Iin),
                    interference_outside_ratio=ints['outside']['value']/ints['dephased_outside_action']['value'],
                    interference_outside_interval=ratio(Iout,outref),
                    interference_inside_ratio=ints['inside']['value']/ints['dephased_inside_action']['value'],
                    interference_inside_interval=ratio(Iin,inref),
                    neutral_edge_vacancy_fraction=ints['neutral_edge_vacancies']['value']/(2*ints['neutral_edges']['value']),
                    neutral_edge_vacancy_fraction_interval=availability,
                    outside_detuning_ge4_fraction_interval=ratio(ints['outside_detuning_ge4']['full_law_interval'],outref),
                    outside_detuning_ge8_fraction_interval=ratio(ints['outside_detuning_ge8']['full_law_interval'],outref),
                    outside_diagonal_detuning_RMS_capped8_interval=[sqrt(delta2[0]),sqrt(delta2[1]) if delta2[1] is not None else None],
                    charged_event_survival_conditional=data['final']['charged_event']['value']/p,
                    charged_event_survival_interval=[max(0.,survival[0]),min(1.,survival[1])] if survival[1] is not None else survival,
                    final_joint_observables=data['final'],lag_state_error=data['lag_state_error'],
                    total_state_error=data['state_error'],integrated_state_error_square_upper=data['integrated_state_error_square_upper'])
                row['probes'].append(probe)
                print(json.dumps(dict(stage='eta',model=name,wait=wait,lag=lag,eta=probe['eta'],
                    interval=probe['eta_interval'],outside=Iout,inside=Iin,
                    interference=probe['interference_outside_interval'],vacancy=availability)),flush=True)
        rows.append(row)
    result=dict(status='PASS',N=n,Q=q,model=name,base_budget=args.budget,
        dimension=H.shape[0],outside_dimension=B.shape[0],H_nnz=H.nnz,
        initial_state='Original equal coherent neutral ray over every labelled fixed-Q graph.',
        formed_state='X=P exp(-iH wait)psi0, with P={c/Q>=1/4,nv/N>=1/4}; normalized psi_formed=X/sqrt(p). The same event and full H are used at every N. This is an explicitly recorded conditional diagnostic, not autonomous state selection.',
        evolution='Full frozen H only. No outside hopping is removed.',
        eta_scope='Ratio of hopping-action-square integrals; homogeneous in X. Not a probability, transport speed, fixed-core metric or a locality certificate.',
        interference_scope='Coherent ||L psi||² divided by <psi,diag_physical(L^dagger L)psi>. The reference is diagonal in the full labelled physical electric/Fock basis, not the orbit basis; no dephasing is applied to the evolved state.',
        Pauli_scope='Fraction of empty endpoints on neutral active edges, weighted by their integrated count. Neutral edges are outside every instantaneous through-flow core. This availability statistic is not a coherent transport rate.',
        detuning_scope='Allowed-hop physical diagonal-energy differences, with the full rewiring escape diagonal. Fractions at |Delta D|>=4 and>=8 and RMS capped at8 are bounded observables; no eigenvalue gap or unbounded energy moment is inferred from norm error.',
        extension=m.get('extension'),rows=rows,dt=args.dt,
        hopping_norm_bound=hopbound,Jacobi_offdiagonal_squares=hop_squares,
        Jacobi_positive_leading_minors=hop_minors,
        sparse_implementation='Real-H LinearOperator with exact adjoint' if n>=6 else 'Explicit shifted sparse generator',
        error_method='Midpoint integration of squared actions through fourth derivative with a sixth-derivative remainder. Full outgoing residual upper sums bound state error and its integrated square; Minkowski controls full-H action integrals. Ordinary floating point plus5e-10 reporting allowance, not outward rounding.',
        seconds=time.monotonic()-start)
    Path(args.out).write_text(json.dumps(include_binary_record_entropy(result),indent=2)+'\n')


if __name__=='__main__':main()
