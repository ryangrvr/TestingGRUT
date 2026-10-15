"""Full CP-12 dynamics and moving-core hopping comparison, original neutral ray."""
from pathlib import Path
import argparse,hashlib,json,sys,time
import numpy as np
from scipy.sparse import vstack
from operators import quotient_model,RF
from integrals import Propagator
sys.path.insert(0,str(RF/'charged_density_dynamics_20261015'))
from diagnostics import arrays,event
from stream import probability_interval

ALLOWANCE=5e-10


def interval(value,radius):return [float(value-radius-ALLOWANCE),float(value+radius+ALLOWANCE)]


def main():
    a=argparse.ArgumentParser();a.add_argument('--n',type=int,required=True)
    a.add_argument('--q',type=int);a.add_argument('--budget',type=int,required=True)
    a.add_argument('--times',default='.25,.5,1,2,3');a.add_argument('--wait',type=float,default=2.)
    a.add_argument('--lags',default='.125,.25,.5');a.add_argument('--dt',type=float,default=1/32)
    a.add_argument('--out',required=True);args=a.parse_args();started=time.monotonic()
    n=args.n;q=args.q or 2*n;budget=args.budget
    m=quotient_model(n,q,budget);H,V,B,BV=(m[k] for k in ['H','V','B','BV'])
    Hc=H-V;Bc=B-BV;psi=m['psi'];diag=arrays(n,q,m['states']);P=event(diag,q,n)
    # n_v changes by one on every hopping transition, and never under R.
    coo=H.tocoo();keep=diag[1][coo.row]!=diag[1][coo.col]
    from scipy.sparse import coo_matrix
    K=coo_matrix((coo.data[keep],(coo.row[keep],coo.col[keep])),shape=H.shape).tocsr()
    times=sorted({float(t) for t in args.times.split(',')}|{args.wait})
    lags=sorted(float(t) for t in args.lags.split(','))
    E=Propagator(H,{'electric':B},q/2,dt=args.dt)
    snapshots=E.run(psi,times);trajectory=[]
    for t,(w,er) in snapshots.items():
        ew=float(er['electric'][0]);p=float(np.linalg.norm(P*w)**2)
        vnorm=float(np.hypot(np.linalg.norm(V@w),np.linalg.norm(BV@w)))
        delta=float(np.linalg.norm(w-psi))
        trajectory.append(dict(time=t,event_probability=p,event_interval=probability_interval(p,ew),
            full_state_error=ew,moving_core_only_event_probability=0.,
            full_vs_moving_core_state_distance=delta,
            state_distance_interval=[max(0.,delta-ew-ALLOWANCE),min(2.,delta+ew+ALLOWANCE)],
            outside_hopping_action_norm=vnorm,outside_hopping_action_norm_interval=[max(0.,vnorm-2*q*ew-ALLOWANCE),vnorm+2*q*ew+ALLOWANCE],
            neutral_origin_exact_Duhamel_upper=t*np.sqrt(2*q)))
    w,ewrap=snapshots[args.wait];ew=float(ewrap['electric'][0]);p=float(np.linalg.norm(P*w)**2)
    columns=np.column_stack([w,P*w]);wait_norms=np.linalg.norm(columns,axis=0)
    Ef=Propagator(H,{'electric':B},q/2,dt=args.dt)
    Ec=Propagator(Hc,{'electric':Bc,'outside_hopping':vstack([V,BV],format='csr'),
        'full_vs_retained_core':vstack([V,B],format='csr')},q/2,dt=args.dt)
    full=Ef.run(columns,lags);core=Ec.run(columns,lags)
    comparisons=[]
    for lag in lags:
        x,ex=full[lag];y,ey=core[lag]
        for col,label in enumerate(['unfiltered_waiting_state','initial_event_branch']):
            ef=float(ex['electric'][col]);ec=float(ey['electric'][col]);norm=float(wait_norms[col])
            dist=float(np.linalg.norm(x[:,col]-y[:,col]));radius=2*ew+ef+ec
            overlap=abs(np.vdot(x[:,col],y[:,col]))
            raydist=float(np.sqrt(max(0.,np.linalg.norm(x[:,col])**2+np.linalg.norm(y[:,col])**2-2*overlap)))
            Iv=float(ey['outside_hopping'][col]);driver=float(ey['full_vs_retained_core'][col])
            # Exact V is bounded by 2Q, from the physical-basis row-sum bound.
            # The initial waiting-state error is common to the two paths.
            pure_duhamel=Iv+2*q*lag*(ew+ec)
            driver_duhamel=driver+ec+2*ew
            direct=[max(0.,dist-radius-ALLOWANCE),dist+radius+ALLOWANCE]
            comparisons.append(dict(preparation=label,lag=lag,initial_norm=norm,
                full_lag_electric_error=ef,core_lag_electric_error=ec,
                state_distance=dist,state_distance_interval=direct,
                phase_optimized_state_distance=raydist,
                phase_optimized_state_distance_interval=[max(0.,raydist-radius-ALLOWANCE),raydist+radius+ALLOWANCE],
                normalized_state_distance_interval=[a/max(norm+ew,1e-300) if i==0 else
                    min(2.,a/max(norm-ew,1e-300)) for i,a in enumerate(direct)],
                outside_hopping_integral_upper_sum=Iv,
                exact_full_vs_exact_core_Duhamel_upper=pure_duhamel,
                exact_full_vs_exact_core_combined_driver_upper=driver_duhamel,
                full_final_event_probability=float(np.linalg.norm(P*x[:,col])**2),
                full_final_event_interval=probability_interval(float(np.linalg.norm(P*x[:,col])**2),ew+ef),
                core_final_event_probability=float(np.linalg.norm(P*y[:,col])**2),
                core_final_event_interval=probability_interval(float(np.linalg.norm(P*y[:,col])**2),ew+ec)))
    action=[]
    for col,label in enumerate(['unfiltered_waiting_state','initial_event_branch']):
        z=columns[:,col];v=float(np.hypot(np.linalg.norm(V@z),np.linalg.norm(BV@z)))
        k=float(np.hypot(np.linalg.norm(K@z),np.linalg.norm(B@z)))
        kc=float(np.hypot(np.linalg.norm((K-V)@z),np.linalg.norm((B-BV)@z)))
        action.append(dict(preparation=label,outside_action_norm=v,
            outside_action_norm_interval=[max(0.,v-2*q*ew-ALLOWANCE),v+2*q*ew+ALLOWANCE],
            total_hopping_action_norm=k,inside_hopping_action_norm=kc,
            outside_action_norm_per_initial_norm=v/wait_norms[col]))
    z=P*w;r=(1-P)*w
    de=-2*float(np.vdot(z,H@r).real)
    bk=np.sqrt(2*q*min(n,q));energy_error=2*bk*(2*ew+ew*ew)+ALLOWANCE
    tv=min(1.,2*ew+ew*ew)
    binary=lambda x:0. if x<=0 or x>=1 else -x*np.log2(x)-(1-x)*np.log2(1-x)
    result=dict(status='PASS',N=n,Q=q,electric_budget=budget,dimension=H.shape[0],
        outgoing_dimension=B.shape[0],H_nonzeros=H.nnz,V_nonzeros=V.nnz,
        construction_seconds=m['construction_seconds'],model_checks=m['checks'],
        operator='chi_e is the directed-cycle indicator. K_inside=sum_e(chi_e K_e+K_e chi_e)/2; V=K-K_inside; H_core=H-V. All mass, electric and graph-rewiring terms are retained.',
        comparison_scope='H_core is a diagnostic intervention, not the frozen law or an emergent effective Hamiltonian. It freezes the neutral ray, so its neutral-origin divergence includes inhibited formation and is not by itself a transport test.',
        initial_state='Original equal coherent neutral ray over ALL labelled Q-edge graphs.',
        primary_event='circulating_edges/Q >= 1/4 and n_v/N >= 1/4',
        trajectory=trajectory,wait=args.wait,wait_state_error=ew,
        initial_event_probability=p,initial_event_interval=probability_interval(p,ew),
        initial_event_record_energy_change=de,initial_event_record_energy_interval=interval(de,energy_error),
        initial_event_record_entropy_bits=binary(p),initial_event_record_entropy_error_bits=binary(min(tv,.5)),
        wait_action_norms=action,post_formation_comparisons=comparisons,
        V_norm_bound=2*q,
        error_method='All outgoing electric rows; analytic midpoint norm upper sums through order five plus sixth-derivative remainder. Duhamel bounds include comparison-law leakage and waiting-state error. Ordinary floating point with separate 5e-10 reporting allowance, not outward rounding.',
        frozen_sources={name:hashlib.sha256((RF/'cp12_ground_connectivity'/name).read_bytes()).hexdigest()
            for name in ['frozen_model.py','exact_destination.py']},
        seconds=time.monotonic()-started)
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(stage='done',out=args.out,N=n,budget=budget,seconds=result['seconds'],
        comparisons=[(x['preparation'],x['lag'],x['state_distance_interval']) for x in comparisons])),flush=True)


if __name__=='__main__':main()
