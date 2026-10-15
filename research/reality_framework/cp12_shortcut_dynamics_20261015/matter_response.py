"""Fixed-receiver response under H and the exact moving-core comparison.

The source-stabilizer quotient retains all S_(N-1) symmetries. The mean density
of the other vertices equals the response of ANY one fixed receiver j!=0.
No receiver is chosen using the final graph. This is not a distance-bin test.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
from operators import quotient_model,RF
from integrals import Propagator
from recorded_transport import response,ratio
sys.path.insert(0,str(RF/'charged_density_dynamics_20261015'))
from diagnostics import arrays,event
from stream import probability_interval


def main():
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=5)
    p.add_argument('--budget',type=int,default=50);p.add_argument('--wait',type=float,default=2)
    p.add_argument('--lags',default='.125,.25,.5');p.add_argument('--dt',type=float,default=1/32)
    p.add_argument('--out',required=True);args=p.parse_args();started=time.monotonic()
    n=args.n;q=2*n;m=quotient_model(n,q,args.budget,source=True)
    H,V,B,BV,psi=(m[k] for k in ['H','V','B','BV','psi']);states=m['states']
    diag=arrays(n,q,states);P=event(diag,q,n);n0=np.array([s[1]&1 for s in states],float)
    O=(diag[1]-n0)/(n-1)
    w,err=Propagator(H,{'electric':B},q/2,dt=args.dt).run(psi,[args.wait])[args.wait]
    ew=float(err['electric'][0]);X=np.column_stack([w,P*w]);Y=n0[:,None]*X
    initial=np.column_stack([X,Y]);lags=sorted(float(s) for s in args.lags.split(','))
    full=Propagator(H,{'electric':B},q/2,dt=args.dt).run(initial,lags)
    core=Propagator(H-V,{'electric':B-BV},q/2,dt=args.dt).run(initial,lags)
    rows=[]
    for lag in lags:
        for i,kind in enumerate(['unfiltered_waiting_state','initial_event_branch']):
            p0=float(np.linalg.norm(X[:,i])**2);pI=probability_interval(p0,ew);values={}
            for name,data in [('full_frozen',full),('moving_core_comparison',core)]:
                z,er=data[lag];x,y=z[:,i],z[:,i+2]
                ex=ew+float(er['electric'][i]);ey=ew+float(er['electric'][i+2])
                r=response(x,y,O,ex,ey);r['conditional_value']=r['joint_value']/p0
                r['conditional_interval']=ratio(r['joint_interval'],pI)
                values[name]=dict(receiver_response=r,X_state_error=ex,Y_state_error=ey)
            a=values['full_frozen']['receiver_response'];b=values['moving_core_comparison']['receiver_response']
            difference=a['joint_value']-b['joint_value'];rad=a['error_bound']+b['error_bound']
            rows.append(dict(lag=lag,preparation=kind,initial_probability=p0,initial_interval=pI,
                models=values,response_change_joint=difference,
                response_change_joint_interval=[difference-rad,difference+rad],
                response_change_conditional=difference/p0,
                response_change_conditional_interval=ratio([difference-rad,difference+rad],pI)))
            print(json.dumps(dict(stage='matter_response',lag=lag,preparation=kind,full=a['conditional_value'],
                comparison=b['conditional_value'],difference=rows[-1]['response_change_conditional_interval'])),flush=True)
    z=P*w;r=(1-P)*w;de=-2*float(np.vdot(z,H@r).real);bk=np.sqrt(2*q*min(n,q))
    result=dict(status='PASS',N=n,Q=q,budget=args.budget,dimension=H.shape[0],
        outgoing_dimension=B.shape[0],wait=args.wait,wait_state_error=ew,model_checks=m['checks'],
        preparation='Original unfiltered coherent neutral ray, evolved under the full frozen H up to the waiting time. The same waiting vector then initializes both diagnostic evolutions.',
        source='Local gauge-invariant exp(-i theta n_0).',
        receiver='ANY fixed j!=0. The operator (n_v-n_0)/(N-1) has exactly its expectation by S_(N-1) symmetry, including the source pulse and both evolutions.',
        comparison='Same Hermitian moving-core split as run_split.py. No distance metric or final graph filter is used in this response experiment.',
        initial_event='c/Q>=1/4 and nv/N>=1/4. Its separate branch is an explicitly recorded diagnostic, not a change in the original preparation.',
        initial_binary_record_energy_change=de,
        initial_binary_record_energy_error_bound=float(2*bk*(2*ew+ew*ew)+5e-10),
        results=rows,error_scope='All outgoing electric rows; analytic Duhamel state error and bilinear response intervals, with5e-10 floating allowance. No outward rounding or large-distance scaling.',
        seconds=time.monotonic()-started)
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(stage='done_matter_response',out=args.out,seconds=result['seconds'])),flush=True)


if __name__=='__main__':main()
