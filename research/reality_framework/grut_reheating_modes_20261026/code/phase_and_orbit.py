"""Labelled phase-noise response and isolated late-quartic diagnostic."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from mode_common import ROOT,WORK,write
from modes import prepared,realvec,matrices,J


def orbit_diagnostic(bg,ks):
    _,_,_,chi,chip=matrices(bg,3.) # chi/chip do not depend on beta here
    init=np.r_[chi[:2],chip[:2]];lam=.001
    def force(t,y):
        q=y[:2];return np.r_[y[2:],-lam*(q@q)*q]
    def event(t,y):return y[:2]@y[2:]
    event.direction=1
    sol=solve_ivp(force,(0,30),init,events=event,dense_output=True,method='DOP853',rtol=2e-12,atol=2e-13)
    assert sol.success and len(sol.t_events[0])>=2
    t0,t1=sol.t_events[0][:2];period=t1-t0
    q0=sol.sol(t0)[:2];q1=sol.sol(t1)[:2]
    angle=np.arctan2(q0[0]*q1[1]-q0[1]*q1[0],q0@q1)
    rotation=np.array([[np.cos(angle),np.sin(angle)],[-np.sin(angle),np.cos(angle)]])
    closing=np.zeros((4,4));closing[:2,:2]=rotation;closing[2:,2:]=rotation
    J4=np.block([[np.zeros((2,2)),np.eye(2)],[-np.eye(2),np.zeros((2,2))]])
    rows=[];sym=0.
    for k in ks:
        def tangent(t,y):
            q=sol.sol(t0+t)[:2];mass=k*k*np.eye(2)+lam*((q@q)*np.eye(2)+2*np.outer(q,q))
            T=y.reshape(4,4);return np.vstack((T[2:],-mass@T[:2])).ravel()
        result=solve_ivp(tangent,(0,period),np.eye(4).ravel(),method='DOP853',rtol=3e-11,atol=3e-13)
        mon=closing@result.y[:,-1].reshape(4,4);eigen=np.linalg.eigvals(mon)
        sym=max(sym,float(np.max(abs(mon.T@J4@mon-J4))))
        mu=max(0.,float(np.log(abs(eigen).max())/period))
        rows.append(dict(k=float(k),mu_conformal=mu,multiplier_max=float(abs(eigen).max())))
    # Circular complex quartic condensate: stable finite-k control.
    radius=float(np.linalg.norm(init[:2]));omega=np.sqrt(lam)*radius;Tc=2*np.pi/omega;circ=[]
    for k in (0.1,.8,2.,4.,6.,10.):
        def tangent(t,y):
            q=radius*np.array([np.cos(omega*t),np.sin(omega*t)])
            M=k*k*np.eye(2)+lam*((q@q)*np.eye(2)+2*np.outer(q,q))
            T=y.reshape(4,4);return np.vstack((T[2:],-M@T[:2])).ravel()
        r=solve_ivp(tangent,(0,Tc),np.eye(4).ravel(),method='DOP853',rtol=2e-12,atol=2e-14)
        circ.append(float(abs(np.linalg.eigvals(r.y[:,-1].reshape(4,4))).max()))
    return dict(initial_chi=init,radial_period=period,rotation_angle=angle,rows=rows,
                symplectic_error=sym,circular_multiplier_error=max(abs(np.array(circ)-1)),
                scope='Isolated quartic diagnostic matched at tau10; six-field, curvature, decay and expansion corrections not evolved in this auxiliary orbit.')


def main():
    saved=np.load(ROOT/'MODE_PROFILES.npz');times=saved['times'];originalks=saved['ks']
    refined=np.load(ROOT/'REFINED_PROFILES.npz');fineks=refined['ks']
    base=json.loads((ROOT/'MODE_RESULTS.json').read_text());profiles=dict(times=times)
    rows=[]
    for entry in base['rows']:
        tag=entry['tag'];b=entry['beta'];mean=abs(entry['late_NX']);delta=1e-5
        vp,_=prepared(b,delta=delta);vm,_=prepared(b,delta=-delta)
        deriv=(vp[2:8]-vm[2:8])/(2*delta)*np.sqrt(2)
        # All physical initial velocities vanish, and gp=H_i g, so p=q
        # for this field-only phase perturbation at a=H/H_i=1.
        seed=np.r_[deriv,deriv];pieces=[]
        for source,ks,prefix in [('coarse',originalks,''),('fine',fineks,'ref_')]:
            for start in range(0,len(ks),8):
                if source=='coarse':
                    d=np.load(WORK/f'{tag}_k{start}.npz');bg=d['bg'];TT=d['trans'];kk=d['ks']
                else:
                    path=WORK/f'ref_{tag}_{start}.bin';kk=ks[start:start+8]
                    raw=np.fromfile(path).reshape(len(times),20+144*len(kk));bg=raw[:,1:20];TT=raw[:,20:].reshape(len(times),len(kk),12,12)
                responses=np.empty((len(times),len(kk),2))
                for n,v in enumerate(bg):
                    _,_,_,chi,chip=matrices(v,b)
                    amp=np.zeros(12);amp[:2]=chi[:2]/(chi[:2]@chi[:2])
                    charge=np.zeros(12);charge[:2]=[-chip[1],chip[0]];charge[6:8]=[chi[1],-chi[0]]
                    perturb=np.einsum('kij,j->ki',TT[n],seed)
                    responses[n,:,0]=perturb@amp;responses[n,:,1]=(perturb@charge)/mean
                pieces.append((kk,responses))
        allks=np.concatenate([r[0] for r in pieces]);allr=np.concatenate([r[1] for r in pieces],axis=1)
        order=np.argsort(allks);allks=allks[order];allr=allr[:,order]
        profiles[tag]=allr;profiles[tag+'_ks']=allks
        intervals=[]
        for end in (3.,10.,40.,80.):
            mask=times<=end;main=allks>=.01;rr=allr[mask][:,main];peak=np.unravel_index(abs(rr[:,:,0]).argmax(),rr[:,:,0].shape)
            intervals.append(dict(end_tau=end,max_amplitude_per_radian=float(abs(rr[:,:,0]).max()),
                max_charge_fraction_per_radian=float(abs(rr[:,:,1]).max()),
                amplitude_at_003rad=float(.03*abs(rr[:,:,0]).max()),charge_fraction_at_003rad=float(.03*abs(rr[:,:,1]).max()),
                peak_k=float(allks[main][peak[1]]),peak_tau=float(times[mask][peak[0]])))
        rows.append(dict(tag=tag,phase_derivative=deriv,intervals=intervals))
        print('PHASE',tag,intervals,flush=True)
    profiles={k:v.astype(np.float32) if np.asarray(v).dtype.kind=='f' and k!='times' else v for k,v in profiles.items()}
    np.savez_compressed(ROOT/'PHASE_PROFILES.npz',**profiles)
    write('PHASE_RESULTS.json',dict(rows=rows,phase_seed_radians=.03,
        scope='Classical initial cosine phase perturbation with locally readjusted conditional amplitudes; no cosmological seed selection or nonlinear evolution. Profile arrays float32, summaries float64.'))
    qbg=saved['b3_g0.5_bg'][np.argmin(abs(times-10))]
    orbit=orbit_diagnostic(qbg,fineks);write('QUARTIC_ORBIT.json',orbit)
    control=orbit['symplectic_error']<1e-7 and orbit['circular_multiplier_error']<1e-7
    print('QUARTIC',max(orbit['rows'],key=lambda r:r['mu_conformal']),orbit['symplectic_error'],orbit['circular_multiplier_error'],flush=True)
    if not control:raise SystemExit('Quartic control failed')

if __name__=='__main__':main()
