"""Supplied tiny homogeneous phi seeds, before nonlinear domain formation.

The background evolves nonlinearly; one charged Fourier perturbation evolves
linearly on it. No finite fragments, mergers or quantum seeds are generated.
"""
from scipy.integrate import solve_ivp
import numpy as np
from common import *
from elliptic import orbit,mode

def initial_growing_X(k):
    bg,t0,P,angle,diag=orbit();_,M,_=mode(k);ids=[0,1,3,4]
    eig,vec=np.linalg.eig(M[np.ix_(ids,ids)]);v=vec[:,np.argmax(abs(eig))]
    v*=np.exp(-1j*np.angle(v[np.argmax(abs(v))])) if np.iscomplexobj(v) else 1.
    v=v.real
    def f(t,y):
        q=bg.sol(t)[:2];H=LX*((q@q)*np.eye(2)+2*np.outer(q,q))+k*k*np.eye(2);T=y.reshape(4,4)
        return np.vstack((T[2:],-H@T[:2])).ravel()
    sol=solve_ivp(f,(0,t0),np.eye(4).ravel(),method='DOP853',rtol=2e-12,atol=2e-14)
    v=np.linalg.solve(sol.y[:,-1].reshape(4,4),v);v/=np.linalg.norm(v[:2])/RC
    return np.r_[v[:2],0.,v[2:],0.]

def run(seed_phi,k,rtol=2e-11):
    bg,t0,P,angle,diag=orbit();cold=diag['initial'];q=np.array([RC,0.,NU*seed_phi]);p=np.array([cold[2],cold[3],0.])
    # At RC, the quadratic phi potential vanishes; match even its tiny
    # quartic energy by the radial momentum, rather than adding energy.
    p[0]=-np.sqrt(max(0.,2*(EX+V0-potential(q))-p[1]**2))
    d=initial_growing_X(k);initial=np.r_[q,p,d]
    def f(t,y):
        q=y[:3];p=y[3:6];d=y[6:];H=hessian(q)+k*k*np.eye(3)
        return np.r_[p,-gradient(q),d[3:],-H@d[:3]]
    def stop(t,y):return abs(y[2])-.1*NU
    stop.direction=1;stop.terminal=True
    sol=solve_ivp(f,(0,120),initial,method='DOP853',rtol=rtol,atol=rtol/100,events=stop,dense_output=True,max_step=.06)
    if not sol.success or not len(sol.t_events[0]):raise RuntimeError('No phi threshold reached')
    ts=np.linspace(0,sol.t[-1],801);yy=sol.sol(ts);qs=yy[:3];ps=yy[3:6];ds=yy[6:]
    energy=np.sum(ps*ps,axis=0)/2+np.array([potential(q) for q in qs.T]);charge=qs[0]*ps[1]-qs[1]*ps[0]
    amplitude=np.sqrt(np.sum(ds[:2]**2,axis=0))/np.sqrt(np.sum(qs[:2]**2,axis=0))
    dq=(ds[0]*ps[1]+qs[0]*ds[4]-ds[1]*ps[0]-qs[1]*ds[3])/Q0
    return dict(seed_phi_over_v=seed_phi,charged_seed_relative=1e-6,tau_to_phi_0p1v=float(sol.t[-1]),
        charge_error=float(max(abs(charge/Q0-1))),energy_error=float(max(abs(energy/(EX+V0)-1))),
        maximum_X_amplitude_for_supplied_seed=float(max(amplitude)*1e-6),
        maximum_X_charge_contrast_for_supplied_seed=float(max(abs(dq))*1e-6),
        times=ts,phi_over_v=qs[2]/NU,X_amplitude_for_supplied_seed=amplitude*1e-6,
        scope='Two supplied homogeneous neutral seeds; cold vacuum, no heat/bias/expansion. Stops at |phi|=.1v. Spatial charged perturbation remains linear; no object population derived.')

def main():
    e=read('ELLIPTIC_RESULTS.json');k=e['peaks']['X']['k'];mux=e['peaks']['X']['mu_X'];mup=e['peaks']['phi']['mu_phi'];ratio=mup/mux
    rows=[run(seed,k) for seed in [1e-6,1e-12]]
    save('ROLL_RESULTS.json',dict(rows=rows,mu_phi_over_mu_X=ratio,
        conditional_seed_competition=[dict(delta_X=d,phi_seed_over_v_needed_for_X_first=d**ratio) for d in [1e-3,1e-6]],
        competition_assumptions='Independent fastest modes and equal unit nonlinear thresholds; envelopes only, not cosmic initial conditions. Actual roll uses .1v threshold and a supplied growing X mode.'))
    print('ROLL',[(r['seed_phi_over_v'],r['tau_to_phi_0p1v'],r['maximum_X_amplitude_for_supplied_seed'],r['maximum_X_charge_contrast_for_supplied_seed']) for r in rows],ratio,flush=True)

if __name__=='__main__':main()
