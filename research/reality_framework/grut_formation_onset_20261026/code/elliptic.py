"""Cold phi=0 local orbit matched to the pinned release energy and net charge.

Relative Floquet spectra of charged X and the neutral binding field are
computed separately. Neither is a free-particle spectrum or a plasma result.
"""
from functools import lru_cache
from scipy.integrate import solve_ivp,quad
from scipy.optimize import minimize_scalar,brentq
import numpy as np
from common import *

J=np.block([[np.zeros((3,3)),np.eye(3)],[-np.eye(3),np.zeros((3,3))]])

@lru_cache(maxsize=4)
def orbit(rtol=2e-12):
    pt=Q0/RC;pr=-np.sqrt(2*(EX-LX*RC**4/4)-pt*pt)
    initial=np.array([RC,0.,pr,pt])
    def f(t,y):
        q=y[:2];return np.r_[y[2:],-LX*(q@q)*q]
    def ev(t,y):return y[:2]@y[2:]
    ev.direction=1
    sol=solve_ivp(f,(0,160),initial,method='DOP853',rtol=rtol,atol=rtol/30,events=ev,dense_output=True,max_step=.15)
    assert sol.success and len(sol.t_events[0])>=5
    t0,t1=sol.t_events[0][:2];P=t1-t0
    q0=sol.sol(t0)[:2];q1=sol.sol(t1)[:2]
    angle=np.arctan2(q0[0]*q1[1]-q0[1]*q1[0],q0@q1)
    ts=np.linspace(t0,t0+4*P,3001);b=sol.sol(ts)
    energy=np.sum(b[2:]**2,axis=0)/2+LX*np.sum(b[:2]**2,axis=0)**2/4
    charge=b[0]*b[3]-b[1]*b[2]
    return sol,t0,float(P),float(angle),dict(initial=initial,period=P,rotation_angle=angle,
        energy_X_m4=EX,charge_m3=Q0,energy_per_net_m=(EX+V0)/Q0,
        max_relative_energy_drift=float(max(abs(energy/EX-1))),max_relative_charge_drift=float(max(abs(charge/Q0-1))),
        radius_min=float(min(np.sqrt(np.sum(b[:2]**2,axis=0)))),radius_max=float(max(np.sqrt(np.sum(b[:2]**2,axis=0)))),
        expansion_H_over_m=K.Background().H_ref/MX)

def closing(angle):
    co,si=np.cos(angle),np.sin(angle);R=np.array([[co,si,0],[-si,co,0],[0,0,1.]])
    C=np.zeros((6,6));C[:3,:3]=R;C[3:,3:]=R;return C

def mode(k,rtol=3e-11,nperiod=1,method='DOP853',profile=False):
    bg,t0,P,angle,diag=orbit()
    def f(t,y):
        q=bg.sol(t0+t)[:2];H=hessian([*q,0.])+k*k*np.eye(3);T=y.reshape(6,6)
        return np.vstack((T[3:],-H@T[:3])).ravel()
    sol=solve_ivp(f,(0,P*nperiod),np.eye(6).ravel(),method=method,rtol=rtol,atol=rtol/100,dense_output=profile,max_step=.15)
    if not sol.success:raise RuntimeError(sol.message)
    M=closing(angle*nperiod)@sol.y[:,-1].reshape(6,6)
    X=M[np.ix_([0,1,3,4],[0,1,3,4])];phi=M[np.ix_([2,5],[2,5])]
    ex=np.linalg.eigvals(X);ep=np.linalg.eigvals(phi)
    mu=lambda eig:max(0.,float(np.log(abs(eig).max())/(P*nperiod)))
    out=dict(k=k,mu_X=mu(ex),mu_phi=mu(ep),multiplier_X=float(max(abs(ex))),multiplier_phi=float(max(abs(ep))),
        symplectic_absolute=float(np.max(abs(M.T@J@M-J))),determinant_error=float(abs(np.linalg.det(M)-1)),nfev=sol.nfev)
    return out,M,sol

def main():
    bg,t0,P,angle,diag=orbit();ks=np.linspace(0,3.,121);rows=[]
    for k in ks:rows.append(mode(float(k))[0])
    peaks={}
    for species in ['X','phi']:
        j=int(np.argmax([r['mu_'+species] for r in rows]));lo=ks[max(j-1,0)];hi=ks[min(j+1,len(ks)-1)]
        opt=minimize_scalar(lambda k:-mode(float(k))[0]['mu_'+species],bounds=(lo,hi),method='bounded',options={'xatol':1e-8})
        candidates=[rows[j],mode(float(opt.x))[0],rows[0]]
        peak=max(candidates,key=lambda r:r['mu_'+species]);peaks[species]=peak
        if abs(peak['mu_'+species]-rows[0]['mu_'+species])<1e-10:
            peaks[species]=rows[0]
    sample=np.linspace(0,4*P,2401);b=bg.sol(t0+sample)
    np.savez_compressed(ROOT/'ELLIPTIC_ORBIT.npz',tau=sample,q=b[:2],p=b[2:],
        neutral_curvature=G*G*np.sum(b[:2]**2,axis=0)-LP*NU*NU)
    save('ELLIPTIC_RESULTS.json',dict(background=diag,rows=rows,peaks=peaks,
        cold_charge_volume_at_X_peak_by_shape={str(f):float(f*Q0*(2*np.pi/peaks['X']['k'])**3) for f in [.1,.3,1.]},
        neutral_linear_charge_perturbation=0.,
        scope='Cold vacuum phi=0 reference with matched charge and energy; no heat/plasma, explicit phi bias, expansion/backreaction, quantum seed, free-mode particles or nonlinear fragments. Neutral phi modes carry zero linear X charge contrast. Finite 0<=k/m<=3 band only.'))
    print('ELLIPTIC',diag,peaks,flush=True)

if __name__=='__main__':main()
