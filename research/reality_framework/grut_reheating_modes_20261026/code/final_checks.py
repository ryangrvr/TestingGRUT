"""Phase tangent and auxiliary orbit controls; aggregate with mode checks."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from mode_common import ROOT,write
from modes import prepared,realvec,rhs,XI,ETA2,build,integrate,matrices

def main():
    b,gam=3.,.5;delta=1e-4;states=[]
    for d in (delta,-delta):
        v,_=prepared(b,delta=d)
        s=solve_ivp(lambda t,y:rhs(t,y,b,gam,XI,ETA2),(0,3),v,method='DOP853',rtol=2e-12,atol=2e-14)
        assert s.success
        _,_,_,q,p=matrices(s.y[:,-1],b);states.append(np.r_[q,p])
    fd=(states[0]-states[1])/(2*delta)
    phase=json.loads((ROOT/'PHASE_RESULTS.json').read_text())['rows'][0]
    seed=np.tile(phase['phase_derivative'],2)
    bg,tt,ctl=integrate(build(),'phase_tangent_control',b,gam,[0.],[0,3],e=0.,rtol=2e-12)
    target=tt[-1,0]@seed;err=float(np.max(abs(fd-target)/(1+abs(target))))
    orbit=json.loads((ROOT/'QUARTIC_ORBIT.json').read_text())
    data=json.loads((ROOT/'VERIFICATION_RESULTS.json').read_text())
    data['groups']['phase_preparation_tangent']=dict(passed=err<1e-5,full_background_finite_difference_error=err,
        scope='k=0 ungauged identity control with gravity differences retained in full backgrounds; finite-k response remains the linear tangent calculation.')
    data['groups']['quartic_orbit_and_circular_control']=dict(passed=orbit['symplectic_error']<1e-7 and orbit['circular_multiplier_error']<1e-7,
        symplectic_error=orbit['symplectic_error'],circular_multiplier_error=orbit['circular_multiplier_error'],
        peak=max(orbit['rows'],key=lambda r:r['mu_conformal']))
    data['passed_groups']=sum(g['passed'] for g in data['groups'].values())
    write('VERIFICATION_RESULTS.json',data)
    print('FINAL CHECKS',data['passed_groups'],'phase',err,flush=True)
    if not all(g['passed'] for g in data['groups'].values()):raise SystemExit(1)

if __name__=='__main__':main()
