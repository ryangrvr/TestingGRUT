"""Check the new vector implementation without repeating the old mode ensemble."""
from timing_common import *

def main():
    exe=build();groups={}
    times=np.array([0.,1.,3.,10.,40.,80.]);ks=np.array([.02,2.13,5.8])
    bg,vec,_=integrate(exe,'verify_vectors',3.,.5,ks,times,rtol=2e-12)
    oldexe=M.build();oldbg,tt,_=M.integrate(oldexe,'timing_matrix_control',3.,.5,ks,times,rtol=2e-12)
    target=np.einsum('tkij,j->tki',tt,phase_seed(3.))
    err=float(np.max(abs(vec-target)/(1+abs(target))))
    berr=float(np.max(abs(bg-oldbg)/(1+abs(oldbg))))
    groups['direct_vector_vs_full_transfer']=dict(passed=max(err,berr)<1e-6,vector_scaled_error=err,background_scaled_error=berr)

    v,_=M.prepared(3.);seed=phase_seed(3.);k=2.13
    def rhs(t,y):
        return np.r_[M.rhs(t,y[:19],3.,.5,M.XI,M.ETA2),M.generator(y[:19],3.,k)@y[19:]]
    sol=solve_ivp(rhs,(0,3),np.r_[v,seed],method='Radau',rtol=2e-9,atol=2e-12)
    assert sol.success
    short,sv,_=integrate(exe,'verify_radau',3.,.5,[k],[0.,3.],rtol=2e-12)
    radau=float(np.max(abs(sol.y[19:,-1]-sv[-1,0])/(1+abs(sv[-1,0]))))
    groups['short_Radau']=dict(passed=radau<1e-6,scaled_state_error=radau,nfev=sol.nfev)

    old=np.load(MODE/'PHASE_PROFILES.npz');times=old['times'];kk=2.122989207841491
    prior_ks=old['b3_g0.5_ks'];idx=int(np.argmin(abs(prior_ks-kk)));assert abs(prior_ks[idx]-kk)<1e-12
    b,vec,_=integrate(exe,'verify_phase_archive',3.,.5,[kk],times,rtol=2e-12)
    amp,dq,nx=responses(3.,b,vec);window=times>=64
    mean=np.trapezoid(nx[window],times[window])/(times[window][-1]-times[window][0])
    archived=old['b3_g0.5'][:,idx,:]
    # Old charge response uses the opposite overall current sign; only its
    # signed convention is converted here, not a particle interpretation.
    ae=float(np.max(abs(amp[:,0]-archived[:,0])/(1+abs(archived[:,0]))))
    ce=float(np.max(abs(dq[:,0]/abs(mean)+archived[:,1])/(1+abs(archived[:,1]))))
    groups['published_phase_archive']=dict(passed=max(ae,ce)<1e-5,amplitude_scaled_error=ae,charge_scaled_error_after_current_convention=ce,
        scope='Archive stores float32 response profiles; compare one selected existing momentum, not the old ensemble.')
    save('RESPONSE_CHECKS.json',dict(groups=groups))
    print((ROOT/'RESPONSE_CHECKS.json').read_text(),flush=True)
    assert all(g['passed'] for g in groups.values())

if __name__=='__main__':main()
