"""Radial action, lab-frequency action and explicit free-field matching controls."""
from pt_common import *
from scipy.optimize import brentq

def radial(E0,Q,lam=LX,mass=0.):
    roots=np.roots([lam/2,mass*mass,-2*E0,Q*Q]);rr=np.sort(roots.real[abs(roots.imag)<1e-8]);neg,lo,hi=rr
    assert neg<0<lo<hi
    def y(u):return lo+(hi-lo)*np.sin(u)**2
    def dt(u):return np.sqrt(2/lam)/np.sqrt(y(u)-neg)
    Tr=2*quad(dt,0,np.pi/2,epsabs=1e-11,epsrel=1e-12)[0]
    Dtheta=2*quad(lambda u:Q/y(u)*dt(u),0,np.pi/2,epsabs=1e-11,epsrel=1e-12)[0]
    action=quad(lambda u:np.sqrt(lam/2)*(hi-lo)**2*np.sin(u)**2*np.cos(u)**2*np.sqrt(y(u)-neg)/y(u),0,np.pi/2,epsabs=1e-9,epsrel=1e-12)[0]/np.pi
    return dict(radial_period=Tr,angular_advance=Dtheta,radial_action=action,rmin=np.sqrt(lo),rmax=np.sqrt(hi))

def main():
    data=np.load(PARENT/'MATCH_L12.npz');y=np.r_[data['light_mean_q'][:2],data['light_mean_p'][:2]]
    from initial import initialize as old_initial
    _,_,original=old_initial();Qnet=original['nQ']
    Q=y[0]*y[3]-y[1]*y[2];E0=.5*np.sum(y[2:]**2)+LX/4*np.sum(y[:2]**2)**2
    orbit=radial(E0,Q);Tr=orbit['radial_period'];Delta=orbit['angular_advance']
    def rhs(t,v):return np.r_[v[2:],-LX*np.sum(v[:2]**2)*v[:2]]
    M=8192;times=np.linspace(0,Tr,M,endpoint=False)
    sol=solve_ivp(rhs,(0,Tr),y,t_eval=times,rtol=2e-13,atol=1e-11,method='DOP853',dense_output=True)
    psi=sol.y[0]+1j*sol.y[1];cn=np.fft.fft(psi*np.exp(-1j*Delta/Tr*times))/M
    frequencies=Delta/Tr+2*np.pi/Tr*np.fft.fftfreq(M)*M
    plus=np.sum((frequencies*abs(cn)**2)[frequencies>0]);minus=-np.sum((frequencies*abs(cn)**2)[frequencies<0])
    phases=np.arange(8)/8;starts=[sol.sol(float(p*Tr)) for p in phases];rows=[];mout=30.
    # These deformations define asymptotic free modes, not a cosmic transition.
    for duration in [0.,.2,2.,20.,80.]:
      for phase,initial in zip(phases,starts):
        if not duration:
            final=initial;work=None;work_error=None;charge_error=None
        else:
            def flow(t,v):
                u=t/duration;f=.5*(1-np.cos(np.pi*u));fp=.5*np.pi/duration*np.sin(np.pi*u)
                m2=mout*mout*f*f;lam=LX*(1-f);r2=v[:2]@v[:2]
                force=-(m2+lam*r2)*v[:2]
                wd=mout*mout*f*fp*r2-LX*fp*r2*r2/4
                return np.r_[v[2:4],force,wd]
            ss=solve_ivp(flow,(0,duration),np.r_[initial,0.],rtol=2e-11,atol=2e-10,method='DOP853')
            final=ss.y[:4,-1];work=ss.y[4,-1]
            Ef=.5*(final[2:]@final[2:]+mout*mout*(final[:2]@final[:2]))
            work_error=abs((Ef-E0-work)/max(Ef,E0));charge_error=abs((final[0]*final[3]-final[1]*final[2])/Q-1)
        Ef=.5*(final[2:]@final[2:]+mout*mout*(final[:2]@final[:2]))
        nm=(Ef/mout-Q)/2
        rows.append(dict(duration=duration,radial_phase=phase,free_nminus_per_charge=nm/Q,free_nminus_per_original_net_charge=nm/Qnet,work=work,
            work_identity_relative=work_error,charge_relative_error=charge_error))
    summary=[]
    for d in sorted(set(x['duration'] for x in rows)):
        a=np.array([r['free_nminus_per_charge'] for r in rows if r['duration']==d])
        summary.append(dict(duration=d,mean=a.mean(),minimum=a.min(),maximum=a.max(),spread=a.ptp() if hasattr(a,'ptp') else np.ptp(a)))
    mass_only=[]
    for phase,initial in zip(phases[::2],starts[::2]):
        duration=80.
        def flow_fixed(t,v):
            u=t/duration;f=.5*(1-np.cos(np.pi*u));fp=.5*np.pi/duration*np.sin(np.pi*u)
            m2=mout*mout*f*f;r2=v[:2]@v[:2]
            return np.r_[v[2:4],-(m2+LX*r2)*v[:2],mout*mout*f*fp*r2]
        ss=solve_ivp(flow_fixed,(0,duration),np.r_[initial,0.],rtol=2e-11,atol=2e-10,method='DOP853')
        f=ss.y[:4,-1];r2=f[:2]@f[:2];Efree=.5*(f[2:]@f[2:]+mout*mout*r2);U=LX*r2*r2/4
        action=radial(Efree+U,Q,lam=LX,mass=mout)['radial_action']
        mass_only.append(dict(radial_phase=phase,duration=duration,quadratic_negative_projection_per_charge=(Efree/mout-Q)/(2*Q),
            final_radial_action_per_charge=action/Q,quartic_fraction_of_final_energy=U/(Efree+U),
            charge_relative_error=abs((f[0]*f[3]-f[1]*f[2])/Q-1),work_identity_relative=abs((Efree+U-E0-ss.y[4,-1])/(Efree+U))))
    result=dict(scope='Homogeneous matching experiment with supplied external mass/quartic ramp, not the binding history or an abundance.',
        initial_charge=Q,original_full_state_net_charge=Qnet,initial_energy=E0,orbit=orbit,radial_action_per_charge=orbit['radial_action']/Q,
        radial_action_per_original_net_charge=orbit['radial_action']/Qnet,lab_negative_action_per_original_net_charge=minus/Qnet,
        lab_negative_frequency_action_per_charge=minus/Q,lab_positive_frequency_action_per_charge=plus/Q,
        lab_charge_reconstruction_relative=abs((plus-minus)/Q-1),matching_mass=mout,rows=rows,summary=summary,
        fixed_coupling_mass_ramp_control=mass_only,
        rule='In an asymptotic free complex scalar, coherent counter-rotation counts in nbar. Slow integrable matching preserves Jr, with nbar_out=Jr. A lab-frequency split of an interacting orbit need not equal Jr. A coherent contribution must remain a field until this matching is justified.',
        physical_nbar=None)
    save('COHERENT_MATCHING.json',result)
    print('coherent lab',minus/Q,'radial action',orbit['radial_action']/Q,'charge',result['lab_charge_reconstruction_relative'],flush=True)
    print(json.dumps(summary,default=lambda x:float(x),indent=2),flush=True)

if __name__=='__main__':main()
