"""Existing derivative inflaton channel: rates, quantum-wave check, receivers.

No imports from prior packet calculations. Every approximation is scoped in
MATHEMATICS.md. A pulse is an external-pump QFT check of the decay coefficient;
the averaged fluid calculation separately includes the conserved receiver.
"""
from pathlib import Path
import json, math, platform, sys, time
import numpy as np
from scipy.integrate import solve_ivp, simpson
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
I=json.loads((ROOT/'INPUTS.json').read_text())
CHECKS=[]


def check(name,ok,value):
    CHECKS.append({'name':name,'passed':bool(ok),'value':value})
    if not ok:raise AssertionError((name,value))


def write(name,data):
    (ROOT/name).write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n')


def calibration(c):
    k=I['kappa_over_c_squared']*c*c;N=I['N_pivot'];M=I['M_GeV']
    xe=math.log1p(math.sqrt(2/k))
    x=brentq(lambda a:k/2*(math.exp(a)-a-math.exp(xe)+xe)-N,xe,20)
    z=math.exp(x);eps=2/(k*(z-1)**2);eta=2*(2-z)/(k*(z-1)**2)
    alpha=math.sqrt(24*math.pi**2*I['A_s']*eps/(-math.expm1(-x))**2)
    V0=I['lambda_H']*I['v_GeV']**4/4
    s0=c*V0/(alpha**2*M**4)
    P=I['b']**2*math.exp(-c*s0)/math.expm1(-s0)**2
    F=P/(k+P);mass=math.sqrt(2/k)*alpha*M
    gamma=F*F*mass**3/(128*math.pi*k*M*M)
    # Finite difference of the literal metric, not of the derived vertex.
    def metric(x):return k+P*math.exp(-x)
    d=1e-3
    slope=(math.log(metric(d))-math.log(metric(-d)))/(2*d*math.sqrt(k))
    slope_error=abs(slope/(-F/math.sqrt(k))-1)
    # Local normal coordinates remove the kinetic cubic; potential restores it.
    # phi_old = phi_RNC - g chi_RNC^2/4 gives |A|=g m^2/2.
    amp_derivative=F*mass*mass/(2*math.sqrt(k)*M)
    amp_RNC=(F/(4*math.sqrt(k)*M))*mass*mass*2
    check('metric vertex c=%g'%c,slope_error<1e-8,slope_error)
    check('normal-coordinate amplitude c=%g'%c,abs(amp_RNC/amp_derivative-1)<1e-14,float(amp_RNC/amp_derivative))
    # Frozen homogeneous sigma must be small not to dominate before this decay.
    s_limit=math.sqrt(3)*gamma/(alpha*M)
    ns=1-6*eps+2*eta;r=16*eps;D=8*(1-ns)-r;E=4*(1-ns)-r
    kobs=8*r/E**2
    mobs=M*math.sqrt(3*math.pi**2*I['A_s']/32)*D
    gobs=F*F*mobs**3/(128*math.pi*kobs*M*M)
    rho_selected=c*c*V0*V0/(alpha*alpha*M**4)
    return dict(c=c,b=I['b'],kappa=k,alpha=alpha,ns=ns,r=r,
        sigma_selected_over_M=s0,pole_fraction_F=F,regular_metric_fraction=k/(k+P),
        inflaton_mass_GeV=mass,partial_width_GeV=gamma,width_over_mass=gamma/mass,
        selected_prepared_rho_GeV4=rho_selected,
        mass_remnant_identity_ratio=rho_selected*mass*mass/(4*V0*V0/(3*M*M)),
        observable_eliminated_mass_GeV=mobs,observable_eliminated_width_GeV=gobs,
        total_density_at_H_equals_partial_width_GeV4=3*M*M*gamma*gamma,
        frozen_sigma_no_domination_scale_over_M=s_limit,
        scope='N=55 single-adiabatic calibration; local post-pole light daughters; density scale assumes this is the total width. The frozen-sigma threshold is a conditional early-domination diagnostic, not a propagated abundance.')


def pump(t,q,T):
    f=q*math.exp(-t*t/(2*T*T));co=math.cos(t);si=math.sin(t)
    z=f*co
    zd=f*(-si-t*co/(T*T))
    zdd=f*((t*t/T**4-1/T**2-1)*co+2*t*si/(T*T))
    return z,zd,zdd


def pulse(q,npoints,refined=False):
    P=I['pulse'];T=P['duration_in_mass_units'];te=P['endpoint_durations']*T
    p=np.linspace(P['p_min'],P['p_max'],npoints)
    z0,zd0,_=pump(-te,q,T)
    v=np.exp(1j*p*te)/np.sqrt(2*p);vp=-1j*p*v
    def rhs(t,y):
        z,zd,zdd=pump(t,q,T)
        return np.concatenate([y[npoints:],-(p*p+zdd/2-zd*zd/4)*y[:npoints]])
    sol=solve_ivp(rhs,(-te,te),np.concatenate([v,vp]),method='DOP853',
                  rtol=2e-12 if refined else 2e-11,atol=2e-14,
                  max_step=.2 if refined else .4)
    if not sol.success:raise RuntimeError(sol.message)
    v,vp=sol.y[:npoints,-1],sol.y[npoints:,-1]
    al=np.sqrt(p/2)*(v-vp/(1j*p))*np.exp(1j*p*te)
    be=np.sqrt(p/2)*(v+vp/(1j*p))*np.exp(-1j*p*te)
    number=abs(be)**2
    # First-Born Fourier transform of a smooth real Gaussian carrier.
    ft=math.sqrt(2*math.pi)*T/2*(np.exp(-.5*T*T*(2*p-1)**2)+np.exp(-.5*T*T*(2*p+1)**2))
    born=(q*p*ft)**2
    rho=float(simpson(p**3*number,x=p)/(2*math.pi**2))
    rhob=float(simpson(p**3*born,x=p)/(2*math.pi**2))
    width_limit=q*q*T/(256*math.sqrt(math.pi))
    corrected=width_limit*(1+5/T**2+15/(4*T**4))
    errw=float(np.max(abs(abs(al)**2-abs(be)**2-1)))
    check('pulse Wronskian q=%g n=%d'%(q,npoints),errw<2e-8,errw)
    check('Born integral vs decay coefficient q=%g n=%d'%(q,npoints),abs(rhob/corrected-1)<2e-8,float(rhob/corrected-1))
    check('quantum production vs Born q=%g n=%d'%(q,npoints),abs(rho/rhob-1)<.01,float(rho/rhob-1))
    name='pulse_q%g_n%d%s.csv'%(q,npoints,'_refined' if refined else '')
    np.savetxt(ROOT/'data'/name,np.column_stack([p,number,born,be.real,be.imag]),delimiter=',',
               header='k_over_m,number,Born_number,beta_real,beta_imag',comments='')
    return {'q':q,'npoints':npoints,'T':T,'rho_over_m4':rho,'Born_rho_over_m4':rhob,
            'infinite_duration_width_limit':width_limit,'finite_envelope_corrected_Born':corrected,
            'relative_quantum_Born':float(rho/rhob-1),'max_Wronskian_error':errw,
            'peak_occupation':float(max(number)),'nfev':sol.nfev,'file':'data/'+name}


def direct_coordinate_control(q=.002):
    T=I['pulse']['duration_in_mass_units'];te=I['pulse']['endpoint_durations']*T
    p=.5
    z,zd,_=pump(-te,q,T);f=math.exp(-z/2)
    v=np.exp(1j*p*te)/math.sqrt(2*p);vp=-1j*p*v
    # Evolve original derivative-coupled variable chi, not the canonical mode.
    def rhs(t,y):
        z,zd,zdd=pump(t,q,T)
        return [y[1],zd*y[1]-p*p*y[0]]
    s=solve_ivp(rhs,(-te,te),[v/f,(vp+zd*v/2)/f],method='DOP853',
                rtol=2e-12,atol=2e-14,max_step=.2)
    z,zd,_=pump(te,q,T);f=math.exp(-z/2)
    chi,chid=s.y[:,-1];vv=f*chi;vvd=f*(chid-zd*chi/2)
    beta=math.sqrt(p/2)*(vv+vvd/(1j*p))*np.exp(-1j*p*te)
    return {'q':q,'p':p,'number':float(abs(beta)**2)}


def zero_mode_control():
    q=.4;H=.08;omega=.9;decay=.015
    def background(t):
        z=q*math.exp(-decay*t)*math.cos(omega*t)
        zd=q*math.exp(-decay*t)*(-decay*math.cos(omega*t)-omega*math.sin(omega*t))
        return z,zd
    def rhs(t,y):
        _,zd=background(t)
        return [y[1],-(3*H-zd)*y[1]]
    initial_v=.2
    s=solve_ivp(rhs,(0,50),[.7,initial_v],method='DOP853',rtol=2e-12,atol=2e-14,dense_output=True)
    tt=np.linspace(0,50,251);y=s.sol(tt)
    moment=np.array([math.exp(3*H*t-background(t)[0])*y[1,j] for j,t in enumerate(tt)])
    e=float(np.max(abs(moment/moment[0]-1)))
    # Every initially stationary patch is exactly stationary, regardless of z(t).
    st=solve_ivp(rhs,(0,50),[1.234,0],rtol=1e-11,atol=1e-13)
    change=float(abs(st.y[0,-1]-1.234))
    check('zero-mode canonical momentum conserved',e<1e-7,e)
    check('stationary homogeneous coordinate not reset',change<1e-14,change)
    return {'momentum_relative_error':e,'stationary_position_change':change,
            'scope':'Exact shift-symmetric kinetic limit; no assertion about the earlier full stochastic transport or non-negligible potential.'}


def receiver(R,ratio,refined=False):
    # tau=Gamma_total t; energies in 3 M^2 Gamma_total^2.
    f=1/(1+ratio)
    def rhs(t,y):
        u,rd,rv,h,N=y
        return [-(3*h+1)*u,-4*h*rd+f*u,-4*h*rv+(1-f)*u,
                -1.5*h*h-.5*(rd+rv),h]
    t0=2/(3*R);end=I['receiver']['end_Gamma_t']
    s=solve_ivp(rhs,(t0,end),[R*R,0,0,R,0],method='DOP853',
            rtol=2e-12 if refined else 2e-10,atol=2e-14,
            max_step=.15 if refined else .3,dense_output=True)
    if not s.success:raise RuntimeError(s.message)
    tt=np.geomspace(t0,end,601);u,rd,rv,h,N=s.sol(tt)
    C=h*h-u-rd-rv
    res=float(np.max(abs(C)/(h*h+abs(u)+abs(rd)+abs(rv))))
    equip=brentq(lambda t:float(s.sol(t)[0]-sum(s.sol(t)[1:3])),t0,end)
    Nt=float(s.sol(equip)[4])
    check('receiver constraint R=%g ratio=%g refined=%s'%(R,ratio,refined),res<2e-7,res)
    actual_fraction=float(rd[-1]/(rd[-1]+rv[-1]))
    check('receiver branching R=%g ratio=%g refined=%s'%(R,ratio,refined),abs(actual_fraction-f)<1e-10,actual_fraction)
    check('receiver positivity R=%g ratio=%g refined=%s'%(R,ratio,refined),min(np.min(u),np.min(rd),np.min(rv))>=-1e-12,float(min(np.min(u),np.min(rd),np.min(rv))))
    data=np.column_stack([tt,N,h,u,rd,rv])
    name='receiver_R%g_vis%g%s.csv'%(R,ratio,'_refined' if refined else '')
    np.savetxt(ROOT/'data'/name,data,delimiter=',',header='Gamma_total_t,N,H_over_Gamma_total,rho_I_scaled,rho_dark_scaled,rho_visible_scaled',comments='')
    return {'R':R,'visible_over_dark_width':ratio,'dark_branch':f,'equality_Gamma_t':equip,
            'equality_N':Nt,'final_N':float(N[-1]),'final_H_over_Gamma':float(h[-1]),
            'final_inflaton_fraction':float(u[-1]/(u[-1]+rd[-1]+rv[-1])),
            'constraint_relative':res,'file':'data/'+name,'nfev':s.nfev}


def receiver_N_control(R):
    # Independent variables: Omega_I and log(H/Gamma), no densities or Raychaudhuri h.
    def rhs(N,y):
        f,l,t=y;return [f*(1-f-math.exp(-l)),-2+f/2,math.exp(-l)]
    def done(N,y):return y[2]-I['receiver']['end_Gamma_t']
    done.terminal=True;done.direction=1
    t0=2/(3*R)
    s=solve_ivp(rhs,(0,40),[1,math.log(R),t0],events=done,method='DOP853',rtol=2e-12,atol=2e-14,max_step=.02)
    if s.status!=1:raise RuntimeError(s.message)
    return {'final_N':float(s.t[-1]),'final_H_over_Gamma':math.exp(float(s.y[1,-1])),
            'final_inflaton_fraction':float(s.y[0,-1])}


def main():
    start=time.perf_counter()
    cal=[calibration(c) for c in I['c_values']]
    check('common-action mass-remnant identity',max(abs(x['mass_remnant_identity_ratio']-1) for x in cal)<2e-14,max(abs(x['mass_remnant_identity_ratio']-1) for x in cal))
    check('inflation-observable rate elimination',max(abs(x['observable_eliminated_width_GeV']/x['partial_width_GeV']-1) for x in cal)<2e-12,max(abs(x['observable_eliminated_width_GeV']/x['partial_width_GeV']-1) for x in cal))
    check('c=1 scalaron-width normalization',abs(1/(128*I['kappa_over_c_squared'])-1/192)<1e-16,1/192)
    modes=[pulse(q,I['pulse']['points']) for q in I['pulse']['q_values']]
    refined=pulse(.001,I['pulse']['refined_points'],True)
    check('pulse grid/time refinement',abs(refined['rho_over_m4']/modes[0]['rho_over_m4']-1)<2e-6,float(refined['rho_over_m4']/modes[0]['rho_over_m4']-1))
    coord=direct_coordinate_control()
    raw=np.loadtxt(ROOT/modes[1]['file'],delimiter=',',skiprows=1)
    peak=raw[np.argmin(abs(raw[:,0]-.5)),1]
    ce=float(abs(coord['number']/peak-1))
    check('canonical versus original-coordinate production',ce<2e-6,ce)
    zero=zero_mode_control()
    R=I['receiver']['H_initial_over_Gamma']
    rec=[receiver(R,x) for x in I['receiver']['comparison_visible_over_dark_width']]
    fine=receiver(R,0,True)
    nc=receiver_N_control(R)
    check('receiver tolerance refinement',abs(fine['final_H_over_Gamma']/rec[0]['final_H_over_Gamma']-1)<2e-7,float(fine['final_H_over_Gamma']/rec[0]['final_H_over_Gamma']-1))
    check('independent receiver formulation',abs(nc['final_H_over_Gamma']/fine['final_H_over_Gamma']-1)<2e-7 and abs(nc['final_N']-fine['final_N'])<2e-7,nc)
    # SM entropy transfer. No contemporary data bound is asserted here.
    ent=4/7*I['gstar_visible_reheating']*(I['gstar_neutrino_decoupling']/I['gs_visible_reheating'])**(4/3)
    dr=[{'visible_over_dark_width':x,'Delta_Neff':ent/x if x else None,
         'scope':'No visible bath exists at ratio zero; other ratios supply an extra visible width and assume standard entropy history.'} for x in I['receiver']['comparison_visible_over_dark_width']]
    # An equivalent direct entropy-ratio calculation, rho_nu,one/rho_vis=7/43 at decoupling.
    direct=(43/7)*(I['gstar_neutrino_decoupling']/I['gstar_visible_reheating'])**(1/3)
    check('entropy conversion identity',abs(ent/direct-1)<1e-14,ent)
    out=dict(calibration=cal,pulse=modes,pulse_refined=refined,original_coordinate_control=coord,
        zero_mode=zero,receivers=rec,receiver_refined=fine,receiver_N_control=nc,
        dark_radiation=dict(coefficient=ent,cases=dr,
            illustrative_visible_over_dark_width_requirement=ent/I['illustrative_Delta_Neff_tolerance'],
            tolerance_is_not_an_observational_limit=True),checks=CHECKS,
        status='PASS' if all(x['passed'] for x in CHECKS) else 'FAIL')
    write('RESULTS.json',out)
    write('ENVIRONMENT.json',dict(python=sys.version,numpy=np.__version__,scipy=__import__('scipy').__version__,
        platform=platform.platform(),elapsed_seconds=time.perf_counter()-start,command='python code/calculate.py'))
    print(json.dumps({'status':out['status'],'checks':len(CHECKS),'c1_width_GeV':cal[1]['partial_width_GeV'],
          'wave_width_relative_errors':[x['relative_quantum_Born'] for x in modes],
          'Delta_Neff_coefficient':ent},indent=2))


if __name__=='__main__':main()
