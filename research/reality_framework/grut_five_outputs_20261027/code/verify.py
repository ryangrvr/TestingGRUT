"""Necessary mathematical/numerical checks, not a physical abundance certification."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import simpson, quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import erf, spherical_jn
from calculate import ROOT, RF, PROFILE, inverse_speed, save


def main():
    read=lambda n:json.loads((ROOT/n).read_text())
    checks=[]
    def check(name,condition,diagnostics):
        checks.append({'name':name,'pass':bool(condition),'diagnostics':diagnostics})
    states=read('STATE_FAMILY.json')['states']
    err=max(max(s['constraint_error'],s['neutrality_error'],abs(s['global_purity']-1)) for s in states)
    check('distinct_legal_states',err<1e-13 and abs(states[0]['X_entropy_bits']-states[1]['X_entropy_bits'])<1e-14
          and states[0]['X_charge']==-states[1]['X_charge'] and states[2]['joint_coherence']*states[3]['joint_coherence']<0,
          {'constraint_or_purity_error':err,'opposite_local_charges':[states[0]['X_charge'],states[1]['X_charge']]})
    rng=np.random.default_rng(2781);err=0.
    for _ in range(50):
        b=.225*np.exp(1j*rng.uniform(0,2*np.pi));j=rng.normal()+1j*rng.normal()
        h=np.array([[1+2*b.real,-2*b.imag],[-2*b.imag,1-2*b.real]])
        force=np.array([j.real,j.imag]);u=-np.linalg.solve(h,force)
        direct=u@h@u+2*force@u
        eft=(-abs(j)**2+2*(b*j*j).real)/(1-4*abs(b)**2)
        err=max(err,abs(direct-eft)/max(1,abs(direct)))
    c=read('COUPLING_RELATIONS.json')
    check('complex_heavy_matching',err<1e-13 and abs(2*c['A_GeV_minus2']/c['B_GeV_minus2']-c['r'])<1e-14,
          {'complex_phase_matrix_error':err,'coefficient_ratio':2*c['A_GeV_minus2']/c['B_GeV_minus2']})
    a=read('ABUNDANCE_LEDGER.json');Y=a['YQ_input'];D=a['target_DM_per_visible_entropy_GeV'];m=a['mX_GeV'];eps=a['object_energy_per_charge_GeV']
    lo,hi=a['core_charge_fraction_interval_at_1pct_total_DM_pairs']
    endpoint=[]
    for eta,p in ((lo,0.),(hi,a['primary_pair_ratio_max'])):
        endpoint.append(abs(Y*(eta*eps+(1-eta)*m+2*m*p)/D-1))
    # Independent two-by-two linear solution at allowance saturation.
    matrix=np.array([[eps-m,2*m],[0.,2*m]])
    rhs=np.array([D/Y-m,.01*D/Y])
    solved=np.linalg.solve(matrix,rhs)
    check('fixed_charge_and_cold_mass_ledger',max(endpoint)<1e-13 and abs(solved[0]-hi)<1e-13
          and a['all_charge_in_reference_objects_mass_over_target']<.99,
          {'endpoint_relative_errors':endpoint,'linear_solve_eta_error':float(abs(solved[0]-hi))})
    v=read('VACUUM_FREEDOM.json');ve=max(abs(x['Friedmann_residual_GeV4'])/v['rho_reference_GeV4_input'] for x in v['examples_with_adjusted_initial_geometry'])
    check('vacuum_constant_and_constraint',ve<1e-14 and len(set(x['H_GeV'] for x in v['examples_with_adjusted_initial_geometry']))==3,
          {'relative_Friedmann_error':ve,'allowed_vacuum_offsets':3})
    # Direct integral of the lab speed density divided by speed, not erf identity.
    vs0,vearth,vesc=220.,232.,544.
    nz=erf(vesc/vs0)-2*(vesc/vs0)*np.exp(-(vesc/vs0)**2)/np.sqrt(np.pi)
    def integrand(speed):
        low=np.exp(-((speed-vearth)/vs0)**2)
        high=np.exp(-((speed+vearth)/vs0)**2) if speed<vesc-vearth else np.exp(-(vesc/vs0)**2)
        return (low-high)/(np.sqrt(np.pi)*vs0*vearth*nz)
    emax=0.
    for minimum in (0.,50.,200.,400.,600.,800.):
        direct=quad(integrand,minimum,vesc+vearth,points=[vesc-vearth] if minimum<vesc-vearth else [],epsabs=1e-14)[0] if minimum<vesc+vearth else 0.
        emax=max(emax,abs(float(inverse_speed(minimum))-direct))
    check('halo_integral',emax<1e-12,{'max_inverse_speed_absolute_error_s_per_km':emax})
    # Independent trapezoidal radial transform and denser recoil quadrature.
    raw=np.loadtxt(PROFILE,delimiter=',',skiprows=1);r,F,Fp,S,Sp=raw.T
    g,vphi=.1,.924089637375362
    def form(q,step=1):
        rr=r[::step];ss=S[::step];w=rr*rr*(1-ss*ss)
        arr=np.atleast_1d(q)
        out=[]
        for chunk in np.array_split(arr,max(1,len(arr)//64)):
            out.extend(simpson(w[None,:]*np.sinc(chunk[:,None]*rr[None,:]/(np.pi*g*vphi)),x=rr,axis=1)/simpson(w,x=rr))
        return np.array(out)
    o=read('OBSERVABLE_PREDICTION.json');node=o['first_source_form_factor_zero_GeV']
    decimated=brentq(lambda q:form(q,2)[0],.06,.09,xtol=1e-14)
    trapped=brentq(lambda q:np.trapezoid(r*r*(1-S*S)*np.sinc(q*r/(np.pi*g*vphi)),x=r),.06,.09,xtol=1e-14)
    mass=o['profile_mass_GeV'];edges=np.linspace(5,40,141);gx,gw=np.polynomial.legendre.leggauss(16)
    ee=((edges[:-1]+edges[1:])[:,None]/2+np.diff(edges)[:,None]/2*gx).ravel()
    data=[]
    for label,A in [('Xe131',131),('Ar40',40),('Ge74',74)]:
        ma=A*.93149410242;q=np.sqrt(2*ma*ee*1e-6)
        rn=np.sqrt((1.23*A**(1/3)-.6)**2+7*np.pi**2*.52**2/3-5*.9**2)
        x=q*rn/.1973269804
        nuc=(3*spherical_jn(1,x)/x*np.exp(-.5*(q*.9/.1973269804)**2))**2*(125.2**2/(125.2**2+q*q))**2
        def prob(dm,fo):
            mu=dm*ma/(dm+ma);vm=299792.458*np.sqrt(ma*ee*1e-6/(2*mu*mu))
            rate=nuc*fo*inverse_speed(vm)
            bins=np.sum(rate.reshape(140,16)*gw[None,:],axis=1)*np.diff(edges)/2
            return bins/bins.sum()
        po=prob(mass,form(q)**2)
        # Freeze ma/nuc/ee inside a closure for subsequent independent optimization.
        data.append((po,ma,nuc))
    def objective(logm):
        dm=np.exp(logm);total=0.
        for po,ma,nuc in data:
            mu=dm*ma/(dm+ma);vm=299792.458*np.sqrt(ma*ee*1e-6/(2*mu*mu))
            rate=nuc*inverse_speed(vm)
            bins=np.sum(rate.reshape(140,16)*gw[None,:],axis=1)*np.diff(edges)/2
            if bins.sum()==0:return 1e5
            qp=bins/bins.sum();mask=po>0
            total+=np.sum(po[mask]*np.log(po[mask]/np.maximum(qp[mask],1e-300)))/3
        return float(total)
    best=minimize_scalar(objective,bounds=(np.log(22.5),np.log(23.5)),method='bounded',options={'xatol':1e-12})
    delta_mass=abs(np.exp(best.x)/o['primary']['common_best_point_mass_GeV']-1)
    delta_kl=abs(best.fun/o['primary']['equal_target_weighted_KL_per_signal_event']-1)
    check('profile_node_and_recoil_quadrature',abs(decimated/node-1)<1e-6 and abs(trapped/node-1)<1e-6
          and delta_mass<2e-4 and delta_kl<2e-3,
          {'decimated_node_relative_change':float(abs(decimated/node-1)),
           'trapezoid_node_relative_change':float(abs(trapped/node-1)),
           'order8_vs16_best_mass_relative_change':float(delta_mass),'order8_vs16_KL_relative_change':float(delta_kl)})
    save('VERIFICATION.json',{'same_author':True,'groups':checks,'all_pass':all(x['pass'] for x in checks)})
    print(json.dumps({'groups':len(checks),'all_pass':all(x['pass'] for x in checks),'failed':[x['name'] for x in checks if not x['pass']]}))
    if not all(x['pass'] for x in checks):raise SystemExit(1)


if __name__=='__main__': main()
