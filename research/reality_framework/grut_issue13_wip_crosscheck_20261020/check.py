#!/usr/bin/env python3
"""Small independent algebra/ODE cross-checks of Claude WIP 559975f.

No Claude computation module is imported. This does not rerun the seven
modules, their pending verifiers, or the new charged-exit lattice evolution.
The constant-g wall/emission ODE is an illustrative consistency diagnostic,
not a corrected Candidate-1 history or an observational prediction.
"""
from pathlib import Path
import hashlib
import itertools
import json
import math
import platform
import sys

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq

P=Path(__file__).resolve().parent
checks={}

def test(name,value):checks[name]=bool(value)

def charge():
    # e(q,S,x)=q^2/(2x)+(S^2-1)^2/4+a*S^2*x/2+b*x^2/4,
    # q=n/(sqrt(lambda_phi)*v^3), x=2|X|^2/v^2.
    a,b,v,lp=1.,.1,1000.,.01
    nc=math.sqrt(b/a**3)
    xs=2*a/(3*(a*a-b)); ns=xs*math.sqrt(a/3)
    def x_broken(q):
        return brentq(lambda x:x*math.sqrt(a-(a*a-b)*x)-q,0,xs)
    def delta(q):
        x=x_broken(q)
        return a*x-.75*(a*a-b)*x*x - (.25+.75*b**(1/3)*q**(4/3))
    nt=brentq(delta,nc,ns)
    nstar=b**(-.25)
    unit=math.sqrt(lp)*v**3
    YB=.02233*1.05375e-5/.938272/2891.2
    entropy_coefficient=2*math.pi**2/45*106.75
    names={'restored_local_spinodal':nc,'homogeneous_global_switch':nt,
           'broken_radial_spinodal':ns,'phase_coexistence_dense_endpoint':nstar}
    table={k:{'n_over_v3':q*math.sqrt(lp),'n_GeV3':q*unit,
              'T_visible_GeV_r1':(q*unit/(YB*entropy_coefficient))**(1/3)} for k,q in names.items()}
    # Integrating out x locally at S=0 gives this coefficient of S^4.
    quartic=.25-a*a/(12*b)
    q=nc
    def local_energy(S):
        x=brentq(lambda x:-q*q/(2*x*x)+a*S*S/2+b*x/2,1e-5,5)
        return q*q/(2*x)+(S*S-1)**2/4+a*S*S*x/2+b*x*x/4
    numerical=[]
    for h in [.02,.01,.005]:
        numerical.append((local_energy(h)-local_energy(0))/h**4)
    H=math.sqrt(math.pi**2*106.75/90)*table['restored_local_spinodal']['T_visible_GeV_r1']**2/2.435323e18
    test('charge_threshold_order',nc<nt<ns<nstar)
    test('homogeneous_switch_equal_energy',abs(delta(nt))<1e-11)
    test('eliminated_quartic_negative',quartic<0)
    test('eliminated_quartic_finite_difference',abs(numerical[-1]-quartic)<2e-4)
    return {'units':'v=1000 GeV, lambda_phi=.01, a=1,b=.1; r=1; same stated baryon conversion as WIP',
            'YB':YB,'thresholds':table,'effective_S4_coefficient_at_nc':quartic,
            'S4_finite_differences':numerical,'H_over_mX_at_nc':H/100,
            'interpretation':'198 TeV is a local spinodal, not a universal forced-restoration threshold or a necessity to supply charge late. Inhomogeneous coexistence and initial state still matter.'}

def correlations():
    # Four cells on a periodic loop. Z2 invariance and one-site marginals
    # do not fix disagreement probabilities or a macroscopic network.
    def stats(states):
        return {'p_plus_each_cell':[sum(p*(s[i]>0) for s,p in states) for i in range(4)],
                'mean_neighbor_disagreement':sum(p*sum(s[i]!=s[(i+1)%4] for i in range(4))/4 for s,p in states)}
    coherent=[((1,1,1,1),.5),((-1,-1,-1,-1),.5)]
    independent=[(s,1/16) for s in itertools.product([-1,1],repeat=4)]
    one,two=stats(coherent),stats(independent)
    test('identical_one_site_sign_marginals',one['p_plus_each_cell']==two['p_plus_each_cell']==[.5]*4)
    test('different_interface_statistics',one['mean_neighbor_disagreement']==0 and two['mean_neighbor_disagreement']==.5)
    return {'globally_correlated_symmetric_mixture':one,'independent_symmetric_signs':two,
            'interpretation':'Both distributions are Z2 invariant. Equal signs at one site do not establish short spatial correlations or percolation. This is a logical counterexample, not a proposal for the actual primordial ensemble.'}

def emission_background(v2,w_emit,rtol):
    # tau=t/t_ref, t_ref=3 M_P^2/(4 A sigma), rho_ref=A sigma/t_ref.
    # rho_w/rho_ref=1/tau, R=(rho_visible/rho_ref)*tau^2,
    # z=rho_emitted/rho_w, Ht=.5 sqrt(R+tau*(1+z)). Constant g*.
    # rho_visible is conserved radiation; products have a chosen constant w.
    def rhs(x,y):
        tau=math.exp(x);R=math.exp(y[0]);z=y[1]
        ht=.5*math.sqrt(R+tau*(1+z))
        q=1-(1+3*v2)*ht
        return [2-4*ht,q+(1-3*(1+w_emit)*ht)*z]
    def event(x,y):return x-y[0]  # rho_wall=rho_visible
    event.terminal=True;event.direction=1
    sol=solve_ivp(rhs,[math.log(1e-10),math.log(2)],[0,0],method='DOP853',
                  rtol=rtol,atol=rtol*1e-2,events=event,dense_output=True)
    if not sol.success or not len(sol.t_events[0]):raise RuntimeError('Toy equality not reached')
    x=float(sol.t_events[0][0]);lnR,z=sol.sol(x);tau=math.exp(x)
    grid=np.linspace(sol.t[0],x,400);ys=sol.sol(grid)
    ht=.5*np.sqrt(np.exp(ys[0])+np.exp(grid)*(1+ys[1]))
    return {'w_emit':w_emit,'tau_wall_equals_visible':tau,
            'rho_emitted_over_rho_wall_at_equality':float(z),
            'T_wall_equality_over_no_backreaction_reference':math.exp((lnR-2*x)/4),
            'minimum_Q_t_over_rho_wall':float(np.min(1-(1+3*v2)*ht)),
            'minimum_emitted_ratio':float(np.min(ys[1]))}

def stress_matching():
    v=.456113137775563  # WIP network_sim headline; a velocity estimator, not a new fit.
    v2=v*v
    w=v2-2/3
    q=(1-3*v2)/2  # Q=q*rho_wall/t in exact radiation domination, if rho_wall=A sigma/t.
    cases=[]
    for we in [0.,1/3]:
        loose=emission_background(v2,we,1e-8)
        tight=emission_background(v2,we,1e-11)
        err=abs(loose['tau_wall_equals_visible']/tight['tau_wall_equals_visible']-1)
        test(f'emission_ODE_refinement_w{we}',err<1e-7)
        test(f'emission_nonnegative_to_equality_w{we}',tight['minimum_Q_t_over_rho_wall']>=-1e-8 and tight['minimum_emitted_ratio']>=-1e-10)
        cases.append(tight)
    test('zero_loss_RD_velocity_control',abs((1-3*(1/math.sqrt(3))**2)/2)<1e-15)
    test('quoted_velocity_requires_loss_if_scaling_energy',q>0)
    return {'input_velocity_estimator':v,'thin_wall_energy_weighted_w_if_velocity_matches':w,
            'isolated_fluid_w_inferred_from_RD_scaling':-1/3,
            'required_Q_t_over_rho_wall':q,'required_Q_over_H_rho_wall':2*q,
            'mono_speed_boost_gamma':1/math.sqrt(1-v2),
            'RD_asymptotic_rho_massive_products_over_rho_wall':2*q,
            'RD_asymptotic_rho_radiation_products_over_rho_wall':q,
            'illustrative_constant_g_source_ODE':cases,
            'interpretation':'Effective dilution w and intrinsic wall pressure differ when the network loses energy. Areal A and energy-density A also require calibration; gamma is a mono-speed control, not a claimed 12% correction to the measured network. ODE cases illustrate bookkeeping only and do not recompute 6.367 keV.'}

def planar_modes():
    # u=z/delta, delta=sqrt(2/lambda_phi)/v. Each small-amplitude mode has
    # [-d_u^2-l(l+1)sech(u)^2] psi = (delta*p_normal)^2 psi.
    # For the real phi fluctuation l=2; complex X l=1; D is noninteger.
    v,lp,g2,h2=1000.,.01,.01,.0001
    delta=math.sqrt(2/lp)/v
    ell_D=(math.sqrt(1+8*h2/lp)-1)/2
    def reflection(ell,k,rtol=1e-10):
        L=22.;start=np.exp(1j*k*L)
        def rhs(x,y):return [y[1],-(k*k+ell*(ell+1)/math.cosh(x)**2)*y[0]]
        sol=solve_ivp(rhs,[L,-L],[start,1j*k*start],method='DOP853',rtol=rtol,atol=rtol*.01)
        if not sol.success:raise RuntimeError('Scattering ODE failed')
        psi,dpsi=sol.y[:,-1];x=-L
        inc=(psi+dpsi/(1j*k))/(2*np.exp(1j*k*x))
        ref=(psi-dpsi/(1j*k))/(2*np.exp(-1j*k*x))
        R=abs(ref/inc)**2;T=1/abs(inc)**2
        exact=math.sin(math.pi*ell)**2/(math.sinh(math.pi*k)**2+math.sin(math.pi*ell)**2)
        return {'ell':ell,'k_delta':k,'R_numeric':R,'T_numeric':T,'R_exact':exact,'flux_error':abs(R+T-1)}
    rows=[]
    for name,ell in [('phi',2.),('X',1.),('D',ell_D)]:
        for k in [.05,.2,1.]:
            row=reflection(ell,k);row['field']=name;rows.append(row)
            test(f'{name}_scattering_flux_k{k}',row['flux_error']<1e-8)
            test(f'{name}_reflection_formula_k{k}',abs(row['R_numeric']-row['R_exact'])<1e-8)
    eig=[]
    for h in [.04,.02]:
        L=16.;u=np.arange(-L+h,L,h)
        X=eigh_tridiagonal(2/h**2-2/np.cosh(u)**2,np.full(len(u)-1,-1/h**2),select='i',select_range=(0,1))[0]
        phi=eigh_tridiagonal(2/h**2-6/np.cosh(u)**2,np.full(len(u)-1,-1/h**2),select='i',select_range=(0,1))[0]
        eig.append({'du':h,'X_eigenvalues_unshifted':X.tolist(),'phi_eigenvalues_unshifted':phi.tolist()})
    test('X_bound_eigenvalue_converges_to_minus1',abs(eig[1]['X_eigenvalues_unshifted'][0]+1)<abs(eig[0]['X_eigenvalues_unshifted'][0]+1)<1e-3)
    test('phi_bound_eigenvalues_converge',all(abs(eig[1]['phi_eigenvalues_unshifted'][i]-target)<abs(eig[0]['phi_eigenvalues_unshifted'][i]-target)<1e-3 for i,target in enumerate([-4,-1])))
    mX_wall=math.sqrt(g2*v*v-1/delta**2)
    test('X_wall_mode_below_bulk_threshold',0<mX_wall<math.sqrt(g2)*v)
    return {'scope':'Leading classical uncharged static planar wall, weak linear probes. This is not a finite charged wall, a capture/production calculation, a lifetime or an ordinary detector signal.',
            'delta_GeV_inverse':delta,'ell_phi':2.,'ell_X':1.,'ell_D':ell_D,
            'X_surface_mode_mass_GeV':mX_wall,'X_bulk_mass_GeV':math.sqrt(g2)*v,
            'phi_translation_mode_mass_squared_GeV2':0.,
            'phi_shape_mode_mass_GeV':math.sqrt(3)/delta,
            'scattering':rows,'eigenvalue_refinement':eig,
            'interpretation':'At the exact tree benchmark, phi/X linear continuum scattering is reflectionless but X has a localized bound mode. This supplies a check for wall-scattering and charge-storage calculations, not a capture rate or a proof of zero friction in a curved interacting network. Running, finite amplitudes, charge and curvature can change it.'}

def main():
    result={'reviewed_commit':'559975fad327672f70221e62f6c11c913f7d30b6',
            'status':'preliminary, retrospective cross-check; pending verification/charged-exit results not incorporated',
            'charge':charge(),'sign_correlations':correlations(),'wall_stress_matching':stress_matching(),
            'planar_wall_modes':planar_modes(),
            'checks':checks,'all_checks_pass':all(checks.values()),
            'environment':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform()},
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (P/'RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))
    if not result['all_checks_pass']:raise SystemExit(1)

if __name__=='__main__':main()
