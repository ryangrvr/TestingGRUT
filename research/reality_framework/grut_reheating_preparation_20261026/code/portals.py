"""Portal matching and source ledgers; not a full relic or detector fit.

The vector source model includes SM inverse decays, A decay back to the SM,
and decay into X pairs, with MB statistics and fixed masses in hot RD.
Inverse X production and hidden rescattering require a separate plasma solver.
"""
import math,json
import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.special import kve
from common import ROOT,write,previous

P=previous();pars=P['normalized_parameters'];MX=P['light_sector']['mX_GeV'];V=P['light_sector']['vPhi_GeV']
H0=pars['H_GeV'];MRED=2.435e18;GS=106.75;MASS=6*1e-7*pars['h_GeV'];GY=.37
HM=math.sqrt(math.pi**2*GS/90)*MASS*MASS/MRED
WIDTH_X=1e-14*MASS/(48*math.pi)*(1-4*MX*MX/MASS**2)**1.5
CEQ=45*3/(4*math.pi**4*GS);KNUM=CEQ*3*math.pi/2
PREV=json.loads((ROOT.parent/'grut_heavy_vector_budget_20261025/BUDGET_RESULTS.json').read_text())['rows'][0]
NX=abs(PREV['N_A_band']/PREV['vectors_per_net_charge'])
EX=json.loads((ROOT.parent/'grut_heavy_vector_budget_20261025/BACKGROUND_RESULTS.json').read_text())['late']['E_X']
RFAC=4*(MRED/H0)**2/EX


def number_energy(eps,rtol=2e-9):
    gsm=41*eps*eps*GY*GY*MASS/(96*math.pi)
    total=gsm+WIDTH_X;bx=WIDTH_X/total
    expected=KNUM*WIDTH_X/HM*(1-bx)
    prompt=5*expected
    if eps==0:return dict(epsilon_BX=eps,Gamma_SM_GeV=0.,branch_X=1.,Y_pair=0.,J=0.,R_thermal=0.,R_primordial=PREV['R_vector_band_only'],R_total=PREV['R_vector_band_only'],number_analytic=0.,number_relative_error=0.,energy_over_prompt_bound=None,nfev=0)
    ymax=45*3/(2*math.pi**4*GS)
    scale=min(ymax,KNUM*gsm/HM)
    def fun(u,y):
        x=math.exp(u);ratio=kve(1,x)/kve(2,x)
        yeq=CEQ*x*x*kve(2,x)*math.exp(-x) if x<745 else 0.
        d=gsm/HM*x*x*ratio;dt=total/HM*x*x*ratio
        return [d*yeq/scale-dt*y[0],WIDTH_X/HM*x*x*ratio*y[0],WIDTH_X/HM*x**3*y[0]]
    xmax=max(100.,30*math.sqrt(HM/total))
    sol=solve_ivp(fun,(math.log(1e-4),math.log(xmax)),[0.,0.,0.],method='Radau',rtol=rtol,atol=rtol*1e-3)
    assert sol.success,sol.message
    ya,yp,j=sol.y[:,-1]*scale
    rprim=PREV['R_vector_band_only']*bx**1.5
    return dict(epsilon_BX=eps,Gamma_SM_GeV=gsm,branch_X=bx,Y_pair=yp,J=j,
                R_thermal=j*RFAC,R_primordial=rprim,R_total=j*RFAC+rprim,
                number_analytic=expected,number_relative_error=abs(yp/expected-1),
                energy_over_prompt_bound=j/prompt,prompt_energy_bound_R=prompt*RFAC,
                f_pair_rest_without_processing=2*MX*yp/P['benchmark']['DM_per_visible_entropy_GeV'],
                remaining_A=ya,nfev=sol.nfev,x_end=xmax,
                negative_A_tail_over_scale=max(0.,float(-sol.y[0].min())),
                hidden_source_energy_over_SM_radiation=4*j/3)


def main():
    epsgrid=[0.,1e-12,3e-12,1e-11,3e-11,1e-10,3e-10,1e-9,3e-9,1e-8,3e-8,1e-7,3e-7,1e-6,3e-6,1e-5,1e-4,1e-3,1e-2]
    rows=[]
    for eps in epsgrid:
        rows.append(number_energy(eps));print('PORTAL_SOURCE',eps,rows[-1]['R_total'],rows[-1]['Y_pair'],flush=True)
    tight=number_energy(1e-7,rtol=2e-10)
    expected_sat=KNUM*WIDTH_X/HM;Rsat=5*expected_sat*RFAC
    q3=quad(lambda x:x**3*kve(1,x)*math.exp(-x),0,700,epsabs=1e-11)[0]
    q4=quad(lambda x:x**4*kve(2,x)*math.exp(-x),0,700,epsabs=1e-11)[0]
    mh,vh,width_h=125.,246.,.0041
    ms=math.sqrt(2*.03)*V
    BMAX=.107;extra_max=width_h*BMAX/(1-BMAX)
    gamma_ss_per_b2=vh*vh/(8*math.pi*mh)*math.sqrt(1-4*ms*ms/mh**2)
    gamma_xx_per_lam2=vh*vh/(16*math.pi*mh)*math.sqrt(1-4*MX*MX/mh**2)
    bmax=math.sqrt(extra_max/gamma_ss_per_b2);lmax=math.sqrt(extra_max/gamma_xx_per_lam2)
    b=5e-11 # equals earlier illustrative kappa_Hphi=1e-10.
    gamma_ss=b*b*gamma_ss_per_b2
    Hh=math.sqrt(math.pi**2*GS/90)*mh*mh/MRED
    # One thermal h -> two s -> two X-pairs; cold-mass/constant-g control.
    kh=135/(8*math.pi**3*GS)
    yp_h=2*kh*gamma_ss/Hh
    higgs=dict(mh_GeV=mh,vh_GeV=vh,SM_width_GeV=width_h,archived_ATLAS_Binv_bound=BMAX,
               experimental_assumptions='SM Higgs production; added decays invisible; archived ATLAS2301.10731, not a global/current exclusion.',
               max_b_Hphi_one_channel=bmax,max_lambda_HX_one_channel=lmax,
               mixing_angle_per_b=2*V*vh/(mh*mh-ms*ms),
               scalar_portal_reference_b=b,reference_h_to_ss_width=gamma_ss,
               reference_vacuum_mass_Higgs_source_Y_pair=yp_h,
               reference_unprocessed_pair_mass_fraction=2*MX*yp_h/P['benchmark']['DM_per_visible_entropy_GeV'],
               source_scope='Vacuum masses, constant g*, thermal h, open s->X Xbar and no erasure. In-medium masses/VEVs can change the source; not a relic exclusion.')
    # Fixed low-energy contact benchmarks, no detector threshold/halo fit.
    eEM=math.sqrt(4*math.pi/137.035999084);mn=.938272;Q=200000
    obj=Q*MX*P['benchmark']['original_E_over_Qm'];conv=.3893793721e-27;av=6.02214076e23
    targets=[]
    for name,A,Z in [('Xe131',131,54),('Ar40',40,18),('Ge74',74,32)]:
        ma=A*.931494;red=obj*ma/(obj+ma)
        G=Q*Z*1e-7*eEM/MASS**2 # epsilon_gamma=1 reference envelope
        sigma=red*red/math.pi*G*G*conv
        flux=.4/obj*8e7;atoms=1e6/A*av;events=flux*atoms*sigma*365.25*86400
        targets.append(dict(target=name,sigma_zero_q_cm2_at_epsilon_gamma1=sigma,
                            events_per_tonne_year_at_fixed_flux_zero_q_envelope=events))
    fn=.30;redp=obj*mn/(obj+mn)
    Cscalar=b*obj*mn*fn/(.03*V*V*mh*mh)
    scalar_sigma=redp*redp/math.pi*Cscalar*Cscalar*conv
    # Reapply included cold rates, no new numerical quantum-rate integration.
    M=pars['M_Y_GeV'];mu=pars['mu_GeV'];h=pars['h_GeV'];k=pars['kappa']
    r=math.sqrt(2)*mu*h/(M*M);alpha=k*k*mu*h/(math.sqrt(2)*(M**4-2*mu*mu*h*h));c6=alpha*MRED*MRED
    cold=json.loads((ROOT.parent/'grut_gauge_point_test_20261025/PAIR_RESULTS.json').read_text())['summaries'][1]
    age=13.8e9*365.25*86400
    leakage=dict(alpha6_GeVm2=alpha,c6=c6,r=r,exact_over_leading=1/(1-r*r),
                 included_cold_core_Gamma_times_age=c6*c6*cold['included_channels_rate_per_second']*age,
                 pair_to_coherent=cold['pair_to_coherent'],scope='Reapplication of previously computed included cold coherent/vacuum-Born pair channels; no new full in-medium lifetime.')
    # Independent mixed derivative of the two-real-X scalar mass trace.
    # V contains b_Hphi phi^2 |H|^2 and lambda_HX |H|^2|X|^2.
    gx2=.01;lam=.002;step=.001
    tr=lambda u,z:2*(.008+gx2*u+lam*z)**2
    cross=(tr(1+step,2+step)-tr(1+step,2-step)-tr(1-step,2+step)+tr(1-step,2-step))/(4*step*step)
    cross_err=abs(cross/(4*gx2*lam)-1)
    # h(x)=x K2/K1 grows with x; its monotone likelihood-ratio proof is
    # recorded in PORTALS.md. Check the displayed integration range too.
    xx=np.geomspace(1e-5,1e5,2000);weight=xx*kve(2,xx)/kve(1,xx)
    groups={
      'inverse_decay_source_moments':dict(passed=max(r['number_relative_error'] for r in rows)<1e-5 and min(r['energy_over_prompt_bound'] for r in rows[1:])>=1-1e-5 and max(r.get('negative_A_tail_over_scale',0) for r in rows)<1e-8,
          max_number_error=max(r['number_relative_error'] for r in rows),min_energy_bound_ratio=min(r['energy_over_prompt_bound'] for r in rows[1:])),
      'source_quadrature_refinement_and_lower_bound':dict(passed=abs(q3/(3*math.pi/2)-1)<1e-9 and abs(q4/(15*math.pi/2)-1)<1e-9 and abs(tight['R_total']/rows[11]['R_total']-1)<1e-5 and Rsat>1.5*PREV['R_vector_band_only'],
          number_Bessel_integral_error=abs(q3/(3*math.pi/2)-1),energy_Bessel_integral_error=abs(q4/(15*math.pi/2)-1),refinement_R_change=abs(tight['R_total']/rows[11]['R_total']-1)),
      'portal_width_and_operator_normalizations':dict(passed=abs(gamma_ss_per_b2/(2*gamma_xx_per_lam2)-1)<1e-4 and cross_err<1e-7 and np.all(np.diff(weight)>0),
          one_loop_mixing_charge_trace=0.,isolated_beta_b_Hphi_over_lambda_HX=.01/(8*math.pi**2),
          mixed_scalar_trace_derivative_error=cross_err,minimum_energy_weight_increment=float(np.diff(weight).min()),
          Higgs_scalar_vs_complex_scalar_width_factor=gamma_ss_per_b2/gamma_xx_per_lam2)
    }
    write('PORTAL_RESULTS.json',dict(source_rows=rows,refined=tight,groups=groups,inputs=dict(gY=GY,gstar=GS,mA_GeV=MASS,Gamma_X_GeV=WIDTH_X,H_at_mA_GeV=HM),
          analytic_bounds=dict(saturated_pair_yield=expected_sat,saturated_prompt_R=Rsat,
               total_R_lower_bound='R_primordial*B_X^(3/2)+R_sat*(1-B_X)',minimum_lower_bound_R=PREV['R_vector_band_only'],
               conclusion='Under the hot-RD/MB source assumptions, opening the SM decay cannot lower the raw hidden vector-energy injection below the original primordial band. This is not a final relic bound.'),
          higgs=higgs,vector_scattering=dict(targets=targets,parameters='rho=.4 GeV/cm3, speed800km/s, pointlike zero-q coherent form factors; epsilon_gamma1 is an upper-reference normalization, rates scale as epsilon_gamma^2.',
                    scope='One heavy-photon contact operator only, not a detector likelihood or all mediator channels.'),
          scalar_zero_q_reference=dict(b=b,nucleon_scalar_fraction=fn,object_mass_GeV=obj,proton_sigma_cm2=scalar_sigma,
                    scope='Prior tree soliton-Higgs identity and unmodified profile parameters, fixed hadronic input. No q spectrum or exclusion fit.'),leakage=leakage,scope=__doc__))
    print('PORTAL_BOUND',Rsat,'HIGGS',higgs,'VECTOR_EVENTS',[r['events_per_tonne_year_at_fixed_flux_zero_q_envelope'] for r in targets],flush=True)
    if not all(g['passed'] for g in groups.values()):raise SystemExit(1)


if __name__=='__main__':main()
