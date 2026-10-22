#!/usr/bin/env python3
"""Fixed-action matter, weak gravity and readable charge labels; no old imports."""
import csv
import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np
import scipy
from scipy.integrate import solve_bvp, solve_ivp, simpson, cumulative_trapezoid
from scipy.sparse import diags, bmat
from scipy.sparse.linalg import eigsh, LinearOperator
from scipy.special import spherical_jn, spherical_yn
from scipy.stats import binom

ROOT=Path(__file__).resolve().parents[1]
I=json.loads((ROOT/'INPUTS.json').read_text())
g=I['g_X']; aa=I['lambda_phi']/g**2; bb=I['lambda_X']/g**2
aD=(I['h_D']/g)**2
Gv2=I['v_over_reduced_M']**2/(8*math.pi)
CHECKS=[]


def check(name,passed,value=None,limit=None):
    CHECKS.append(dict(name=name,passed=bool(passed),value=value,limit=limit))


class Seed:
    def __init__(self,path):
        d=np.load(path)
        self.x=d['rho']; self.y=d['y']; self.p=np.array([float(d['Omega'])])
        self.integral=cumulative_trapezoid(self.x**2*self.y[0]**2,self.x,initial=0)

    def sol(self,r):
        rows=[np.interp(r,self.x,row) for row in self.y]
        rows.append(np.interp(r,self.x,self.integral))
        return np.array(rows)


def solve_profile(Q,L,tol,seed):
    def rhs(r,y,p):
        F,Fp,S,Sp,J=y; om=p[0]
        return np.array([Fp,(S*S+bb*F*F-om*om)*F,Sp,
                         (aa*(S*S-1)+F*F)*S,r*r*F*F])

    def boundary(left,right,p):
        return [left[1],left[3],right[0],right[2]-1,left[4],4*math.pi*p[0]*right[4]-g*g*Q]

    r=np.linspace(0,L,701)
    guess=seed.sol(r)
    outside=r>seed.x[-1]
    if np.any(outside):
        guess[:,outside]=np.array([0,0,1,0,float(seed.sol(np.array([seed.x[-1]]))[4,0])])[:,None]
    sol=solve_bvp(rhs,boundary,r,guess,p=seed.p,S=np.diag([0.,-2.,0.,-2.,0.]),
                  tol=tol,bc_tol=tol,max_nodes=30000)
    if not sol.success:
        raise RuntimeError(f'Q={Q},L={L}: '+sol.message)
    test=sol.sol(np.linspace(0,L,4001))
    if not (0<sol.p[0]<1 and test[0].min()>-1e-7 and test[2].min()>0 and test[2].max()<1.000001):
        raise RuntimeError(f'Non-admissible prescribed branch at Q={Q}')
    return sol


def continuation(Q0,Q1,seed,L,tol):
    if Q0==Q1:
        return solve_profile(Q1,L,tol,seed)
    count=max(1,math.ceil(abs(math.log(Q1/Q0))/.1))
    for Q in np.geomspace(Q0,Q1,count+1)[1:]:
        seed=solve_profile(float(Q),L,tol,seed)
    return seed


def profile_quantities(sol,Q):
    r=np.linspace(0,sol.x[-1],24001)
    F,Fp,S,Sp,_=sol.sol(r); om=float(sol.p[0])
    U=aa/4*(S*S-1)**2+.5*S*S*F*F+bb/4*F**4
    grad=.5*(Fp*Fp+Sp*Sp); time=.5*om*om*F*F
    energy=grad+time+U; pr=grad+time-U; pt=-grad+time-U
    integrate=lambda x:float(4*math.pi*simpson(r*r*x,x=r))
    Ebar=integrate(energy); charge=om*integrate(F*F)/g**2
    mass_curve=4*math.pi*cumulative_trapezoid(r*r*energy,r,initial=0)
    active=energy+pr+2*pt
    R90=float(np.interp(.9*mass_curve[-1],mass_curve,r))
    V=aD*(S*S-1)
    phi=-Gv2*4*math.pi*simpson(r*active,x=r)
    virial=integrate(grad)+3*(integrate(U)-integrate(time))
    # A globally bounded source interpolation is used only beyond the actual BVP endpoint.
    def potential(rr):
        rr=np.asarray(rr)
        ss=sol.sol(np.minimum(rr,sol.x[-1]))[2]
        return np.where(rr<=sol.x[-1],aD*(ss*ss-1),0.)

    def zero_rhs(rr,y):
        return [y[1],float(potential(rr))*y[0]]

    zero=solve_ivp(zero_rhs,(0,sol.x[-1]),[0.,1.],method='DOP853',
                   rtol=1e-12,atol=1e-13,max_step=.06)
    if not zero.success:
        raise RuntimeError(zero.message)
    scattering_length=float(sol.x[-1]-zero.y[0,-1]/zero.y[1,-1])
    born_length=float(simpson(r*r*V,x=r))
    no_bound_integral=float(simpson(r*np.maximum(0.,-V),x=r))
    depletion=integrate(1-S*S)
    F_equation=integrate(Fp*Fp+S*S*F*F+bb*F**4-om*om*F*F)
    S_equation=integrate(Sp*Sp+aa*S*S*(S*S-1)+S*S*F*F)
    return dict(Q=Q,computed_charge=charge,Omega=om,F0=float(F[0]),S0=float(S[0]),
                Ebar=Ebar,mass_over_gv=Ebar/g**2,energy_per_charge_over_gv=Ebar/(g*g*Q),
                R90_over_inverse_gv=R90,central_clock_shift=float(-phi),
                schwarzschild_radius_over_R90=2*Gv2*Ebar/R90,
                gravitational_Coulomb_parameter_at_min_k=Gv2*Ebar*(I['m_D0_over_gv']**2+aD+2*I['probe_k_min']**2)/I['probe_k_min'],
                relative_virial=float(virial/Ebar),active_mass_over_energy=integrate(active)/Ebar,
                scattering_length_over_inverse_gv=scattering_length,
                born_scattering_length_over_inverse_gv=born_length,
                no_reporter_bound_state_sufficient_integral=no_bound_integral,
                born_mass_relation_ratio=-4*math.pi*I['lambda_phi']*born_length/(I['h_D']**2*Ebar),
                exact_mass_relation_ratio=-4*math.pi*I['lambda_phi']*scattering_length/(I['h_D']**2*Ebar),
                mass_depletion_relative_residual=(Ebar-aa*depletion)/Ebar,
                integrated_F_equation_relative_residual=F_equation/Ebar,
                integrated_S_equation_relative_residual=S_equation/Ebar,
                exact_scattering_length_lower_bound=born_length/(1-no_bound_integral),
                exact_scattering_length_upper_bound=born_length,
                max_BVP_rms_residual=float(np.max(sol.rms_residuals)))


def phases(sol,ks,refined=False,zero=False):
    ks=np.asarray(ks); ll=np.array(I['probe_partial_waves'])[:,None]
    kk=np.broadcast_to(ks,(len(ll),len(ks)))
    ell=np.broadcast_to(ll,kk.shape); size=kk.size; r0=1e-3

    def potential(r):
        if zero or r>sol.x[-1]:
            return 0.
        S=float(sol.sol(r)[2]); return aD*(S*S-1)

    V0=potential(0)
    c2=(V0-kk*kk)/(4*ell+6)
    f0=1+c2*r0*r0; fp0=2*c2*r0
    dfactor=np.where(ell==0,1.,9.)
    delta0=-V0*kk**(2*ell+1)*r0**(2*ell+3)/((2*ell+3)*dfactor)

    def rhs(r,y):
        f=y[:size].reshape(kk.shape); fp=y[size:2*size].reshape(kk.shape)
        d=y[2*size:].reshape(kk.shape); V=potential(r)
        z=kk*r
        J=z*spherical_jn(ell,z); Y=z*spherical_yn(ell,z)
        phase=-V/kk*(np.cos(d)*J-np.sin(d)*Y)**2
        return np.concatenate([fp.ravel(),(-2*(ell+1)/r*fp+(V-kk*kk)*f).ravel(),phase.ravel()])

    ini=np.concatenate([f0.ravel(),fp0.ravel(),delta0.ravel()])
    tol=1e-12 if refined else 2e-10
    out=solve_ivp(rhs,(r0,sol.x[-1]),ini,method='DOP853',rtol=tol,atol=tol*.01,
                  max_step=.06 if refined else .15)
    if not out.success:
        raise RuntimeError(out.message)
    f=out.y[:size,-1].reshape(kk.shape); fp=out.y[size:2*size,-1].reshape(kk.shape)
    delta=out.y[2*size:,-1].reshape(kk.shape)
    rr=sol.x[-1]; z=kk*rr
    J=z*spherical_jn(ell,z); Y=z*spherical_yn(ell,z)
    Jp=kk*(spherical_jn(ell,z)+z*spherical_jn(ell,z,derivative=True))
    Yp=kk*(spherical_yn(ell,z)+z*spherical_yn(ell,z,derivative=True))
    u=rr**(ell+1)*f; up=rr**ell*((ell+1)*f+rr*fp)
    wave=np.arctan2(u*Jp-up*J,u*Yp-up*Y)
    difference=(wave-delta+math.pi/2)%math.pi-math.pi/2
    return delta,dict(max_wave_phase_difference=float(np.max(abs(difference))),
                      evaluation_count=out.nfev)


def hessian(sol,N):
    L=sol.x[-1]; h=L/(N+1); r=h*np.arange(1,N+1)
    F,_,S,_,_=sol.sol(r); om=float(sol.p[0])
    off=-np.ones(N-1)/h**2
    A=diags([off,2/h**2+S*S+3*bb*F*F-om*om,off],[-1,0,1])
    B=diags([off,2/h**2+aa*(3*S*S-1)+F*F,off],[-1,0,1])
    C=diags(2*F*S)
    H=bmat([[A,C],[C,B]],format='csr')
    rg=np.linspace(0,L,24001); fg=sol.sol(rg)[0]
    integral=4*math.pi*simpson(rg*rg*fg*fg,x=rg)
    q=np.concatenate([math.sqrt(4*math.pi)*r*F,np.zeros(N)])
    factor=4*om*om/integral
    op=LinearOperator(H.shape,matvec=lambda v:H@v+factor*h*q*np.dot(q,v),dtype=float)
    lowest=float(eigsh(op,k=1,which='SA',tol=1e-10,maxiter=30000,
                      v0=np.linspace(.5,1,2*N),return_eigenvectors=False)[0])
    return dict(N=N,h=h,lowest_fixed_charge_radial_amplitude_eigenvalue=lowest)


def binomial_success(N,p0,p1):
    if p0<p1:
        p0,p1=p1,p0
    cut=math.ceil(N*math.log((1-p1)/(1-p0))/(math.log(p0/p1)+math.log((1-p1)/(1-p0))))
    return float(1-.5*(binom.cdf(cut-1,N,p0)+binom.sf(cut-1,N,p1)))


def readout(ks,d0,d1,smallmass):
    w=np.exp(-.5*((ks-I['packet_k_center'])/I['packet_k_sigma'])**2)
    w/=simpson(w,x=ks)
    overlap=simpson(w*.5*np.sum(np.exp(2j*(d1-d0)),axis=0),x=ks)
    c2=float(abs(overlap)**2)
    # Stable formulas when the outgoing states are close.
    logc2=math.log(c2)
    success=lambda n: .5*(1+math.sqrt(-math.expm1(n*logc2)))
    goal=I['readout_success_target']
    Nopt=math.ceil(math.log(1-(2*goal-1)**2)/logc2)
    theta0=d0[0]-d0[1]; theta1=d1[0]-d1[1]
    mid0=float(np.interp(I['packet_k_center'],ks,theta0))
    mid1=float(np.interp(I['packet_k_center'],ks,theta1))
    beta=mid0+mid1+math.pi/2
    p0=float(simpson(w*.5*(1+np.cos(2*theta0-beta)),x=ks))
    p1=float(simpson(w*.5*(1+np.cos(2*theta1-beta)),x=ks))
    lo,hi=0,1
    while binomial_success(hi,p0,p1)<goal:
        hi*=2
        if hi>10**9:
            raise RuntimeError('Declared readout probe-count bound exceeded')
    while hi-lo>1:
        mid=(lo+hi)//2
        if binomial_success(mid,p0,p1)>=goal:
            hi=mid
        else:
            lo=mid
    mD=math.sqrt(I['m_D0_over_gv']**2+aD)
    kinetic=float(simpson(w*(np.sqrt(ks*ks+mD*mD)-mD),x=ks))
    curve=[]
    for n in sorted(set([0,1,10,100,1000,10000,Nopt,hi])):
        visibility=math.exp(.5*n*logc2)
        pg=success(n)
        curve.append(dict(probes=n,optimal_success=pg,register_visibility=visibility,
                          information_dephasing_identity=(2*pg-1)**2+visibility**2,
                          fixed_measurement_success=binomial_success(n,p0,p1) if n else .5))
    return dict(overlap_real=float(overlap.real),overlap_imag=float(overlap.imag),
                overlap_squared=c2,single_probe_optimal_success=success(1),
                ideal_independent_probe_count_optimal_95=Nopt,
                optimal_success_at_N=success(Nopt),optimal_success_before_N=success(Nopt-1),
                measurement_beta=beta,measurement_plus_probabilities=[p0,p1],
                single_probe_fixed_measurement_success=.5*(1+abs(p0-p1)),
                fixed_measurement_probe_count_95=hi,
                fixed_measurement_success_at_N=binomial_success(hi,p0,p1),
                fixed_measurement_success_before_N=binomial_success(hi-1,p0,p1),
                mean_probe_kinetic_energy_over_gv=kinetic,
                incident_kinetic_budget_at_optimal_count_over_smaller_mass=Nopt*kinetic/smallmass,
                incident_kinetic_budget_at_fixed_count_over_smaller_mass=hi*kinetic/smallmass,
                no_relative_phase_monochromatic_control_success=.5,curve=curve)


def write_csv(path,columns,rows):
    with path.open('w',newline='') as file:
        writer=csv.writer(file,lineterminator='\n');writer.writerow(columns);writer.writerows(rows)


def main():
    source=ROOT.parent/I['source_profile_relative_to_parent']
    assert hashlib.sha256(source.read_bytes()).hexdigest()==I['source_profile_sha256']
    data=ROOT/'data';data.mkdir(exist_ok=True)
    seed=Seed(source)
    central=solve_profile(39944,I['radial_domain'],I['profile_tolerance'],seed)
    sols={39944:central}
    for Q in I['charges']:
        if Q!=39944:
            sols[Q]=continuation(39944,Q,central,I['radial_domain'],I['profile_tolerance'])
    stats={Q:profile_quantities(sols[Q],Q) for Q in I['charges']}
    ks=np.linspace(I['probe_k_min'],I['probe_k_max'],I['probe_k_points'])
    phase={}; refined_sols={}; refinements=[]; spectra=[]
    for Q in I['main_charges']:
        delta,ctrl=phases(sols[Q],ks)
        fine=solve_profile(Q,I['refined_radial_domain'],I['refined_profile_tolerance'],sols[Q])
        refined_sols[Q]=fine
        fstats=profile_quantities(fine,Q); fdelta,fctrl=phases(fine,ks,refined=True)
        # Published readout uses the refined profile and spectrum throughout.
        phase[Q]=fdelta
        grid=np.linspace(0,fine.x[-1],I['raw_profile_points'])
        raw=fine.sol(grid)
        write_csv(data/f'profile_Q{Q}.csv',['rho','F','Fprime','S','Sprime'],np.column_stack([grid,raw[:4].T]))
        write_csv(data/f'phases_Q{Q}.csv',['k_over_gv','delta_l0','delta_l1'],np.column_stack([ks,fdelta.T]))
        rec=dict(Q=Q,max_phase_refinement_difference=float(np.max(abs(delta-fdelta))),
                 relative_mass_difference=abs(stats[Q]['mass_over_gv']/fstats['mass_over_gv']-1),
                 scattering_length_refinement_difference=abs(stats[Q]['scattering_length_over_inverse_gv']-
                                                             fstats['scattering_length_over_inverse_gv']),
                 base_phase_control=ctrl,refined_phase_control=fctrl)
        refinements.append(rec)
        stats[Q]=fstats
        spectra.append(dict(Q=Q,radial_amplitude_runs=[hessian(fine,n) for n in I['radial_hessian_grid_points']]))
        check(f'Q{Q}_phase_methods',max(ctrl['max_wave_phase_difference'],fctrl['max_wave_phase_difference'])<2e-8,
              max(ctrl['max_wave_phase_difference'],fctrl['max_wave_phase_difference']),2e-8)
        check(f'Q{Q}_phase_refinement',rec['max_phase_refinement_difference']<2e-8,rec['max_phase_refinement_difference'],2e-8)
        check(f'Q{Q}_mass_refinement',rec['relative_mass_difference']<2e-8,rec['relative_mass_difference'],2e-8)
        check(f'Q{Q}_no_linear_D_storage_mode',fstats['no_reporter_bound_state_sufficient_integral']<1,
              fstats['no_reporter_bound_state_sufficient_integral'],1)
        check(f'Q{Q}_finite_radial_minimum',min(x['lowest_fixed_charge_radial_amplitude_eigenvalue'] for x in spectra[-1]['radial_amplitude_runs'])>0,
              min(x['lowest_fixed_charge_radial_amplitude_eigenvalue'] for x in spectra[-1]['radial_amplitude_runs']),0)
    deriv=(stats[39945]['mass_over_gv']-stats[39943]['mass_over_gv'])/2
    derivative_error=abs(deriv-stats[39944]['Omega'])
    for Q,row in stats.items():
        check(f'Q{Q}_charge_constraint',abs(row['computed_charge']/Q-1)<2e-8,abs(row['computed_charge']/Q-1),2e-8)
        check(f'Q{Q}_classical_free_carrier_binding',row['energy_per_charge_over_gv']<1,row['energy_per_charge_over_gv'],1)
        check(f'Q{Q}_virial',abs(row['relative_virial'])<2e-7,abs(row['relative_virial']),2e-7)
        check(f'Q{Q}_active_mass',abs(row['active_mass_over_energy']-1)<2e-7,abs(row['active_mass_over_energy']-1),2e-7)
    check('chemical_potential_identity',derivative_error<2e-6,derivative_error,2e-6)
    sum_error=max(abs(x['mass_depletion_relative_residual']) for x in stats.values())
    field_error=max(max(abs(x['integrated_F_equation_relative_residual']),
                        abs(x['integrated_S_equation_relative_residual'])) for x in stats.values())
    check('finite_charge_mass_readout_sum_rule',sum_error<2e-7,sum_error,2e-7)
    check('integrated_field_equations',field_error<2e-7,field_error,2e-7)
    check('exact_zero_energy_lengths_in_controlled_intervals',all(
          x['exact_scattering_length_lower_bound']<=x['scattering_length_over_inverse_gv']<=x['exact_scattering_length_upper_bound']
          for x in stats.values()))
    check('finite_momentum_gravity_hierarchy',max(x['gravitational_Coulomb_parameter_at_min_k'] for x in stats.values())<2e-7,
          max(x['gravitational_Coulomb_parameter_at_min_k'] for x in stats.values()),2e-7)
    no_contrast,zctrl=phases(sols[39944],ks,refined=True,zero=True)
    check('no_contrast_no_phase',float(np.max(abs(no_contrast)))<1e-14,float(np.max(abs(no_contrast))),1e-14)
    check('free_wave_matching',zctrl['max_wave_phase_difference']<2e-8,zctrl['max_wave_phase_difference'],2e-8)
    lowks=np.array([.001,.002,.004])
    low=[]
    for Q in I['main_charges']:
        # Low-energy fit independently checks the zero-energy radial integration.
        sol=solve_profile(Q,I['refined_radial_domain'],I['refined_profile_tolerance'],sols[Q])
        dd,_=phases(sol,lowks,refined=True)
        aa0=stats[Q]['scattering_length_over_inverse_gv']
        error=abs(-dd[0,0]/lowks[0]-aa0)
        low.append(dict(Q=Q,k=lowks.tolist(),delta_l0=dd[0].tolist(),zero_energy_length=aa0,
                        finite_k_length_error=float(error)))
        check(f'Q{Q}_zero_energy_scattering',error<5e-5,float(error),5e-5)
    Q0,Q1=I['readout_charges']
    info=readout(ks,phase[Q0],phase[Q1],stats[Q0]['mass_over_gv'])
    finer_ks=np.linspace(I['probe_k_min'],I['probe_k_max'],2*I['probe_k_points']-1)
    pair={}
    for Q in [Q0,Q1]:
        pair[Q],_=phases(refined_sols[Q],finer_ks,refined=True)
        write_csv(data/f'phases_packet_refined_Q{Q}.csv',['k_over_gv','delta_l0','delta_l1'],np.column_stack([finer_ks,pair[Q].T]))
    finer_info=readout(finer_ks,pair[Q0],pair[Q1],stats[Q0]['mass_over_gv'])
    packet_control=dict(points=len(finer_ks),overlap_squared=finer_info['overlap_squared'],
                        absolute_overlap_squared_difference=abs(finer_info['overlap_squared']-info['overlap_squared']),
                        optimal_count_95=finer_info['ideal_independent_probe_count_optimal_95'],
                        fixed_count_95=finer_info['fixed_measurement_probe_count_95'])
    check('packet_quadrature_refinement',packet_control['absolute_overlap_squared_difference']<2e-10,
          packet_control['absolute_overlap_squared_difference'],2e-10)
    check('packet_refinement_retains_counts',packet_control['optimal_count_95']==info['ideal_independent_probe_count_optimal_95']
          and packet_control['fixed_count_95']==info['fixed_measurement_probe_count_95'])
    check('probe_states_distinguishable',info['single_probe_optimal_success']>.5,info['single_probe_optimal_success'],.5)
    check('fixed_measurement_below_quantum_optimum',info['single_probe_fixed_measurement_success']<=info['single_probe_optimal_success']+1e-12,
          info['single_probe_optimal_success']-info['single_probe_fixed_measurement_success'],0)
    check('fixed_measurement_probe_count_consistent',info['fixed_measurement_probe_count_95']>=info['ideal_independent_probe_count_optimal_95'],
          info['fixed_measurement_probe_count_95'],info['ideal_independent_probe_count_optimal_95'])
    check('information_dephasing_identity',max(abs(x['information_dephasing_identity']-1) for x in info['curve'])<2e-14,
          max(abs(x['information_dephasing_identity']-1) for x in info['curve']),2e-14)
    result=dict(classification='One fixed covariant EFT gives linked matter/gravity/probe outputs and an ideal readable prepared charge label; no complete universe or observer derived.',
                profiles=[stats[Q] for Q in I['charges']],refinements=refinements,radial_spectra=spectra,
                chemical_potential=dict(finite_difference=deriv,Omega=stats[39944]['Omega'],absolute_error=derivative_error),
                low_energy_controls=low,readout=info,packet_refinement=packet_control,checks=CHECKS,check_count=len(CHECKS),
                failed_checks=sum(not x['passed'] for x in CHECKS))
    (ROOT/'RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    (ROOT/'ENVIRONMENT.json').write_text(json.dumps(dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                  platform=platform.platform(),input_sha256=hashlib.sha256((ROOT/'INPUTS.json').read_bytes()).hexdigest()),indent=2)+'\n')
    print(json.dumps(dict(checks=len(CHECKS),failed=result['failed_checks'],profiles=result['profiles'],readout=info,
                          chemical_potential=result['chemical_potential'],refinements=refinements),indent=2))
    if result['failed_checks']:
        raise SystemExit(1)


if __name__=='__main__':
    main()
