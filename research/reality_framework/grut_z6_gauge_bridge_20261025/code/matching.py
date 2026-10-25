#!/usr/bin/env python3
"""Exact quadratic heavy-Y matching and controls of a new bosonic EFT."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import root

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT.parent/'grut_cstar_audit_20261024'
PARAM=dict(lambda_X=.001,lambda_Y=1.,lambda_Sigma=1.,kappa=.001,
           Higgs_amplitude_over_M=.1,mu_over_M=.5,M_GeV=1e12,
           gauge_coupling=1e-7)


def matched(x,s,kappa,mu,m2=1.):
    denom=m2*m2-4*abs(mu*s)**2
    assert denom>0
    j=kappa*x**3
    y=-(m2*j-2*np.conj(mu)*s*np.conj(j))/denom
    # mu*s.conjugate() multiplies Y^2 in the potential.
    value=-(m2*abs(j)**2-2*np.real(mu*np.conj(s)*j*j))/denom
    return y,float(value)


def vy(y,x,s,kappa,mu,ly=0.,m2=1.):
    return float(m2*abs(y)**2+ly*abs(y)**4+
                 2*np.real(kappa*np.conj(y)*x**3+mu*np.conj(s)*y*y))


def main():
    p=PARAM;h=p['Higgs_amplitude_over_M'];s0=h/np.sqrt(2)
    k=p['kappa'];mu=p['mu_over_M'];den=1-4*(mu*s0)**2
    rng=np.random.default_rng(61925)
    errors=[];root_errors=[];quartic_changes=[];ward=[];mass=[]
    for amplitude in [.001,.01,.1,.3]:
        for _ in range(12):
            # Also test complex couplings and a non-unitary gauge choice.
            x=amplitude*np.exp(1j*rng.uniform(-np.pi,np.pi))
            s=s0*np.exp(1j*rng.uniform(-np.pi,np.pi))
            kap=k*np.exp(1j*rng.uniform(-np.pi,np.pi))
            muv=mu*np.exp(1j*rng.uniform(-np.pi,np.pi))
            y,value=matched(x,s,kap,muv)
            b=2*np.conj(muv)*s
            mat=np.array([[1+b.real,b.imag],[b.imag,1-b.real]])
            j=kap*x**3
            direct=np.linalg.solve(mat,-np.array([j.real,j.imag]))
            yd=direct[0]+1j*direct[1]
            scale=max(abs(value),1e-100)
            errors.append(max(abs(yd-y)/max(abs(y),1e-100),abs(vy(y,x,s,kap,muv)-value)/scale))
            ys=max(abs(y),1e-100)
            def residual(z):
                yy=ys*(z[0]+1j*z[1])
                res=(yy+2*np.conj(muv)*s*np.conj(yy)+kap*x**3+2*p['lambda_Y']*abs(yy)**2*yy)/ys
                return [res.real,res.imag]
            sol=root(residual,[y.real/ys,y.imag/ys],tol=1e-11)
            res=np.linalg.norm(residual(sol.x));root_errors.append(float(res))
            assert res<1e-10
            yy=ys*(sol.x[0]+1j*sol.x[1])
            quartic_changes.append(dict(X_over_M=amplitude,relative_energy_change=abs(vy(yy,x,s,kap,muv,p['lambda_Y'])-value)/scale))
            phase=rng.uniform(-2,2)
            moved=vy(y*np.exp(3j*phase),x*np.exp(1j*phase),s*np.exp(6j*phase),kap,muv)
            ward.append(abs(moved-vy(y,x,s,kap,muv))/scale)
            mass.extend(np.linalg.eigvalsh(mat).tolist())
    alpha=mu*k*k*s0/den # in units M^-2, multiplies X^6+h.c.
    beta=-k*k/den       # in units M^-2, multiplies |X|^6
    angles=np.arange(512)*2*np.pi/512
    energies=[];xf=.07
    for angle in angles:
        xx=xf*np.exp(1j*angle);jj=k*xx**3
        direct=np.linalg.solve(np.diag([1+2*mu*s0,1-2*mu*s0]),-np.array([jj.real,jj.imag]))
        yy=direct[0]+1j*direct[1]
        energies.append(vy(yy,xx,s0,k,mu))
    beta_fft=float(np.mean(energies)/xf**6)
    alpha_fft=complex(np.mean(np.array(energies)*np.exp(-6j*angles))/xf**6)
    # Exact angular Gauss reduction, with canonical real amplitudes f,h.
    gauss=[];kinetic=[];gauge=[]
    for _ in range(100):
        f,hg=np.exp(rng.uniform(-3,3,2));td,sd,ec=rng.normal(size=3)
        ec=abs(ec)+.01
        av=(f*f*td+6*hg*hg*sd)/(ec*(f*f+36*hg*hg))
        vx,vs=td-ec*av,sd-6*ec*av
        inertia=f*f*hg*hg/(f*f+36*hg*hg)
        delta=6*td-sd
        scale=max(abs(f*f*vx),abs(6*hg*hg*vs),1.)
        gauss.append(abs(f*f*vx+6*hg*hg*vs)/scale)
        kin=.5*f*f*vx*vx+.5*hg*hg*vs*vs
        kinetic.append(abs(kin-.5*inertia*delta*delta)/max(kin,1e-100))
        shift=rng.normal()
        av2=(f*f*(td+ec*shift)+6*hg*hg*(sd+6*ec*shift))/(ec*(f*f+36*hg*hg))
        gauge.append(abs((av2-av)-shift)/max(1.,abs(shift)))
    rows=json.loads((PARENT/'CHARGE_LEAKAGE_RESULTS.json').read_text())['rows']
    base={d:next(r for r in rows if r['power']==d and r['radius']==60. and r['potential_scale']==1.) for d in [5,6]}
    mred=base[6]['reference_reduced_Planck_GeV'];M=p['M_GeV']
    ac=alpha/M**2;eff=ac*mred*mred
    rate=base[6]['initial_fractional_core_charge_loss_rate_per_second_for_c1']*eff*eff
    age=base[6]['reference_age_seconds']
    point=json.loads((ROOT.parent/'grut_portal_debris_window_20261023/FIXED_YIELD_RESULTS.json').read_text())['refinement']
    v=.924089637375362
    f0=v*point['F0'];absx0=f0/np.sqrt(2)
    proxy=mred/np.sqrt(base[5]['initial_fractional_core_charge_loss_rate_per_second_for_c1']/base[6]['initial_fractional_core_charge_loss_rate_per_second_for_c1'])
    out=dict(scope='New tree-level gauge/mediator EFT matching and phase kinematics; no UV-QG completion, full quantum matching, abundance or formation prediction.',parameters=p,
        quadratic_matching=dict(samples=len(errors),max_relative_error=max(errors),max_nonlinear_root_residual=max(root_errors),min_mediator_mass2_over_M2=min(mass),max_gauge_potential_invariance_error=max(ward),full_quartic_energy_corrections=quartic_changes),
        stability=dict(sufficient_reserved_quartic_kappa_max=(p['lambda_X']**3*p['lambda_Y']/27)**.25,benchmark_kappa=k,
                       proof='Young bound reserves lambda_X/2 and lambda_Y/2. For |Sigma|/M<1 the remaining Y mass term is nonnegative. Above 1, the lower bound is (s^2-.005)^2-(s-1)^2/2 >= .490025*s^2 >0. Phi-X cross term is positive.'),
        matching=dict(D_over_M4=den,alpha6_times_M2=alpha,beta6_times_M2=beta,alpha6_GeV_minus2=ac,beta6_GeV_minus2=beta/M**2,
                      Fourier_relative_errors=dict(alpha=abs(alpha_fft-alpha)/abs(alpha),beta=abs(beta_fft-beta)/abs(beta)),
                      abs_torque_to_radial_sextic_ratio=2*abs(alpha/beta),effective_c6_at_reduced_Planck=eff,
                      initial_parent_response_rate_per_second=rate,rate_times_reference_age=rate*age,
                      reference_target_fraction=.01,gauge_vector_mass_GeV=6*p['gauge_coupling']*h*M,
                      reference_core_inertia_relative_correction=f0*f0/(f0*f0+36*(h*M)**2),
                      reference_core_radial_sextic_to_quartic_potential=abs(beta/M**2)*absx0**2/p['lambda_X'],
                      Higgs_radial_mass_GeV=np.sqrt(2*p['lambda_Sigma'])*h*M,
                      mediator_squared_masses_GeV2=[(1-2*mu*s0)*M*M,(1+2*mu*s0)*M*M]),
        compact_phase=dict(samples=len(gauss),max_Gauss_residual=max(gauss),max_kinetic_reduction_relative_error=max(kinetic),max_A0_gauge_shift_error=max(gauge),
                           inertia_fractions=[dict(f_over_h=r,relative_phase_inertia_to_single_field=36/(36+r*r)) for r in [.01,1.,6.,30.]]),
        amplitude_ratio_caution=dict(rate_ratio_overlap_scale_GeV=proxy,reference_canonical_f0_GeV=f0,reference_abs_X0_GeV=absx0,
                                     inference='The ratio defines an overlap scale; different harmonic Green functions and source/profile powers prevent identifying it with a measured central amplitude.'),
        instanton_caution=dict(conditional_S_min_if_abs_c5_equals_exp_minusS=float(-np.log(base[5]['max_abs_coefficient_for_one_percent_initial_rate_age_target'])),Wilson_coefficient_or_instanton_predicted=False))
    (ROOT/'MATCHING_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('MATCHING',len(errors),'controls; alpha6=',ac,'Gamma6=',rate)


if __name__=='__main__':main()
