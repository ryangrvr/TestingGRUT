"""Complementary algebra/linear-response controls, not a heavy-point replay."""
from pathlib import Path
import json
import math
import numpy as np
from scipy.integrate import quad, solve_ivp

ROOT = Path(__file__).resolve().parents[1]
MRED, GSTAR = 2.435e18, 106.75


def potential(z, p):
    x, y, sigma = z
    return (p['lx']*abs(x)**4+p['ly']*abs(y)**4
            +p['ls']*(abs(sigma)**2-p['h']**2/2)**2+p['M']**2*abs(y)**2
            +2*np.real(p['k']*np.conj(y)*x**3+p['mu']*np.conj(sigma)*y**2)
            -p['cH']*p['H']**2*abs(x)**2)


def force(z, p):
    x, y, sigma = z
    return np.array([
        2*p['lx']*abs(x)**2*x+3*np.conj(p['k'])*y*np.conj(x)**2-p['cH']*p['H']**2*x,
        p['M']**2*y+2*p['ly']*abs(y)**2*y+p['k']*x**3+2*np.conj(p['mu'])*sigma*np.conj(y),
        2*p['ls']*(abs(sigma)**2-p['h']**2/2)*sigma+p['mu']*y*y])


def dilate(p, factor):
    return {key:(value*factor if key in ['M','h','mu','H'] else value) for key,value in p.items()}


def matching(p):
    sigma = p['h']/math.sqrt(2)
    determinant = p['M']**4-4*abs(p['mu']*sigma)**2
    return (p['mu']*p['k']**2*np.conj(sigma)/determinant,
            -p['M']**2*abs(p['k'])**2/determinant)


def scalar_homogeneity():
    rng = np.random.default_rng(2026102502)
    p = dict(lx=.001,ly=1.,ls=1.,M=20.,h=25.,mu=3.5,H=1.,cH=10.,k=.002)
    maxima = dict(potential=0.,force=0.,source=0.,number=0.)
    for _ in range(40):
        z = 5*(rng.normal(size=3)+1j*rng.normal(size=3))
        velocity = rng.normal(size=3)+1j*rng.normal(size=3)
        base_v = potential(z,p); base_f = force(z,p)
        base_source = 2*np.imag(np.conj(z)*base_f)
        base_n = -2*np.imag(np.conj(z)*velocity)
        for factor in [.07,.3,2.,11.]:
            scaled = dilate(p,factor)
            for key, got, wanted in [
                ('potential',potential(factor*z,scaled),factor**4*base_v),
                ('force',force(factor*z,scaled),factor**3*base_f),
                ('source',2*np.imag(np.conj(factor*z)*force(factor*z,scaled)),factor**4*base_source),
                ('number',-2*np.imag(np.conj(factor*z)*(factor**2*velocity)),factor**3*base_n)]:
                denominator = max(float(np.max(np.abs(wanted))),1e-30)
                maxima[key] = max(maxima[key],float(np.max(np.abs(got-wanted)))/denominator)
    return dict(passed=max(maxima.values())<1e-11,samples=40,scale_factors=[.07,.3,2.,11.],
                maximum_relative_residuals=maxima,
                scope='Polynomial homogeneity on a supplied background; no trajectory or candidate yield.')


def entropy(hubble):
    temperature = (90*MRED*MRED*hubble*hubble/(math.pi**2*GSTAR))**.25
    return 2*math.pi**2/45*GSTAR*temperature**3


def normalization_scaling():
    p = dict(lx=.001,ly=1.,ls=1.,M=1e12,h=1e11,mu=5e11,H=1e11,cH=10.,k=.001)
    alpha, beta = matching(p)
    normalized_number = 300.
    normalized_debris = 5.
    reference_yield = normalized_number*p['H']**3/entropy(p['H'])
    rows=[]
    max_residual=0.
    for scale in [.1,.3,1.,3.,10.]:
        q=dilate(p,scale)
        anew,bnew=matching(q)
        actual_yield=normalized_number*q['H']**3/entropy(q['H'])
        ratios=dict(yield_ratio=actual_yield/reference_yield,
                    alpha_ratio=float(abs(anew/alpha)),beta_ratio=float(abs(bnew/beta)),
                    included_late_rate_ratio=float(abs(anew/alpha)**2),
                    string_tension_ratio=(q['h']/p['h'])**2,
                    spectator_fraction_ratio=scale**2,
                    conformal_nonrelativistic_deposit_ratio=(normalized_debris*q['M']/q['H'])/(normalized_debris*p['M']/p['H']))
        expected=dict(yield_ratio=scale**1.5,alpha_ratio=scale**-2,beta_ratio=scale**-2,
                      included_late_rate_ratio=scale**-4,string_tension_ratio=scale**2,
                      spectator_fraction_ratio=scale**2,conformal_nonrelativistic_deposit_ratio=1.)
        residual=max(abs(ratios[k]/expected[k]-1) for k in expected)
        max_residual=max(max_residual,residual)
        rows.append(dict(scale=scale,**ratios))
    return dict(passed=max_residual<1e-12,maximum_residual=max_residual,
                scope='Fixed low-energy core, fixed reduced Planck mass, fixed couplings/g* and heavy-daughter limit.',
                rows=rows,
                included_rate_reference='Only the matched X6 channels on the unchanged core; not a complete lifetime.')


def bounds_and_quartic():
    lambda_x,lambda_y,epsilon=.001,1.,.146
    kappa_max=(16*lambda_x**3*lambda_y/27)**.25
    amplitude_ratio_at_min=(kappa_max/(2*lambda_y))**(1/3)
    plants=[]
    for fraction in [.9,1.,1.1]:
        kk=kappa_max*fraction
        yy=(kk/(2*lambda_y))**(1/3)
        minimum=lambda_x+lambda_y*yy**4-2*kk*yy
        formula=lambda_x-3*kk**(4/3)/(2**(4/3)*lambda_y**(1/3))
        plants.append(dict(kappa_over_limit=fraction,minimum_quartic_at_X1=minimum,
                           analytic_minimum=formula))
    rows=[];max_match=0.
    for cH in [1.,10.]:
        envelope=math.sqrt(cH/(2*lambda_x))*(1+2/cH)**.25
        nstar=2*math.sqrt(2*lambda_x)/3**.75*envelope**3
        eta_required=epsilon*nstar/(12*envelope**6)
        for ratio_r in [.3,.5,.8]:
            kappa_over_u=math.sqrt(2*(1-ratio_r**2)*eta_required/ratio_r)
            upper_u=math.sqrt(3**.25*(cH*(cH+2))**.75*math.sqrt(lambda_y/lambda_x)
                              *ratio_r/((1-ratio_r**2)*epsilon))
            max_match=max(max_match,abs(upper_u-kappa_max/kappa_over_u)/upper_u)
            rows.append(dict(cH=cH,r=ratio_r,target_epsilon=epsilon,
                             conformal_quartic_envelope=envelope,eta_minimum=eta_required,
                             minimum_kappa_at_u13=13*kappa_over_u,
                             minimum_kappa_at_u43=43*kappa_over_u,
                             maximum_u_from_quartic_and_assumed_envelope=upper_u))
    # Independent dimensionful matching versus ratio form.
    matching_errors=[]
    for u in [13.,20.,43.]:
        for ratio_r in [.07071067811865475,.3,.8]:
            p=dict(M=u,h=25.,mu=ratio_r*u*u/(math.sqrt(2)*25),k=.002,H=1.)
            alpha,beta=matching(p)
            ratio_alpha=ratio_r*p['k']**2/(2*(1-ratio_r**2)*u*u)
            ratio_beta=-p['k']**2/((1-ratio_r**2)*u*u)
            matching_errors.append(max(abs(alpha/ratio_alpha-1),abs(beta/ratio_beta-1)))
    passed=(max_match<1e-12 and max(matching_errors)<1e-12
            and plants[0]['minimum_quartic_at_X1']>0 and plants[2]['minimum_quartic_at_X1']<0
            and abs(plants[1]['minimum_quartic_at_X1'])<1e-16)
    return dict(passed=passed,kappa_max_lambdaX_001_lambdaY_1=kappa_max,
                quartic_minimizing_absY_over_absX=amplitude_ratio_at_min,
                quartic_plants=plants,maximum_matching_formula_residual=max(matching_errors),
                maximum_combined_bound_formula_residual=max_match,conditional_envelope_rows=rows,
                scope='Quartic positivity exact. Epsilon bounds conditional on constant Sigma and the declared leading quartic envelope.')


def ramp_filter(x):
    if abs(x)<1e-3:
        zz=x/2
        return np.exp(1j*zz)*(1-zz*zz/10+zz**4/280)
    zz=x/2
    return np.exp(1j*zz)*3*(math.sin(zz)-zz*math.cos(zz))/zz**3


def preparation_response():
    rows=[];max_integral_error=0.;max_solve_error=0.
    for frequency in [.3,1.,3.,10.,30.,100.]:
        exact=ramp_filter(frequency)
        numerical=quad(lambda u:6*u*(1-u)*math.cos(frequency*u),0,1,epsabs=1e-12)[0]+1j*quad(
            lambda u:6*u*(1-u)*math.sin(frequency*u),0,1,epsabs=1e-12)[0]
        max_integral_error=max(max_integral_error,abs(exact-numerical))
        energies={}
        for prepared in [True,False]:
            def rhs(time,state):
                target=1-3*time*time+2*time**3
                return [state[1],-frequency**2*(state[0]-target)]
            solution=solve_ivp(rhs,[0,1],[1. if prepared else 0.,0.],method='DOP853',
                               rtol=2e-12,atol=2e-14)
            assert solution.success
            q,velocity=solution.y[:,-1]
            energy_ratio=q*q+(velocity/frequency)**2
            predicted=abs(exact-(0 if prepared else 1))**2
            max_solve_error=max(max_solve_error,abs(energy_ratio-predicted))
            energies['relaxed' if prepared else 'displaced']=dict(
                normalized_free_energy=energy_ratio,analytic=predicted)
        rows.append(dict(m_times_release_duration=frequency,**energies))
    return dict(passed=max_integral_error<1e-11 and max_solve_error<1e-9,
                maximum_filter_integral_error=max_integral_error,
                maximum_energy_solve_error=max_solve_error,rows=rows,
                scope='Canonical driven real oscillator only; neither the full gauge model nor the heavy spot point.')


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    xx=np.geomspace(.03,100,1800)
    ff=np.array([ramp_filter(x) for x in xx])
    fig,ax=plt.subplots(figsize=(8,4.8))
    ax.loglog(xx,abs(ff)**2,label='Initially follows equilibrium',color='#12616c',lw=2)
    ax.loglog(xx,abs(ff-1)**2,label='Initially displaced from equilibrium',color='#ac4b32',lw=2)
    ax.set_ylim(1e-8,4)
    ax.set_xlabel(r'Heavy-mode mass $\times$ release duration, $m\tau$')
    ax.set_ylabel(r'Final free energy / $(m^2q_0^2/2)$')
    ax.set_title('Preparation controls real heavy-mode debris',loc='left',fontsize=12)
    ax.grid(True,which='major',alpha=.2)
    ax.legend(loc='lower left',frameon=False)
    fig.text(.13,.015,'Exact smooth-ramp oscillator illustration; not a simulation of the GRUT spot point.',fontsize=9)
    fig.tight_layout(rect=[0,.04,1,1])
    fig.savefig(ROOT/'preparation_response.png',dpi=180)
    plt.close(fig)


def main():
    groups=dict(scalar_homogeneity=scalar_homogeneity(),normalization_scaling=normalization_scaling(),
                source_bound_and_quartic=bounds_and_quartic(),preparation_response=preparation_response())
    result=dict(author='Codex / AInstein; same-author controls',
                groups_passed=int(sum(g['passed'] for g in groups.values())),groups_total=len(groups),groups=groups,
                physical_status=dict(exact_heavy_spot_inputs_received=False,
                                     heavy_spot_independently_reproduced=False,
                                     full_kinetic_or_vector_release_computed=False,
                                     new_successful_cosmology_demonstrated=False))
    def native(value):
        if isinstance(value,np.generic):return value.item()
        raise TypeError(f'Unsupported result type: {type(value).__name__}')
    (ROOT/'RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=native)+'\n')
    plot()
    print(json.dumps(dict(groups_passed=result['groups_passed'],groups_total=result['groups_total'])))
    if result['groups_passed']!=result['groups_total']:raise SystemExit(1)


if __name__=='__main__':main()
