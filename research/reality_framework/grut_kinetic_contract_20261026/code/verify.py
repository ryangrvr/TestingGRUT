"""Mathematical checks of the handoff, not a production Boltzmann solution."""
import json
import math
import numpy as np
from scipy.integrate import quad,solve_ivp
from scipy.special import kve
from audit import ROOT,MX,MS,UC,bose_moments,save


def bessel_moments(m,T,mu=0.):
    # Independent Bessel-series reference for the direct quadrature in audit.
    k=np.arange(1,int(42*T/(m-mu))+80,dtype=float)
    z=k*m/T
    e=np.exp(-k*(m-mu)/T)
    n=m*m*T/(2*math.pi**2)*np.sum(e*kve(2,z)/k)
    rho=m*m*T/(2*math.pi**2)*np.sum(e*(3*T*kve(2,z)/(k*k)+m*kve(1,z)/k))
    return n,rho


def main():
    sources=json.loads((ROOT/'DECAY_SOURCE.json').read_text())['rows']
    groups={}
    pdf=lambda x:x*math.exp(-x*x/2)
    moments=[quad(lambda x:x**j*pdf(x),0,np.inf,epsabs=2e-13,epsrel=2e-13)[0] for j in (0,1,2)]
    analytic=[1.,math.sqrt(math.pi/2),2.]
    # An independently integrated time-domain decay: u=Gamma*t,
    # a/a_unit=sqrt(2u), daughter radiation work integral accumulates a*dN.
    def rhs(u,y):return [-y[0],y[0],math.sqrt(2*u)*y[0]]
    sol=solve_ivp(rhs,(0,40),[1.,0.,0.],method='DOP853',rtol=2e-12,atol=2e-14)
    n,d,er=sol.y[:,-1]
    moment_error=max(abs(x/y-1) for x,y in zip(moments,analytic))
    ode_error=max(abs(n+d-1),abs(er/math.sqrt(math.pi/2)-1))
    budget_error=max(abs(r['R_from_energy_moment']/r['R_parent_mass_approx']-1) for r in sources)
    early=max(r['estimated_decay_fraction_before_parent_p_over_m01'] for r in sources)
    groups['decay_number_energy_and_time'] = dict(passed=moment_error<1e-11 and ode_error<1e-9 and budget_error<1e-12 and early<1e-10,
        moment_quadrature_relative_error=moment_error,time_ODE_error=ode_error,
        parent_budget_relative_error=budget_error,ODE_evaluations=sol.nfev,
        estimated_early_decay_fraction_before_parent_p_over_m01=early)

    # Species order (X, anti-X, s). Columns: s decay, ss->X anti-X,
    # 2X+2anti-X->X+anti-X, 3s->2s. Charge is always conserved;
    # pair inventory n_minus+n_s is conserved by inverse decay alone.
    stoich=np.array([[1,1,-1,0],[1,1,-1,0],[-1,-2,0,-1]],int)
    charge=np.array([1,-1,0]);pair=np.array([0,1,1])
    q_change=charge@stoich;p_change=pair@stoich
    # Differentiate P=n_minus+K n_plus*n_minus at fixed T and net charge.
    derivative_errors=[]
    for K,nm,nq in ((.1,.2,.4),(3.,2.,.1),(30.,.001,1.)):
        P=lambda x:x+K*x*(x+nq)
        h=1e-5*nm
        fd=(P(nm+h)-P(nm-h))/(2*h)
        derivative_errors.append(abs(fd/(1+K*(2*nm+nq))-1))
    groups['reaction_inventory_and_fast_s_denominator']=dict(
        passed=bool(np.all(q_change==0) and p_change.tolist()==[0,-1,-1,-1] and max(derivative_errors)<1e-9),
        stoichiometry=stoich.tolist(),charge_change=q_change.tolist(),pair_inventory_change=p_change.tolist(),
        fixed_T_inventory_derivative_error=max(derivative_errors),
        reviewer_term='Existing rate_two s_channel already contains ss->X anti-X and this fixed-T inventory denominator. It is not inverse decay treated as a sink.')

    # Bose collision products for explicitly on-shell kinematics. These are
    # equilibrium tests, not thermal reaction coefficients.
    rng=np.random.default_rng(781024)
    db=[];energy=[]
    for _ in range(80):
        T=MX*10**rng.uniform(-.25,1.2)
        mu=MX*rng.uniform(-.7,.7)
        boost=rng.uniform(0,.75)
        gm=1/math.sqrt(1-boost*boost)
        p=math.sqrt(MS*MS/4-MX*MX)
        es=gm*MS;e1=gm*(MS/2+boost*p);e2=gm*(MS/2-boost*p)
        cases=[([es],[0.],[e1,e2],[mu,-mu]),
               ([MS,MS],[0.,0.],[MS,MS],[mu,-mu]),
               ([MX]*4,[mu,mu,-mu,-mu],[2*MX,2*MX],[mu,-mu])]
        for ein,uin,eout,uout in cases:
            fi=[1/math.expm1((e-u)/T) for e,u in zip(ein,uin)]
            fo=[1/math.expm1((e-u)/T) for e,u in zip(eout,uout)]
            forward=math.prod(fi)*math.prod(1+f for f in fo)
            backward=math.prod(fo)*math.prod(1+f for f in fi)
            db.append(abs(forward-backward)/(forward+backward))
            energy.append(abs(sum(ein)-sum(eout))/sum(ein))
    groups['on_shell_Bose_detailed_balance']=dict(passed=max(db)<1e-12 and max(energy)<1e-14,
        samples=len(db),max_forward_reverse_residual=max(db),max_energy_residual=max(energy))

    # Check the quadrature for both species, including charged chemical
    # potentials, against a different representation (Bessel series).
    errors=[]
    for m in (MX,MS):
        for t in (.3,1.,12.,30.):
            for mu in (0.,.5*m):
                direct=bose_moments(m,t*MX,mu)
                bn,br=bessel_moments(m,t*MX,mu)
                errors += [abs(direct['n']/bn-1),abs(direct['rho']/br-1)]
    th=json.loads((ROOT/'THERMODYNAMIC_CONTROL.json').read_text())
    target=th['onset_comparators']['rho_GeV4']
    energy_errors=[abs(s['rho']/target-1) for s in th['onset_comparators']['states']]
    groups['thermodynamic_reservoir_quadrature']=dict(passed=max(errors)<1e-9 and max(energy_errors)<1e-9,
        species_states=16,max_Bessel_comparison_error=max(errors),max_fixed_energy_match_error=max(energy_errors))

    # Exact finite-Q work identity tested with an independent central
    # directional derivative. This differentiable test energy is NOT a fit
    # to GRUT objects or a proposed population equation of state.
    fd_errors=[]
    for _ in range(40):
        Q=10**rng.uniform(2,6);nl=10**rng.uniform(-6,1)
        a,b=rng.uniform(.2,2,2);dn,dq=rng.uniform(-.3,.3,2)
        E=lambda q:a*q+b*q**(2/3)
        omega=a+(2/3)*b*Q**(-1/3)
        h=1e-5
        f=lambda x:nl*(1+x*dn)*E(Q*(1+x*dq))
        fd=(f(h)-f(-h))/(2*h)
        exact=omega*(Q*nl*(dn+dq))+(E(Q)-Q*omega)*nl*dn
        fd_errors.append(abs(fd-exact)/max(E(Q)*nl,1e-200))
    groups['finite_Q_object_number_work']=dict(passed=max(fd_errors)<1e-9,
        derivative_samples=len(fd_errors),max_normalized_derivative_error=max(fd_errors),
        test_energy='Arbitrary smooth aQ+bQ^(2/3) only; not a GRUT fit.')

    result=dict(groups=groups,passed=sum(bool(g['passed']) for g in groups.values()),total=len(groups),
                full_Boltzmann_run=False,new_threshold=False,physical_formation_kernel=False,
                scope='Mathematical checks of supplied sources and conservation structure; no physical viability certification.')
    save('VERIFICATION_RESULTS.json',result)
    print((ROOT/'VERIFICATION_RESULTS.json').read_text())
    if result['passed']!=result['total']:raise SystemExit(1)


if __name__=='__main__':main()
