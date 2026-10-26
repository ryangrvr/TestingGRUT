"""Same-author derivation/numerics checks for the limited onset problem."""
from scipy.integrate import solve_ivp,quad
from scipy.optimize import linear_sum_assignment,linprog
import numpy as np
from common import *
import circular as C
import elliptic as E
import roll

def main():
    groups={};rng=np.random.default_rng(3499)
    he=[]
    for _ in range(80):
        q=rng.normal(size=3)*12;h=1e-5*max(1.,np.linalg.norm(q));H=np.column_stack([(gradient(q+np.eye(3)[j]*h)-gradient(q-np.eye(3)[j]*h))/(2*h) for j in range(3)])
        he.append(np.max(abs(H-hessian(q)))/(1+np.max(abs(H))))
    spectra=[]
    for _ in range(300):
        s=C.state(rng.uniform(.01,1.2));k=rng.uniform(.001,3.)
        direct=np.linalg.eigvals(C.matrix(s,k))**2;poly=np.repeat(np.roots(C.polynomial(s,k)),2)
        i,j=linear_sum_assignment(abs(direct[:,None]-poly));spectra.append(max(abs(direct[i]-poly[j]))/(1+max(abs(poly))))
    groups['potential_and_dispersion']=dict(passed=max(he)<1e-7 and max(spectra)<1e-9,hessian=max(he),cubic_vs_matrix=max(spectra))

    end=C.row(1.);exact=max(C.growth(C.state(1.),k) for k in np.linspace(0,3,31))
    small=C.row(.02);relative=abs(small['full_growth']/small['adiabatic_growth']-1)
    groups['adiabatic_limit_and_endpoint']=dict(passed=exact<1e-8 and end['full_growth']==0 and end['adiabatic_growth']>.6 and relative<.005,
        exact_endpoint_growth_roundoff=exact,adiabatic_endpoint_growth=end['adiabatic_growth'],small_amplitude_growth_relative=relative,
        scope='Endpoint is marginal in neutral phi, not nonlinearly stable. No wavelength assigned without positive exponential growth.')

    bg,t0,P,angle,diag=E.orbit();roots=np.roots([LX/4,0,-EX,Q0*Q0/2]);z=sorted(roots.real);rn,rp=np.sqrt(z[1]),np.sqrt(z[2]);zn=z[0]
    def dt(theta):
        r=rn+(rp-rn)*np.sin(theta)**2;zz=r*r
        speed=np.sqrt(LX/(2*zz)*(rp*rp-zz)*(zz-rn*rn)*(zz-zn))
        return 2*(rp-rn)*np.sin(theta)*np.cos(theta)/speed
    period=2*quad(dt,0,np.pi/2,epsabs=0,epsrel=1e-10)[0]
    ang=2*quad(lambda th:dt(th)*Q0/(rn+(rp-rn)*np.sin(th)**2)**2,0,np.pi/2,epsabs=0,epsrel=1e-10)[0]
    groups['matched_orbit']=dict(passed=diag['max_relative_charge_drift']<1e-9 and diag['max_relative_energy_drift']<1e-9 and abs(period/P-1)<1e-9 and abs(ang-angle)<1e-9,
        energy=diag['max_relative_energy_drift'],charge=diag['max_relative_charge_drift'],period_quadrature=abs(period/P-1),angle_quadrature=abs(ang-angle))

    result=read('ELLIPTIC_RESULTS.json');sy=max(r['symplectic_absolute'] for r in result['rows']);de=max(r['determinant_error'] for r in result['rows'])
    _,M,_=E.mode(0.);b=bg.sol(t0);q,p=b[:2],b[2:];phase=np.r_[-q[1],q[0],0.,-p[1],p[0],0.]
    lq=np.array([p[1],-p[0],0.,-q[1],q[0],0.]);pe=np.linalg.norm(M@phase-phase)/np.linalg.norm(phase);qe=np.linalg.norm(lq@M-lq)/np.linalg.norm(lq)
    groups['canonical_evolution_and_charge']=dict(passed=sy<1e-8 and de<1e-8 and pe<1e-9 and qe<1e-9,
        symplectic=sy,determinant=de,phase_tangent=pe,charge_left_vector=qe,
        scope='The tiny k=0 charged exponent from splitting a unit Jordan multiplier is numerical, not a claimed instability.')

    comparisons=[];duration=[]
    for species,ids in [('X',[0,1,3,4]),('phi',[2,5])]:
        k=result['peaks'][species]['k'];r,M,_=E.mode(k);tight,Mt,_=E.mode(k,rtol=3e-13);rad,Mr,_=E.mode(k,rtol=2e-11,method='Radau')
        comparisons.append(dict(species=species,tolerance=abs(tight['mu_'+species]/r['mu_'+species]-1),
            Radau=abs(rad['mu_'+species]/r['mu_'+species]-1)))
        r4,M4,_=E.mode(k,rtol=3e-12,nperiod=4)
        sub=M[np.ix_(ids,ids)];direct=M4[np.ix_(ids,ids)]
        duration.append(float(np.linalg.norm(direct-np.linalg.matrix_power(sub,4))/np.linalg.norm(direct)))
    groups['method_tolerance_and_direct_growth']=dict(passed=max(max(r['tolerance'],r['Radau']) for r in comparisons)<1e-8 and max(duration)<1e-8,
        comparisons=comparisons,four_period_transfer=duration,scope='Distinct methods by one author, not external independent replication.')

    inv=read('INVENTORY_RESULTS.json');eta=inv['eta'];lp=[]
    for e in [1.,eta,2.]:
        for pp in [0.,.01,.1,1.,100.]:
            fit=linprog(-np.array([1,1,K.RS,e,e]),A_eq=np.array([[1,-1,0,1,-1],[0,1,1,0,1]]),b_eq=[1,pp],bounds=(0,None),method='highs')
            assert fit.success
            lp.append(abs(-fit.fun-(e+max(K.RS,2*e)*pp))/(1+abs(fit.fun)))
    groups['extended_inventory_inequality']=dict(passed=max(lp)<1e-10 and inv['random_rest_bound_max_residual']<=0,
        independent_linear_program=max(lp),random_mixtures_max_residual=inv['random_rest_bound_max_residual'],scope='Supplied upper object E/Q; no arbitrary transient configurations.')

    rr=read('ROLL_RESULTS.json');ref=roll.run(1e-6,result['peaks']['X']['k'],rtol=2e-13);base=rr['rows'][0]
    rollerr=abs(ref['tau_to_phi_0p1v']/base['tau_to_phi_0p1v']-1)
    groups['supplied_seed_roll']=dict(passed=rollerr<1e-8 and max(r['energy_error'] for r in rr['rows'])<1e-9 and max(r['charge_error'] for r in rr['rows'])<1e-9,
        time_tolerance=rollerr,energy=max(r['energy_error'] for r in rr['rows']),charge=max(r['charge_error'] for r in rr['rows']),
        scope='Homogeneous neutral seeds and one linear charged mode, not nonlinear fragmentation.')
    save('VERIFICATION_RESULTS.json',dict(groups=groups,physical_initial_spectrum_derived=False,physical_formation_kernel_derived=False,physical_f_sym=None))
    print(groups,flush=True)
    if not all(r['passed'] for r in groups.values()):raise SystemExit(1)

if __name__=='__main__':main()
