"""Fixed-potential coexistence and global classical energy-charge lower bound."""
from common import *
from scipy.optimize import minimize_scalar,root

def main():
    a=LX/G2;b=LP/G2;c=EPS/G2
    u=lambda s:b/4*(s*s-1)**2-c*(s-1)**3
    us=lambda s:b*(s*s-1)*s-3*c*(s-1)**2
    f=lambda s:s*s+2*np.sqrt(a*u(s))
    opt=minimize_scalar(f,bounds=(0,1),method='bounded',options={'xatol':1e-14})
    # For these frozen a,b,c the squared stationary equation is
    # 2.1 S^4+3.76 S^3+2.864 S^2-.288 S-.036=0.
    # Exactly one positive root by its single sign change (descending order).
    # The unsquared equation has positive sides on (0,1), so no spurious root.
    quartic=np.array([2.1,3.76,2.864,-.288,-.036])
    S=float(brentq(lambda s:np.polyval(quartic,s),0,1,xtol=1e-15))
    F2=2*np.sqrt(u(S)/a);F=float(np.sqrt(F2));om=float(np.sqrt(f(S)))
    def eq(z):
        f,s,omega=z
        return [omega*omega-s*s-a*f*f,us(s)+f*f*s,u(s)-a*f**4/4]
    alternative=root(eq,[F,S,om],tol=1e-12)
    if not alternative.success and np.linalg.norm(eq(alternative.x))>1e-11:raise RuntimeError(alternative.message)
    rng=np.random.default_rng(20261027);samples=rng.uniform(-4,4,(100000,2));s=samples[:,0];f2=samples[:,1]**2
    U=u(s)+.5*s*s*f2+a/4*f2*f2
    margin=U-.5*om*om*f2
    rows=[]
    for density in (.1,.5,1.,1.5,2.,2.63):
        s=brentq(lambda s:us(s)+density*density*s,0,1)
        omega=np.sqrt(s*s+a*density*density)
        E=.5*omega*omega*density*density+u(s)+.5*s*s*density*density+a/4*density**4
        Q=omega*density*density
        rows.append({'F':density,'S':s,'Omega':omega,'E_over_Q':float(E/Q),
                     'density_over_bulk':float(Q/(om*F2)),
                     'minimum_Q_weighted_energy_fraction_left_outside_ground_bulk':float(1-om/(E/Q))})
    out={'scope':__doc__,'bulk':{'S':S,'F':F,'Omega':om,'energy_per_charge_over_mX':om,
            'canonical_charge_density_over_mX3':om*F2/G2},
         'stationary_quartic_descending':quartic.tolist(),
         'stationary_quartic_roots':[[float(z.real),float(z.imag)] for z in np.roots(quartic)],
         'bounded_minimizer_S_difference':float(S-opt.x),
         'direct_three_equation_difference':float(np.max(abs(alternative.x-[F,S,om]))),
         'direct_three_equation_residual':float(np.max(np.abs(eq(alternative.x)))),
         'sampled_potential_bound_min_margin':float(margin.min()),
         'bound':'E >= Omega_bulk*m_X*abs(Q) for the classical M1 scalar Hamiltonian',
         'proof':'U(-S)>U(S) for S>0; S>1 gives Omega^2>=1; minimize F^2 analytically then S in [0,1]. Complete kinetic squares and add nonnegative gradients.',
         'constant_density_phi_minimum':'unique positive S in (0,1] at every F>=0: F^2=b(1-S^2)+3c(S-2+1/S) decreases strictly on (0,1). No negative root because base V_phi derivative and F^2*S are both negative there.',
         'uniform_controls':rows,'formed_mass_distribution':None}
    save('BULK_RESULTS.json',out);print(json.dumps(out,indent=2))

if __name__=='__main__':main()
