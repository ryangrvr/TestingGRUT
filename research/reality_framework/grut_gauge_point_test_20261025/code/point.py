#!/usr/bin/env python3
"""Six-condition diagnostics at the unchanged gauge point, with explicit variants."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize,root

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT.parent
MRED=2.435e18
HBAR=6.582119569e-25
HINF=4.7e13
CH=1.
DELTA=1e-22 # new, proposed wall-bias coefficient; absent from 7737487.


def parameters():
    a=json.loads((BASE/'grut_z6_gauge_bridge_20261025/MATCHING_RESULTS.json').read_text())
    p=a['parameters'];M=p['M_GeV']
    return p,a,M,p['Higgs_amplitude_over_M']*M,p['mu_over_M']*M


def inflation_functions(p,M,h,mu):
    lx,ly,ls,k=p['lambda_X'],p['lambda_Y'],p['lambda_Sigma'],p['kappa']
    mh=M/HINF;hh=h/HINF;u=mu/HINF
    # Positive amplitudes, at the phase valley making both interaction terms negative.
    def potential(z):
        x,y,s=z
        return lx*x**4+ly*y**4+ls*(s*s-hh*hh/2)**2+mh*mh*y*y-2*k*y*x**3-2*u*s*y*y-CH*x*x
    def grad(z):
        x,y,s=z
        return np.array([4*lx*x**3-6*k*y*x*x-2*CH*x,
            4*ly*y**3+2*mh*mh*y-2*k*x**3-4*u*s*y,
            4*ls*s*(s*s-hh*hh/2)-2*u*y*y])
    def hessian(z):
        x,y,s=z
        return np.array([[12*lx*x*x-12*k*y*x-2*CH,-6*k*x*x,0.],
            [-6*k*x*x,12*ly*y*y+2*mh*mh-4*u*s,-4*u*y],
            [0.,-4*u*y,12*ls*s*s-2*ls*hh*hh]])
    return potential,grad,hessian


def main():
    p,old,M,h,mu=parameters();potential,grad,hessian=inflation_functions(p,M,h,mu)
    candidates=[]
    for x in [15.,22.,30.]:
        for y,s in [(.1,.01),(1.,.1),(3.,1.)]:
            z=minimize(potential,[x,y,s],jac=grad,bounds=[(1e-12,None)]*3,method='L-BFGS-B',
                       options=dict(maxiter=3000,ftol=1e-15,gtol=1e-10))
            polished=root(grad,z.x,jac=hessian,tol=1e-12)
            zz=polished.x if np.all(polished.x>0) and np.linalg.norm(grad(polished.x))<1e-7 else z.x
            candidates.append(dict(amplitudes_over_H=zz.tolist(),potential_over_H4=float(potential(zz)),gradient_norm=float(np.linalg.norm(grad(zz)))))
    best=min(candidates,key=lambda x:x['potential_over_H4']);z=np.array(best['amplitudes_over_H']);x,y,s=z
    v1=np.array([3.,-1.,0.]);v2=np.array([0.,2.,-1.]);q=np.array([1.,3.,6.])
    angular=2*p['kappa']*y*x**3*np.outer(v1,v1)+2*(mu/HINF)*s*y*y*np.outer(v2,v2)
    phase_base=np.array([x,-y,-s],dtype=complex)
    def phase_potential(t):
        xx,yy,ss=phase_base*np.exp(-1j*np.asarray(t))
        return 2*np.real(p['kappa']*np.conj(yy)*xx**3+(mu/HINF)*np.conj(ss)*yy*yy)
    eps=.001;fd=np.zeros((3,3));zero=np.zeros(3)
    for i in range(3):
        ei=zero.copy();ei[i]=eps
        fd[i,i]=(phase_potential(ei)+phase_potential(-ei)-2*phase_potential(zero))/eps**2
        for j in range(i):
            ej=zero.copy();ej[j]=eps
            fd[i,j]=fd[j,i]=(phase_potential(ei+ej)+phase_potential(-ei-ej)-phase_potential(ei-ej)-phase_potential(-ei+ej))/(4*eps**2)
    inv=np.diag(1/(np.sqrt(2)*z))
    eig=np.linalg.eigvalsh(inv@angular@inv)
    r_single=np.sqrt(CH/(2*p['lambda_X']))*HINF
    d=old['matching']['D_over_M4']
    Y_quadratic=p['kappa']*r_single**3/(M*M*(1-np.sqrt(2)*mu*h/M**2))
    inflation=dict(assumed_H_inf_GeV=HINF,assumed_c_H=CH,
        preparation='Added Hubble mass -c_H H^2 |X|^2. Not in 7737487. A proposed covariant realization is -xi_X R|X|^2 with xi_X=-c_H/12 in de Sitter; its reheating dynamics are not solved.',
        full_radial_minimum=best,multiple_start_controls=candidates,
        full_radial_amplitudes_GeV=(z*HINF).tolist(),radial_mass2_over_H2=(np.linalg.eigvalsh(hessian(z))/2).tolist(),
        angular_mass2_over_H2=eig.tolist(),angular_gauge_null_residual=float(np.linalg.norm(angular@q)),
        angular_Hessian_finite_difference_relative_error=float(np.linalg.norm(fd-angular)/np.linalg.norm(angular)),
        gauge_vector_mass_over_H_at_min=float(np.sqrt(2)*p['gauge_coupling']*np.sqrt(x*x+9*y*y+36*s*s)),
        vacuum_mediator_mass_over_H=M/HINF,vacuum_Higgs_mass_over_H=old['matching']['Higgs_radial_mass_GeV']/HINF,
        vacuum_vector_mass_over_H=old['matching']['gauge_vector_mass_GeV']/HINF,
        single_field_reference_abs_X_GeV=r_single,
        initial_single_field_f_over_6h=np.sqrt(2)*r_single/(6*h),
        single_field_estimate_abs_kappa_X2_over_M2=p['kappa']*(r_single/M)**2,
        quadratic_Y_estimate_over_M=Y_quadratic/M,
        nonlinear_to_quadratic_Y_force_at_that_estimate=2*p['lambda_Y']*(Y_quadratic/M)**2,
        full_Higgs_amplitude_over_vacuum=float(s*HINF/(h/np.sqrt(2))),
        quoted_one_field_c6_window=[137.,1600.],matched_effective_c6=old['matching']['effective_c6_at_reduced_Planck'],
        full_potential_analytic_lower_bound_over_H4=-CH**2/(2*p['lambda_X']),
        minimal_nonminimal_coupling_variant=dict(xi_X=-CH/12,effective_Planck_mass2_over_Mred2=1+CH*(x*HINF/MRED)**2/6),
        hidden_energy_fraction_after_release_estimate=float((potential(z)+CH*x*x)*HINF**2/(3*MRED**2)),
        applicability_of_review_one_field_window=False)
    book=json.loads((BASE/'grut_portal_debris_window_20261023/FIXED_YIELD_RESULTS.json').read_text())
    row=next(a for a in book['rows'] if a['Q']==200000 and a['lambda_phi']==.03)
    v=row['v_phi_GeV'];lp=.03;mx=.1*v;ms=np.sqrt(2*lp)*v
    sigma=2*np.sqrt(2*lp)*v**3/3
    # Polynomial tilt keeps both +/-v stationary, with only +v at zero energy.
    wall_rows=[]
    for A in [.5,.8,2.]:
        for C in [1.,2.,10.]:
            for gs in [10.75,61.75]:
                tann=C*A*sigma/DELTA;H=1/(2*tann)
                T=np.sqrt(H*MRED/np.sqrt(np.pi*np.pi*gs/90))
                rad=3*MRED*MRED*H*H;entropy=2*np.pi**2*gs*T**3/45
                last=A*sigma/tann+DELTA
                wall_rows.append(dict(area_parameter=A,annihilation_coefficient=C,gstar=gs,
                    time_seconds=tann*HBAR,temperature_GeV=T,terminal_wall_plus_bias_fraction=last/rad,
                    terminal_energy_per_entropy_GeV=last/entropy))
    walls=dict(status='Absent in the baseline; separately proposed bias variant. Timing uses a scaling-network estimate, not a simulation.',
        DeltaV_GeV4=DELTA,sigma_GeV3=sigma,sigma_cube_root_GeV=sigma**(1/3),
        potential='V0 + DeltaV/4 [2-3 phi/v+(phi/v)^3]',
        factorization='(phi/v-1)^2/4 * [lambda_phi v^4 (phi/v+1)^2+DeltaV(phi/v+2)]',
        bias_over_old_barrier=DELTA/(lp*v**4/4),
        bracket_global_minimum_GeV4=DELTA-DELTA**2/(4*lp*v**4),
        false_vacuum_mass2_GeV2=2*lp*v*v-3*DELTA/(2*v*v),
        plus_vacuum_mass2_shift_GeV2=3*DELTA/(2*v*v),rows=wall_rows,
        minimum_estimated_annihilation_temperature_GeV=min(z['temperature_GeV'] for z in wall_rows),
        maximum_terminal_energy_to_calibrated_DM=max(z['terminal_energy_per_entropy_GeV'] for z in wall_rows)/book['calibrated_total_DM_GeV_per_entropy'],
        vacuum_zero='The added constant sets the +v reference to zero by construction; no radiative vacuum-energy selection is supplied.')
    target=book['calibrated_total_DM_GeV_per_entropy'];Y=book['original_charge_yield']
    limits=[]
    for allowed in [.01,.1,1.]:
        yp=allowed*target/(2*mx)
        limits.append(dict(allowed_extra_DM_fraction=allowed,extra_cold_pair_yield_max=yp,
            extra_pairs_per_net_charge_max=yp/Y,
            epsilon_min_in_quoted_remainder_model=(1+allowed/.9)**(-.75)))
    quench=[]
    for T in [40.,69.,86.]:
        for gs in [86.25,106.75]:
            ent=2*np.pi**2*gs*T**3/45
            value=2*mx/ms*(lp*v**4/4)/ent
            quench.append(dict(T_GeV=T,gstar_entropy=gs,stable_X_rest_energy_per_entropy_if_zero_momentum_scalar_decay_GeV=value,
                               ratio_to_calibrated_DM=value/target))
    remainder=dict(status='Unknown dynamics, not erased or fitted away.',fixed_net_yield=Y,calibrated_DM_GeV_per_entropy=target,
        extra_pair_budgets=limits,quoted_remainder_law='f_extra≈0.9(epsilon^(-4/3)-1), from the reviewer estimate; not independently derived here.',
        quench_debris_diagnostic=quench,
        quench_scope='Conditional release of lambda_phi v^4/4 into zero-momentum scalar quanta decaying to X anti-X. Not a lower bound, complete transition simulation, or extra energy to add if already present in a source ledger.')
    out=dict(scope='Stated-point tests and tagged proposed variants; no joint cosmological viability claim.',parameters=p,
        inflation=inflation,wall_bias_variant=walls,symmetric_remainder=remainder,
        kinetic_mixing=dict(one_loop_sum_hypercharge_times_dark_charge=0,
            assumption='chi=0 at matching, no bi-charged fields, and real scalar couplings with dark charge-conjugation symmetry. The bias and xi|X|^2 preserve it. A new CP/phase-fixing sector needs a fresh audit.',
            predicted_by_quantum_gravity=False))
    (ROOT/'POINT_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('FULL_INFLATION_AMPLITUDES/H',z,'GRAD',best['gradient_norm'],'PHASE_MASS2/H2',eig)
    print('EARLY_THRESHOLDS',inflation['vacuum_mediator_mass_over_H'],inflation['single_field_estimate_abs_kappa_X2_over_M2'])
    print('WALL_VARIANT_T_MIN_GeV',walls['minimum_estimated_annihilation_temperature_GeV'])


if __name__=='__main__':main()
