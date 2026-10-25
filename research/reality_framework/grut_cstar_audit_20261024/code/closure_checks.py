#!/usr/bin/env python3
"""Mathematical identities and counterexamples; not a cosmological data fit."""
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def monomials(max_dimension=6,modulus=None):
    rows=[]
    for a in range(max_dimension+1):
        for b in range(a):
            for phi in range(0,max_dimension+1,2):
                for dark in range(0,max_dimension+1,2):
                    for higgs in range(max_dimension//2+1):
                        dimension=a+b+phi+dark+2*higgs
                        if not 1<=dimension<=max_dimension:continue
                        charge=a-b
                        if modulus and charge%modulus:continue
                        rows.append(dict(X=a,Xstar=b,phi=phi,D=dark,HdagH=higgs,dimension=dimension,charge_violation=charge))
    return sorted(rows,key=lambda r:(r['dimension'],r['charge_violation'],r['X'],r['phi'],r['D']))


def main():
    M=2.435e18;hbar=6.582119569e-25
    H0=67.4*1000/3.085677581491367e22;omegaL=.685
    HL=H0*np.sqrt(omegaL);rho=3*M*M*(hbar*HL)**2
    entropy=8*np.pi**2*M*M/(hbar*HL)**2
    z=np.diag([1.,-1.]);I=np.eye(2);H=np.kron(z,I)-np.kron(I,z)
    Q=np.kron(I,z);psi=np.array([1,0,0,1],complex)/np.sqrt(2)
    eig=np.linalg.eigvalsh(H);physical_dimension=int(sum(abs(eig)<1e-12))
    clock=[]
    for theta in [0.,np.pi/2,np.pi]:
        bra=np.array([1.,np.exp(-1j*theta)])/np.sqrt(2)
        conditional=bra@psi.reshape(2,2);conditional/=np.linalg.norm(conditional)
        sigma_x=np.array([[0,1],[1,0]])
        clock.append(dict(clock_phase=float(theta),conditional_sigma_x=float(np.vdot(conditional,sigma_x@conditional).real)))
    gamma=-np.log(.99)/(13.8e9*365.25*86400)
    times=np.linspace(0,3,31)
    closed_dS_error=float(np.max(abs(np.tanh(times)**2+1/np.cosh(times)**2-1)))
    ward=[]
    for power in [5,6]:
        theta=.37;radius=.81;coupling=.2+.17j;step=1e-6
        potential=lambda angle:2*np.real(coupling*(radius*np.exp(1j*angle))**power)
        numerical=(potential(theta+step)-potential(theta-step))/(2*step)
        exact=-2*power*np.imag(coupling*(radius*np.exp(1j*theta))**power)
        ward.append(dict(power=power,torque_finite_difference_relative_error=float(abs(numerical-exact)/abs(exact))))
    out=dict(scope='Exact elementary identities/counterexamples and operator enumeration, not empirical verification of C* or a quantum-gravity completion.',
        deSitter=dict(convention='Reduced Planck mass; entropy in nats; reference H0 and OmegaLambda are adopted, not inferred.',M_reduced_GeV=M,H0_km_s_Mpc=67.4,OmegaLambda=omegaL,rho_GeV4=float(rho),entropy_nats=float(entropy),entropy_bits=float(entropy/np.log(2)),identity_relative_error=float(abs(rho*entropy/(24*np.pi**2*M**4)-1)),fixed_entropy_implies_fixed_positive_density=True),
        finite_closed_clock=dict(physical_H_constraint_dimension=physical_dimension,H_state_residual=float(np.linalg.norm(H@psi)),exact_U1_commutator_norm=float(np.linalg.norm(H@Q-Q@H)),relational_clock_values=clock,
            meaning='A finite closed constrained quantum toy has a nontrivial physical state space and exact U(1). It is not a gravitational UV completion.'),
        fading_counterexample=dict(w=-2/3,power_law_a_t=2.,rho_power_t=-2.,H_power_t=-1.,future_event_horizon_radius_over_t=1.,meaning='rho tends to zero and w>-1 while an accelerating power-law FLRW future retains an event horizon; Minkowski global fate is not implied.'),
        compact_GR_counterexample=dict(spatial_topology='S3',scale_factor='cosh(H t)/H',positive_density='3 M_reduced^2 H^2',pressure='-rho',Friedmann_identity_max_error=closed_dS_error,meaning='A compact boundaryless classical GR solution satisfies its Hamiltonian/Friedmann constraint with positive vacuum energy and nonzero pressure. Closure alone does not force the q-theory equilibrium condition.'),
        vacuum_exchange=dict(intrinsic_w=-1.,effective_w_formula='-1 + Qtransfer/(3 H rho_vac)',meaning='A decaying vacuum transferring energy can retain intrinsic p=-rho. w>-1 requires a separately conserved/effective interpretation.'),
        one_dimensional_algebra=dict(matrix_algebra_dimension=1,maximum_nonzero_orthogonal_projectors=1,meaning='Ordinary operators on one-dimensional space cannot realize two distinct orthogonal records. Observer/code rules are additional structure.'),
        charge_retention_benchmark=dict(age_seconds=13.8e9*365.25*86400,fractional_loss_target=.01,constant_exponential_fractional_rate_target_per_second=float(gamma),energy_width_target_GeV=float(hbar*gamma),meaning='Engineering reference target, not an observational dark-matter lifetime bound.'),
        explicit_breaking_Ward_controls=ward,
        operator_scan=dict(scope='Non-derivative pure scalar/Higgs-invariant monomials through dimension six, preserving the two declared scalar sign symmetries. One representative of each conjugate pair; not the full fermionic/derivative EFT.',without_X_protection=monomials(),residual_Z5=monomials(modulus=5),residual_Z6=monomials(modulus=6)))
    (ROOT/'CLOSURE_CHECKS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('de Sitter entropy',entropy,'rho',rho,'clock',clock,flush=True)
    print('Z5 breaking operators',out['operator_scan']['residual_Z5'],'Z6',out['operator_scan']['residual_Z6'],flush=True)


if __name__=='__main__':main()
