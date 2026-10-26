"""Inputs and conservation tests for a kinetic handoff, not a relic solver.

Uses the already published finite vector band. No new release or mode run.
All thermodynamic examples use fixed vacuum masses and a stated equilibrium
comparator; they are not thermal mass calculations or formation predictions.
"""
from pathlib import Path
import json
import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
PREVIOUS = ROOT.parent / 'grut_heavy_vector_budget_20261025'
MX = .09240896373753621
GX, LX, LP, V = .1, .001, .03, .924089637375362
MS = math.sqrt(2 * LP) * V
UC = .8611951308217225
YQ, DM = 4.8294009076988175e-9, 4.4e-10


def save(name, value):
    def convert(x):
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    (ROOT / name).write_text(json.dumps(value, indent=2, allow_nan=False,default=convert) + '\n')


def bose_moments(mass, temperature, chemical=0.):
    """Independent direct momentum quadrature; one real bosonic degree."""
    z, u = mass / temperature, chemical / temperature

    def integrand(x, moment):
        e = math.sqrt(x*x + z*z)
        v = e - u
        f = 0. if v > 700 else 1. / math.expm1(v)
        return x*x * f * (1. if moment == 0 else e if moment == 1 else x*x/(3*e))

    n = temperature**3 / (2*math.pi**2) * quad(integrand, 0, np.inf, args=(0,), epsabs=2e-11, epsrel=2e-11)[0]
    rho = temperature**4 / (2*math.pi**2) * quad(integrand, 0, np.inf, args=(1,), epsabs=2e-11, epsrel=2e-11)[0]
    pressure = temperature**4 / (2*math.pi**2) * quad(integrand, 0, np.inf, args=(2,), epsabs=2e-11, epsrel=2e-11)[0]
    return dict(n=n, rho=rho, P=pressure)


def bose_charge(temperature, chemical):
    """Charge quadrature without cancellation of the two large densities."""
    z, u = MX / temperature, chemical / temperature

    def fun(x):
        e = math.sqrt(x*x + z*z)
        return x*x * 2*math.exp(-e)*math.sinh(u) / ((-math.expm1(-(e-u)))*(-math.expm1(-(e+u))))

    return temperature**3 / (2*math.pi**2) * quad(fun, 0, np.inf, epsabs=2e-11, epsrel=2e-11)[0]


def eq_state(temperature, net_charge, include_s):
    cap = UC*MX
    if bose_charge(temperature, cap) < net_charge:
        mu = cap
        nc = net_charge - bose_charge(temperature, mu)
    else:
        mu = brentq(lambda u: bose_charge(temperature, u)-net_charge, 0., cap, xtol=2e-15)
        nc = 0.
    plus = bose_moments(MX, temperature, mu)
    minus = bose_moments(MX, temperature, -mu)
    neutral = bose_moments(MS, temperature) if include_s else dict(n=0.,rho=0.,P=0.)
    rho = plus['rho']+minus['rho']+neutral['rho']+cap*nc
    return dict(T_GeV=temperature,T_over_m=temperature/MX,mu_over_m=mu/MX,
                condensed_charge=nc,plus=plus,minus=minus,neutral=neutral,rho=rho,
                pair_inventory=minus['n']+neutral['n'],
                scalar_energy_fraction=neutral['rho']/rho)


def thermodynamics(R, epsilon):
    xc = math.sqrt(LP*V*V/(2*GX*GX))
    nc = 2*math.sqrt(2*LX)*xc**3
    rhoc = 3*LX*xc**4
    nq = epsilon*nc
    energy = (1+R)*rhoc+LP*V**4/4
    states=[]
    for include_s in (False, True):
        def resid(lt):
            return eq_state(math.exp(lt),nq,include_s)['rho']-energy
        T = math.exp(brentq(resid,math.log(.03*MX),math.log(100*MX),xtol=2e-12))
        states.append(dict(include_s=include_s,**eq_state(T,nq,include_s)))
    return dict(R=R,nQ_GeV3=nq,rho_GeV4=energy,critical_amplitude_GeV=xc,
                n_circular_critical_GeV3=nc,rho_critical_GeV4=rhoc,states=states,
                scope='Fixed vacuum-mass chemical-equilibrium control, linear lump energy omega*Q. Not an initial-state derivation or a new relic threshold.')


def source(row, background):
    H = background['parameters']['H_GeV']
    nx, ex = abs(background['late']['N_X']), background['late']['E_X']
    mA, width, number = row['mA_GeV'],row['Gamma_A_GeV'],row['N_A_band']
    gamma = width/H
    pstar = math.sqrt(mA*mA/4-MX*MX)
    aunit = 1/math.sqrt(gamma)
    q0 = pstar*aunit # dimensional comoving momentum; a_end=1
    # Coordinate only: match the conserved classical charge to the case-L
    # density definition. This does NOT assign a gas or lump population.
    xc = math.sqrt(LP*V*V/(2*GX*GX))
    nq = abs(background['late']['epsilon'])*2*math.sqrt(2*LX)*xc**3
    aref = H*(nx/nq)**(1/3)
    aenergy = H*(ex/(3*LX*xc**4))**.25
    mean_a = aunit*math.sqrt(math.pi/2)
    mean_q = q0*math.sqrt(math.pi/2)
    log_survival=-gamma*(aref*aref-1)/2
    kmax_over_H=5*math.sqrt(mA/H)
    a_parent_cold=10*kmax_over_H/(mA/H)
    early_decay=-math.expm1(-gamma*(a_parent_cold*a_parent_cold-1)/2)
    x = np.geomspace(1e-4,10.,161)
    # dN/(H_i^3 d ln q) for EACH daughter, q=q0*x.
    pdf_log = number*x*x*np.exp(-x*x/2)
    # Integral of massless daughter conformal energy / initial X energy,
    # accounting for both species. This is valid for the hard part here.
    e_log = 2*(q0/H)*x*pdf_log/ex
    quantiles=[]
    for f in (.01,.1,.5,.9,.99):
        xx=math.sqrt(-2*math.log1p(-f))
        quantiles.append(dict(number_fraction=f,a_decay=xx*aunit,q_over_H=xx*q0/H,
                              p_at_reference_over_m=xx*q0/(aref*MX)))
    return dict(name=row['name'],e=row['e'],H_GeV=H,mA_GeV=mA,Gamma_A_GeV=width,
                gamma=gamma,N_A_over_H3=number,N_charge_over_H3=nx,E_X_over_H4=ex,
                daughter_number_each_over_H3=number,pairs_per_net_charge=number/nx,
                f_rest_if_unprocessed=2*MX*YQ*(number/nx)/DM,
                a_unit=aunit,a_mean_decay=mean_a,p_star_GeV=pstar,q0_over_H=q0/H,
                mean_q_over_H=mean_q/H,R_from_energy_moment=2*number*mean_q/(H*ex),
                R_parent_mass_approx=row['R_vector_band_only'],
                reference_a_from_charge=aref,reference_a_from_energy=aenergy,
                reference_scale_ratio=aref/aenergy,
                mean_p_at_reference_GeV=mean_q/aref,
                mean_p_at_reference_over_m=mean_q/(aref*MX),
                parent_p_over_m01_scale=a_parent_cold,
                estimated_decay_fraction_before_parent_p_over_m01=early_decay,
                parent_log_survival_at_reference=log_survival,
                surviving_parent_fraction_at_reference=math.exp(log_survival) if log_survival>-745 else 0.,
                quantiles=quantiles,
                spectrum=[dict(x=float(xx),q_over_H=float(xx*q0/H),
                               dN_each_over_H3_dlnq=float(nn),dR_dlnq=float(ee))
                          for xx,nn,ee in zip(x,pdf_log,e_log)],
                domain='Computed primordial finite band only; nonrelativistic parent at decay; prompt constant-g radiation during decay; vacuum daughter kinematics; no subsequent collisions. q=a*p, a_end=1. Reference a is a density coordinate, not a thermalization claim.')


def reservoir_examples():
    rows=[]
    for t in (.3,1.,3.,12.,30.):
        T=t*MX
        X=bose_moments(MX,T)
        s=bose_moments(MS,T)
        rows.append(dict(T_over_m=t,mu_plus=0.,mu_minus=0.,mu_s=0.,n_s_over_nminus=s['n']/X['n'],
                         scalar_fraction_of_gas_energy=s['rho']/(2*X['rho']+s['rho'])))
    coupling=2*GX*GX*V
    width=coupling**2/(16*math.pi*MS)*math.sqrt(1-4*MX*MX/(MS*MS))
    return dict(rows=rows,ms_GeV=MS,trilinear_s_X_Xbar_GeV=coupling,
                vacuum_Gamma_s_to_X_Xbar_GeV=width,
                scope='Vacuum masses, one neutral scalar and two charged degrees in chemical equilibrium. The vacuum width is not an in-medium rate.')


def main():
    previous=json.loads((PREVIOUS/'BUDGET_RESULTS.json').read_text())
    background=json.loads((PREVIOUS/'BACKGROUND_RESULTS.json').read_text())
    rows=[source(r,background) for r in previous['rows']]
    save('DECAY_SOURCE.json',dict(rows=rows,source_commit='225028b56fbc6f11174b8213394874ba6fadad4e',
                                source_band=[.2,5.],new_quantum_modes_run=False,
                                number_convention='One degree per daughter: n=int dp p^2 f/(2*pi^2). Comoving numbers are divided by H_i^3; q is a*p and a_end=1.'))
    thermal=thermodynamics(rows[0]['R_from_energy_moment'],abs(background['late']['epsilon']))
    save('THERMODYNAMIC_CONTROL.json',dict(reservoir=reservoir_examples(),onset_comparators=thermal,
                                         new_relic_threshold_computed=False))
    for row in rows:
        print('SOURCE',row['e'],row['R_from_energy_moment'],row['pairs_per_net_charge'],row['mean_p_at_reference_over_m'])
    print('THERMO',[(s['T_over_m'],s['scalar_energy_fraction'],s['condensed_charge']) for s in thermal['states']])


if __name__=='__main__':main()
