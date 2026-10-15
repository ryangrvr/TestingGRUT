"""All-time sector bounds and a charged circulation energy window.

No Hamiltonian/coupling changes. Scalar functions evaluate analytically
derived inequalities; their floating-point values are not interval proofs.
"""
from pathlib import Path
import json,math

P=Path(__file__).resolve().parent
EPS=.25;MASS=HOP=1.


def residual_probability(delta,sigma,b):
    if delta<=sigma:return 1.
    x=(delta*sigma+b*math.sqrt(delta*delta+b*b-sigma*sigma))/(delta*delta+b*b)
    return min(1.,x*x)


def vacuum_bound(n,q,required_edges):
    b=abs(HOP)*math.sqrt(2*q)
    return residual_probability(EPS*required_edges,b,b)


def low_matter_bound(n,q,required_edges,k):
    sigma=abs(HOP)*math.sqrt(2*q)
    b=sigma*(math.sqrt(k+1)+math.sqrt(k))
    return min(1.,((sigma+b)/(EPS*required_edges))**2)


def mean_matter_lower(p,nu=2.,s=2.):
    a=2*MASS+EPS/2;z=2*nu*HOP**2-a*EPS*s*p
    if z<0:return None
    return ((abs(HOP)*math.sqrt(2*nu)-math.sqrt(z))/a)**2


def dressed_moments(eta,cosphase=1.):
    a=2*MASS+3*EPS/2
    bright=math.sqrt(2*eta*(1-eta))
    local_mean=a*eta-2*HOP*bright*cosphase
    local_second=2*HOP**2+a*a*eta-2*a*HOP*bright*cosphase
    escape=2*HOP**2*(1-1.5*eta+eta*eta)
    return dict(matter_density=eta,electric_square_per_vertex=4+3*eta,
        mean_energy_per_vertex=1+2*EPS+local_mean,
        variance_per_vertex=local_second-local_mean**2+escape,
        local_band_variance=local_second-local_mean**2,
        outside_band_residual_per_vertex=escape,
        full_circulation_loss_T2_per_vertex=2*HOP**2*(1-eta/2)*(1-eta),
        full_protection_loss_T2_per_vertex=0.)


def run():
    a=2*MASS+3*EPS/2
    eta_ground=(1-a/math.sqrt(a*a+8*HOP**2))/2
    eta_neutral=(60-32*math.sqrt(2))/291
    grid=[]
    for n in [32,64,128,256,512,1024,4096,16384,65536,1048576]:
        grid.append(dict(N=n,Q=2*n,vacuum_spanning=vacuum_bound(n,2*n,n),
            vacuum_fully_circulating=vacuum_bound(n,2*n,2*n),
            full_circulation_nv_le_1=low_matter_bound(n,2*n,2*n,1),
            full_circulation_nv_le_sqrt_N=low_matter_bound(n,2*n,2*n,math.isqrt(n))))
    out=dict(frozen_couplings=dict(m=1,t=1,epsilon_E=.25,kappa=.5,Gamma=.1),
        all_time_bounds=grid,
        energy_probability_tradeoff=[dict(probability=p,necessary_mean_matter_density=mean_matter_lower(p)) for p in [.1,.5,.9,1.]],
        charged_comparators=dict(ground_band=dressed_moments(eta_ground),
            neutral_mean_matched=dressed_moments(eta_neutral),
            eta_neutral_exact='(60 - 32 sqrt(2))/291',
            eta_neutral_polynomial='873 eta^2 - 360 eta + 16 = 0'),
        mean_and_variance_matched_supplied_pairs=[dict(N=n,Q=2*n,matter_density=.25,
            h=math.sqrt(33/(16*(n-1))),
            phase_cosines=[(35/32-sign*math.sqrt(33/(16*(n-1))))/math.sqrt(3/2) for sign in [1,-1]],
            full_mean_energy=n,full_variance=4*n,full_circulation_probability=1,
            construction='Pure superposition of opposite Euler-flow branches on a supplied graph. Not claimed dynamically reachable from the neutral ray or in its same pure symmetry sector.') for n in [128,1024,11664]],
        scope='Probability bounds hold uniformly in evolution time for the specified neutral input energy and variance and full unbounded frozen CP12. They constrain low-matter circulation sectors, not unrestricted charged phases. Charged comparators are supplied states, not autonomously formed geometries. Tables are floating evaluations of analytic formulas.')
    (P/'ANALYTIC_BOUNDS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out['charged_comparators']),flush=True)
    return out


if __name__=='__main__':run()
