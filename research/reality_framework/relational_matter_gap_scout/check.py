#!/usr/bin/env python3
"""Independent SCOUT-RG01: local matter-defined neighbors, 1D phonons, Gibbs gaps.
Requires numpy, scipy. DOES NOT simulate a full quantum dynamical circumference.
"""
import json
from pathlib import Path
from math import pi, sqrt, log
import numpy as np
from scipy.special import kv, gammaln, jv
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

G,P,M,HBAR=1.,1.,100.,1.
a0=sqrt(G/P)
c=sqrt(2*G/(M*a0))
wmax=sqrt(8*G/(M*a0**3))

def hessian(n):
    k=2*G/a0**3
    H=np.zeros((n,n))
    for i in range(n):
        j=(i+1)%n
        H[i,i]+=k;H[j,j]+=k
        H[i,j]-=k;H[j,i]-=k
    return H

def dispersion(n):
    return wmax*abs(np.sin(pi*np.arange(n)/n))

def force_check():
    n=8
    eq=np.arange(n)*a0
    def gradient(u):
        x=eq+u
        gaps=np.diff(np.r_[x,x[0]+n*a0])
        return -G/np.roll(gaps,1)**2+G/gaps**2
    f=np.zeros(n);step=1e-5;jac=np.empty((n,n))
    for j in range(n):
        e=np.zeros(n);e[j]=step
        jac[:,j]=(gradient(f+e)-gradient(f-e))/(2*step)
    err=float(np.max(abs(jac-hessian(n))))
    assert err<1e-8
    return err

def gibbs(beta):
    z=2*beta*sqrt(G*P)
    # Gauge-fixed isothermal-isobaric configurational gap measure.
    Z=2*sqrt(G/P)*kv(1,z)
    direct=quad(lambda l:np.exp(-beta*(G/l+P*l)),0,np.inf,
                epsabs=1e-11,epsrel=1e-11)[0]
    assert abs(Z-direct)/Z<1e-10
    return dict(beta=beta,log_gap_partition=log(Z),
        mean_gap=float(a0*kv(2,z)/kv(1,z)),
        mean_inverse_gap=float(sqrt(P/G)*kv(0,z)/kv(1,z)))

def main():
    rng=np.random.default_rng(922337)
    samples=[]
    for n in (8,16,32,64,128):
        theta=np.arange(n)*a0
        L=n*a0
        def energy(gaps):
            return float(G*np.sum(1/np.array(gaps))+P*np.sum(gaps))
        exact_min=2*n*sqrt(G*P)
        assert abs(energy(np.ones(n)*a0)-exact_min)<1e-11
        for _ in range(25):
            gaps=np.exp(rng.normal(0,0.2,n));gaps*=L/sum(gaps)
            assert energy(gaps)>=exact_min-1e-9
        modes=dispersion(n)
        computed=np.sqrt(np.maximum(np.linalg.eigvalsh(hessian(n)/M),0))
        assert max(abs(np.sort(modes)-computed))<1e-8
        zpe=HBAR*modes.sum()/2
        zpe_formula=HBAR*wmax/(2*np.tan(pi/(2*n)))
        assert abs(zpe-zpe_formula)<1e-10
        gap_var=2*HBAR/(M*n*wmax*np.tan(pi/(2*n)))
        samples.append(dict(N=n,L=L,classical_minimum=exact_min,
            first_phonon_frequency=float(modes[1]),
            predicted_continuum_first_frequency=float(2*pi*c/L),
            relative_gap_error=float(modes[1]/(2*pi*c/L)-1),
            harmonic_zero_point_energy=float(zpe),
            zero_point_energy_per_site=float(zpe/n),
            quantum_neighbor_gap_rms=float(sqrt(gap_var))))
    assert samples[-1]['relative_gap_error']>-0.00011
    exact_zpe_density=HBAR*wmax/pi
    assert abs(samples[-1]['zero_point_energy_per_site']/exact_zpe_density-1)<1e-4
    # The one-loop harmonic quantum correction to the volume equilibrium.
    A=HBAR/pi*sqrt(8*G/M)
    e_density=lambda a:G/a+P*a+A*a**(-1.5)
    opt=minimize_scalar(e_density,bounds=(0.5,2.),method='bounded',
        options={'xatol':1e-12})
    assert opt.success and 1.06<opt.x<1.07
    # Exact impulse-velocity Green function: in the N->infinity ring limit,
    # v_r(t) = (p_impulse/m) J_{2r}(omega_max*t). Cyclic images vanish
    # for a sufficiently large ring at finite time.
    Nimp=128
    q=np.arange(Nimp)
    finite_bessel=[]
    for t in (2.,10.,25.):
        for distance in (0,1,2,5,10):
            finite=float(np.mean(np.cos(dispersion(Nimp)*t)*
                       np.cos(2*pi*q*distance/Nimp)))
            analytic=float(jv(2*distance,wmax*t))
            assert abs(finite-analytic)<1e-11
            finite_bessel.append(dict(t=t,relational_hops=distance,
                full_mode_velocity_amplitude=finite,
                infinite_ring_bessel=analytic))
    # Exact discrete config-count comparison with q=2N physical phase bins.
    counts=[]
    for n in (16,32,64,128,256):
        q=2*n;nm=n*(n-1)//2
        logphase=float(gammaln(q+1)-gammaln(n+1)-gammaln(q-n+1))
        loggraph=float(gammaln(nm+1)-gammaln(n+1)-gammaln(nm-n+1))
        counts.append(dict(N=n,log_phase_states=logphase,
                           log_fixed_N_edge_graph_states=loggraph))
        assert logphase<2*n*np.log(2)
    # Only neighbor springs; offsite order t^(2d+1) in momentum-kick response.
    k=hessian(64)
    assert abs(k[10,0])<1e-12
    # Pairwise logarithmic repulsion, by contrast, has a nonzero offsite
    # Hessian between every pair on a fixed phase circle.
    i,j,n=0,10,64
    pairwise_offsite=-(1/(4*n*np.sin(pi*(j-i)/n)**2))
    assert pairwise_offsite!=0
    out=dict(
       statement="Constrained S1 phase topology supplied; connections NOT independent.",
       parameters=dict(g=G,P=P,m=M,hbar=HBAR),
       analytic=dict(circumference="N sqrt(g/P)",gap="sqrt(g/P)",
        minimum_energy="2 N sqrt(gP)", dispersion="omega_k^2=(8g/m a^3)sin^2(pi k/N)",
        sound_speed="sqrt(2g/(m a))",
        one_gap_partition="2 sqrt(g/P) K1(2 beta sqrt(gP))",
        harmonic_zpe="hbar omega_max/2 cot(pi/(2N))"),
       sound_speed=c, force_hessian_error=force_check(),
       spectra=samples, finite_temperature=[gibbs(x) for x in (.5,1.,3.,10.)],
       harmonic_corrected_equilibrium_gap=float(opt.x),
       harmonic_corrected_energy_per_site=float(opt.fun),
       phase_graph_entropy_comparison=counts,
       impulse_velocity_bessel_controls=finite_bessel,
       pairwise_log_nonlocal_hessian_example=float(pairwise_offsite),
       limits=["S1 topology and cyclic interaction supplied",
          "quantum modes only harmonic; no exact full quantum thermal state",
          "volume equilibrium and phonon dynamics separately evaluated",
          "not Lorentz/diffeomorphism invariant; no gravity or empirical test"])
    Path(__file__).with_name("results.json").write_text(
        json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))
    print("SCOUT_RG01_RELATIONAL_GAP_CHECK_PASS")

if __name__=="__main__":
    main()
