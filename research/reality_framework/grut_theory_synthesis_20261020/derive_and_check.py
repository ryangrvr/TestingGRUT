"""Check the classical scalar-sector derivations and resolve a planar interface.

These are author calculation checks, not laboratory tests or preregistered
collision outcomes. The same a=1,b=.1 action is retained. NumPy/SciPy required.
"""
from pathlib import Path
import hashlib
import json
import platform

import numpy as np
import scipy
from scipy.integrate import simpson, solve_bvp
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent
A, B = 1.0, 0.1


def potential(f, s, a=A, b=B):
    return .25*a*(s*s-1)**2 + .5*s*s*f*f + .25*b*f**4


def envelope(f, a, b):
    y = f*f
    return np.where(y <= a, .5*y+.25*(b-1/a)*y*y, .25*a+.25*b*y*y)


def mixed_matrix(k, f, s, omega, a=A, b=B):
    """Exact full two-field tangent generator around a uniform rotating state."""
    matrix = np.zeros((6,6))
    matrix[:3,3:] = np.eye(3)
    matrix[3,0] = -k*k-2*b*f*f
    matrix[3,2] = -2*s*f
    matrix[3,4] = -2*omega
    matrix[4,1] = -k*k
    matrix[4,3] = 2*omega
    matrix[5,0] = -2*s*f
    matrix[5,2] = -k*k-2*a*s*s
    return matrix


def wall(length, points, tol):
    omega = (A*B)**.25
    f0 = (A/B)**.25
    x = np.linspace(-length, length, points)
    th = np.tanh(.5*x)
    derivative = .5/np.cosh(.5*x)**2
    guess = np.array([f0*(1-th)/2, -f0*derivative/2,
                      (1+th)/2, derivative/2])

    def ode(x, y):
        f, df, s, ds = y
        return np.array([df, (s*s+B*f*f-omega*omega)*f,
                         ds, (A*(s*s-1)+f*f)*s])

    def bc(left, right):
        return np.array([left[0]-f0, left[2], right[0], right[2]-1])

    solution = solve_bvp(ode, bc, x, guess, tol=tol, max_nodes=50000)
    if not solution.success:
        raise RuntimeError(solution.message)
    xx = np.linspace(-length, length, 20001)
    f, df, s, ds = solution.sol(xx)
    kinetic = .5*(df*df+ds*ds)
    grand = potential(f, s)-.5*omega**2*f*f
    tension = simpson(kinetic+grand, x=xx)
    # Translational first integral is exactly zero on the infinite-domain wall.
    first_integral = kinetic-grand
    row = dict(length=length, initial_points=points, tolerance=tol,
               final_nodes=solution.x.size, tension=float(tension),
               max_first_integral_error=float(np.abs(first_integral).max()),
               endpoint_slope=float(max(abs(df[0]), abs(ds[0]), abs(df[-1]), abs(ds[-1]))),
               minimum_F=float(f.min()), minimum_S=float(s.min()),
               max_rms_residual=float(solution.rms_residuals.max()))
    return row, dict(x=xx, F=f, Fprime=df, S=s, Sprime=ds, grand_potential=grand)


def main():
    checks = []
    def check(name, passed, detail):
        checks.append(dict(name=name, passed=bool(passed), detail=detail))
    rng = np.random.default_rng(20261020)
    worst_ratio = worst_envelope = worst_bound = 0.0
    parameter_checks = []
    # Both sides of ab=1, and asymmetric values. The proofs are in THEORY.md;
    # these checks only guard normalization and implementation of the formulas.
    for a, b in [(1,.1), (.3,.7), (3,.02), (1,1), (1,1.2), (.2,8), (4,.4)]:
        mu = min(1.0, (a*b)**.25)
        f = np.geomspace(1e-5, 20*max(1,(a/b)**.25), 2000)
        smin = np.sqrt(np.maximum(0,1-f*f/a))
        err = float(np.max(abs(potential(f,smin,a,b)-envelope(f,a,b))))
        worst_envelope = max(worst_envelope, err)
        # Independently optimize the original two-field potential ratio.
        candidates = []
        for start in [[(a/b)**.25,0], [.02,1], [1,.5]]:
            result = minimize(lambda z: 2*potential(z[0],z[1],a,b)/z[0]**2,
                              start, bounds=[(1e-5,100),(0,20)],
                              method='L-BFGS-B', options={'ftol':1e-14,'gtol':1e-10})
            candidates.append(float(result.fun))
        numeric = min(candidates)
        ratio_error = abs(numeric-mu*mu)
        worst_ratio = max(worst_ratio, ratio_error)
        u,w,pu,pw = rng.normal(size=(4,10000))
        s = rng.uniform(-3,3,10000)
        f = np.hypot(u,w)
        q = w*pu-u*pw
        epsilon = .5*(pu*pu+pw*pw)+potential(f,s,a,b)
        gap = float(np.min(epsilon-mu*np.abs(q)))
        worst_bound = min(worst_bound, gap)
        parameter_checks.append(dict(a=a,b=b,mu=mu,optimized_ratio=numeric,
                                     ratio_error=ratio_error,minimum_sample_energy_gap=gap))
    check('Local envelope agrees with original potential at its S minimum', worst_envelope < 1e-8, worst_envelope)
    check('Independent optimization recovers sharp potential-ratio infimum', worst_ratio < 2e-7, worst_ratio)
    check('Sampled pointwise charge-energy bound has no violation', worst_bound >= -1e-11, worst_bound)

    omega = (A*B)**.25
    f0 = (A/B)**.25
    n0 = omega*f0*f0
    epsilon0 = potential(f0,0)+.5*omega**2*f0*f0
    pressure0 = .5*omega**2*f0*f0-potential(f0,0)
    check('Bulk coexistence solves amplitude equation and zero-pressure condition',
          abs(B*f0*f0-omega**2)<1e-12 and abs(pressure0)<1e-12 and abs(epsilon0-A)<1e-12,
          dict(amplitude_equation=B*f0*f0-omega**2,pressure=pressure0,energy_density=epsilon0))

    # Independent eigenvalues of the real-time first-order linearized equations.
    # Variables (eta,xi,eta_dot,xi_dot) in the rotating frame.
    maximum_spectrum_error = 0.
    spectra = []
    for k in np.concatenate(([0.],np.geomspace(1e-3,5,80))):
        matrix = np.array([[0,0,1,0],[0,0,0,1],[-k*k-2*omega**2,0,0,-2*omega],
                           [0,-k*k,2*omega,0]],float)
        numerical = np.sort(np.abs(np.linalg.eigvals(matrix).imag))
        lower = k*k+3*omega**2-np.sqrt(9*omega**4+4*omega**2*k*k)
        upper = k*k+3*omega**2+np.sqrt(9*omega**4+4*omega**2*k*k)
        expected = np.repeat(np.sqrt(np.maximum(0,[lower,upper])),2)
        maximum_spectrum_error=max(maximum_spectrum_error,float(abs(numerical-expected).max()))
        spectra.append(dict(k=float(k),nu_minus_squared=float(lower),nu_plus_squared=float(upper),
                            nu_S_squared=float(k*k+f0*f0-A)))
    check('Bulk dispersion agrees with real-time matrix eigenvalues', maximum_spectrum_error<1e-10, maximum_spectrum_error)
    k=1e-3
    low=k*k+3*omega**2-np.sqrt(9*omega**4+4*omega**2*k*k)
    check('Long-wavelength bulk sound speed approaches 1/sqrt(3)', abs(low/(k*k)-1/3)<1e-6, float(low/(k*k)))

    # A separate exact homogeneous solution on the mixed S!=0 branch.
    # Its unstable long-wavelength mode supplies a conditional fragmentation test,
    # not an assertion about the nonstationary finite Gaussian preparations.
    fm=.5
    sm=np.sqrt(1-fm*fm/A)
    om=np.sqrt(1+(B-1/A)*fm*fm)
    mixed_cs2=(1-om*om)/(1-3*om*om)
    wave_numbers=np.linspace(0,1,2001)
    growth=[float(np.max(np.linalg.eigvals(mixed_matrix(k,fm,sm,om)).real)) for k in wave_numbers]
    small_k=1e-4
    numerical_slope=float(np.max(np.linalg.eigvals(mixed_matrix(small_k,fm,sm,om)).real))/small_k
    check('Dilute mixed-branch instability has the predicted long-wave growth',
          mixed_cs2<0 and abs(numerical_slope-np.sqrt(-mixed_cs2))<2e-7,
          dict(predicted=np.sqrt(-mixed_cs2),computed=numerical_slope))
    cutoff_squared=-(B*fm*fm+A*sm*sm)+np.sqrt((B*fm*fm-A*sm*sm)**2+4*sm*sm*fm*fm)
    cutoff=np.sqrt(cutoff_squared)
    check('Mixed-branch growth ends at its zero-frequency amplitude boundary',
          max(np.linalg.eigvals(mixed_matrix(.99*cutoff,fm,sm,om)).real)>1e-3 and
          max(np.linalg.eigvals(mixed_matrix(1.01*cutoff,fm,sm,om)).real)<1e-10,
          dict(cutoff=float(cutoff)))
    fastest=int(np.argmax(growth))
    mixed=dict(F=fm,S=float(sm),omega=float(om),long_wave_cs_squared=float(mixed_cs2),
               instability_cutoff=float(cutoff),sampled_fastest_k=float(wave_numbers[fastest]),
               sampled_largest_growth_rate=growth[fastest],scan_spacing=.0005,
               scope='Uniform coherent rotating solution in infinite flat space; finite Gaussian clouds and cosmology require their own dynamics.')
    mixed['growth_scan']=[dict(k=float(k),growth=float(g)) for k,g in zip(wave_numbers[::4],growth[::4])]

    wall_cases=[]
    for length,points,tolerance in [(16,401,1e-6),(24,801,1e-8)]:
        row,profile=wall(length,points,tolerance)
        wall_cases.append(row)
    (ROOT/'wall_profile.json').write_text(json.dumps({key:value[::40].tolist() for key,value in profile.items()},indent=2)+'\n')
    check('Planar interface tension survives domain/tolerance refinement',
          abs(wall_cases[0]['tension']-wall_cases[1]['tension'])<2e-6,
          abs(wall_cases[0]['tension']-wall_cases[1]['tension']))
    check('Refined planar interface first integral and tail slopes are small',
          wall_cases[-1]['max_first_integral_error']<2e-7 and wall_cases[-1]['endpoint_slope']<2e-6,
          wall_cases[-1])
    sigma=wall_cases[-1]['tension']
    coefficient=(36*np.pi)**(1/3)*sigma*n0**(-2/3)
    output=dict(status='Author analytic derivations with numerical implementation checks; no new collision or laboratory experiment.',
                a=A,b=B,benchmark=dict(omega0=omega,F0=f0,n0=n0,epsilon0=epsilon0,
                                      sound_speed=1/np.sqrt(3),amplitude_gap=np.sqrt(6)*omega,
                                      mediator_gap=np.sqrt(f0*f0-A)),
                parameter_checks=parameter_checks,dispersion=spectra,wall_cases=wall_cases,
                mixed_branch_prediction=mixed,
                thin_wall=dict(surface_tension=sigma,E_surface_Q_two_thirds_coefficient=coefficient,
                               rough_binding_crossover_Qbar=(coefficient/(1-omega))**3,
                               warning='Large-Q asymptotic only. This crossover is not a computed critical charge or a fit to the existing Qbar399 object.'),
                checks=checks,all_checks_pass=all(c['passed'] for c in checks),
                runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (ROOT/'DERIVATION_RESULTS.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:output[k] for k in ['benchmark','wall_cases','thin_wall','all_checks_pass']},indent=2))
    if not output['all_checks_pass']:
        raise SystemExit('Derivation implementation check failed; inspect results')


if __name__=='__main__':
    main()
