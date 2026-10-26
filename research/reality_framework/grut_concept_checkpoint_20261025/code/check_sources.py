"""Exact charge identities and the source/free-mode distinction.

Algebraic checks only. No evolved preparation or cosmological exclusion.
All dimensionful examples use M_Y=1; h/M_Y=.1 and mu/M_Y=.5.
"""
from pathlib import Path
import json
import math
import numpy as np
import sympy as sp
from scipy.optimize import root, brentq
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]


def symbolic_currents():
    ax, ay, ass, phase_k, phase_m = sp.symbols('ax ay ass phase_k phase_m', real=True)
    amp_k, amp_m = sp.symbols('amp_k amp_m', positive=True)
    phase_1 = 3*ax-ay+phase_k
    phase_2 = 2*ay-ass+phase_m
    potential = 2*amp_k*sp.cos(phase_1)+2*amp_m*sp.cos(phase_2)
    derived = sp.Matrix([sp.diff(potential, q) for q in [ax, ay, ass]])
    predicted = sp.Matrix([-6*amp_k*sp.sin(phase_1),
                           2*amp_k*sp.sin(phase_1)-4*amp_m*sp.sin(phase_2),
                           2*amp_m*sp.sin(phase_2)])
    exact = [sp.simplify(x) == 0 for x in derived-predicted]
    gauge = sp.simplify(derived[0]+3*derived[1]+6*derived[2])
    combined = sp.simplify(derived[0]+3*derived[1]+12*amp_m*sp.sin(phase_2))
    return dict(pass_group=all(exact) and gauge == 0 and combined == 0,
                exact_source_matches=exact, gauge_source=str(gauge),
                combined_source_residual=str(combined),
                convention='n_i=-2 Im(z_i* D_t z_i); phase coordinates z_i=r_i exp(i alpha_i)')


def complex_sources():
    rng = np.random.default_rng(2026102501)
    max_force_error, max_angle_error, max_gauge = 0., 0., 0.
    for _ in range(48):
        z = (rng.normal(size=3)+1j*rng.normal(size=3))*.4
        coupling_k = .07*np.exp(1j*rng.uniform(-math.pi, math.pi))
        coupling_m = .13*np.exp(1j*rng.uniform(-math.pi, math.pi))
        x, y, sigma = z
        forces = np.array([
            2*.001*abs(x)**2*x+3*np.conj(coupling_k)*y*np.conj(x)**2,
            y+2*abs(y)**2*y+coupling_k*x**3+2*np.conj(coupling_m)*sigma*np.conj(y),
            2*(abs(sigma)**2-.005)*sigma+coupling_m*y*y])
        source = 2*np.imag(np.conj(z)*forces)
        term_k = coupling_k*np.conj(y)*x**3
        term_m = coupling_m*np.conj(sigma)*y*y
        wanted = np.array([-6*term_k.imag,
                           2*term_k.imag-4*term_m.imag,
                           2*term_m.imag])
        scale = max(np.max(np.abs(wanted)), 1e-15)
        max_force_error = max(max_force_error, float(np.max(np.abs(source-wanted))/scale))
        max_gauge = max(max_gauge, abs(float(source@[1,3,6]))/scale)

        def potential(zz):
            xx, yy, ss = zz
            return (.001*abs(xx)**4+abs(yy)**4+(abs(ss)**2-.005)**2+abs(yy)**2
                    +2*np.real(coupling_k*np.conj(yy)*xx**3+coupling_m*np.conj(ss)*yy**2))

        for idx in range(3):
            shift = 1e-4
            vals = []
            for integer in [2,1,-1,-2]:
                trial = z.copy()
                trial[idx] *= np.exp(1j*integer*shift)
                vals.append(potential(trial))
            finite_diff = (-vals[0]+8*vals[1]-8*vals[2]+vals[3])/(12*shift)
            max_angle_error = max(max_angle_error, abs(float(finite_diff)-wanted[idx])/scale)
    return dict(pass_group=bool(max_force_error < 1e-12 and max_angle_error < 1e-8 and max_gauge < 1e-12),
                samples=48, complex_couplings=True,
                max_force_source_relative_error=max_force_error,
                max_angular_difference_relative_error=float(max_angle_error),
                max_gauge_source_normalized=max_gauge)


def forced_response():
    mediator_mass, sigma_h, coupling_m, coupling_k = 1., .1, .5, .001
    sigma = sigma_h/math.sqrt(2)
    matrix_a, matrix_b = mediator_mass**2, 2*coupling_m*sigma
    determinant = matrix_a**2-abs(matrix_b)**2
    real_matrix = np.array([[matrix_a+matrix_b.real, matrix_b.imag],
                            [matrix_b.imag, matrix_a-matrix_b.real]])
    alpha6 = coupling_m*coupling_k**2*np.conj(sigma)/determinant
    rows = []
    for amplitude in [.2,.1,.05,.025]:
        x = amplitude*np.exp(1j*math.pi/12)
        forcing = coupling_k*x**3
        ystar = (-matrix_a*forcing+matrix_b*np.conj(forcing))/determinant
        real_solution = np.linalg.solve(real_matrix, -np.array([forcing.real,forcing.imag]))
        y_matrix = real_solution[0]+1j*real_solution[1]
        normalized_residual = abs(matrix_a*ystar+matrix_b*np.conj(ystar)+forcing)/abs(forcing)
        source_x = -6*np.imag(coupling_k*np.conj(ystar)*x**3)
        source_d = -12*np.imag(coupling_m*np.conj(sigma)*ystar**2)
        source_effective = -12*np.imag(alpha6*x**6)

        def nonlinear(v):
            yy = v[0]+1j*v[1]
            ff = matrix_a*yy+matrix_b*np.conj(yy)+forcing+2*abs(yy)**2*yy
            return np.array([ff.real, ff.imag])/abs(forcing)

        solved = root(nonlinear, real_solution, tol=1e-11)
        # Hybrid root can report stalled relative progress below machine precision;
        # the independently evaluated equation residual determines acceptance.
        nonlinear_residual = float(np.linalg.norm(nonlinear(solved.x)))
        ynonlinear = solved.x[0]+1j*solved.x[1]
        nonlinear_source = -12*np.imag(coupling_m*np.conj(sigma)*ynonlinear**2)
        rows.append(dict(X_over_M=amplitude, abs_Ystar_over_X=float(abs(ystar)/abs(x)),
            matching_parameter_kappa_X2_over_M2=coupling_k*abs(x)**2,
            quartic_relative_force=2*abs(ystar)**2/mediator_mass**2,
            quadratic_stationarity_residual=float(normalized_residual),
            matrix_relative_error=float(abs(ystar-y_matrix)/abs(ystar)),
            source_nX=float(source_x), source_nD=float(source_d),
            source_effective=float(source_effective),
            source_match_relative_error=float(max(abs(source_x-source_d),
                                                  abs(source_x-source_effective))/abs(source_x)),
            instantaneous_free_displacement_energy=0.,
            full_nonlinear_Y_source=float(nonlinear_source),
            full_nonlinear_Y_residual=nonlinear_residual,
            nonlinear_relative_source_shift=float(abs(nonlinear_source/source_d-1))))
    passed = all(r['quadratic_stationarity_residual'] < 1e-12
                 and r['matrix_relative_error'] < 1e-12
                 and r['source_match_relative_error'] < 1e-12
                 and r['full_nonlinear_Y_residual'] < 1e-10
                 and r['nonlinear_relative_source_shift'] < 1e-8
                 and abs(r['source_nD']) > 0 for r in rows)
    return dict(pass_group=passed, gap_squared_over_M2=float(matrix_a-abs(matrix_b)),
                matched_alpha6_in_M_units=float(alpha6), examples=rows,
                scope='Instantaneous states; no evolved adiabatic orbit or viable cosmology demonstrated.')


def adiabatic_yield_bound():
    """Conditional leading-EFT radiation-envelope bound, not a new release.

    Zero physical velocities and the leading inflationary radial relation
    give H_i=sqrt(2 lambda_X/c_H) r_i. In radiation the quartic conformal
    energy sets R=r_i*(1+2/c_H)**(1/4). Assume |X|<=R/a and constant Sigma.
    Triangle inequality then gives |a^3 n_X|<=12 |alpha6| R^6/H_i.
    """
    mediator_mass, sigma_h, coupling_m, coupling_k = 1e12, 1e11, 5e11, .001
    reduced_planck, lambda_x, degrees, needed_yield = 2.435e18, .001, 106.75, 4.8294009076988175e-9
    sigma = sigma_h/math.sqrt(2)
    gap_squared = mediator_mass**2-2*coupling_m*sigma
    sigma_mass = math.sqrt(2)*sigma_h
    alpha6 = coupling_m*coupling_k**2*sigma/(mediator_mass**4-4*coupling_m**2*sigma**2)

    def evaluate(ratio, hubble_coefficient):
        initial_r = ratio*mediator_mass
        inflation_hubble = math.sqrt(2*lambda_x/hubble_coefficient)*initial_r
        envelope = initial_r*(1+2/hubble_coefficient)**.25
        temperature = (90*reduced_planck**2*inflation_hubble**2/(math.pi**2*degrees))**.25
        entropy = 2*math.pi**2/45*degrees*temperature**3
        number_bound = 12*abs(alpha6)*envelope**6/inflation_hubble
        forced_y_bound = abs(coupling_k)*envelope**3/gap_squared
        return dict(cH=hubble_coefficient, initial_X_over_M=ratio,
                    H_i_GeV=inflation_hubble, conformal_envelope_over_M=envelope/mediator_mass,
                    yield_upper_bound=number_bound/entropy,
                    upper_bound_over_required_YQ=number_bound/entropy/needed_yield,
                    kappa_X2_over_M2=coupling_k*envelope**2/mediator_mass**2,
                    X_frequency_over_Y_gap=math.sqrt(2*lambda_x)*envelope/math.sqrt(gap_squared),
                    X_frequency_over_Sigma_mass=math.sqrt(2*lambda_x)*envelope/sigma_mass,
                    Y_quartic_relative_force_bound=2*forced_y_bound**2/gap_squared,
                    Sigma_relative_forced_shift_bound=coupling_m*forced_y_bound**2/(sigma_mass**2*sigma))

    rows = [evaluate(ratio, ch) for ch in [1.,10.] for ratio in [.1,1.,2.]]
    extrapolations = []
    for ch in [1.,10.]:
        root_ratio = brentq(lambda q: math.log(evaluate(q,ch)['upper_bound_over_required_YQ']), .01, 10000.)
        trial = evaluate(root_ratio, ch)
        trial['outside_matching_regime'] = bool(trial['Y_quartic_relative_force_bound'] > 1
                                              or trial['Sigma_relative_forced_shift_bound'] > 1
                                              or trial['X_frequency_over_Y_gap'] > 1)
        extrapolations.append(trial)
    integral, error = quad(lambda scale: scale**-2, 1., np.inf, epsabs=1e-12)
    measured_power = math.log(evaluate(.2,1.)['yield_upper_bound']/evaluate(.1,1.)['yield_upper_bound'])/math.log(2)
    passed = abs(integral-1) < 1e-12 and abs(measured_power-3.5) < 1e-12 and all(
        x['outside_matching_regime'] for x in extrapolations)
    return dict(pass_group=bool(passed), alpha6_GeV_minus2=alpha6,
                triangle_integral=float(integral), quadrature_error=float(error),
                yield_power_in_initial_amplitude=measured_power,
                assumed_envelope='|X| <= r_i*(1+2/c_H)^(1/4)/a; constant broken Sigma; radiation',
                examples=rows,
                invalid_regime_extrapolation_to_required_yield=extrapolations,
                scope='Conditional leading-EFT bound only. No viable low-debris history or full-action exclusion.')


def main():
    groups = dict(symbolic_currents=symbolic_currents(),
                  complex_force_and_phase_checks=complex_sources(),
                  forced_response_distinction=forced_response(),
                  conditional_adiabatic_yield_bound=adiabatic_yield_bound())
    result = dict(author='Codex / AInstein; same-author mathematical controls',
                  groups_passed=int(sum(x['pass_group'] for x in groups.values())),
                  groups_total=len(groups), groups=groups,
                  physical_scope=dict(new_release_runs=0,new_parameter_scan_points=0,
                      successful_genesis_demonstrated=False,
                      all_preparations_excluded=False,
                      structural_debris_bound_proved=False,
                      additional_closure_prediction_derived=False,
                      pending_debris_review_replicated=False))
    (ROOT/'CHECKS.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(groups_passed=result['groups_passed'],
                          groups_total=result['groups_total']), indent=2))
    if result['groups_passed'] != result['groups_total']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
