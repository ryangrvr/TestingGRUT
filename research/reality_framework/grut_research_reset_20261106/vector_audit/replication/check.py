"""Compare frozen measurements and run separate stress and perturbation checks."""
import json

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.interpolate import make_interp_spline
from scipy.optimize import brentq
from scipy.special import expit

from replicate import (ROOT, UNITS, FROZEN, Expansion, MASS, PLANCK, G_RHO, G_S,
                       K_HARD_MIN, K_HARD_MAX, MPC, S_TODAY, TARGET, high_q_tail)


def omitted_hard_integral(bg, samples):
    """Instantaneous-frequency mixing and two integrations by parts."""
    knots = [-18., -5., -bg.width, bg.width, 1., 5., 8., bg.n_reh-2.,
             bg.n_reh+2., bg.n_reh+8.]
    abs_u = 0.
    integrals = np.zeros(4)
    for left, right in zip(knots[:-1], knots[1:]):
        n = np.linspace(left, right, samples)
        logh, epsilon = np.array([bg.state(x) for x in n]).T
        b = np.exp(n+logh)
        u = b*b*(2-epsilon)
        abs_u = max(abs_u, max(abs(u)))
        spline = make_interp_spline(n, u, k=5)
        e_prime = make_interp_spline(n, epsilon, k=5)(n, 1)
        n1, n2, n3 = (spline(n, order) for order in (1, 2, 3))
        d1 = b*n1
        d2 = b*b*(n2+(1-epsilon)*n1)
        d3 = b**3*(n3+3*(1-epsilon)*n2+
                   (2*(1-epsilon)**2-e_prime)*n1)
        integrals += [simpson(value/b, x=n) for value in
                      (abs(d1), abs(d3), abs(d1*d2), abs(d1)**3)]

    # Include margin for quadrature and the negligible outer adiabatic ends.
    u1, u3, u12, u111 = 2*integrals
    k = K_HARD_MAX
    w = np.sqrt(1-abs_u/k**2)
    j = u1/(4*w*w)
    if j/k**2 >= 1:
        raise ValueError("The ultraviolet band has large instantaneous mixing")
    c4 = u3/(16*w**4)
    c6 = 5*u12/(16*w**6)+j**3
    c8 = 9*u111/(32*w**8)
    # For all K >= k, |beta| <= (c4+c6/k²+c8/k⁴)/K⁴.
    return float((c4+c6/k**2+c8/k**4)**2/(4*k**4))


def radiation_check(bg, measured, k, occupation):
    uv_estimates = [omitted_hard_integral(bg, size) for size in (481, 961)]
    uv_upper = max(uv_estimates)
    band = measured["hard_K4_occupation_integral"]
    hard_unit = (bg.h_end/PLANCK)**2*np.exp((2*bg.p-4)*bg.n_reh)/(6*np.pi**2)
    hard_each = hard_unit*(band+uv_upper)

    xr = np.exp(bg.n_reh-bg.n_star)
    h_r = bg.h_reh/MASS
    q_edge = K_HARD_MIN*bg.k_end/bg.k_star
    exponent, coefficient = high_q_tail(bg)
    pt_cap = 1.1*UNITS["A_s_pivot"]*measured["r_at_frozen_As"]
    gap_each = (4*pt_cap/(48*xr**4*h_r*h_r)*coefficient
                *q_edge**(1-exponent)/(1-exponent))
    x_bbn = xr*bg.case["nominal_reheat_temperature_GeV"]/.001*(G_S/10.75)**(1/3)
    h_bbn = h_r*(.001/bg.case["nominal_reheat_temperature_GeV"])**2*np.sqrt(10.75/G_RHO)
    extra_gw = 10*pt_cap/24*np.log(xr*h_r/(x_bbn*h_bbn))
    factor = (8/7)*(11/4)**(4/3)*(G_RHO/2)*(UNITS["g_entropy_today"]/G_S)**(4/3)
    hard_number = simpson(k**3*occupation, x=np.log(k))+uv_upper/K_HARD_MAX
    hard_rest_ratio = MASS*(bg.k_end/MPC)**3*hard_number/(2*np.pi**2*S_TODAY*TARGET)
    return dict(uv_hard_integral_estimates=uv_estimates,
                omitted_to_resolved_hard_upper=uv_upper/band,
                resolved_R_one_polarization=float(hard_unit*band),
                hard_R_bound_one_polarization=float(hard_each),
                hard_vector_rest_density_over_target_upper=float(hard_rest_ratio),
                growing_mode_gap_R_bound_one_polarization=float(gap_each),
                Delta_N_eff_upper=float(factor*(3*(hard_each+gap_each)+extra_gw
                                               +measured["max_resolved_soft_energy_fraction"])))


def perturbative_closure(bg):
    """Separate-universe test of a specified local barotropic-fluid closure."""
    rate = bg.rate
    contrast = 2*bg.p-4
    initial, final = -1.5, 2.
    step = 1e-4

    def background(n):
        return (-2*bg.p*n+contrast/rate*np.logaddexp(0., rate*n), rate*n)

    def trajectory(shift, seed):
        log_rho, z = background(initial+shift)

        def conservation(_n, y):
            f = expit(y[1])
            # Q / (H rho f(1-f)) and the differential redshift term.
            q_scaled = rate-contrast
            return -2*bg.p+contrast*f, q_scaled+contrast

        sol = solve_ivp(conservation, (initial, final+1.), (log_rho, z+seed),
                        dense_output=True, rtol=2e-12, atol=2e-13,
                        max_step=.025)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.sol

    plus, minus = trajectory(step, 0.), trajectory(-step, 0.)
    seed_plus, seed_minus = trajectory(0., step), trajectory(0., -step)

    def crossing(solution, target):
        return brentq(lambda n: solution(n)[0]-target, initial+.1, final+.1)

    error = []
    for fixed in (-.75, 0., 1.):
        target = background(fixed)[0]
        response = (crossing(plus, target)-crossing(minus, target))/(2*step)
        error.append(abs(response+1))
    target = background(final)[0]
    response = (crossing(seed_plus, target)-crossing(seed_minus, target))/(2*step)
    expected = contrast/(4*rate)*(1-expit(rate*initial))
    f_blend = -np.logaddexp(0., -rate*(bg.width-bg.n_reh))/np.log(10.)
    return dict(adiabatic_uniform_density_time_shift_max_error=float(max(error)),
                independent_seed_to_delta_N=float(response),
                analytic_seed_to_delta_N=float(expected),
                log10_radiation_seed_at_blend=float(f_blend))


def main():
    directory = ROOT / "measured"
    measured = {name: {grid: json.loads((directory / f"{name}_{grid}.json").read_text())
                       for grid in ("fine", "coarse")} for name in FROZEN["histories"]}
    audit = {}
    for name in ("narrow", "wide"):
        model = Expansion(name, True)
        with np.load(directory / f"{name}_fine.npz") as packet:
            radiation = radiation_check(model, measured[name]["fine"],
                                        packet["hard_K"], packet["hard_occupation"])
        audit[name] = dict(radiation=radiation, perturbations=perturbative_closure(model))
    (ROOT / "independent_checks.json").write_text(json.dumps(audit, indent=2)+"\n")

    # Open the reference only after all fixed-input measurements exist.
    archive = json.loads((ROOT / "archived_reference.json").read_text())
    archived = archive["histories"]
    def fit(rows, case, grid, key):
        return rows[case][grid][key]["ratio"]

    def fractional(x, y):
        return x/y-1

    comparisons = {}
    baseline = measured["prompt"]["fine"]
    control = measured["early_matter"]["fine"]
    for name, values in measured.items():
        fine, coarse = values["fine"], values["coarse"]
        ref = archived[name]
        ratios = {}
        for comparator in ("prompt", "early_matter"):
            obtained = fine["fitted"]["ratio"]/fit(measured, comparator, "fine", "fitted")
            previous = coarse["fitted"]["ratio"]/fit(measured, comparator, "coarse", "fitted")
            reference = ref["density_over_target"]/archived[comparator]["density_over_target"]
            ratios[comparator] = dict(measured=obtained, reference=reference,
                                      fractional_error=fractional(obtained, reference),
                                      grid_shift=fractional(obtained, previous))
        tests = dict(
            fitted_density_error=fractional(fine["fitted"]["ratio"], ref["density_over_target"]),
            r_error=fractional(fine["r_at_frozen_As"], ref["r"]),
            pivot_efolds_difference=fine["pivot_efolds"]-ref["pivot_efolds"],
            entropy_error=fractional(fine["entropy_per_a_star_cubed_GeV3"],
                                     ref["entropy_per_a_star_cubed_GeV3"]))
        gates = dict(abundance=abs(tests["fitted_density_error"]) < .005 and
                     max(abs(x["fractional_error"]) for x in ratios.values()) < .005,
                     abundance_refinement=max(abs(x["grid_shift"])
                                              for x in ratios.values()) < .002,
                     tensor_normalization=abs(tests["r_error"]) < .002,
                     pivot=abs(tests["pivot_efolds_difference"]) < .001,
                     entropy=abs(tests["entropy_error"]) < .0002,
                     longitudinal_stress_identity=fine["longitudinal_conservation_relative_residual"] < 1e-8)
        if name in audit:
            key = "hard_K4_occupation_integral"
            tests.update(hard_resolved_error=fractional(fine[key], ref[key]),
                         hard_grid_shift=fractional(fine[key], coarse[key]))
            gates.update(hard_resolved=abs(tests["hard_resolved_error"]) < .01,
                         hard_grid=abs(tests["hard_grid_shift"]) < .005,
                         hard_tail=audit[name]["radiation"]["omitted_to_resolved_hard_upper"] < .05,
                         radiation=audit[name]["radiation"]["Delta_N_eff_upper"] < .1,
                         long_wavelength_adiabatic_mode=audit[name]["perturbations"][
                             "adiabatic_uniform_density_time_shift_max_error"] < 1e-6)
        comparisons[name] = dict(fitted_abundance_ratios=ratios, differences=tests,
                                 matched_spectrum_to_prompt=fine["matched_stiff_power"]["ratio"]/
                                 baseline["matched_stiff_power"]["ratio"],
                                 matched_spectrum_to_early_matter=fine["matched_stiff_power"]["ratio"]/
                                 control["matched_stiff_power"]["ratio"],
                                 flat_spectrum_to_prompt=fine["common_flat_power"]["ratio"]/
                                 baseline["common_flat_power"]["ratio"],
                                 gates={k: bool(v) for k, v in gates.items()})
    document = dict(archive=archive["archived_commit"], comparisons=comparisons,
                    passed=all(all(row["gates"].values()) for row in comparisons.values()))
    (ROOT / "comparison.json").write_text(json.dumps(document, indent=2, allow_nan=False)+"\n")
    print(json.dumps(document, indent=2), flush=True)


if __name__ == "__main__":
    main()
