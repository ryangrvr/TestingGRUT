"""Reduce the stiff consistency audit to stress, radiation and cold-matter gates."""
import json

import numpy as np
from scipy.integrate import cumulative_trapezoid, simpson
from scipy.optimize import brentq

from stiff_consistency import (H_END, H_R, K_MAX, K_MIN, LATE, NAMES, N_R, P_STIFF,
                               Q_END, TENSOR_PIVOT, Smooth, fit_frozen_power,
                               potential_derivatives, uv_envelope)
from transfer import DM_PER_ENTROPY, G, G_S0, M, MPL, RESULTS, T_CMB, write

RHO_UNIT = H_END**2/(24*np.pi**2*MPL**2)
X_R = LATE.x_reheat
X_EQ = X_R*LATE.reheat_temperature/T_CMB*(G/G_S0)**(1/3)/3401
B_R = np.exp((1-P_STIFF)*N_R)
TENSOR_CAP = 1.1*TENSOR_PIVOT
G_BBN = 10.75
T_BBN = 1e-3
X_BBN = X_R*LATE.reheat_temperature/T_BBN*(G/G_BBN)**(1/3)
Q_BBN = X_BBN*LATE.h_reheat*(T_BBN/LATE.reheat_temperature)**2*np.sqrt(G_BBN/G)
DELTA_N_FACTOR = (8/7)*(11/4)**(4/3)*(G/2)*(G_S0/G)**(4/3)


def moments(model, row, k, occ, energy, coherent):
    """Envelope the part of the stiff era after the sampled ending transient."""
    e_fold = np.linspace(3., 8., 121)
    integral = []
    for time in e_fold:
        b = model.horizon(time)
        field, slope = model.basis(k, time)
        reference = abs(field)**2+abs((slope-field)*b/k)**2
        excitation = 2*(occ+np.sqrt(occ*(1+occ)))*reference
        uv_tail = 8*(energy+coherent)
        integral.append(simpson(k**4*excitation, x=np.log(k))+uv_tail)
    early_phase_bound = max(RHO_UNIT*np.exp(.8*e_fold)*integral)
    # At N>=8 the whole hard band is inside the horizon. A 0.1 H^-1
    # average suppresses the fast coherence by at least the 20 aH/k bound.
    coherence_inverse_moment = (simpson(k**3*np.sqrt(occ*(1+occ)), x=np.log(k))
                               +coherent/K_MAX)
    late_coherence = (8*RHO_UNIT*np.exp(.8*8.)*20*model.horizon(8.)
                      *coherence_inverse_moment)
    hard_radiation = 4*RHO_UNIT*np.exp(.8*N_R)*(row["hard_integral"]+energy)
    # The factor two allows for finite-wavelength stress during entry.
    late_phase_bound = 2*hard_radiation+late_coherence
    return hard_radiation, float(early_phase_bound), float(late_phase_bound)


def local_scale(model):
    n = np.unique(np.r_[np.linspace(-model.width, model.width, 1801),
                        np.linspace(model.width, 3., 401)])
    b, u, u1, u2, _, _ = potential_derivatives(model, n)
    a, a_n, c = -u/2, -u1/(2*b), (u2-u*u)/8
    rho = a*a+b*b*(a_n-a)
    pressure = 4*c/3-a*a/3+b*b*(a_n-a)
    h = b*np.exp(-n)
    # Logarithmic size of local fourth-order terms, in units of rho_bg.
    # The finite GR counterterms are fixed by the spectator EFT convention.
    return float(max(100*RHO_UNIT*np.exp(-4*n)/h**2
                     *(b**4+np.maximum(abs(rho), abs(pressure)))))


def gap_and_cold(k, occ, uv_energy):
    exponent, coeff = LATE.uv_power
    q_hard = K_MIN*Q_END
    integral = coeff*q_hard**(1-exponent)/(1-exponent)
    # The same massless growing mode applies to the tensor and hard longitudinal
    # gap. The factor four covers entry and the low-q transition region.
    gap_relativistic = 4*TENSOR_CAP*integral/(48*X_R**4*LATE.h_reheat**2)
    radiation_entrants = 10*TENSOR_CAP/24*np.log(X_R*LATE.h_reheat/Q_BBN)
    hard_density_unit = M*(M*Q_END)**3/(2*np.pi**2*LATE.entropy_reference*DM_PER_ENTROPY)
    hard_rest_density = hard_density_unit*(simpson(k**3*occ, x=np.log(k))
                                            +uv_energy/K_MAX)

    def free_streaming_upper(q):
        x = max(q, X_R)
        before_radiation = (X_R**(P_STIFF-1)/((P_STIFF-1)*LATE.k_star)
                            +2/(LATE.k_star*Q_END))
        after_radiation = ((x-X_R)+q*np.log(X_EQ/x))/(LATE.k_star*LATE.h_reheat*X_R**2)
        return before_radiation+after_radiation

    q_stream = brentq(lambda q: free_streaming_upper(q)-.01, 128., X_EQ/10)
    q_velocity = .01*X_EQ/np.sqrt(1-.01**2)
    q_pass = min(q_stream, q_velocity)
    density_unit = M*M*MPL*MPL/(16*LATE.entropy_reference*DM_PER_ENTROPY)
    q_tail_density = density_unit*TENSOR_CAP*coeff*q_pass**(-exponent)/exponent
    return dict(gap_relativistic_upper=float(gap_relativistic),
                additional_RD_GW_upper=float(radiation_entrants),
                hard_rest_density_over_target_upper=float(hard_rest_density),
                free_stream_q_limit=float(q_stream), velocity_q_limit=float(q_velocity),
                late_rejected_rest_density_fraction_upper=float(q_tail_density+hard_rest_density),
                soft_hard_boundary_q=float(q_hard))


def summarize(name):
    model = Smooth(name)
    fine = json.loads((RESULTS/f"consistency_{name}_fine.json").read_text())
    coarse = json.loads((RESULTS/f"consistency_{name}_coarse.json").read_text())
    initial_uv = json.loads((RESULTS/"consistency_uv_attempt.json").read_text())[name]
    with np.load(RESULTS/f"consistency_{name}_fine.npz") as packet:
        k = packet["k_end_units"]
        occ = packet["occupation"]
    envelopes = [uv_envelope(model, K_MAX, mesh_scale=scale) for scale in (1, 2)]
    uv_energy = min(max(row[key] for row in envelopes)
                    for key in ("absolute_energy", "oscillatory_energy"))
    coherent = max(row["coherent_pressure_integral"] for row in envelopes)
    hard_rel, early_stiff, late_stiff = moments(model, fine, k, occ, uv_energy, coherent)
    band = gap_and_cold(k, occ, uv_energy)
    gap = band["gap_relativistic_upper"]
    rd_entrants = band["additional_RD_GW_upper"]
    local = local_scale(model)
    # Before the end blend the adiabatic local stress is of order H^4.
    # During the blend the absolute unaveraged hard tail has this larger cap.
    after_start = np.linspace(-model.width, 3., 801)
    c_max = max(RHO_UNIT*np.exp(-4*n-2*model.background(n)[0]) for n in after_start)
    ending = max(fine["max_early_hard_energy_abs_fraction"],
                 fine["max_early_hard_pressure_abs_fraction"])+8*c_max*(uv_energy+coherent)+local
    stiff_upper = max(ending, early_stiff, late_stiff)+gap+fine["max_soft_energy_fraction"]
    vector_radiation = hard_rel+gap+fine["max_soft_energy_fraction"]
    gw_radiation = 2*(hard_rel+gap)+rd_entrants
    delta_neff = DELTA_N_FACTOR*(vector_radiation+gw_radiation)
    # At the end of the spectator window total standard matter/radiation is 1%.
    # The frozen cosmology has a cold matter share below 0.85 of that total.
    early_radiation_vector_upper = .0085+vector_radiation
    grid_shift = abs(fine["hard_integral"]/coarse["hard_integral"]-1)
    uv_fraction = uv_energy/fine["hard_integral"]
    uv_previous = initial_uv["ultraviolet_hard_integral_upper"]
    cold_miss = band["late_rejected_rest_density_fraction_upper"]
    density = fine["vector_rest_density_over_frozen_target"]
    overlap = k <= .0017
    q_overlap = k[overlap]*Q_END
    exponent, coefficient = LATE.uv_power
    expected_occ = (np.pi**2*MPL**2/(8*M*M*Q_END**3)*coefficient
                    *q_overlap**(-exponent)*fit_frozen_power(q_overlap, model)/k[overlap]**3)
    overlap_error = max(abs(occ[overlap]/expected_occ-1))
    seams = np.r_[np.linspace(-.3, .3, 601), np.linspace(N_R-.9, N_R+.9, 601)]
    dn = 3e-5
    continuity_error = max(abs((model.background(n+dn)[0]-model.background(n-dn)[0])
                               /(2*dn)+model.background(n)[1]) for n in seams)
    k_cdf = cumulative_trapezoid(k**4*occ, x=np.log(k), initial=0.)
    k_quantiles = np.exp(np.interp([.05, .5, .95], k_cdf/k_cdf[-1], np.log(k)))
    checks = dict(stress=early_radiation_vector_upper < .01 and
                  2*stiff_upper+rd_entrants < .01,
                  dark_radiation=delta_neff < .1, mesh_and_start=grid_shift < .05,
                  ultraviolet=uv_fraction < .05, abundance=abs(density-1) < .01,
                  cold=1-cold_miss/density >= .99,
                  effective_background=fine["epsilon_end_blend_min"] > 0 and
                  fine["epsilon_end_blend_max"] < 3 and model.rate > .8 and
                  abs(fine["radiation_asymptote_H_ratio"]-1) < 1e-10)
    return dict(checks={k: bool(v) for k, v in checks.items()},
                vector_rest_density_over_target=density,
                hard_to_growing_mode_overlap_max_relative_error=float(overlap_error),
                H_epsilon_derivative_check_absolute_error=float(continuity_error),
                hard_momentum_K_quantiles_05_50_95=k_quantiles.tolist(),
                hard_integral_mesh_start_shift=grid_shift,
                fraction_hard_GW_energy_above_original_end=fine["hard_fraction_above_original_end"],
                initial_numerical_derivative_UV_bound_fraction=uv_previous/fine["hard_integral"],
                ODE_derivative_UV_bounds=envelopes, omitted_hard_energy_fraction_upper=uv_fraction,
                local_fourth_order_stress_scale=local,
                max_ending_hard_stress_fraction_upper=float(ending),
                max_stiff_phase_hard_stress_fraction_upper=max(early_stiff, late_stiff),
                vector_stress_fraction_upper_through_early_radiation=early_radiation_vector_upper,
                vector_stress_fraction_upper_through_reheating=stiff_upper,
                tensor_stress_fraction_upper_through_early_radiation=2*stiff_upper+rd_entrants,
                vector_extra_radiation_fraction_upper=vector_radiation,
                gravitational_wave_extra_radiation_fraction_upper=gw_radiation,
                delta_N_eff_upper=float(delta_neff), **band)


def main():
    rows = {name: summarize(name) for name in NAMES}
    result = dict(sharp_background_passed=False,
                  sharp_failure="Step in a''/a produces logarithmically divergent excitation energy.",
                  H_end_gev=H_END, H_reheat_reference_gev=H_R, horizon_k_BBN_over_k_star=Q_BBN,
                  k_reheat_over_k_end=B_R, x_equality=X_EQ, relativistic_to_delta_N_eff=DELTA_N_FACTOR,
                  effective_smooth_gates_pass=all(all(row["checks"].values()) for row in rows.values()),
                  rows=rows)
    write(RESULTS/"stiff_consistency.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
