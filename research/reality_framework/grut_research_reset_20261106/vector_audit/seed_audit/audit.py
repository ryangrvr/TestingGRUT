"""Causal seed inventory and perturbation audit at the frozen vector fit."""
import json
import os
from pathlib import Path

import numpy as np
from scipy.integrate import quad, simpson
from scipy.optimize import brentq
from scipy.special import zeta

HERE = Path(__file__).resolve().parent
PACKET = Path(os.environ.get("GRUT_REPLICATION_PACKET", HERE.parent / "replication"))
INPUTS = json.loads((PACKET / "frozen_inputs.json").read_text())
ARCHIVE = json.loads((PACKET / "archived_reference.json").read_text())
UNITS = INPUTS["units"]
PLANCK = UNITS["Mbar_GeV"]
LN10 = np.log(10.)
HBAR_C_GEV_METRE = 1.973269804e-16
GEV_JOULE = 1.602176634e-10
SPHERE = 4*np.pi/3
HARD_BIN = (8., 16.)
SCATTER_START, SCATTER_STOP = 4., 4.1


def log_fraction(logit):
    return -np.logaddexp(0., -logit)


def seed_inventory(hubble, p, sigma, start, midpoint, a_end):
    log_f = log_fraction(sigma*(start-midpoint))
    rho_total = 3*(PLANCK*hubble)**2
    log_rho = np.log(rho_total)+log_f
    log_temperature = (log_rho+np.log(30/(np.pi**2*UNITS["g_energy"])))/4
    log_number_density = 3*log_temperature+np.log(UNITS["g_energy"]*zeta(3)/np.pi**2)
    log_volume = np.log(SPHERE)-3*np.log(hubble)
    log_count = log_rho+log_volume-np.log(hubble)
    horizon_metre = HBAR_C_GEV_METRE/hubble
    causal_ratio = np.expm1((p-1)*(midpoint-start))/(p-1)
    pivot_ratio = (a_end*np.exp(start)*hubble*UNITS["GeV_to_Mpc_inverse"]
                   / UNITS["k_pivot_Mpc_inverse"])
    one_quantum_fraction = hubble**2/(4*np.pi*PLANCK**2)
    return dict(N_initial=start, N_half_radiation=midpoint, H_initial_GeV=hubble,
                total_density_GeV4=rho_total, log10_seed_fraction=log_f/LN10,
                log10_seed_density_GeV4=log_rho/LN10,
                log10_seed_density_J_m3=(log_rho+np.log(GEV_JOULE)-3*np.log(HBAR_C_GEV_METRE))/LN10,
                hubble_radius_m=horizon_metre, hubble_volume_m3=SPHERE*horizon_metre**3,
                log10_seed_energy_in_hubble_sphere_GeV=(log_rho+log_volume)/LN10,
                log10_seed_quanta_per_hubble_sphere_upper=log_count/LN10,
                log10_seed_quanta_in_nominal_conversion_past_cone_upper=(
                    log_count+3*np.log(causal_ratio))/LN10,
                log10_seed_quanta_in_pivot_sphere_upper=(log_count+3*np.log(pivot_ratio))/LN10,
                initial_radius_of_nominal_conversion_past_cone_upper_m=horizon_metre*causal_ratio,
                initial_pivot_sphere_radius_m=horizon_metre*pivot_ratio,
                formal_temperature_GeV=np.exp(log_temperature),
                formal_temperature_over_H=np.exp(log_temperature)/hubble,
                log10_formal_thermal_quanta_per_hubble_sphere=(log_number_density+log_volume)/LN10,
                one_horizon_quantum_fraction=one_quantum_fraction,
                homogeneous_midpoint_with_one_quantum_per_initial_hubble_sphere=(
                    start+(np.log1p(-one_quantum_fraction)-np.log(one_quantum_fraction))/sigma))


def hard_moments(name, resolution, h_end):
    with np.load(PACKET / "measured" / f"{name}_{resolution}.npz") as data:
        k, occupation = data["hard_K"], data["hard_occupation"]
    inside = (k > HARD_BIN[0]) & (k < HARD_BIN[1])
    band = np.r_[HARD_BIN[0], k[inside], HARD_BIN[1]]
    occ = np.interp(np.log(band), np.log(k), occupation)
    momentum = band*h_end*np.exp(-SCATTER_START)
    return tuple(float(simpson(occ*momentum**power, x=np.log(band))/(2*np.pi**2))
                 for power in (3, 4, 5))


def gauge_injection(name, resolution, h_end, hubble, sigma, midpoint):
    number, rho, second = hard_moments(name, resolution, h_end)
    hi = hubble(SCATTER_START)
    h_stop = hubble(SCATTER_STOP)
    dt = SCATTER_STOP-SCATTER_START
    # Identical initial pairs count once. The cross section for one
    # Maxwell species already includes the identical final-state factor.
    # Angle selection is a restriction of the isotropic pair measure.
    angular = 13/24
    event_rate = angular*rho*rho/(480*np.pi*PLANCK**4)
    power = 2*angular*rho*second/(480*np.pi*PLANCK**4)
    photon_density = power*quad(
        lambda t: np.exp(-9*t-4*(dt-t))/hubble(SCATTER_START+t),
        0., dt, epsabs=1e-28, epsrel=1e-11)[0]
    fraction = photon_density/(3*(PLANCK*h_stop)**2)
    log_target = log_fraction(sigma*(SCATTER_STOP-midpoint))
    # The following midpoint is only the homogeneous fluid thought experiment.
    ideal_fluid_midpoint = SCATTER_STOP+(np.log1p(-fraction)-np.log(fraction))/sigma
    event_count = event_rate*quad(
        lambda t: np.exp(-8*t+3*t)/hubble(SCATTER_START+t),
        0., dt, epsabs=1e-28, epsrel=1e-11)[0]*SPHERE/hi**3
    return dict(N_interval=[SCATTER_START, SCATTER_STOP], K_bin=list(HARD_BIN),
                restricted_mu=[-.5, .5],
                longitudinal_number_in_initial_hubble_sphere=SPHERE*number/hi**3,
                longitudinal_energy_fraction_at_N4=rho/(3*(PLANCK*hi)**2),
                tensor_energy_fraction_at_N4_two_polarizations=2*rho/(3*(PLANCK*hi)**2),
                band_min_p_over_H_at_N4=HARD_BIN[0]*h_end*np.exp(-SCATTER_START)/hi,
                band_max_p_over_Mbar_at_N4=HARD_BIN[1]*h_end*np.exp(-SCATTER_START)/PLANCK,
                visible_injection_energy_per_H_over_background_at_N4=power/(hi*3*(PLANCK*hi)**2),
                visible_fraction_at_N4p1_before_amplification=fraction,
                log10_target_fraction_at_N4p1=log_target/LN10,
                log10_injection_to_frozen_seed_ratio=(np.log(fraction)-log_target)/LN10,
                expected_pair_events_in_initial_comoving_hubble_sphere_during_interval=event_count,
                midpoint_if_direct_injection_immediately_forms_homogeneous_Q_fluid=ideal_fluid_midpoint,
                advance_in_that_fluid_thought_experiment_efolds=midpoint-ideal_fluid_midpoint)


def independent_scalar_response(p, sigma, start, midpoint, r_frozen):
    z_initial = sigma*(start-midpoint)
    f_initial = np.exp(log_fraction(z_initial))
    contrast = 2*p-4
    duration = midpoint+3-start

    def log_density(z, n):
        return (-2*p*n+contrast/sigma*(np.logaddexp(0., z+sigma*n)
                                      -np.logaddexp(0., z)))

    target = log_density(z_initial, duration)

    def shift(dz):
        return brentq(lambda n: log_density(z_initial+dz, n)-target,
                      duration-1, duration+1, xtol=1e-13)-duration

    analytic = contrast/(4*sigma)*(1-f_initial)
    step = .001
    derivative = (shift(step)-shift(-step))/(2*step)
    assert abs(derivative/analytic-1) < 1e-8
    a_s = UNITS["A_s_pivot"]
    examples = []
    for rms in (1e-4, 1e-3, 1e-2):
        curvature_power = (derivative*rms)**2
        power_ratio = curvature_power/a_s
        examples.append(dict(seed_logit_root_power=rms, linear_delta_N_root_power=derivative*rms,
                             positive_finite_shift=shift(rms), negative_finite_shift=shift(-rms),
                             added_scalar_power=curvature_power, added_over_frozen_scalar_power=power_ratio,
                             r_at_fixed_tensors_and_changed_scalar_power=r_frozen/(1+power_ratio)))
    return dict(logit_to_delta_N=derivative, analytic_logit_to_delta_N=analytic,
                r_frozen=r_frozen,
                max_seed_logit_root_power_for_10_percent_of_frozen_scalar_power=np.sqrt(.1*a_s)/derivative,
                examples=examples)


def audit_case(name):
    case = INPUTS["histories"][name]
    reference = ARCHIVE["histories"][name]
    h_end, a_end = reference["H_end_GeV"], reference["a_end"]
    p, sigma = case["pre_radiation_H_power"], case["reheat_rate_per_efold"]
    start = case["inflation_end_blend_halfwidth"]
    nominal = case["nominal_reheat_temperature_GeV"]
    h_r = np.sqrt(np.pi**2*UNITS["g_energy"]/90)*nominal**2/PLANCK
    midpoint = np.log(h_end/h_r)/p

    def hubble(n):
        return h_end*np.exp(-p*n+(p-2)/sigma*np.logaddexp(0., sigma*(n-midpoint)))

    seed = seed_inventory(hubble(start), p, sigma, start, midpoint, a_end)
    production = {grid: gauge_injection(name, grid, h_end, hubble, sigma, midpoint)
                  for grid in ("coarse", "fine")}
    return dict(seed=seed, hard_band_and_one_Maxwell_example=production,
                scalar_response=independent_scalar_response(p, sigma, start, midpoint, reference["r"]))


def main():
    results = dict(replication_archive="7294e556154aa4c27ffab782ae30cd281ebe7d5d",
                   cases={name: audit_case(name) for name in ("narrow", "wide")})
    target = HERE / "results.json"
    target.write_text(json.dumps(results, indent=2, allow_nan=False)+"\n")
    print(json.dumps(results, indent=2), flush=True)


if __name__ == "__main__":
    main()
