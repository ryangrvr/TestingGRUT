"""Reproduce the C3 reporter phase and calculate its static angular response."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import cumulative_simpson, solve_ivp
from scipy.interpolate import BPoly, CubicSpline, PPoly
from scipy.special import eval_legendre, spherical_jn, spherical_yn


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REF = ROOT / "grut_candidate_1_20261018/code/results_profile.npz"
C1 = ROOT / "grut_c1_quantum_channels_20261127"
C2 = ROOT / "grut_c2_continuum_stability_20261127"


def load_json(path):
    return json.loads(path.read_text())


def archived_reference():
    with np.load(REF) as data:
        r, y = data["rho"], data["y"]
    spline = CubicSpline(r, y[2])
    return lambda x: np.where(x <= r[-1], spline(np.minimum(x, r[-1])), 1.)


def parent_profiles(charge):
    with np.load(C1 / "PROFILES.npz") as data:
        r = data[f"fine_{charge}_knots"]
        coeff = data[f"fine_{charge}_coefficients"]
    f, _, s, ds = PPoly(coeff, r)(r).T.copy()
    second = (s*s-1+f*f)*s
    second[0] /= 3
    second[1:] -= 2*ds[1:]/r[1:]
    # Match the C2 stationary Hermite predictor, including its vacuum taper.
    ds[0] = 0
    rho = np.r_[r, 60.125]
    values = list(zip(np.r_[s, 1.], np.r_[ds, 0.], np.r_[second, 0.]))
    field = BPoly.from_derivatives(rho, values)
    return r, PPoly(coeff, r), lambda x: np.where(x <= rho[-1], field(np.minimum(x, rho[-1])), 1.)


def swave(k, field, controls, stage):
    radius = controls["radii"][stage]

    def flow(r, phase):
        potential = 0.01*(field(r)**2-1)
        return -potential*np.sin(k*r+phase)**2/k

    result = solve_ivp(flow, (0., radius), [0.], method="DOP853",
                       rtol=controls["tolerances"][stage],
                       atol=controls["absolute_tolerances"][stage],
                       max_step=controls["step_caps"][stage])
    if not result.success:
        raise ArithmeticError(result.message)
    return float(result.y[0, -1])


def outward_match(k, r, potential, ell, h, base_origin):
    """Numerov with regular inner seed and exterior Riccati-Bessel matching."""
    # Keep h*sqrt(ell*(ell+1))/rho small at the first Numerov step.
    starts = np.maximum(0, 2*(ell+1)-round(base_origin/h)).astype(int)
    core_k = math.sqrt(k*k-potential[0])
    arguments = core_k*r[starts]
    next_arguments = core_k*r[starts+1]
    previous = np.ones_like(ell, dtype=float)
    current = next_arguments*spherical_jn(ell, next_arguments)/(arguments*spherical_jn(ell, arguments))
    if not np.all(np.isfinite(current)):
        raise ArithmeticError("Regular inner seed underflow")
    centrifugal = np.outer(1/r**2, ell*(ell+1))
    a = h*h*(centrifugal+potential[:, None]-k*k)/12
    for j in range(1, len(r)-1):
        active = starts < j
        recurrence = (2*(1+5*a[j])*current-(1-a[j-1])*previous)/np.where(active, 1-a[j+1], 1)
        previous = np.where(active, current, previous)
        current = np.where(active, recurrence, current)
        if j % 128 == 0:
            norm = np.maximum(np.maximum(1, abs(current)), abs(previous))
            previous /= norm
            current /= norm

    z, z0 = k*r[-1], k*r[-2]
    jf, j0 = z*spherical_jn(ell, z), z0*spherical_jn(ell, z0)
    nf, n0 = z*spherical_yn(ell, z), z0*spherical_yn(ell, z0)
    phase = np.arctan((current*j0-previous*jf)/(current*n0-previous*nf))
    if not np.all(np.isfinite(phase)):
        raise ArithmeticError("Nonfinite exterior phase")
    return phase


def partial_waves(k, field, method):
    h = method["step"]
    origin = method["origin"]
    r = origin + np.arange(round((method["radius"]-origin)/h)+1)*h
    ell = np.arange(method["ell_max"]+1)
    potential = 0.01*(field(r)**2-1)
    raw = outward_match(k, r, potential, ell, h, origin)
    vacuum = outward_match(k, r, np.zeros_like(r), ell, h, origin)
    return raw-vacuum, vacuum


def cross_sections(k, phases, angles):
    ell = np.arange(len(phases))
    partials = 4*np.pi*(2*ell+1)*np.sin(phases)**2/(k*k)
    amplitude = ((2*ell+1)*np.exp(1j*phases)*np.sin(phases))/k
    differential = []
    for a in angles:
        legendre = eval_legendre(ell, np.cos(np.deg2rad(a)))
        differential.append(float(abs(np.dot(amplitude, legendre))**2))
    return partials, differential


def match_observables(settings, bound, reference, field):
    rows = []
    for k in settings["momenta_in_g_X_v"]:
        inherited = swave(k, reference, settings["continuum_phase"], 1)
        target0, target = [swave(k, field, settings["continuum_phase"], i) for i in (0, 1)]
        coarse, coarse_vac = partial_waves(k, field, settings["numerov"]["coarse"])
        fine, fine_vac = partial_waves(k, field, settings["numerov"]["fine"])
        angles = settings["angles_degrees"]
        fine_cross, fine_diff = cross_sections(k, fine, angles)
        coarse_cross, coarse_diff = cross_sections(k, coarse, angles)
        # Bound the exact-vs-predictor potential by the C2 ball in L1.
        potential_bound = 0.04*13*bound
        s_error = potential_bound/k * math.exp((0.03*60.125+potential_bound)/k)
        swave_error = abs(fine[0]-target)
        rows.append({
            "k": k, "reference_delta0": inherited, "certified_target_predictor_delta0": target,
            "swave_predictor_numerical_change": abs(target-target0),
            "exact_profile_delta0_analytic_bound": s_error,
            "sigma0": 4*np.pi/k**2*math.sin(target)**2,
            "sigma_total_retained": float(sum(fine_cross)),
            "sigma_total_resolution_change": float(abs(sum(coarse_cross)-sum(fine_cross))),
            "sigma_total_last_16_wave_increment": float(sum(fine_cross[len(coarse):])),
            "s_wave_numerov_error_from_phase": swave_error,
            "max_phase_change_on_common_waves": float(np.max(abs(fine[:len(coarse)]-coarse))),
            "maximum_raw_vacuum_phase": float(max(np.max(abs(fine_vac)), np.max(abs(coarse_vac)))),
            "sigma_ell": fine_cross.tolist(),
            "delta_ell": fine.tolist(),
            "angle_degrees": angles, "dsigma_domega_at_angles": fine_diff,
            "dsigma_domega_resolution_change": [abs(a-b) for a,b in zip(fine_diff, coarse_diff)],
            "swave_agreement": bool(swave_error < settings["numerov"]["required_swave_agreement"]),
        })
        print(f"k={k} reference phase={inherited:.12g} target phase={target:.12g} "
              f"sigma_total={sum(fine_cross):.12g}", flush=True)
    return rows


def physical_scales(settings, rows, c1, knots, parent, field):
    q = settings["target_charge"]
    constants = settings["constants_si"]
    c, h, G, energy_unit = (constants[key] for key in ("c", "h", "G", "GeV_joules"))
    hbar = h/(2*np.pi)
    m_planck_gev = math.sqrt(hbar*c**5/(8*np.pi*G))/energy_unit
    unit = 1e-6*m_planck_gev
    length = hbar*c/(unit*energy_unit)
    E = c1["exact_charge_spectra"][str(q)]["E0plus1_finite"]
    pole = c1["direct_emission"]["pole"]
    r = np.linspace(0, knots[-1], 12001)
    f = parent(r)[:, 0]
    integral = cumulative_simpson(r*r*f*f, x=r, initial=0)
    r90 = float(np.interp(0.9*integral[-1], integral, r))
    k = settings["kinematic_shift"]["fixed_k"]
    row = next(item for item in rows if item["k"] == k)
    m_d = math.sqrt(1.01)
    m_d_pole = settings["kinematic_shift"]["m_D_pole_over_v"] / 0.1
    probe_energy = math.sqrt(k*k+1.01)
    at_same_energy = math.sqrt(probe_energy**2-m_d_pole**2)-k
    cm_shift = k*E/math.sqrt(E*E+m_d*m_d+2*E*probe_energy)-k
    step = settings["kinematic_shift"]["dphase_step"]
    derivative = (swave(k+step, field, settings["continuum_phase"], 1)-
                  swave(k-step, field, settings["continuum_phase"], 1))/(2*step)
    return {
        "reduced_planck_mass_gev": m_planck_gev, "g_X_v_gev": unit, "rho_unit_m": length,
        "target_rest_energy_gev": E*unit, "target_rest_energy_j": E*unit*energy_unit,
        "target_mass_kg": E*unit*energy_unit/c**2, "target_charge_R90_rho": r90,
        "target_charge_R90_m": r90*length, "free_carrier_rest_energy_j": q*pole*unit*energy_unit,
        "binding_release_for_cold_free_total_j": (q*pole-E)*unit*energy_unit,
        "cold_single_capture_release_over_g_X_v": -c1["direct_emission"]["value"],
        "binding_scalar_mass_tree_over_g_X_v": math.sqrt(2),
        "reporter_pair_threshold_j": 2*m_d_pole*unit*energy_unit,
        "probe_total_energy_j_at_k5": probe_energy*unit*energy_unit,
        "probe_p_gev_at_k5": k*unit, "probe_speed_over_c_at_k5": k/probe_energy,
        "sigma0_m2_at_k5": row["sigma0"]*length**2,
        "sigma_total_retained_m2_at_k5": row["sigma_total_retained"]*length**2,
        "ideal_signal_one_fluence_m_minus_2_for_sigma0": 1/(row["sigma0"]*length**2),
        "ideal_signal_one_fluence_m_minus_2_for_total": 1/(row["sigma_total_retained"]*length**2),
        "s_wave_elastic_unitarity_cap_m2_at_k5": 4*np.pi*length**2/k**2,
        "max_backscatter_recoil_estimate_over_g_X_v": 2*k*k/E,
        "one_loop_reporter_vacuum_pole_fractional_shift": m_d_pole/m_d-1,
        "phase_slope_wrt_k_at_k5": derivative,
        "kinematic_phase_shift_pole_at_fixed_probe_energy": derivative*at_same_energy,
        "kinematic_phase_shift_cm_recoil_estimate": derivative*cm_shift,
        "alpha_X_D_scalar_beta_at_matching": 2*0.01*0.0001/(16*np.pi**2),
        "alpha_X_D_leading_log_at_5g_X_v": 2*0.01*0.0001*math.log(k)/(16*np.pi**2),
    }


def main(destination):
    protocol = HERE / "PROTOCOL.json"
    settings = load_json(protocol)
    c1, c2 = load_json(C1 / "RESULTS.json"), load_json(C2 / "RESULTS.json")
    q = settings["target_charge"]
    target = c2["charges"][str(q)]
    assert target["linear_stability"]
    knots, parent, field = parent_profiles(q)
    rows = match_observables(settings, target["profile_norm_radius"], archived_reference(), field)
    archive_delta = next(item["reference_delta0"] for item in rows if item["k"] == 5)
    ref = settings["reference"]
    reference_matches = abs(archive_delta-ref["archived_phase_at_k5"]) < ref["maximum_absolute_discrepancy"]
    if not reference_matches or not all(row["swave_agreement"] for row in rows):
        raise ArithmeticError("Expected reference or two-method S-wave agreement failed")
    result = {
        "protocol_sha256": hashlib.sha256(protocol.read_bytes()).hexdigest(), "scipy": scipy.__version__,
        "target_charge": q, "archive_delta0_absolute_discrepancy": abs(archive_delta-ref["archived_phase_at_k5"]),
        "scattering_reference_pass": reference_matches, "rows": rows,
        "physical_scales": physical_scales(settings, rows, c1, knots, parent, field),
    }
    with destination.open("x") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, required=True)
    main(cli.parse_args().output)
