"""Run C4 stress checks and enclose the static scalar tidal source."""
import argparse
import hashlib
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import cumulative_simpson, simpson
from scipy.interpolate import BPoly, PPoly


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
C1 = ROOT / "grut_c1_quantum_channels_20261127"
C2 = ROOT / "grut_c2_continuum_stability_20261127"
C3 = ROOT / "grut_c3_empirical_feasibility_20261127"
sys.path.insert(0, str(C2))
from certify import pi_interval
from intervals import I, ZERO, add, hermite, hull, integral, multiply, scale


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def as_interval(rational):
    value = Fraction(rational)
    return I(value.numerator) / value.denominator


def power(value, exponent):
    result = I(1)
    for _ in range(exponent):
        result = result * value
    return result


def symmetric(upper):
    return I(-upper.hi, upper.hi)


def reconstruct(cfg):
    q = cfg["target_charge"]
    with np.load(C1 / "PROFILES.npz") as archive:
        r = archive[f"fine_{q}_knots"]
        omega = float(archive[f"fine_{q}_omega"])
        f, df, s, ds = PPoly(archive[f"fine_{q}_coefficients"], r)(r).T.copy()
    df[0] = ds[0] = 0
    second = np.array([(s*s+0.1*f*f-omega*omega)*f, (s*s-1+f*f)*s])
    second[:, 0] /= 3
    second[:, 1:] -= 2*np.array([df, ds])[:, 1:]/r[1:]
    r = np.r_[r, cfg["profile"]["c2_taper_rho"]]
    data = (np.c_[np.r_[f, 0.], np.r_[df, 0.], np.r_[second[0], 0.]],
            np.c_[np.r_[s, 1.], np.r_[ds, 0.], np.r_[second[1], 0.]])
    return r, data, omega


def stress(fields, r, omega):
    f, s = (field(r) for field in fields)
    df, ds = (field(r, nu=1) for field in fields)
    kt, kr = 0.5*omega**2*f*f, 0.5*(df*df+ds*ds)
    w = 0.25*(s*s-1)**2+0.5*s*s*f*f+0.025*f**4
    return f, s, kt+kr+w, kt+kr-w, kt-kr-w, 4*kt-2*w


def numerical_stress(fields, r, omega, alpha):
    f, s, e, pr, pt, active = stress(fields, r, omega)
    r2 = r*r
    cumulative = cumulative_simpson(r2*e, x=r, initial=0.)
    charge = cumulative_simpson(r2*f*f, x=r, initial=0.)
    areal_center = -alpha*simpson(r*(e+pr), x=r)
    active_center = -alpha*simpson(r*active, x=r)
    df, ds = (field(r[1:], nu=1) for field in fields)
    ddf, dds = (field(r[1:], nu=2) for field in fields)
    wf = s[1:]**2*f[1:]+0.1*f[1:]**3
    ws = (s[1:]**2-1+f[1:]**2)*s[1:]
    pressure_prime = omega**2*f[1:]*df+df*ddf+ds*dds-wf*df-ws*ds
    return {
        "dimensionless_mass_integral": float(cumulative[-1]),
        "energy_R90_rho": float(np.interp(.9*cumulative[-1], cumulative, r)),
        "charge_R90_rho": float(np.interp(.9*charge[-1], charge, r)),
        "max_twice_local_compactness": float(max(2*alpha*cumulative[1:]/r[1:])),
        "center_lapse_areal": float(areal_center),
        "center_lapse_active_poisson": float(active_center),
        "center_lapse_sensitivity_to_accumulation": float(abs(areal_center-active_center)),
        "integrated_stress_trace_relative": float(simpson(r2*(pr+2*pt), x=r)/cumulative[-1]),
        "max_flat_stress_conservation_residual": float(
            max(abs(pressure_prime+2*(pr[1:]-pt[1:])/r[1:]))),
    }


def stationary_energy(nodes, h, polynomials, omega, eta):
    rho = [I(nodes[:-1]), h]
    f, s = polynomials
    u = multiply(rho, f)
    w = multiply(rho, add(s, [I(-1)]))
    f2, s2 = multiply(f, f), multiply(s, s)
    sdev = add(s2, [I(-1)])
    potential = add(add(scale(multiply(sdev, sdev), I(1)/4),
                        scale(multiply(s2, f2), I(1)/2)), scale(multiply(f2, f2), I(1)/40))
    pot_total = integral(multiply(multiply(rho, rho), potential), h)
    u2, w2 = integral(multiply(u, u), h), integral(multiply(w, w), h)
    du, dw = [[i*y[i]/h for i in range(1, len(y))] for y in (u, w)]
    gradient = integral(add(multiply(du, du), multiply(dw, dw)), h).sqrt()
    return (omega.square()*(u2.sqrt()+symmetric(eta)).square()/2
            + (gradient+symmetric(4*eta)).square()/2
            + pot_total+symmetric(10*((u2+w2).sqrt()+eta)*eta))


def tidal_force_norm(nodes, h, polynomials, omega, eta, tail_a, kappa):
    denominator = I(np.where(nodes[:-1] == 0, 1, nodes[:-1]))
    graph = (12*eta).hi
    field_error = np.minimum(graph, (3*eta/denominator).hi)
    field_error[0] = graph
    shift = I(-field_error, field_error)
    f, s = [hull(poly)+shift for poly in polynomials]
    f2, s2 = f.square(), s.square()
    j_f = 4*omega.square()*f-2*s2*f-(I(1)/5)*f2*f
    j_s = -2*s*(s2-1+f2)
    inside = (h*power(I(nodes[1:]), 6)*(j_f.square()+j_s.square())).sum_cells()
    radius = I(nodes[-1])
    tail_moment = sum((math.comb(4, j)*math.factorial(j)*power(radius, 4-j)
                       /power(2*kappa, j+1) for j in range(5)), ZERO)
    return I((inside+484*tail_a.square()*tail_moment).hi)


def bounds_from_c2(nodes, hermite_data, cert, alpha, cfg, physical):
    h = I(nodes[1:])-I(nodes[:-1])
    fields = [hermite(*data.T, h) for data in hermite_data]
    omega = I(*cert["omega_enclosure"])
    eta = I(cert["profile_norm_radius"])
    mass = stationary_energy(nodes, h, fields, omega, eta)
    mass_lower = I(mass.lo)
    R = I(nodes[-1])
    kappa = as_interval(cfg["profile"]["tail_decay_minimum"])
    a = I(cert["tail_at_R_upper"])
    b = I(cert["tail_derivative_at_R_upper"])
    B = (b+a/R).square()/2+3*a.square()
    exterior_mass = B/(2*kappa)
    pressure_lapse = alpha*exterior_mass*(I(1)/R+I(1)/(2*kappa*R.square()))
    active_lapse = 4*alpha*B/((2*kappa).square()*R.square())
    relative_gradient = 4*exterior_mass/mass_lower

    tidal = cfg["single_response"]
    exterior_floor = as_interval(tidal["exterior_floor"])
    gap = 2*exterior_floor/(tidal["potential_norm_cap"]*R.square())
    force_sq = tidal_force_norm(nodes, h, fields, omega, eta, a, kappa)
    susceptibility = force_sq/I(gap.lo)
    mass_ref = 4*pi_interval()*mass_lower
    size_ref = I(tidal["profile_radius_reference_rho"])
    l0, c = I(physical["rho_unit_m"]), I(physical["c"])
    t = I(tidal["example_tide_s_minus_2"])*(l0/c).square()
    rho_readout = I(tidal["example_readout_m"])/l0
    quadrupole_readout = alpha * susceptibility*t/(5*power(rho_readout, 3))
    return {
        "mass_integral_interval": mass.display(),
        "tail_energy_integral_upper_at_R": exterior_mass.display()[1],
        "tail_mass_fraction_upper_at_R": (exterior_mass/mass_lower).display()[1],
        "tail_twice_compactness_difference_upper_at_R": (2*alpha*exterior_mass/R).display()[1],
        "tail_lapse_pressure_separated_upper_at_R": pressure_lapse.display()[1],
        "tail_lapse_improved_upper_at_R": active_lapse.display()[1],
        "tail_lapse_fraction_upper_at_R": (active_lapse/(alpha*mass_lower/R)).display()[1],
        "tail_lapse_gradient_fraction_upper_at_R": relative_gradient.display()[1],
        "quadrupole_stiffness_gap_lower": gap.display()[0],
        "normalized_tidal_force_L2_squared_upper": force_sq.display()[1],
        "material_susceptibility_upper_per_dimensionless_tide": susceptibility.display()[1],
        "illustrative_dimensionless_tide": t.display(),
        "illustrative_material_moment_over_target_E_Rref_squared_upper": (
            susceptibility*t/(mass_ref*size_ref.square())).display()[1],
        "illustrative_material_lapse_upper_at_readout": quadrupole_readout.display()[1],
        "illustrative_upper_material_lapse_ratio_to_monopole": (
            susceptibility*t/(5*mass_lower*rho_readout.square())).display()[1],
    }


def main(output):
    if np.finfo(np.longdouble).nmant < 63:
        raise RuntimeError("Extended precision interval arithmetic is required")
    protocol = HERE / "PROTOCOL.json"
    cfg = read(protocol)
    q = cfg["target_charge"]
    c1, c2, c3 = (read(d / "RESULTS.json") for d in (C1, C2, C3))
    cert = c2["charges"][str(q)]
    assert cert["stationary_profile_certified"] and cert["all_angular_sectors"]
    assert digest(C1 / "PROFILES.npz") == c2["profiles_sha256"]
    assert read(C2 / "PROTOCOL.json")["domain"]["taper_radius"] == cfg["profile"]["c2_taper_rho"]
    c1_params = read(C1 / "PROTOCOL.json")
    g = float(Fraction(cfg["couplings"]["g_X"]))
    vm = float(Fraction(cfg["couplings"]["v_over_M_P"]))
    assert g == c1_params["g_X"] and vm == c1_params["v_over_reduced_Planck_mass"]
    for key in ("lambda_phi", "lambda_X"):
        assert float(Fraction(cfg["couplings"][key])) == c1_params[key]
    alpha = as_interval(cfg["couplings"]["v_over_M_P"]).square()/2
    nodes, data, omega = reconstruct(cfg)
    assert len(nodes)-1 == cert["cells"] <= cfg["budget"]["maximum_flat_profile_cells"]
    assert cert["omega_predictor"][0] < omega < cert["omega_predictor"][1]
    fields = [BPoly.from_derivatives(nodes, values) for values in data]
    n = cfg["gravitation_checks"]
    rows = [numerical_stress(fields, np.linspace(0, nodes[-1], intervals+1), omega, vm**2/2)
            for intervals in n["radial_subdivisions"]]
    fine = rows[-1]
    stress_points = np.array(n["stress_samples_rho"])
    f, s, e, pr, pt, active = stress(fields, stress_points, omega)
    stress_samples = [dict(zip(("rho", "F", "S", "e", "p_r", "p_t", "J"), row))
                      for row in np.stack((stress_points, f, s, e, pr, pt, active), axis=1).tolist()]
    tail_samples = {}
    for rho in n["tail_samples_rho"]:
        rr = float(rho)
        mesh = rr + np.linspace(0., nodes[-1]-rr, n["radial_subdivisions"][-1]+1)
        _, _, energy, _, _, _ = stress(fields, mesh, omega)
        tail_samples[str(rho)] = float(simpson(mesh**2*energy, x=mesh)/fine["dimensionless_mass_integral"])

    e0 = c1["exact_charge_spectra"][str(q)]["E0"]
    e1 = c1["exact_charge_spectra"][str(q)]["E0plus1_finite"]
    energy = 4*np.pi*fine["dimensionless_mass_integral"]/g**2
    norm = c3["physical_scales"]
    si = read(C3 / "PROTOCOL.json")["constants_si"]
    scales = {"rho_unit_m": norm["rho_unit_m"], "c": si["c"]}
    energy_unit_j = norm["g_X_v_gev"]*si["GeV_joules"]
    mass = energy*energy_unit_j/si["c"]**2
    mass_one_loop = e1*energy_unit_j/si["c"]**2
    radius = fine["energy_R90_rho"]*scales["rho_unit_m"]
    free_e = q*c1["direct_emission"]["pole"]
    result = {
        "protocol_sha256": digest(protocol),
        "inputs": {"c1_profiles": digest(C1 / "PROFILES.npz"), "c2": digest(C2 / "RESULTS.json"),
                   "c1": digest(C1 / "RESULTS.json"), "c3": digest(C3 / "RESULTS.json")},
        "numpy": np.__version__, "scipy": scipy.__version__, "target_charge": q,
        "classical_stress_diagnostics": rows, "stress_samples": stress_samples,
        "estimated_fractions_of_energy_outside": tail_samples,
        "analytical_bounds": bounds_from_c2(nodes, data, cert, alpha, cfg, scales),
        "physical_scales": {
            "rho_unit_m": scales["rho_unit_m"], "g_X_v_gev": norm["g_X_v_gev"],
            "tree_rest_energy_in_g_X_v": energy, "classical_mass_kg": mass,
            "scalar_one_loop_rest_mass_estimate_kg": mass_one_loop,
            "scalar_one_loop_rest_energy_difference_to_free_g_X_v": e1-free_e,
            "scalar_one_loop_mass_defect_kg": (free_e-e1)*energy_unit_j/si["c"]**2,
            "fractional_scalar_one_loop_mass_defect": (free_e-e1)/free_e,
            "energy_R90_m": radius,
            "charge_R90_m": fine["charge_R90_rho"]*scales["rho_unit_m"],
            "total_C_GM_over_energy_R90_c_squared": si["G"]*mass/(radius*si["c"]**2),
            "schwarzschild_radius_from_classical_mass_m": 2*si["G"]*mass/si["c"]**2,
            "self_gravity_dimensional_scale_j": si["G"]*mass**2/radius,
            "monopole_acceleration_at_example_readout_m_s2": (
                si["G"]*mass/cfg["single_response"]["example_readout_m"]**2),
            "interior_resolution_momentum_c_gev_at_energy_R90": norm["g_X_v_gev"]/fine["energy_R90_rho"],
        },
        "checks": {
            "energy_agreement_C1": abs(energy-e0)<n["energy_reference_absolute_tolerance_in_g_X_v"],
            "virial_trace": abs(fine["integrated_stress_trace_relative"])<n["integrated_stress_trace_relative_tolerance"],
            "stress_conservation": fine["max_flat_stress_conservation_residual"]<n["stress_conservation_absolute_tolerance"],
            "lapse_routes_agree": fine["center_lapse_sensitivity_to_accumulation"]<n["two_lapse_routes_absolute_tolerance"],
        },
    }
    with output.open("x") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print("C4 numerical gates", result["checks"], flush=True)
    if not all(result["checks"].values()):
        raise ArithmeticError("Frozen stress verification gate failed")


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", required=True, type=Path)
    main(cli.parse_args().output)
