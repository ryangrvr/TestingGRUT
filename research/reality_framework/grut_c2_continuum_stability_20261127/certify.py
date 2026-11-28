"""Run the frozen C2 half-line certificate on the archived C1 profiles."""
import argparse
import hashlib
import json
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_bvp
from scipy.interpolate import PPoly

import spectrum
from intervals import I, ZERO, add, hermite, hessian, hull, integral, multiply, residual_norm, scale


HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "PROTOCOL.json"
SOURCE = HERE.parent / "grut_c1_quantum_channels_20261127/PROFILES.npz"


def pi_interval():
    """Machin identity with alternating rational series enclosures."""
    def arctangent_bounds(denominator):
        terms = [Fraction((-1)**n, (2*n+1)*denominator**(2*n+1)) for n in range(33)]
        return sum(terms[:32]), sum(terms)

    a, b = arctangent_bounds(5), arctangent_bounds(239)
    lo, hi = 16*a[0]-4*b[1], 16*a[1]-4*b[0]
    # These rational endpoints are exact in the comparison with Machin's bounds.
    assert Fraction("3.14159265358979323846264338327950") < lo
    assert hi < Fraction("3.14159265358979323846264338327951")
    return I(np.nextafter(np.longdouble("3.14159265358979323846264338327950"), -np.inf),
             np.nextafter(np.longdouble("3.14159265358979323846264338327951"), np.inf))


def predictor(archive, q, settings):
    nodes = archive[f"fine_{q}_knots"]
    omega = float(archive[f"fine_{q}_omega"])
    matrix = archive[f"fine_{q}_coefficients"]
    spline = PPoly(matrix, nodes)

    def tangent_ode(x, y):
        f, _, s, _ = spline(x).T
        p, dp, t, dt = y
        return np.array([dp, (s*s+0.3*f*f-omega*omega)*p+2*f*s*t-2*omega*f,
                         dt, (3*s*s-1+f*f)*t+2*f*s*p])

    def bc(a, b):
        return np.array([a[1], a[3], b[0], b[2]])

    tangent = solve_bvp(tangent_ode, bc, nodes, np.zeros((4, len(nodes))),
                        S=np.diag([0., -2., 0., -2.]),
                        tol=settings["profile"]["tangent_tolerance"],
                        max_nodes=settings["profile"]["maximum_tangent_nodes"])
    if not tangent.success:
        raise ArithmeticError("Tangent predictor failed: " + tangent.message)
    yr = spline(nodes).T.copy()
    vr = tangent.sol(nodes)
    yr[[1, 3], 0] = 0
    vr[[1, 3], 0] = 0
    f, df, s, ds = yr
    nonlinear = np.array([(s*s+0.1*f*f-omega*omega)*f, (s*s-1+f*f)*s])
    tangent_source = tangent_ode(nodes, vr)[[1, 3]]

    def second(rhs, first):
        result = rhs.copy()
        result[:, 0] /= 3
        result[:, 1:] -= 2*first[:, 1:]/nodes[1:]
        return result

    ydd = second(nonlinear, yr[[1, 3]])
    vdd = second(tangent_source, vr[[1, 3]])
    end = settings["domain"]["taper_radius"]
    nodes = np.r_[nodes, end]
    if len(nodes)-1 > settings["budget"]["maximum_cells_per_profile"]:
        raise ArithmeticError("Polynomial cell budget exceeded")
    h = I(nodes[1:])-I(nodes[:-1])
    fpoly = hermite(np.r_[f, 0.], np.r_[df, 0.], np.r_[ydd[0], 0.], h)
    spoly = hermite(np.r_[s, 1.], np.r_[ds, 0.], np.r_[ydd[1], 0.], h)
    ppoly = hermite(np.r_[vr[0], 0.], np.r_[vr[1], 0.], np.r_[vdd[0], 0.], h)
    tpoly = hermite(np.r_[vr[2], 0.], np.r_[vr[3], 0.], np.r_[vdd[1], 0.], h)
    trial = [nodes*np.r_[df, 0.], nodes*np.r_[ds, 0.]]
    return nodes, h, I(omega), (fpoly, spoly, ppoly, tpoly), trial


def defect_data(nodes, h, omega, polys):
    f, s, p, t = polys
    rho = [I(nodes[:-1]), h]
    u, w = multiply(rho, f), multiply(rho, add(s, [I(-1)]))
    v1, v2 = multiply(rho, p), multiply(rho, t)
    f2, s2 = multiply(f, f), multiply(s, s)
    phase = add(add(s2, scale(f2, I(1)/10)), [-omega.square()])
    va = add(phase, scale(f2, I(1)/5))
    vb = scale(multiply(f, s), 2)
    vc = add(add(scale(s2, 3), f2), [I(-1)])
    stat = (add(scale(hessian(u, h), -1), multiply(phase, u)),
            add(scale(hessian(w, h), -1), multiply(multiply(rho, s), add(add(s2, f2), [I(-1)]))))
    tan = (add(add(scale(hessian(v1, h), -1), multiply(va, v1)),
               add(multiply(vb, v2), scale(u, -2*omega))),
           add(add(scale(hessian(v2, h), -1), multiply(vb, v1)), multiply(vc, v2)))
    norm_u_squared = integral(multiply(u, u), h)
    norm_v_squared = integral(add(multiply(v1, v1), multiply(v2, v2)), h)
    pairing = integral(multiply(u, v1), h)
    return (residual_norm(stat, h), residual_norm(tan, h), norm_u_squared,
            norm_v_squared, pairing, [hull(a) for a in (va, vb, vc)], [hull(phase)])


def upper(values):
    return I(np.maximum.reduce([I.of(x).hi for x in values]))


def both_signs(bound):
    return I(-bound.hi, bound.hi)


def matrix_cap(va, vb, vc):
    return upper([2, np.max((va.magnitude()+vb.magnitude()).hi),
                  np.max((vc.magnitude()+vb.magnitude()).hi)])


def prove_charge(archive, q, settings, result, pi):
    nodes, h, omega, polys, trial = predictor(archive, q, settings)
    result["cells"] = len(nodes)-1
    result["omega_predictor"] = omega.display()
    stat, tan, ii, vv, pair, amplitude, phase = defect_data(nodes, h, omega, polys)
    target = I(q)/(400*pi)
    charge_defect = omega*ii-target
    result["stationary_residual_l2"] = stat.display()
    result["tangent_residual_l2"] = tan.display()
    result["charge_defect"] = charge_defect.display()
    ext = 1-omega.square()
    radial_trials = []
    result["radial_trials"] = radial_trials
    gap = None
    for ratio in settings["spectrum"]["radial_gap_trials"]:
        test_gap = float(Fraction(ratio))
        upper_trial = spectrum.trial_form(trial, h, amplitude, -test_gap)
        attempt = {"gap": test_gap, "negative_trial_shifted_form": upper_trial.display()}
        radial_trials.append(attempt)
        if upper_trial.hi >= 0:
            continue
        try:
            attempt["positive_shift_upper_count"] = spectrum.count(nodes, h, amplitude, 0, test_gap, ext)
        except ArithmeticError as err:
            attempt["enclosure_failure"] = str(err)
            continue
        if attempt["positive_shift_upper_count"] == 1:
            gap = I(test_gap)
            break
    if gap is None:
        raise ArithmeticError("Unable to bound the radial amplitude index and inverse gap")

    norm_u, norm_v = upper([ii.sqrt()]), upper([vv.sqrt()])
    v_error = tan/gap
    norm_v_exact = upper([norm_v+v_error])
    s0 = ii+2*omega*pair+both_signs(2*omega.magnitude()*norm_u*v_error)
    result["predictor_slope_Qbar"] = s0.display()
    if s0.hi >= 0:
        raise ArithmeticError("Predictor charge-slope denominator is not bounded below zero")
    d = (1+norm_v_exact)/I(-s0.hi)
    k = I(1)/gap+norm_v_exact*d
    vcap = matrix_cap(*amplitude)
    gamma = upper([d, k, (vcap*k+2*omega.magnitude()*norm_u*d+1)/16])
    eps = upper([stat, charge_defect.magnitude()])
    radius = upper([2*gamma*eps])
    contraction = upper([360*radius*gamma])
    result["border_inverse_norm_upper"] = gamma.display()[1]
    result["profile_norm_radius"] = radius.display()[1]
    result["contraction_upper"] = contraction.display()[1]
    if contraction.hi >= 0.5:
        raise ArithmeticError("Bordered Newton contraction exceeds one half")
    f, s, _, _ = polys
    fcap = I(np.max(hull(f).magnitude().hi))
    scap = I(max(1, np.max(hull(s).magnitude().hi)))
    if any(v.hi >= lim for v, lim in [(fcap+12*radius, 3), (scap+12*radius, 2), (omega+radius, 1),
                                      (norm_u+radius, 8), (radius, 1)]):
        raise ArithmeticError("Newton ball exceeds the frozen polynomial derivative caps")

    exact_omega = omega+both_signs(radius)
    potential_error = 336*radius
    exact_gap = gap-potential_error
    result["potential_error_upper"] = potential_error.display()[1]
    result["omega_enclosure"] = exact_omega.display()
    result["radial_gap_lower"] = exact_gap.display()[0]
    if exact_gap.lo <= 0 or exact_omega.lo <= 0:
        raise ArithmeticError("Exact-saddle error reaches zero frequency or the radial spectral gap")

    # The predictor vanishes in y variables outside R.
    tail_fields = both_signs(12*radius)
    s_out = 1+tail_fields
    a = s_out.square()+(I(1)/10)*tail_fields.square()-exact_omega.square()
    b = 2+3*tail_fields+tail_fields.square()
    c = (tail_fields*s_out).magnitude()
    ext_stationary = min((a-c/2).lo, (b-c/2).lo)
    exterior_amplitude = min((s_out.square()+(I(3)/10)*tail_fields.square()-exact_omega.square()-2*c).lo,
                             (3*s_out.square()-1+tail_fields.square()-2*c).lo)
    if ext_stationary <= (I(1)/9).hi or (a.magnitude()+b.magnitude()+c).hi >= 4:
        raise ArithmeticError("Exact exterior stationary decay inequalities are unresolved")
    if min(exterior_amplitude, a.lo) <= 0.125:
        raise ArithmeticError("Exact exterior spectral floor is below one eighth")
    result["tail_at_R_upper"] = (3*radius).display()[1]
    result["tail_derivative_at_R_upper"] = (36*radius).display()[1]
    result["exterior_spectral_floor"] = I(min(exterior_amplitude, a.lo)).display()[0]

    true_tan_defect = tan+potential_error*norm_v+18*radius
    slope = ii+2*omega*pair+both_signs(17*radius+18*radius*norm_v+16*true_tan_defect/exact_gap)
    result["exact_slope_Qbar"] = slope.display()
    result["exact_slope_Q"] = (400*pi*slope).display()
    if slope.hi >= 0:
        raise ArithmeticError("Exact continuum charge-frequency derivative is unresolved")
    result["fixed_charge_radial_positive"] = True
    result["stationary_profile_certified"] = True

    symmetry = float(Fraction(settings["spectrum"]["symmetry_gap"]))
    gap_symmetry = I(symmetry)-potential_error
    if gap_symmetry.lo <= 0:
        raise ArithmeticError("Certified profile error reaches symmetry-complement gap")
    for key, potential, ell in [("phase_ell0", phase, 0), ("amplitude_ell1", amplitude, 1)]:
        n = spectrum.count(nodes, h, potential, ell, symmetry, ext)
        result[key] = {"upper_count_below_shift": n, "complement_gap_lower": gap_symmetry.display()[0]}
        if n != 1:
            raise ArithmeticError(f"Exact zero mode plus one-mode upper count unresolved in {key}")

    result["all_angular_sectors"] = True
    result["reporter_floor"] = 1
    result["linear_stability"] = True


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(destination):
    if np.finfo(np.longdouble).nmant < 63:
        raise RuntimeError("Certificate requires extended precision long double")
    t0 = time.monotonic()
    cfg = json.loads(PROTOCOL.read_text())
    output = {"protocol_sha256": sha(PROTOCOL), "profiles_sha256": sha(SOURCE),
              "numpy": np.__version__, "scipy": scipy.__version__, "charges": {}}
    pi = pi_interval()
    with destination.open("x") as stream:
        json.dump(output, stream, indent=2)
    with np.load(SOURCE) as archive:
        for charge in cfg["charges"]:
            data = {}
            output["charges"][str(charge)] = data
            try:
                if time.monotonic()-t0 > cfg["budget"]["wall_seconds"]:
                    raise ArithmeticError("Frozen wall budget exhausted")
                prove_charge(archive, charge, cfg, data, pi)
            except ArithmeticError as error:
                data["unresolved"] = str(error)
            print(charge, json.dumps(data), flush=True)
            destination.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
    passed = sum(row.get("linear_stability", False) for row in output["charges"].values())
    partial = any(row.get("stationary_profile_certified", False) for row in output["charges"].values())
    output["verdict"] = "A" if passed == len(cfg["charges"]) else "B" if partial else "D"
    output["complete_targets"] = passed
    destination.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
    print("Verdict", output["verdict"], flush=True)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, required=True)
    main(cli.parse_args().output)
