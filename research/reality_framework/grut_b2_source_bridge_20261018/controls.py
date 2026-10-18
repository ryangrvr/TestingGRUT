#!/usr/bin/env python3
"""Frozen B2-SB rational identities and small diagnostic equality controls.

No source Hamiltonian trajectory, source eigenvalue computation or search.
The atomic witness is a countermodel to scalar inference, not the B2 law.
"""
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent


def encode(x):
    if isinstance(x, F):
        return {"numerator": x.numerator, "denominator": x.denominator,
                "decimal_diagnostic": float(x)}
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [encode(v) for v in x]
    return x


def main():
    checked = []

    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        checked.append(name)

    jp, jm, delta = F(1, 64), F(1, 96), F(5, 4)
    A, J = F(13, 46080), F(1, 3840)
    omega_energy = 2 - 2 * A
    b0, bp = jp + jm, jp - jm
    q0 = (2 * b0**2 / delta) / (1 + 2 * b0**2 / delta**2)
    d0_J = (q0 - 2 * A) / J
    dp_J = (2 * A - 2 * bp**2 / (delta - 2 * bp)) / J
    check("lower_band_margin_exact", d0_J == F(13811, 6918))
    check("upper_band_margin_exact", dp_J == F(1427, 714))
    check("witness_strictly_inside_reference_band", min(d0_J, dp_J) > F(7, 4))

    alpha, beta = F(768, 791), F(23, 791)
    grid = list(range(-112, 113))
    weights = {n: alpha / 225 + (beta if n == 0 else 0) for n in grid}
    mass = sum(weights.values(), F(0))
    mean = sum((weights[n] * F(n, 64) for n in grid), F(0))
    variance = sum((weights[n] * F(n, 64)**2 for n in grid), F(0))
    central_mass = sum((weights[n] for n in grid if abs(n) <= 96), F(0))
    check("positive_atomic_measure", min(weights.values()) > 0)
    check("mass_exactly_one", mass == 1)
    check("dimensionless_mean_exactly_zero", mean == 0)
    check("dimensionless_variance_exactly_one", variance == 1)
    check("central_mass_exact", central_mass == F(51133, 59325))
    check("central_mass_exceeds_five_ninths", central_mass > F(5, 9))
    check("high_measure_zero", omega_energy + J * F(7, 4) < F(5, 2))
    sine_lower = F(1, 64) - F(1, 6 * 64**3)
    early_upper = beta + alpha / (225 * sine_lower)
    check("positive_sine_lower_bound", sine_lower > 0)
    check("continuous_early_upper_below_031", early_upper < F(31, 100))
    check("continuous_early_upper_below_065", early_upper < F(13, 20))

    L, theta_pi_coefficient = 3072, 128
    check("declared_length_allowed", L >= 2016 and L % 24 == 0)
    check("recurrence_phase_integer_turns", all(F(theta_pi_coefficient * n, 64) == 2 * n for n in grid))
    check("pi_lower_puts_recurrence_above_window_start", theta_pi_coefficient * 3 == L // 8)
    check("pi_upper_puts_recurrence_below_window_end", theta_pi_coefficient * 4 == L // 6)

    central_lower = F(5, 9)
    eta = F(24, 25)
    check("central_capture_first_radical", central_lower * eta == F(8, 15))
    check("central_capture_second_radical", (1 - central_lower) * (1 - eta) == F(4, 225))
    check("positive_rational_radical_lower", F(73, 100)**2 < F(8, 15))
    amplitude_lower = F(73, 100) - F(2, 15)
    arrival_lower = amplitude_lower**2
    check("capture_rational_amplitude", amplitude_lower == F(179, 300))
    check("capture_rational_probability", arrival_lower == F(32041, 90000))
    check("capture_closes_original_source_threshold", arrival_lower > F(17, 50))
    joint_lower = (F(17, 50) - F(1, 2000) - F(1, 1000))**2
    check("correct_joint_lower", joint_lower == F(458329, 4000000))
    check("joint_lower_exceeds_one_tenth", joint_lower > F(1, 10))
    sigma_upper = F(6, 27613)
    check("twice_high_amplitude_below_00005", 2 * sigma_upper < F(1, 2000))
    check("high_weight_below_1e6", sigma_upper**2 < F(1, 10**6))
    retained_amplitude_lower = 1 - F(2, 10**6) - (F(2, 3) + F(1, 10**6))**2 - F(1, 1000)
    check("retention_closure_exceeds_one_quarter", retained_amplitude_lower**2 > F(1, 4))
    check("left_ray_margin", F(11, 10) / 8 - F(1, 8) == F(1, 80))
    check("right_ray_margin", F(3, 8) - F(11, 5) / 6 == F(1, 120))

    # Independent elementary 2D arithmetic check of the sharp lemma.
    # Exact proof is in DECOMPOSITION.md, not inferred from these examples.
    equality_controls = []
    for w, e in [(F(5, 9), F(24, 25)), (F(4, 5), F(7, 10)), (F(1, 2), F(4, 5))]:
        wf, ef = float(w), float(e)
        chi = (math.sqrt(ef), math.sqrt(1 - ef))
        zeta = (-math.sqrt(1 - ef), math.sqrt(ef))
        psi = tuple(math.sqrt(wf) * chi[n] + math.sqrt(1 - wf) * zeta[n] for n in (0, 1))
        exact_formula_diagnostic = (math.sqrt(wf * ef) - math.sqrt((1 - wf) * (1 - ef)))**2
        probability_diagnostic = psi[0]**2
        cross = -2 * math.sqrt(wf * ef * (1 - wf) * (1 - ef))
        check("capture_equality_" + str(w), abs(probability_diagnostic - exact_formula_diagnostic) < 1e-14)
        check("orthogonality_" + str(w), abs(sum(chi[n] * zeta[n] for n in (0, 1))) < 1e-14)
        check("normalization_" + str(w), abs(sum(x*x for x in psi) - 1) < 1e-14)
        equality_controls.append({"w": w, "eta": e, "probability_diagnostic": probability_diagnostic,
                                  "coherent_cross_term_diagnostic": cross,
                                  "equality_defect_diagnostic": abs(probability_diagnostic - exact_formula_diagnostic)})

    # Only the single predeclared recurrence point; proof uses integer phases.
    theta = theta_pi_coefficient * math.pi
    return_diagnostic = sum(float(weights[n]) * complex(math.cos(theta*n/64), -math.sin(theta*n/64)) for n in grid)
    check("declared_recurrence_diagnostic", abs(abs(return_diagnostic) - 1) < 1e-13)

    payload = {
        "status": "all_frozen_controls_pass",
        "scope": "exact scalar countermodel and conditional bridge checks; no physical source violation",
        "arithmetic": "Fraction expressions exact; trigonometric and 2D floats diagnostic only",
        "checks": checked,
        "witness": {"grid_denominator": 64, "grid_maximum": 112, "atoms": 225,
                    "alpha": alpha, "beta": beta, "mass": mass, "mean": mean,
                    "variance": variance, "central_mass": central_mass,
                    "early_continuous_upper": early_upper,
                    "L": L, "theta_recurrence": "128*pi",
                    "theta_window": [L//8, L//6], "exact_return_modulus": 1,
                    "return_modulus_diagnostic": abs(return_diagnostic)},
        "band_margins_in_J": [d0_J, dp_J],
        "bridge": {"conditional_central_capture": eta, "one_source_lower_if_capture": arrival_lower,
                   "unchanged_joint_sufficient_lower": joint_lower,
                   "joint_retention_lower_if_fuel_bound": retained_amplitude_lower**2},
        "capture_equality_controls": equality_controls,
    }
    output = ROOT / "CONTROLS.json"
    output.write_text(json.dumps(encode(payload), indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": payload["status"], "checks": len(checked),
                      "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest()}))


if __name__ == "__main__":
    main()
