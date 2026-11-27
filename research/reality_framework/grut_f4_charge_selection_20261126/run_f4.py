"""Evaluate the frozen F4 fixtures with exact arithmetic."""
import argparse
import json
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent


def fraction_sum(terms):
    return sum(terms, Fraction(0))


def check_fixture(k, inputs):
    fields = inputs["fields"]
    fermions = [field for field in fields if field["label"] != "H"]
    weights = {field["label"]: field["y0"] + 6 * k * field["D"] for field in fields}
    q, u, d, l, e, n, h = [weights[field["label"]] for field in fields]
    anomalies = {
        "color_squared_y": 2 * q + u + d,
        "weak_squared_y": 3 * q + l,
        "gravity_gravity_y": sum(
            f["color"] * f["weak"] * weights[f["label"]] for f in fermions
        ),
        "y_cubed": sum(
            f["color"] * f["weak"] * weights[f["label"]] ** 3 for f in fermions
        ),
        "color_cubed": 2 - 1 - 1,
    }
    mixed_D = {
        "color_squared_D": 2 * fields[0]["D"] + fields[1]["D"] + fields[2]["D"],
        "weak_squared_D": 3 * fields[0]["D"] + fields[3]["D"],
        "gravity_gravity_D": sum(f["color"] * f["weak"] * f["D"] for f in fermions),
        "y_squared_D": sum(
            f["color"] * f["weak"] * weights[f["label"]] ** 2 * f["D"]
            for f in fermions
        ),
    }
    center_residues = {
        f["label"]: (weights[f["label"]] + 2 * f["triality"] + 3 * f["weak_center"])
        % inputs["quotient_order"]
        for f in fields
    }
    # Two Higgs fields plus one derivative in C. The Lorentz scalar C^2
    # conjugate has weight -4h, and the singlet bilinear has weight 2n.
    qqql_weight = 3 * q + l
    dressed_weight = 2 * n - 3 * k * 4 * h
    dimension_dressing = 3 + 3 * abs(k) * 2 * 3
    slope_key = "k_negative" if k < 0 else "k_nonnegative"
    slope = Fraction(inputs["triangle_slopes"][slope_key])
    dimension_lower = 2 * abs(6 * k + slope)

    up, down = Fraction(1, 2) + Fraction(q, 6), Fraction(-1, 2) + Fraction(q, 6)
    nu, lepton = Fraction(1, 2) + Fraction(l, 6), Fraction(-1, 2) + Fraction(l, 6)
    proton, neutron = 2 * up + down, up + 2 * down
    neutron_ratio = neutron / abs(lepton)
    weinberg_allowed = 2 * (l + h) == 0
    b_y = Fraction(2, 3) * inputs["families"] * fraction_sum(
        f["color"] * f["weak"] * Fraction(weights[f["label"]], 6) ** 2
        for f in fermions
    ) + Fraction(2, 3) * Fraction(h, 6) ** 2
    ratio_threshold = Fraction(inputs["minimum_nonzero_abs_neutron_to_electron_charge_ratio"])
    checks = {
        "fixed_higgs_and_yukawas": h == inputs["higgs_weight"]
        and [q + u + h, q + d - h, l + e - h, l + n + h] == [0] * 4,
        "compact_group_residues": all(value == 0 for value in center_residues.values()),
        "perturbative_gauge_anomalies": all(value == 0 for value in anomalies.values()),
        "weak_doublet_parity": inputs["families"] * (3 + 1) % 2 == 0,
        "renormalizable_continuous_B_minus_L_anomaly_free": all(
            value == 0 for value in mixed_D.values()
        ),
        "registered_breaking_operator_weights_zero": (qqql_weight, dressed_weight) == (0, 0),
        "classical_breaking_directions_independent": 1 * (-2) - 1 * 0 != 0,
        "dimension_bounds_ordered": dimension_lower <= dimension_dressing,
        "unbroken_fermion_charges_vectorlike": [up, down, lepton, nu]
        == [Fraction(-weight, 6) for weight in (u, d, e, n)],
        "neutral_hydrogen_control": proton + lepton == 0,
        "weinberg_and_majorana_only_at_zero": weinberg_allowed == (k == 0)
        and (2 * n == 0) == (k == 0),
        "physical_neutron_charge_gap": neutron_ratio == 0
        if k == 0
        else abs(neutron_ratio) >= ratio_threshold,
        "one_loop_running_polynomial": b_y == Fraction(41 + 192 * k + 576 * k**2, 6),
    }
    return {
        "k": k,
        "weights": weights,
        "charges": {
            label: str(charge)
            for label, charge in {
                "u": up, "d": down, "e": lepton, "nu": nu, "p": proton, "n": neutron,
                "hydrogen_atom": proton + lepton, "neutron_over_abs_e": neutron_ratio,
            }.items()
        },
        "gauge_anomalies": anomalies,
        "mixed_D_anomalies": mixed_D,
        "center_residues": center_residues,
        "breaking_operator_weights": {"QQQL": qqql_weight, "dressed": dressed_weight},
        "lower_dimension": int(dimension_lower),
        "dressed_dimension": dimension_dressing,
        "weinberg_allowed": weinberg_allowed,
        "bY": str(b_y),
        "checks": checks,
    }


def evaluate(inputs):
    fields = inputs["fields"]
    coefficients = {}
    for slope_label, encoded_slope in inputs["triangle_slopes"].items():
        slope = Fraction(encoded_slope)
        coefficients[slope_label] = {
            field["label"]: str(
                abs(field["y0"] - slope * field["D"]) / Fraction(field["dimension"])
            )
            for field in fields
        }
    exp = inputs["neutron_experiment"]
    center = Fraction(exp["charge_ratio_center"])
    uncertainty = Fraction(exp["charge_ratio_reported_uncertainty"])
    scale = abs(center) + int(exp["comparison_multiplier"]) * uncertainty
    threshold = Fraction(inputs["minimum_nonzero_abs_neutron_to_electron_charge_ratio"])
    fixtures = [check_fixture(k, inputs) for k in inputs["fixtures_k"]]
    common_checks = {
        "triangle_constants": all(
            Fraction(value) <= 3 for row in coefficients.values() for value in row.values()
        ),
        "neutron_benchmark_selected_only_zero": scale < threshold,
        "registered_ratio_bound_attained": any(
            abs(Fraction(row["charges"]["neutron_over_abs_e"])) == threshold
            for row in fixtures
        ),
        "zero_member_dimension_exact": any(
            row["k"] == 0 and row["lower_dimension"] == row["dressed_dimension"] == 3
            for row in fixtures
        ),
    }
    check_total = len(common_checks) + sum(len(row["checks"]) for row in fixtures)
    checks_passed = sum(common_checks.values()) + sum(
        sum(row["checks"].values()) for row in fixtures
    )
    return {
        "checks": common_checks,
        "checks_passed": checks_passed,
        "checks_total": check_total,
        "fixtures": fixtures,
        "neutron_benchmark": {
            "absolute_comparison_scale": str(scale),
            "data_center": str(center),
            "data_reported_uncertainty": str(uncertainty),
            "effect_gap": str(threshold),
            "gap_divided_by_comparison_scale": str(threshold / scale),
            "doi": exp["doi"],
        },
        "triangle_ratios": coefficients,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    inputs = json.loads((HERE / "PARAMETERS.json").read_text())
    results = evaluate(inputs)
    content = json.dumps(results, sort_keys=True, indent=2) + "\n"
    if args.output:
        args.output.write_text(content)
    else:
        print(content, end="")
    raise SystemExit(int(results["checks_total"] != results["checks_passed"]))


if __name__ == "__main__":
    main()
