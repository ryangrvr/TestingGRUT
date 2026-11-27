"""Independent exact checks of the preregistered GRUT-F4 fixtures."""

import json
from fractions import Fraction
from pathlib import Path


def weights_from_quark(quark, higgs):
    """Solve the four Yukawa equations and the weak mixed anomaly."""
    lepton = -3 * quark
    return {
        "Q": quark,
        "u^c": -quark - higgs,
        "d^c": -quark + higgs,
        "L": lepton,
        "e^c": -lepton + higgs,
        "n^c": -lepton - higgs,
        "H": higgs,
    }


def weighted_trace(values, dimensions):
    return sum(values[name] * dimensions[name] for name in values)


def anomaly_coefficients(fermions, y):
    """Drop common nonzero Dynkin indices and family multiplicities."""
    dimensions = {name: row["color"] * row["weak"] for name, row in fermions.items()}
    d = {name: row["D"] for name, row in fermions.items()}
    color_names = [name for name, row in fermions.items() if row["color"] == 3]
    weak_names = [name for name, row in fermions.items() if row["weak"] == 2]
    weak_multiplicity = sum(fermions[name]["color"] for name in weak_names)
    return {
        "color_cubic": sum(
            fermions[name]["weak"] * fermions[name]["triality"] for name in color_names
        ),
        "color_color_y": sum(fermions[name]["weak"] * y[name] for name in color_names),
        "weak_weak_y": sum(fermions[name]["color"] * y[name] for name in weak_names),
        "gravity_y": weighted_trace(y, dimensions),
        "y_cubic": weighted_trace({name: value**3 for name, value in y.items()}, dimensions),
        "color_color_D": sum(fermions[name]["weak"] * d[name] for name in color_names),
        "weak_weak_D": sum(fermions[name]["color"] * d[name] for name in weak_names),
        "gravity_D": weighted_trace(d, dimensions),
        "D_y_y": weighted_trace({name: d[name] * y[name]**2 for name in d}, dimensions),
        "D_D_y": weighted_trace({name: d[name]**2 * y[name] for name in d}, dimensions),
        "D_cubic": weighted_trace({name: value**3 for name, value in d.items()}, dimensions),
        "lifted_weak_doublet_parity": weak_multiplicity % 2,
    }


def fixture(k, inputs):
    fields = {row["label"]: row for row in inputs["fields"]}
    fermions = {name: row for name, row in fields.items() if name != "H"}
    order = inputs["quotient_order"]
    charges = weights_from_quark(1 + order * k, inputs["higgs_weight"])
    assert weights_from_quark(1, inputs["higgs_weight"]) == {
        name: row["y0"] for name, row in fields.items()
    }
    assert all(
        charges[name] - row["y0"] == order * k * row["D"] for name, row in fields.items()
    )

    # The Higgs restricts a kernel element to theta=pi*m/3, weak center (-1)^m.
    # Its action on Q must then be a multiple of the color center.
    kernel_centers = []
    for m in range(order):
        for triality_center in range(3):
            compatible = all(
                (2 * triality_center * row["triality"] + m * (3 * row["weak_center"] + charges[name]))
                % order == 0
                for name, row in fields.items()
            )
            if compatible:
                kernel_centers.append([m, triality_center])
    assert kernel_centers == [[m, m % 3] for m in range(order)]

    anomaly = anomaly_coefficients(fermions, {name: charges[name] for name in fermions})
    assert not any(anomaly.values())

    # In coordinates (B,L) modulo the gauged generator, QQQL and the dressed
    # n^c n^c may be assigned charges (1,1) and (0,-2). Higgs has B=L=0.
    qqql_phases = (1, 1)
    dressed_phases = (0, -2)
    phase_determinant = qqql_phases[0] * dressed_phases[1] - qqql_phases[1] * dressed_phases[0]
    assert 3 * charges["Q"] + charges["L"] == 0
    assert 2 * charges["n^c"] - 4 * charges["H"] * 3 * k == 0
    assert phase_determinant != 0
    slope_checks = {}
    for slope in inputs["triangle_slopes"].values():
        slope = Fraction(slope)
        ratios = [
            abs(row["y0"] - slope * row["D"]) / Fraction(row["dimension"])
            for row in fields.values()
        ]
        assert all(ratio <= 3 for ratio in ratios)
        slope_checks[str(slope)] = str(max(ratios))
    slope = Fraction(inputs["triangle_slopes"]["k_negative" if k < 0 else "k_nonnegative"])
    lower = abs(order * k + slope) * 6 / 3
    bilinear_dimension = 2 * Fraction(fields["n^c"]["dimension"])
    dressing_dimension = 2 * (2 * Fraction(fields["H"]["dimension"]) + 1)
    upper = bilinear_dimension + 3 * abs(k) * dressing_dimension
    assert lower.denominator == upper.denominator == 1
    assert lower <= upper

    # Right-handed particles are conjugates of the named singlets. The Dirac
    # Yukawa constraints give the same electric charges for both chiralities.
    u, d, e, nu = (-Fraction(charges[name], order) for name in ("u^c", "d^c", "e^c", "n^c"))
    assert u - d == nu - e == 1
    proton, neutron = 2 * u + d, u + 2 * d
    neutron_ratio = neutron / abs(e)
    electric = {
        "u": u, "d": d, "e": e, "nu": nu, "p": proton, "n": neutron,
        "hydrogen_atom": proton + e, "neutron_over_abs_e": neutron_ratio,
    }

    # Above all retained thresholds, Y = y/6 and the index for a U(1)
    # representation is multiplicity times Y squared.
    y2 = {name: Fraction(charge, order)**2 for name, charge in charges.items()}
    fermion_index = sum(row["color"] * row["weak"] * y2[name] for name, row in fermions.items())
    scalar_index = fields["H"]["weak"] * y2["H"]
    b_y = Fraction(2, 3) * inputs["families"] * fermion_index + Fraction(1, 3) * scalar_index
    return {
        "anomaly_residuals": anomaly,
        "bY": str(b_y),
        "charges": {name: str(value) for name, value in electric.items()},
        "dressed_dimension": int(upper),
        "k": k,
        "kernel_cover_centers": kernel_centers,
        "lower_dimension": int(lower),
        "phase_breaking_determinant": phase_determinant,
        "triangle_max_ratios": slope_checks,
        "weights": charges,
        "weinberg_allowed": charges["L"] + charges["H"] == 0,
    }


def main():
    params = Path(__file__).with_name("PARAMETERS.json")
    inputs = json.loads(params.read_text(encoding="utf-8"))
    rows = [fixture(k, inputs) for k in inputs["fixtures_k"]]
    assert len({row["charges"]["neutron_over_abs_e"] for row in rows}) == len(rows)
    assert all(row["charges"]["hydrogen_atom"] == "0" for row in rows)
    for row in rows:
        k = row["k"]
        assert row["weinberg_allowed"] == (k == 0)
        assert row["bY"] == str(Fraction(41 + 192 * k + 576 * k**2, 6))

    benchmark = inputs["neutron_experiment"]
    uncertainty = Fraction(benchmark["charge_ratio_reported_uncertainty"])
    experimental_scale = abs(Fraction(benchmark["charge_ratio_center"])) + (
        benchmark["comparison_multiplier"] * uncertainty
    )
    minimum_nonzero = min(
        abs(Fraction(row["charges"]["neutron_over_abs_e"])) for row in rows if row["k"]
    )
    assert minimum_nonzero == Fraction(inputs["minimum_nonzero_abs_neutron_to_electron_charge_ratio"])
    assert experimental_scale < minimum_nonzero
    print(json.dumps({"fixtures": rows, "neutron_comparison_scale": str(experimental_scale)}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
