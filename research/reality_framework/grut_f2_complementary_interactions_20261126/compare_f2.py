"""Compare exact values from the two frozen F2 implementations."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import json
from pathlib import Path


def compare(production, critic):
    counts = Counter()

    def verify(group, left, right):
        if isinstance(left, dict):
            if left.keys() != right.keys():
                raise AssertionError(f"Coordinate mismatch in {group}")
            for key in left:
                verify(group, left[key], right[key])
        elif isinstance(left, list):
            if len(left) != len(right):
                raise AssertionError(f"Array size mismatch in {group}")
            for a, b in zip(left, right):
                verify(group, a, b)
        else:
            if F(left) != F(right):
                raise AssertionError(f"{group}: {left} differs from {right}")
            counts[group] += 1

    spring = production["f1_pair"]
    for row in spring["rows"]:
        for name, primary in row["networks"].items():
            independent = critic["spring"][name]["by_shift"][str(row["z"])]
            for a, b in (
                ("base", "baseline"), ("root_pin", "root_pin"), ("uniform_mass", "uniform_mass"),
                ("single_masses", "masses"), ("collective_channel", "aux_collective"),
                ("local_channel", "aux_local"),
            ):
                verify("f1_linear", primary[a], independent[b])
            verify("f1_linear", {i: values[0] for i, values in primary["single_pins"].items()},
                   independent["pins"])
            if "reconstructed_stiffness" in primary:
                independent = critic["spring"][name]["reconstruction"]
                verify("f1_reconstruction", primary["reconstructed_stiffness"], independent["stiffness"])
                verify("f1_reconstruction", primary["single_pins"], independent["singletons"])
                verify("f1_reconstruction", primary["pair_pins"], independent["pairs"])

    for name, primary in spring["quartic"].items():
        independent = critic["spring"][name]["quartic"]
        verify("quartic", primary["cubic_time_power"], independent["quartic_time_order"])
        verify("quartic", primary["cubic_coefficient"], independent["quartic_coefficient"])
        verify("quartic", primary["travel_power"], independent["root_to_tag_order"])

    for row in production["singleton_obstruction"]["rows"]:
        for name, data in row["networks"].items():
            independent = critic["stress"]["graphs"][name]["by_shift"][str(row["z"])]
            verify("singleton_obstruction", data["base"], independent["baseline"])
            verify("singleton_obstruction", data["pair_pin"], independent["pair"])
            verify("singleton_obstruction", data["single_pins"], [
                independent["singletons"][str(i)] for i in range(1, 17)
            ])

    theory = production["continuum_certificates"]
    independent = critic["continuum"]
    weak = independent["channels"]["weak"]
    strong = independent["channels"]["strong"]
    for left, right in (
        (theory["box_derivative_bound"], independent["box_cap"]),
        (theory["volume_over_pi"], independent["volume_over_pi"]),
        (theory["core_volume_over_pi"], independent["plateau_volume_over_pi"]),
        (theory["free_massive_wave_bound"], independent["free_wave_bounds"]),
        (theory["massive_neumann_norm"], strong["mass_kernel_contraction"]),
        (theory["massive_distorted_wave_floor"], strong["mass_wave_floor"]),
        (theory["massive_signal_lower"]["strong"], strong["mass_gap_lower"]),
        (theory["massive_signal_lower"]["weak"], weak["mass_gap_lower"]),
        (theory["weak_minimal_neumann_norm"], weak["minimal_contraction"]),
        (theory["weak_minimal_remainder_bound"], weak["minimal_remainder_upper"]),
        (theory["weak_minimal_remainder_ratio"], weak["minimal_remainder_to_linear_floor"]),
        (theory["weak_minimal_signal_lower"], weak["minimal_gap_lower"]),
        (theory["strong_minimal_norm_certificate"], strong["minimal_contraction"]),
    ):
        verify("continuum", left, right)

    primary = next(row for row in spring["rows"] if row["z"] == 1)
    verify("f1_linear", primary["contrasts"]["tagged_mass_gap"], critic["finite_time"]["untruncated_gap"])
    outer_lo, outer_hi = map(F, independent["minimal_to_massive_linearized_ratio_bounds"])
    narrow_lo, narrow_hi = map(F, theory["linearized_signal_ratio_interval"])
    ratio_consistent = outer_lo <= 1 / narrow_hi <= 1 / narrow_lo <= outer_hi
    tail_checks = (
        F(primary["contrasts"]["tagged_mass_gap"]) - F(spring["finite_tail_cap"]) > F("0.0066242")
        and F(critic["finite_time"]["gap_lower"]) > F("0.0066242")
    )
    if not all((production["passed"], ratio_consistent, tail_checks, independent["early_causal_blind"],
                independent["future_causal_blind"])):
        raise AssertionError("Registered or independent interval check failed")
    return {
        "exact_scalar_comparisons": sum(counts.values()),
        "comparisons_by_group": dict(counts),
        "exact_agreement": True,
        "both_time_tail_bounds_separate": tail_checks,
        "linearized_ratio_bounds_compatible": ratio_consistent,
        "independent_causal_blind_controls": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("builder", type=Path)
    parser.add_argument("critical", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = compare(json.loads(args.builder.read_text()), json.loads(args.critical.read_text()))
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"{result['exact_scalar_comparisons']} exact cross-comparisons passed")


if __name__ == "__main__":
    main()
