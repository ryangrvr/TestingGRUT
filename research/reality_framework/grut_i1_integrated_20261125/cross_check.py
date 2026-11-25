#!/usr/bin/env python3
"""Compare the frozen scalar reduction to the separate density implementation."""

import argparse
import itertools
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text())


def key(row):
    return tuple(row[name] for name in ("regime", "assay", "phase"))


def compare():
    cfg = read("PARAMETERS.json")
    reduced = read("RESULTS.json")["bridge"]
    direct = read("REVIEW_CHECK.json")
    n = cfg["direct_density_slices"]
    expected = set(itertools.product(cfg["regimes"], cfg["assays"], cfg["phase_over_pi"]))
    physical_rows = direct["rows"]
    if len(physical_rows) != len(expected) or {key(r) for r in physical_rows} != expected:
        raise AssertionError("Direct grid differs from the freeze")
    if any(r["slices"] != n for r in physical_rows):
        raise AssertionError("Direct runs used different contact counts")
    lookup = {key(r): r for r in reduced if r["slices"] == n}
    errors = []

    def check(ref, matrix, current_cont, memory_contrast):
        for u, v in (
            (ref["fidelity_assisted"], matrix["assisted"]),
            (ref["fidelity_blind"], matrix["blind"]),
            (ref["thermal"]["excited_probability"], matrix["excited_probability"]),
            (ref["memory_initial_contrast"], matrix["memory_initial_contrast"]),
            (current_cont, matrix["memory_current_contrast"]),
        ):
            errors.append(abs(u - v))
        for a, b in zip(ref["jitter"], matrix["jitter"], strict=True):
            if a["epsilon"] != b["epsilon"]:
                raise AssertionError("Jitter grid differs")
            errors.extend((abs(a["assisted"] - b["assisted"]), abs(a["blind"] - b["blind"])))
        z = complex(*ref["phase_factor"])
        coherence = math.sqrt(ref["thermal"]["attenuation"])
        branch = [complex(*v) for v in matrix["raw_conditional_offdiag"]]
        expected_zero = coherence * (z.real + 1j * memory_contrast * z.imag) / 4
        errors.extend((abs(branch[0] - expected_zero),
                       abs(branch[1] - expected_zero.conjugate()),
                       abs(complex(*matrix["raw_offdiag"]) - coherence * z.real / 2)))

    for actual in physical_rows:
        ref = lookup[key(actual)]
        d = ref["memory_initial_contrast"]
        check(ref, actual, ref["memory_current_sign_contrast"], d)
        limit = ref["continuum"]
        d_limit = limit["memory_initial_contrast"]
        rate = cfg["regimes"][ref["regime"]]["source_rate"]
        check(limit, actual["continuum"], d_limit * math.exp(-2 * rate * cfg["episode_time"]), d_limit)
        for side in ("hot", "cold"):
            errors.append(abs(ref["thermal"][side + "_heat_into_cells"] - actual[side + "_heat"]))

    thresholds = cfg["thresholds"]
    if max(errors) > thresholds["density_crosscheck_abs"]:
        raise AssertionError(max(errors))
    if direct["min_eigenvalue"] < -thresholds["density_crosscheck_abs"]:
        raise AssertionError("Direct evolution lost positivity")
    if direct["max_energy_error"] > thresholds["energy_abs"]:
        raise AssertionError("Direct energy check failed")
    if direct["min_entropy_production"] < thresholds["entropy_floor"]:
        raise AssertionError("Direct entropy check failed")
    if direct["max_entropy_identity_error"] > thresholds["density_crosscheck_abs"]:
        raise AssertionError("Direct entropy identity check failed")
    return {"rows": len(physical_rows), "scalar_comparisons": len(errors),
            "max_abs_discrepancy": max(errors),
            "max_direct_energy_residual": direct["max_energy_error"],
            "max_direct_entropy_identity_residual": direct["max_entropy_identity_error"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    result = compare()
    serialized = json.dumps(result, indent=2) + "\n"
    args = parser.parse_args()
    if args.output:
        args.output.write_text(serialized)
    print(serialized, end="")


if __name__ == "__main__":
    main()
