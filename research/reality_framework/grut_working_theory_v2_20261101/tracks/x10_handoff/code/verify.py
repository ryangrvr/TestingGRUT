#!/usr/bin/env python3
"""Verify the X10 finite-control evidence and the physical H55 scope."""

import argparse
import hashlib
import json
import math
from pathlib import Path

TRACK = Path(__file__).resolve().parents[1]
RESULTS = TRACK / "results"
ARCHIVES = ("controls", "h55_seed", "h55_probe", "h55_live", "handoff", "retention", "hybrid")
TIMES = ("1.25", "2.25", "3.25")
RULES = ("energy", "energy_plus_one", "energy_plus_two")


def read(folder, name):
    return json.loads((folder / f"{name}.json").read_text())


def matches(archive, replay):
    if isinstance(archive, dict):
        return archive.keys() == replay.keys() and all(matches(v, replay[k]) for k, v in archive.items())
    if isinstance(archive, list):
        return len(archive) == len(replay) and all(matches(x, y) for x, y in zip(archive, replay, strict=True))
    if isinstance(archive, float):
        return math.isclose(archive, replay, rel_tol=2e-8, abs_tol=2e-10)
    return archive == replay


def jump(row):
    return max(map(abs, row["jump_over_reference_density"].values()))


def endpoint(row):
    return max(row["terminal_relative_errors"].values())


def budget(row):
    return (row["energy_pressure_error"] < 1e-5 and row["reciprocal_quantum_work_error"] < 1e-5
            and row["momentum_error"] < 1e-9)


def retained(row):
    j = row["jump_over_reference_density"]
    bound = row["normalized_energy_jump_bound"]
    return (0 <= bound <= 0.5e-6 and abs(j["rho"]) <= bound + 1e-12
            and abs(j["P"] - j["rho"] / 3) < 1e-12
            and max(abs(j["momentum"]), abs(j["work"])) < 1e-12
            and row["orders_min"] >= 1 and row["pauli_roundoff_clip"] < 1e-10
            and 0 < row["stored_covariance_fraction"] <= 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    replay = parser.add_mutually_exclusive_group()
    replay.add_argument("--replay", type=Path)
    replay.add_argument("--spot-replay", type=Path)
    args = parser.parse_args()
    archive = {name: read(RESULTS, name) for name in ARCHIVES}
    controls, seed, probe, live, handoff, retention, hybrid = (archive[name] for name in ARCHIVES)
    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)

    inputs = read(RESULTS, "inputs")
    for path, expected in inputs["sha256"].items():
        actual = hashlib.sha256((TRACK.parent / path).read_bytes()).hexdigest()
        check(f"input {path}", actual == expected)
    homogeneous = controls["homogeneous"]
    error = [row["max_state_error"] for row in homogeneous["steps"]]
    check("homogeneous axial rotation", homogeneous["axial_rotation_error"] < 1e-12)
    check("homogeneous split converges against dense solution", error[-1] < 5e-6 and error[0] / error[-1] > 15)
    spatial = controls["spatial_ode"]
    error = [row["max_error"] for row in spatial["comparisons"]]
    check("spatial split converges against independent ODE", error[-1] < 2e-6 and error[0] / error[-1] > 15)
    check("fixed-background reversibility", spatial["reversal_error"] < 1e-10)
    action = controls["action"]
    check("field variation and volume derivative", action["field_variation_error"] < 1e-9
          and action["pressure_relative_error"] < 1e-7)
    check("matched vacuum derivative interpolation", max(row["loga_derivative"] for row in action["vacuum_table"]) < 1e-7)
    static, massless = (controls["vacuum"][name] for name in ("static", "massless"))
    check("static coupled vacuum subtracts to roundoff", max(static[k] for k in
          ("max_quantum_energy", "max_quantum_pressure", "max_quantum_work")) < 1e-10)
    check("massless drive gives zero quantum work", massless["max_quantum_work"] == 0 and budget(massless["ledger"]))

    for name, case in handoff["cases"].items():
        check(f"{name} slab energy pressure and momentum ledger", budget(case["ledger"]))
        check(f"{name} finite run ends before discrete recurrence", case["terminal"]["t"] < case["recurrence"])
    ref = handoff["cases"]["reference"]
    bare = ref["static_projection"]["bare"]
    check("bare occupations inject energy into static coupled vacuum", bare["rho_change"] > 0.01
          and abs(bare["count_change"]) < 1e-12)
    dressed = [rows["dressed"] for rows in ref["handoffs"].values()]
    check("interacting energy occupations preserve energy but lose work and pressure", all(
          abs(row["jump_over_reference_density"]["rho"]) < 1e-12 and
          min(abs(row["jump_over_reference_density"][k]) for k in ("P", "work")) > 1e-6 and
          endpoint(row) > 0.01 for row in dressed))
    separated = [rows["separated"] for rows in ref["handoffs"].values()]
    check("separated blocks keep scalar work and lose interaction energy", all(
          abs(row["jump_over_reference_density"]["work"]) < 1e-12 and
          abs(row["jump_over_reference_density"]["rho"]) > 1e-6 for row in separated))
    fine = retention["cases"]["bath_spacing"]
    check("fixed seven-direction energy miss survives bath spacing test", ref["passes_all_boundaries"]["r7"]
          and not fine["passes_all_boundaries"]["r7"] and
          fine["boundaries"]["3.25"]["r7"]["jump_over_reference_density"]["rho"] > 1e-6)

    for name, case in retention["cases"].items():
        for method in RULES[:2]:
            check(f"{name} state splice {method}", tuple(case["boundaries"]) == TIMES
                  and all(retained(row[method]) and budget(row[method]) and jump(row[method]) < 1e-6
                          and endpoint(row[method]) < 0.01 for row in case["boundaries"].values())
                  and max(case["boundary_spreads"][method].values()) < 0.01)
    for order, examples in retention["cases"]["reference"]["derivative_hierarchy"].items():
        check(f"frozen Hamiltonian derivative hierarchy through order {order}", all(
              len(err) == int(order) + 2 and max(err[:-1]) < 1e-12 and err[-1] > 1e-10
              for err in examples.values()))

    for name, case in hybrid["cases"].items():
        boundary = case["boundaries"]
        check(f"{name} hybrid generator switch preserves boundary observables", tuple(boundary) == TIMES and all(
              max(map(abs, row["model_switch_difference"].values())) < 1e-12
              for rows in boundary.values() for row in rows.values()))
        check(f"{name} energy-selected hybrid retains the early endpoint miss", all(
              jump(rows["energy"]) < 1e-6 for rows in boundary.values())
              and endpoint(boundary[TIMES[0]]["energy"]) > 0.01)
        expected = name not in ("bath_cutoff", "bath_holdout")
        check(f"{name} plus-one decision preserves cutoff misses", case["passes_all_boundaries"][RULES[1]] == expected
              and (endpoint(boundary[TIMES[0]][RULES[1]]) < 0.01) == expected)
        check(f"{name} two extra directions satisfy moved boundaries",
              case["passes_all_boundaries"][RULES[2]] and all(
              retained(rows[RULES[2]]) and budget(rows[RULES[2]]) and jump(rows[RULES[2]]) < 1e-6
              and endpoint(rows[RULES[2]]) < 0.01
              and max(rows[RULES[2]]["downstream_normalized_stress_error"].values()) < 1e-6
              for rows in boundary.values()) and max(case["boundary_spreads"][RULES[2]].values()) < 0.01)
    resolutions = handoff["terminal_resolution_differences"]
    check("bath spacing cutoff and spatial reference shifts remain visible",
          all(abs(resolutions[k]["daughter_total"]) > 0.01 for k in
              ("bath_spacing", "bath_cutoff", "longitudinal", "transverse"))
          and abs(resolutions["half_step"]["daughter_total"]) < 1e-4)
    work = handoff["instantaneous_work_resolution_differences"]
    check("phase-sensitive work is timestep stable and spatially underresolved",
          work["half_step"] < 2e-5 and min(work["longitudinal"], work["transverse"]) > 0.1)
    base = hybrid["cases"]["reference"]["reference_terminal"]
    wide = hybrid["cases"]["bath_holdout"]["reference_terminal"]
    check("holdout production changes exceed the splice error", wide["daughter_total"] / base["daughter_total"] - 1 > 0.03
          and wide["daughter_difference"] / base["daughter_difference"] - 1 > 0.15)

    errors = seed["diagnostics"]
    check("H55 scalar seed agrees with archived initial field and gradient diagnostics",
          max(abs(v["relative_difference"]) for k, v in errors.items() if k != "xdot_std") < 2e-6)
    check("H55 seed velocity mismatch is measured", 3e-4 < abs(errors["xdot_std"]["relative_difference"]) < 5e-4)
    check("planar scalar gradients capture less than nine ten-thousandths of the seed",
          all(0 < row["first_harmonic_gradient_fraction"] < row["gradient_fraction"] < 9e-4
              for name in ("sigma", "chi") for row in seed["axis_projections"][name]))
    check("physical probes pin the same early seed and trace", seed["diagnostic_sha256"] ==
          probe["diagnostic_sha256"] == live["diagnostic_sha256"] ==
          inputs["sha256"]["x2_strip/results/lower_H_diag.txt.gz"] and
          probe["source_trace_sha256"] == live["source_trace_sha256"])
    small, half, uniform = (probe["cases"][k] for k in ("spatial", "spatial_half", "homogeneous"))
    field = "spatial_sterile_bath_norm"
    check("physical axial spatial correlations survive timestep refinement", small["terminal"][field] > 2e-5
          and abs(half["terminal"][field] / small["terminal"][field] - 1) < 1e-4
          and uniform["terminal"][field] < 1e-9)
    check("physical prescribed probes are unitary before recurrence", all(
          row["conformal_duration"] < row["minimum_discrete_recurrence"] and row["orthonormality_error"] < 1e-9
          for row in probe["cases"].values()))
    check("prescribed energy work discrepancy shrinks with timestep",
          half["max_relative_energy_work_error"] < 4e-7 and
          small["max_relative_energy_work_error"] / half["max_relative_energy_work_error"] > 3.5)
    l0, l1, denser, off = (live["cases"][k] for k in ("reference", "half_step", "bath_spacing", "no_bath"))
    check("physical axial reciprocal ledger", all(row["max_total_energy_work_error"] < 1e-5 for row in live["cases"].values())
          and l0["max_total_energy_work_error"] / l1["max_total_energy_work_error"] > 3.5
          and abs(l1["integrated_axial_work"] / l0["integrated_axial_work"] - 1) < 0.01)
    check("axial sterile–bath correlations require the coupled bath", denser["terminal_correlations"][field] > 2e-5
          and off["terminal_correlations"][field] < 1e-9)
    z0, z1, zoff = (complex(*r["sigma_mode_terminal"]) for r in (l0, l1, off))
    check("extra bath feedback on axial amplitude below timestep change", abs(z0 - zoff) < abs(z0 - z1))
    check("archive preserves the restricted scope", all(d["scientific_acceptance"] is False for d in archive.values()))

    if args.replay:
        for name, expected in archive.items():
            check(f"full replay {name}", matches(expected, read(args.replay, name)))
    if args.spot_replay:
        for name in ("controls", "h55_seed", "h55_probe", "h55_live"):
            check(f"fresh replay {name}", matches(archive[name], read(args.spot_replay, name)))
        check("fresh reference handoff replay", matches(ref, read(args.spot_replay / "handoff", "reference")))
        replayed = read(args.spot_replay, "hybrid")["cases"]
        for name in ("reference", "bath_holdout"):
            check(f"fresh hybrid replay {name}", matches(hybrid["cases"][name], replayed[name]))

    passed = sum(checks.values())
    result = dict(checks=checks, passed=passed, total=len(checks),
                  replay_mode="full" if args.replay else "spot" if args.spot_replay else "archive",
                  finite_hybrid_boundaries=3 * len(hybrid["cases"]),
                  H55_particle_handoff_certified=False, inherited_strict_acceptance="0/135",
                  scientific_acceptance=False)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"X10 {result['replay_mode']} evidence checks {passed}/{len(checks)}")
    for name, ok in checks.items():
        if not ok:
            print(f"FAIL {name}")
    raise SystemExit(passed != len(checks))


if __name__ == "__main__":
    main()
