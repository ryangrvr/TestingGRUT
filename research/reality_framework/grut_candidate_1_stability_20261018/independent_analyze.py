#!/usr/bin/env python3
"""Describe the frozen local minimization; do not rerun or alter its results."""
import hashlib
import json
from pathlib import Path

import numpy as np

here = Path(__file__).resolve().parent
results_path = here / "independent_results.json"
branch_path = here / "branch_results.json"
r = json.loads(results_path.read_text())
b = json.loads(branch_path.read_text())
comparisons = []
for fraction in r["protocol"]["charge_fractions"]:
    target = b["targets"][str(fraction)]
    reference = next(x["fine"] for x in target["found"] if x["segment"] == "lower_Omega")
    row = {"fraction": fraction, "Qbar": fraction * r["protocol"]["target_Qbar"],
           "reference_Ebar": reference["Ebar"], "reference_Omega": reference["Omega"], "meshes": {}}
    for n in r["protocol"]["cell_counts"]:
        runs = [x for x in r["runs"] if x["fraction"] == fraction and x["N"] == n and not x.get("skipped")]
        energies = [x["Ebar"] for x in runs]
        omegas = [x["Omega"] for x in runs]
        e = float(np.mean(energies))
        omega = float(np.mean(omegas))
        row["meshes"][str(n)] = {
            "mean_Ebar": e, "mean_Omega": omega,
            "Ebar_error_vs_BVP": e - reference["Ebar"],
            "Omega_error_vs_BVP": omega - reference["Omega"],
            "Ebar_start_spread": max(energies) - min(energies),
            "Omega_start_spread": max(omegas) - min(omegas),
            "max_projected_gradient_inf": max(x["projected_gradient_inf"] for x in runs),
        }
    coarse = row["meshes"]["800"]
    fine = row["meshes"]["1600"]
    row["energy_error_reduction_factor"] = coarse["Ebar_error_vs_BVP"] / fine["Ebar_error_vs_BVP"]
    row["second_order_Richardson_Ebar"] = (4 * fine["mean_Ebar"] - coarse["mean_Ebar"]) / 3
    row["Richardson_Ebar_error_vs_BVP"] = row["second_order_Richardson_Ebar"] - reference["Ebar"]
    row["second_order_Richardson_Omega"] = (4 * fine["mean_Omega"] - coarse["mean_Omega"]) / 3
    row["Richardson_Omega_error_vs_BVP"] = row["second_order_Richardson_Omega"] - reference["Omega"]
    comparisons.append(row)

equal_half = {}
for n in (800, 1600):
    ref = next(x for x in comparisons if x["fraction"] == 1.0)["meshes"][str(n)]
    half = next(x for x in comparisons if x["fraction"] == 0.5)["meshes"][str(n)]
    equal_half[str(n)] = 2 * half["mean_Ebar"] - ref["mean_Ebar"]
equal_half["second_order_Richardson"] = (4 * equal_half["1600"] - equal_half["800"]) / 3
equal_half["BVP_reference"] = (
    2 * next(x for x in comparisons if x["fraction"] == 0.5)["reference_Ebar"]
    - next(x for x in comparisons if x["fraction"] == 1.0)["reference_Ebar"]
)
out = {
    "independent_results_sha256": hashlib.sha256(results_path.read_bytes()).hexdigest(),
    "branch_results_sha256": hashlib.sha256(branch_path.read_bytes()).hexdigest(),
    "averaging_rule": "Use arithmetic means over both frozen starts; report their spreads. No result was selected by best agreement.",
    "comparison": comparisons,
    "equal_half_branch_fission_difference": equal_half,
    "gradient_preflight_max_relative_discrepancy": max(x["relative_discrepancy"] for x in r["gradient_preflight"]),
    "number_requested": 12,
    "number_optimizer_success": sum(bool(x.get("optimizer_success")) for x in r["runs"]),
    "number_requested_projected_gradient_met": sum(bool(x.get("projected_gradient_meets_requested_gtol")) for x in r["runs"]),
    "limits": [
        "All optimizer stops use relative-energy ftol; none satisfies the frozen projected-gradient threshold.",
        "Two cell counts at fixed L40 are mesh refinement, not domain refinement or interval bounds.",
        "Richardson estimates assume leading second-order discretization and are cross-checks, not certified error bounds.",
        "Fields are restricted to the nodeless positive F>=0,0<=S<=1 radial class; no global or nonradial minimization is claimed.",
        "Both source-BVP and explicit Gaussian starts share the same action and implementation context; this is an independent method, not external independent validation.",
    ],
}
(here / "independent_comparison.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
