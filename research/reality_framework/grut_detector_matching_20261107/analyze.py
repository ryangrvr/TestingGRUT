"""Frozen M0-01 analysis. Reports missed and unresolved gates without fitting."""
import argparse
import json
from pathlib import Path

from continuum import tv


def analyze(directory, card):
    ref = json.loads((directory / "reference.json").read_text())
    target = ref["histories"]
    controls = ref["controls"]
    verify = ref["oscillator_verification"]
    thresholds = card["gates"]
    gates = {
        "reference_implementation": max(verify["tv_gaussian_vs_high"], verify["tv_cutoff_change"])
                                     < thresholds["oscillator_tv"],
        "continuum_reference_control": max(ref["maximum_history_tail_bound"],
                                            controls["mode_doubling_tv"],
                                            controls["quadrature_tightening_tv"])
                                        < thresholds["reference_tv"],
    }
    rows = []
    points = {}
    for path in sorted(directory.glob("point_*.json")):
        d = json.loads(path.read_text())
        points[path.stem] = d
        distances = {k: tv(d["histories"][k], p) for k, p in target.items()}
        rows.append({"file": path.name, "N": d["N"],
                     "bare_zero_control": d["bare_zero_control"],
                     "solver_refinement": d["solver_refinement"],
                     "mean": d["mean"], "mean_difference": d["mean"] - ref["spectrum"]["mean"],
                     "gap": d["gap"], "gap_difference": d["gap"] - ref["spectrum"]["first_gap"],
                     "variance": d["covariance"]["0"][0],
                     "variance_difference": d["covariance"]["0"][0] - ref["spectrum"]["variance"],
                     "tv_by_protocol": distances, "max_tv": max(distances.values()),
                     "numerical_diagnostics_pass": max(d["eigen_residuals"]) < thresholds["eigen_residual"]
                         and max(d["free_norm_error"],d["probability_norm_error"]) < thresholds["probability_norm"]})
    primary_names = [f"point_N{n}_primary" for n in card["primary_N"]]
    all_primary = all(name in points for name in primary_names)
    refine_name = f"point_N{card['refined_N']}_refined"
    normal_name = f"point_N{card['refined_N']}_primary"
    solver_tv = None
    if refine_name in points and normal_name in points:
        solver_tv = max(tv(points[refine_name]["histories"][k],points[normal_name]["histories"][k])
                        for k in target)
    gates["lattice_numerical_diagnostics"] = (all_primary and solver_tv is not None and
        solver_tv < thresholds["solver_refinement_tv"] and all(r["numerical_diagnostics_pass"] for r in rows))
    comparisons = {}
    if all_primary:
        last = points[primary_names[-1]]["histories"]
        penultimate = points[primary_names[-2]]["histories"]
        for k, p in target.items():
            errors = [tv(points[name]["histories"][k],p) for name in primary_names]
            step = tv(last[k],penultimate[k])
            comparisons[k] = {"primary_tv": errors, "finest_step_tv": step,
                              "fine_three_nonincreasing": errors[1]>=errors[2]>=errors[3]}
        gates["finite_resolution_matching"] = all(v["primary_tv"][-1]<=thresholds["matching_tv"]
                                                   for v in comparisons.values())
        gates["resolution_diagnostic"] = all(v["fine_three_nonincreasing"] and
                    v["finest_step_tv"]<=thresholds["finest_step_tv"] for v in comparisons.values())
    else:
        gates["finite_resolution_matching"] = None
        gates["resolution_diagnostic"] = None
    if not gates["reference_implementation"]:
        verdict = "implementation_verification_failed"
    elif not all_primary or solver_tv is None:
        verdict = "incomplete_or_resource_limited"
    elif not gates["lattice_numerical_diagnostics"] or not gates["continuum_reference_control"]:
        verdict = "numerical_control_failed"
    elif not gates["finite_resolution_matching"] or not gates["resolution_diagnostic"]:
        verdict = "finite_resolution_matching_gate_missed"
    else:
        verdict = "finite_protocol_numerical_matching_passes_not_a_uniform_physical_certificate"
    return {"protocol": card["id"], "verdict": verdict, "gates": gates,
            "solver_refinement_max_tv": solver_tv, "reference_mean":ref["spectrum"]["mean"],
            "reference_variance":ref["spectrum"]["variance"], "reference_gap":ref["spectrum"]["first_gap"],
            "comparison": comparisons, "points": rows,
            "limits": ["Finite-interval mass matching remains a specified prescription.",
                       "Instantaneous extended instruments are idealized; causal realization error is not bounded.",
                       "No massive-history error bound or evaluator-reliability certification follows."]}


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("directory",type=Path)
    args=parser.parse_args()
    card=json.loads((Path(__file__).parent/"protocol_m0_01.json").read_text())
    print(json.dumps(analyze(args.directory,card),indent=2,allow_nan=False))
