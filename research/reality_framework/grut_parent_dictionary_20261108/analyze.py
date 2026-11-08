"""Frozen P1 analysis: necessary gates, trends and explicit unresolved errors."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent/"grut_detector_matching_20261107/raw_m1_01"


def tv(p,q):
    return .5*sum(abs(p.get(k,0)-q.get(k,0)) for k in set(p)|set(q))


def comparisons(p,q):
    if set(p) != set(q):
        raise ValueError("Unmatched history protocols")
    return {k:tv(p[k],q[k]) for k in p}


def analyze(out):
    out = Path(out)
    card = json.loads((HERE/"protocol_p1_01.json").read_text())
    missing = [p["label"] for p in card["points"] if not (out/(p["label"]+".json")).exists()]
    if missing:
        return {"verdict":"incomplete_resource_or_implementation", "complete":False, "missing":missing,
                "scope":"No physical matching verdict from incomplete prerequisites"}
    points = {p["label"]:json.loads((out/(p["label"]+".json")).read_text()) for p in card["points"]}
    old10 = json.loads((ARCHIVE/"boson_m05_E10_q128.json").read_text())
    old12 = json.loads((ARCHIVE/"boson_m05_E12_q128.json").read_text())
    lattice = json.loads((ARCHIVE/"lattice_N20.json").read_text())
    new12, new14, quad = [points[k] for k in ("E12_q128","E14_q128","E14_q256")]
    replay = comparisons(old12["histories"],new12["histories"])
    previous = comparisons(old10["histories"],old12["histories"])
    current = comparisons(new12["histories"],new14["histories"])
    quadrature = comparisons(new14["histories"],quad["histories"])
    matched = comparisons(lattice["histories"],new14["histories"])
    old_matched = comparisons(lattice["histories"],old12["histories"])
    overflow = {k:p["overflow"] for k,p in new14["histories"].items()}
    nested, g = new14["nested_E12_checks"], card["gates"]
    diagnostics = all(p["hermiticity_residual"]<g["hermiticity"]
        and p["eigen_residual"]<g["eigen_residual"]
        and p["minimum_defect_eigenvalue"]>=-g["negative_defect"]
        and p["probability_norm_error"]<g["probability_norm"] for p in points.values())
    gates = {
        "archive_replay":max(replay.values())<=g["archive_replay_tv"],
        "numerical_diagnostics":diagnostics,
        "nested_hamiltonian":nested["hamiltonian_max_entry"]<=g["nested_hamiltonian_max_entry"],
        "nested_weyl":nested["weyl_max_entry"]<=g["nested_weyl_max_entry"],
        "nested_defect_identity":nested["defect_decomposition_max_entry"]<=g["nested_defect_max_entry"],
        "reference_cutoff_stability":max(current.values())<=g["reference_cutoff_tv"],
        "reference_quadrature_stability":max(quadrature.values())<=g["reference_quadrature_tv"],
        "parent_instrument_overflow":max(overflow.values())<=g["overflow"],
        "finite_lattice_matching":max(matched.values())<=g["finite_matching_tv"],
        "overflow_lower_bound_consistency":all(matched[k]+1e-12>=overflow[k] for k in matched)
    }
    return {
        "verdict":"necessary_finite_gates_passed_not_continuum_certificate" if all(gates.values()) else "physical_matching_unresolved_prerequisite_or_comparison_miss",
        "complete":True, "gates":gates, "archive_replay_tv":replay,
        "cutoff_previous_TV_E10_E12":previous, "cutoff_current_TV_E12_E14":current,
        "cutoff_step_decreased_every_history":all(current[k]<previous[k] for k in current),
        "cutoff_max_step_ratio":max(current.values())/max(previous.values()),
        "quadrature_TV_E14":quadrature, "comparison_TV_N20_E12":old_matched,
        "comparison_TV_N20_E14":matched, "overflow_E14":overflow,
        "comparison_above_overflow_floor":{k:matched[k]-overflow[k] for k in matched},
        "nested_E12_checks":nested,
        "cost_seconds":{k:p["seconds"] for k,p in points.items()},
        "scope":"Reused solver, fixed-physics cutoff trend and necessary gates; no absolute continuum, lattice-spacing, boundary-mass or causal-detector error bound; M0/M1/E1/E2 verdicts unchanged"
    }
