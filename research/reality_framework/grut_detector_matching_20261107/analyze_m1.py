"""Frozen analysis of M1-01; no extrapolation or data-dependent fit."""
import json
from pathlib import Path


def distance(p,q):
    return .5*sum(abs(p.get(s,0)-q.get(s,0)) for s in set(p)|set(q))


def max_distance(p,q):
    return max(distance(p[k],q[k]) for k in p)


def analyze(out,card):
    out=Path(out)
    def read(name):
        return json.loads((out/name).read_text()) if (out/name).exists() else None
    verification=read("verification.json")
    bosons={e:read(f"boson_m05_E{e}_q128.json") for e in card["energy_cutoffs"]}
    free=read("boson_m0_E10_q128.json")
    free6=read("boson_m0_E6_q128.json")
    quad=read("boson_m05_E10_q256.json")
    lattice={n:read(f"lattice_N{n}.json") for n in card["lattice_N"]}
    refined=read("lattice_N14_refined.json")
    complete=all(x is not None for x in [verification,free,free6,quad,refined,*bosons.values(),*lattice.values()])
    if not complete:
        return {"verdict":"incomplete_or_resource_unresolved","complete":False,
                "missing_boson_cutoffs":[e for e,b in bosons.items() if b is None],
                "missing_lattice_sites":[n for n,b in lattice.items() if b is None]}
    g=card["gates"]
    cutoff=max_distance(bosons[10]["histories"],bosons[12]["histories"])
    qchange=max_distance(bosons[10]["histories"],quad["histories"])
    overflow=max(p["overflow"] for p in bosons[12]["histories"].values())
    step=max_distance(lattice[18]["histories"],lattice[20]["histories"])
    solver=max_distance(lattice[14]["histories"],refined["histories"])
    comparisons={str(n):{k:distance(p,bosons[12]["histories"][k])
                          for k,p in b["histories"].items()} for n,b in lattice.items()}
    all_bosons=[*bosons.values(),free,free6,quad]
    gates={
        "algebraic_verification":verification["vertex_max_error"]<g["vertex_error"] and max(
            verification["normal_factor_images_vs_integral"],verification["normal_factor_image_doubling"])<g["normal_factor_error"],
        "free_reference":max(free["free_reference_tv"].values())<=g["free_reference_tv"],
        "boson_numerical_diagnostics":all(b["hermiticity_residual"]<g["hermiticity"] and b["eigen_residual"]<g["eigen_residual"] and b["minimum_defect_eigenvalue"]>=-g["negative_defect"] and b["probability_norm_error"]<g["probability_norm"] for b in all_bosons),
        "reference_cutoff_stability":cutoff<=g["reference_cutoff_tv"],
        "reference_quadrature_stability":qchange<=g["reference_quadrature_tv"],
        "instrument_overflow":overflow<=g["overflow"],
        "lattice_numerical_diagnostics":all(max(b["eigen_residuals"])<g["eigen_residual"] and b["free_norm_error"]<g["probability_norm"] and b["probability_norm_error"]<g["probability_norm"] for b in [*lattice.values(),refined]),
        "lattice_solver_refinement":solver<g["solver_refinement_tv"],
        "lattice_refinement":step<=g["lattice_finest_step_tv"],
        "finite_matching":max(comparisons["20"].values())<=g["matching_tv"]}
    return {"verdict":"finite_numerical_gates_passed" if all(gates.values()) else "matching_unresolved_prerequisite_or_comparison_miss",
            "complete":True,"gates":gates,"cutoff_max_tv":cutoff,"quadrature_max_tv":qchange,
            "overflow_max":overflow,"lattice_step_max_tv":step,"solver_max_tv":solver,
            "comparison_tv":comparisons,
            "boson_summary":[{k:b[k] for k in ["cutoff","dimension","gap","energies_low","seconds"]} for b in bosons.values()],
            "scope":"Finite numerical gates only; no absolute interacting continuum or causal-instrument error theorem."}
