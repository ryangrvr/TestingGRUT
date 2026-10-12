"""Fast adversarial verification of this executed CP-12 follow-up packet.

Recomputes outward certificates and exact rational comparisons; checks
recorded complete-Hamiltonian identities, cutoff comparisons and bounds.
This is not an empirical validation or a second many-body eigensolver.
"""
from pathlib import Path
from decimal import Decimal
from fractions import Fraction
import ast
import hashlib
import json
import math
import platform

import numpy as np
import scipy
import networkx as nx

from variational_exclusion import finite_certificate
from near_regular_exclusion import certificate, asymptotic_defect_corollary

P = Path(__file__).resolve().parent


def read(name):
    return json.loads((P/name).read_text())


def main():
    groups = []
    for file in sorted(P.glob("*.py")):
        ast.parse(file.read_text(), filename=str(file))
    sources = read("parent_sources.json")
    model = P.parents[1]/"constructive_physics/attempt12/model.py"
    assert hashlib.sha256(model.read_bytes()).hexdigest() == sources["frozen_model_sha256"]
    for dep in sources["cluster_dependencies"]:
        raw = (P.parents[2]/dep["path"]).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == dep["sha256"]
        assert hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest() == dep["sha"]
    parent = P.parent/"cp12_even_width_cluster_control/check.py"
    account = read("even_width_execution_account.json")
    assert hashlib.sha256(parent.read_bytes()).hexdigest() == account["original_sha256"]
    try:
        ast.parse(parent.read_text())
    except SyntaxError:
        pass
    else:
        raise AssertionError("Pinned original was expected to have its recorded syntax error")
    for row in read("results.json")["results"]:
        a = row["exact_embedding_weight_method"]
        b = row["independent_energy_telescoping_method"]
        assert a > 0 and abs(a-b) < 1e-13
        assert [(r["vertices"], r["embedding_density_difference"])
                for r in row["contributing_clusters"]] == [(8, .125), (9, -1.)]
    groups.append("pinned sources, execution failure and repaired even-width control")

    for file in sorted(P.glob("regular_bound_*.json")):
        x = read(file.name)
        b = .1/x["M"]*math.sqrt((x["d"]+1)*x["Q"]*(x["M"]-x["Q"]))
        assert abs(b-x["universal_physical_bound"]) < 1e-13
        assert x["exact_neutral_coupling_norm"] <= b+1e-10
        assert x["max_nonregular_predecessors"] <= x["d"]+1
    groups.append("complete regular graph transition counts and operator norms")

    for file in sorted(P.glob("mixing_N*.json")):
        x = read(file.name)
        assert x["projected_formula_error"] < 1e-10 and x["residual"] < 1e-9
        assert x["full_energy"] <= min(x["best_graph_block_energy"], x["uniform_fiber_variational_energy"])+1e-9
        assert x["full_energy"] >= x["best_graph_block_energy"]-x["offdiag_norm_exact"]-1e-9
        assert x["rewiring_energy"] >= x["analytic_gap_lower_bound"]*(1-x["uniform_fiber_projection"])-1e-9
    assert abs(read("mixing_N4_Q3_L3.json")["full_energy"]-read("mixing_N4_Q3_L2.json")["full_energy"]) < 1e-8
    groups.append("exact Johnson fibers and independent projected hopping matrix")

    for file in sorted(P.glob("regular_full_*.json")):
        x = read(file.name)
        assert x["two_sector_energy_lower"] <= x["full_energy"] <= x["two_sector_energy_upper"]+1e-9
        assert x["coupling_norm"] <= x["universal_coupling_bound"]+1e-10
        assert x["residual"] < 1e-9
        for row in x.get("coherent_density_response", []):
            assert row["response_error"] <= row["proven_error_bound"]+1e-10
            assert row["departure"] <= row["departure_bound"]+1e-10
    for n in [4, 5]:
        assert abs(read(f"regular_full_N{n}_d2_L3.json")["full_energy"]-
                   read(f"regular_full_N{n}_d2_L2.json")["full_energy"]) < 1e-8
    x = read("regular_full_N4_d2_L3.json")
    assert x["regular_sector_energy"] < x["best_nonregular_fixed_block_energy"]
    assert x["nonregular_sector_energy"] < x["regular_sector_energy"]
    groups.append("full-law sector reversal and coherent prepared response")

    for host in read("star_bound_results.json")["controls"]:
        assert host["decomposition_error"] < 1e-10
        for x in host["stars"]:
            assert x["minimum_diagonal_slack"] >= -1e-10
            assert x["hopping_block_norm"] <= x["hopping_norm_bound"]+1e-10
            assert x["actual_star_minimum"] >= x["analytical_star_lower"]-1e-10
    groups.append("Gauss-law star lower bound and complete block decomposition")

    v = read("variational_exclusion_results.json")
    for x in v["small_hamiltonian_controls"]:
        assert x["projected_matrix_error"] < 1e-10
        assert x["rewiring_annihilation_error"] < 1e-10
        assert x["full_energy"] <= x["trial_energy"]+1e-9
    for saved in v["certificates"]:
        fresh = finite_certificate(saved["N"])
        for key in ["rayleigh_upper_total", "regular_lower_density",
                    "regular_minus_trial_gap_lower", "regular_ground_probability_upper"]:
            assert fresh[key] == saved[key], (saved["N"], key)
        assert Decimal(saved["regular_minus_trial_gap_lower"]) > 0
    groups.append("coherent matching trial and exact-regular outward certificates")

    near = read("near_regular_results.json")
    assert asymptotic_defect_corollary() == near["asymptotic_defect_corollary"]
    for saved in near["certificates"]:
        assert certificate(saved["N"], saved["K"]) == saved
    for row in near["full_graph_matrix_control"]["rows"]:
        assert row["max_incoming_crossings"] <= row["proven_incoming_bound"]
        for w in row["weight_checks"]:
            assert w["actual_max_weighted_row"] <= w["proven_row_bound"]+1e-8
    assert Fraction(near["asymptotic_defect_corollary"]["mean_absolute_degree_defect_density_liminf_lower"]) == Fraction(39, 400000)
    groups.append("nearly regular weighted mixing bounds and positive defect-density corollary")

    hosts = read("finite_host_tail_results.json")
    for h in hosts:
        assert abs(h["exact_relative_per_vertex"]-h["reconstructed_per_vertex"]) < 1e-11
        assert h["residual"] < 1e-9
        for row in h["partials"]:
            assert abs(row["partial_relative_per_vertex"]+row["exact_remainder"]-h["exact_relative_per_vertex"]) < 1e-13
    for n in [4, 6, 8]:
        values = [h["exact_relative_per_vertex"] for h in hosts if h["N"] == n]
        assert max(values)-min(values) < 1e-11
    groups.append("independent complete periodic hosts and finite remainders")

    schur = read("schur_N4_Q3_L3.json")
    assert abs(schur["full_energy"]-schur["exact_schur_root"]) < 1e-10
    for row in schur["rows"]:
        lo, hi = schur["complement_spectral_bracket"]
        a, b = lo-row["point"], hi-row["point"]
        r = (math.sqrt(b)-math.sqrt(a))/(math.sqrt(b)+math.sqrt(a))
        tail = 2/math.sqrt(a*b)*r**(row["degree"]+1)/(1-r)
        assert math.isclose(tail, row["resolvent_uniform_tail_bound"], rel_tol=1e-12)
        if row["below_complement_spectrum"]:
            assert row["energy_lower"] <= schur["full_energy"] <= row["energy_upper"]
    assert schur["rows"][-1]["energy_halfwidth"] < 1e-8
    groups.append("finite full-geometry Schur root and conditional polynomial remainder")

    out = dict(status="CP12_REMAINDER_AND_MIXING_CHECKS_PASS", group_count=len(groups),
               groups=groups, dependencies=dict(python=platform.python_version(),
               numpy=np.__version__, scipy=scipy.__version__, networkx=nx.__version__),
               scope="Executed model controls and analytical certificate checks; not empirical validation or a completed theory of spacetime.")
    (P/"verification.json").write_text(json.dumps(out, indent=2)+"\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
