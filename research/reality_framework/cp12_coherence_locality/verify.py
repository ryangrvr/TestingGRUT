"""Verify the executed coherence/locality packet and outward certificates.

This is model-level scientific verification, not independent experimental
validation. Heavy complete-H scripts are provided separately for reruns.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal
import ast
import hashlib
import json
import math
import platform

import numpy as np
import scipy
import networkx as nx

from locality_pool import certificate
from light_cone_control import chain_tail

P=Path(__file__).resolve().parent
ROOT=P.parents[2]


def read(name):return json.loads((P/name).read_text())


def main():
    groups=[]
    source=read("parent_sources.json")
    for entry in [source["frozen_model"]]+source["inherited_dependencies"]:
        assert hashlib.sha256((ROOT/entry["path"]).read_bytes()).hexdigest()==entry["sha256"]
    for file in P.glob("*.py"):ast.parse(file.read_text())
    groups.append("pinned frozen law and inherited analytical dependency hashes")

    ground=[read(file.name) for file in sorted(P.glob("full_ground_*.json"))]
    for x in ground:
        assert x["residual"]<1e-9
        assert x["coherence_energy"]>0
        assert math.isclose(sum(x["edge_record_costs"]),2*x["coherence_energy"],abs_tol=1e-10)
        assert max(x["edge_record_costs"])-min(x["edge_record_costs"])<1e-9
        for r in x["record_controls"]:
            expected=math.sin(r["theta"])**2*x["coherence_energy"]
            assert abs(r["actual_system_energy_increase"]-expected)<1e-10
        assert abs(x["completely_graph_dephased_energy"]-x["full_energy"]-x["coherence_energy"])<1e-12
    groups.append("joint system-reference energy and exact single-edge record sum rule")

    for x in ground:
        rows=x["graph_probabilities"]
        assert abs(sum(r["probability"] for r in rows)-1)<1e-10
        total=math.comb(x["N"]*(x["N"]-1)//2,x["Q"])
        assert total==len(rows)
        uniform=sum(r["connected"] for r in rows)/total
        if x["N"]==4:assert uniform==.8
        if x["N"]==5:assert abs(uniform-37/42)<1e-15
        weighted=sum(r["probability"]*r["connected"] for r in rows)
        assert weighted>uniform
        assert abs(weighted-x["geometry_statistics"]["connected_probability"])<1e-12
    a=read("full_ground_N5_Q5_L2.json");b=read("full_ground_N5_Q5_L3.json")
    assert abs(a["coherence_energy"]-b["coherence_energy"])<2e-9
    for k in a["geometry_statistics"]:
        assert abs(a["geometry_statistics"][k]-b["geometry_statistics"][k])<2e-8
    groups.append("complete graph distributions, exact uniform baseline and electric cutoff controls")

    pools=read("locality_pool_results.json")
    for x in pools["neutral_norm_controls"]:
        for r in x["rows"]:
            assert abs(r["internal_adjacency_norm"]-r["quotient_norm"])<1e-9
            assert abs(r["crossing_adjacency_norm_square"]-r["exact_boundary_square"])<1e-8
    phys=pools["full_physical_norm_and_dynamics"]
    for r in phys["norm_controls"]:
        assert abs(r["actual_internal_norm"]-r["expected_internal_norm"])<1e-10
        assert abs(r["actual_crossing_norm_square"]-r["expected_crossing_norm_square"])<1e-10
    groups.append("exact neighborhood quotient and complete physical mixing norms")

    assert abs(phys["prepared_mean_energy"]-phys["predicted_mean_energy"])<1e-10
    assert abs(phys["prepared_energy_variance"]-phys["predicted_energy_variance"])<1e-10
    for r in phys["time_controls"]:
        assert r["mean_long_bonds"]<=r["proven_mean_bound"]+1e-10
    groups.append("full-law coherent preparation moments and outside-bond time bound")

    for name in ["cycle_pool","branched_pool"]:
        x=read(f"prepared_response_{name}.json")
        assert x["source_prepared_sector_norm"]<=1+1e-10
        for r in x["rows"]:
            assert r["full_pool_error"]<=r["pool_error_bound"]+1e-10
            assert r["full_fixed_error"]<=r["total_error_bound"]+1e-10
            if "impulse_error" in r:assert r["impulse_error"]<1e-7
        for r in x["connectivity_response_controls"]:
            assert abs(r["exact_block_connectivity_response"])<1e-11
            assert abs(r["full_connectivity_response"])<=r["proven_response_bound"]+1e-10
    groups.append("gauge-invariant source impulse, detector response and suppressed connectivity response")

    light=read("light_cone_results.json")
    for x in light["controls"]:
        D=x["maximum_incidence_degree"]
        for r in x["rows"]:
            bound=2*D/(D+1)*chain_tail(2*(D+1)*r["T"],r["incidence_distance"])
            assert abs(min(2.,bound)-r["chain_bound"])<1e-12
            assert r["actual_commutator_norm"]<=bound+1e-11
        assert all(r["order"]<r["distance"] for r in x["integer_nested_zero_controls"])
    groups.append("conditional fermionic commutator cone and exact integer support-order zeros")

    host=read("three_dimensional_host_results.json")
    for x in host["host_controls"]:
        assert x["pool_maximum_degree"]==5
        assert abs(x["sparse_gap"]-x["independently_computed_bloch_gap"])<1e-9
    for x in host["analytical_response_and_preparation_bounds"]:
        n=x["N"];M=n*(n-1)//2;q=x["Q"];S=x["S"];r=x["vacancies"];T=x["observation_time"]
        b2=F(q*(M-S)*(r+1),100*M*M);vin=F(q*r,10*M)
        assert F(x["exact_response_error_bound"])==4*b2*T*T+4*vin*T
        assert F(x["prepared_energy_variance"])==2*q+b2
        assert F(x["prepared_mean_energy"])==F(q,2)+F(q*(M-S),10*M)
    groups.append("independent supplied 3D spectra and exact macroscopic response/energy fractions")

    for saved in pools["certificates"]:
        assert certificate(saved["N"],F(saved["threshold_long_density"]))==saved
        assert Decimal(saved["total_graph_record_energy_lower"])>0
    groups.append("outward ground coherence and fixed-calibration probability certificates")

    out=dict(status="CP12_COHERENCE_LOCALITY_CHECKS_PASS",group_count=len(groups),groups=groups,
             dependencies=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,networkx=nx.__version__),
             scope="Executed frozen-model controls, supplied-host locality and analytical certificates; no host formation, Lorentz invariance or gravity established.")
    (P/"verification.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))


if __name__=="__main__":main()
