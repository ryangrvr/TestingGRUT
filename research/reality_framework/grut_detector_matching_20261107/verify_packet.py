"""Verify source provenance and reconstruct the archived decisions; no production rerun."""
import argparse
import hashlib
import json
import math
from pathlib import Path

from analyze import analyze as analyze_m0
from analyze_m1 import analyze as analyze_m1, distance

HERE=Path(__file__).resolve().parent


def verify():
    checks={}
    for label,cardname in [("m0_01","protocol_m0_01.json"),("m0_02","protocol_m0_02.json"),("m1_01","protocol_m1_01.json")]:
        directory=HERE/("raw_"+label)
        card=json.loads((HERE/cardname).read_text())
        actual=(analyze_m1 if label.startswith("m1") else analyze_m0)(directory,card)
        expected=json.loads((directory/"analysis.json").read_text())
        checks[label+"_analysis_exact"]=actual==expected
        provenance=json.loads((directory/"provenance.json").read_text())
        hashes=provenance["source_sha256"]
        checks[label+"_frozen_sources"]=all(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==value for name,value in hashes.items())
        checks[label+"_miss_preserved"]="miss" in actual["verdict"]
    d=json.loads((HERE/"raw_d1_01/result.json").read_text())
    checks["d1_frozen_sources"]=all(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==value for name,value in d["source_sha256"].items())
    checks["d1_distances_and_bounds"]=True
    for record in d["records"].values():
        for reduction in record["reduced"].values():
            p=reduction["probabilities"]
            measured=distance(p,record["parent"]["probabilities"])
            bound=min(1.,sum(2*math.sqrt(p[s])*v+v*v for s,v in reduction["branch_vector_bounds"].items()))
            checks["d1_distances_and_bounds"] &= (abs(measured-reduction["tv_from_parent"])<1e-12 and
                abs(bound-reduction["tv_bound"])<1e-12 and measured<=bound+1e-10)
    checks["d1_miss_preserved"]=d["verdict"]=="finite_parent_reduction_gate_miss" and not d["gates"]["K32_accuracy"] and not d["gates"]["K32_bound"]
    checks["probabilities_normalized"]=True
    count=0
    def walk(value):
        nonlocal count
        if isinstance(value,dict):
            if value and all(k=="overflow" or (k and set(k)<=set("+-")) for k in value) and all(isinstance(v,(int,float)) for v in value.values()):
                count+=1
                checks["probabilities_normalized"] &= min(value.values())>=-1e-10 and abs(sum(value.values())-1)<1e-9
            else:
                for item in value.values():walk(item)
        elif isinstance(value,list):
            for item in value:walk(item)
    for directory in ["raw_m0_01","raw_m0_02","raw_m1_01"]:
        for p in (HERE/directory).glob("*.json"):
            value=json.loads(p.read_text())
            if "histories" in value:walk(value["histories"])
    for record in d["records"].values():
        walk(record["parent"]["probabilities"]);walk(record["dephased"])
        for value in record["reduced"].values():walk(value["probabilities"])
    return {"checks":checks,"probability_distributions_checked":count,
            "verification_pass":all(checks.values()),
            "scientific_status":"Matching and reduction misses remain; this verifies preservation and analysis, not physical acceptance."}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=verify()
    if args.output:
        with args.output.open("x") as f:json.dump(result,f,indent=2);f.write("\n")
    print(json.dumps(result,indent=2))
    raise SystemExit(0 if result["verification_pass"] else 1)
