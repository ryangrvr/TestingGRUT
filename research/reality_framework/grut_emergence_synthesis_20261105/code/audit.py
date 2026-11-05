"""Source integrity, saved-result reaggregation and package consistency.

This does not rerun any archived physical solver and is not independent review.
"""
import argparse
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
checks=[]


def record(name,passed,detail=None):
    checks.append(dict(name=name,pass_=bool(passed),detail=detail))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-root",type=Path,
                    default=ROOT.parent/"grut_source_snapshot")
    args=ap.parse_args()
    m=json.loads((ROOT/"SOURCE_MANIFEST.json").read_text())
    sources={r["id"]:r for r in m["local"]}
    for r in m["local"]:
        b=(args.source_root/r["path"]).read_bytes()
        sha=hashlib.sha256(b).hexdigest()
        oid=hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
        record(r["id"]+" exact source",sha==r["sha256"] and oid==r["git_blob"])
    w=json.loads((ROOT/"results/witnesses.json").read_text())
    record("witness source identity",
           hashlib.sha256((ROOT/"code/witnesses.py").read_bytes()).hexdigest()==w["source_sha256"])
    record("all finite checks pass",w["all_checks_pass"] and all(r["pass"] for r in w["checks"]))
    record("finite check count",w["check_count"]==len(w["checks"])==48)
    d=json.loads((args.source_root/sources["R09"]["path"]).read_text())
    bracketed={k:r["primary"]["free X"] for k,r in d["roots"].items()
               if r["primary"]["free X"]["self_consistent"] and not r["primary"]["free X"]["extrapolated"]}
    numbers={k:r["e_crit"] for k,r in bracketed.items() if isinstance(r["e_crit"],(int,float))}
    no_cross=sum(r["e_crit"] is None for r in bracketed.values())
    single=sum(isinstance(r["e_crit"],str) for r in bracketed.values())
    band=d["ecrit_band"]
    strict=[k for k,v in numbers.items() if band[0]<=v<=band[1]]
    qkeys=[k for k,r in d["roots"].items() if r["strict_M2_Q"]]
    diagnostics=d["computed_pair_extrapolation_diagnostic_only"]
    record("135 keys at each simulated node",
           len(d["roots"])==135 and all(n["grid_points"]==135 for n in d["nodes"].values()))
    record("37 bracketed self-consistent free-X roots",len(bracketed)==37)
    record("32 absent, 3 single-node, 2 numeric thresholds",(no_cross,single,len(numbers))==(32,3,2))
    record("registered band unchanged",band==[1.03e-6,4.12e-6])
    record("zero strict M2-X matches archive",strict==d["strict_M2_X_keys"]==[])
    record("two isolated Q roots",qkeys==d["strict_M2_Q_keys"]==["1.380:+0","1.395:+2"])
    record("15 diagnostic extrapolations only",len(diagnostics)==15 and all(
        d["roots"][k]["primary"]["free X"]["extrapolated"] for k in diagnostics))
    threshold_ratio=min(numbers.values())/band[1]
    record("lower threshold 5.30 times upper bound",abs(threshold_ratio-5.299898135292915)<1e-8)
    low=json.loads((args.source_root/sources["R10"]["path"]).read_text())
    fr=low["comparison"]["f_rad"]
    f55=fr["H55_fresh_scalar_lattice"]["3.0"]
    f84=fr["H84_reused_lattice"]["3.0"]
    record("hand-over radiation fractions",abs(f55-.6695852)<1e-12 and abs(f84-.8321434)<1e-12)
    record("H55 not promoted to new kinetic/Dirac result",not any(low[k] for k in [
        "full_3D_Dirac_evolved","new_kinetic_splice_computed","used_in_registered_verdict"]))
    # Independent polynomial assembly of dV/dphi in units lambda_phi*v^3.
    from numpy.polynomial import Polynomial as P
    s=P([0,1])
    direct=s*(s*s-1)-3*(.002/.03)*(s-1)**2
    factored=(s-1)*(s*s+.8*s+.2)
    record("tree potential derivative factorization",np.max(np.abs(direct.coef-factored.coef))<1e-14)
    record("quadratic has no real stationary root",abs(.8**2-4*.2+.16)<1e-14)
    claims=json.loads((ROOT/"CLAIMS.json").read_text())
    cid={c["id"] for c in claims}
    sid={r["id"] for r in m["local"]+m["external"]}|{"W01"}
    allowed={"ADOPTED","STANDARD","CONDITIONAL","FAILED-IN-SCOPE","OPEN","PROPOSED"}
    record("unique claim identifiers",len(cid)==len(claims))
    record("typed sourced scoped claims",all(c["status"] in allowed and c["premises"] and c["challenge"]
           and set(c["sources"])<=sid for c in claims))
    record("no independent-review promotion",all(not c["independent_review_in_this_packet"] for c in claims))
    graph={c["id"]:c["depends_on"] for c in claims}
    def visit(k,path):
        if k in path or k not in graph: return False
        return all(visit(j,path|{k}) for j in graph[k])
    record("claim dependencies valid and acyclic",all(visit(k,set()) for k in graph))
    em=json.loads((ROOT/"EMERGENCE_MAP.json").read_text())
    nodes={n["id"] for n in em["nodes"]}
    record("map references valid",all(set(n["sources"])<=sid for n in em["nodes"])
           and all(e["source"] in nodes and e["target"] in nodes for e in em["edges"]))
    ET.parse(ROOT/"EMERGENCE_MAP.svg")
    record("map SVG parses",True)
    # Make the report target real before checking its incoming links. Never expose
    # a stale passing report if the remaining checks fail or this run is interrupted.
    (ROOT/"results/audit.json").write_text(json.dumps({"status":"running","all_checks_pass":False})+"\n")
    missing=[]
    for f in ROOT.glob("*.md"):
        for target in re.findall(r"\]\(([^)]+)\)",f.read_text()):
            if not re.match(r"https?://|#",target):
                if not (f.parent/target.split("#")[0]).exists():missing.append([f.name,target])
    record("local Markdown links resolve",not missing,missing)
    txt=(ROOT/"GRUT_WORKING_THEORY.md").read_text()
    refs=set(re.findall(r"\b(?:R|S|W)\d\d\b",txt))
    record("manuscript citation identifiers known",refs<=sid,sorted(refs-sid))
    out={"scope":"same-author audit; no physical solver rerun; no independent review",
         "base_commit":m["base_commit"],"check_count":len(checks),
         "all_checks_pass":all(c["pass_"] for c in checks),
         "checks":checks,
         "recomputed":{"strict_M2_X":len(strict),"bracketed_roots":len(bracketed),
           "numeric_thresholds":numbers,"minimum_threshold_to_band_ceiling":threshold_ratio,
           "Q_roots":qkeys,"diagnostic_extrapolations":len(diagnostics),
           "H55_minus_H84_radiation_percentage_points":100*(f55-f84)}}
    (ROOT/"results/audit.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:out[k] for k in ["check_count","all_checks_pass","recomputed"]},indent=2))
    if not out["all_checks_pass"]:
        print(json.dumps([c for c in checks if not c["pass_"]],indent=2));raise SystemExit(1)


if __name__=="__main__": main()
