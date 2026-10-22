"""Internal algebra/data/replay checks; not independent replication."""
import argparse
import ast
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from scipy.integrate import simpson

ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT.parent / "grut_unified_working_theory_20261022"
checks = []


def check(name, condition):
    checks.append((name, bool(condition)))
    if not condition:
        print("FAIL:", name)


def leaves(x, prefix=""):
    if isinstance(x, dict):
        return {p: v for k, a in sorted(x.items()) for p, v in leaves(a, prefix+"/"+k).items()}
    if isinstance(x, list):
        return {p: v for k, a in enumerate(x) for p, v in leaves(a, prefix+"/"+str(k)).items()}
    return {prefix: x}


def main(replay=False):
    result = json.loads((ROOT / "RESULTS.json").read_text())
    check("three profiles and twelve explicit variants", len(result["profiles"]) == 3 and len(result["examples"]) == 12)
    check("halo probability normalization", abs(result["halo"]["normalization_quad"]-1) < 1e-12)
    check("mean inverse speed direct quadrature", all(math.isclose(x["analytic"], x["quad"], rel_tol=1e-8, abs_tol=1e-16) for x in result["halo"]["eta_checks"]))
    check("contact recoil integral normalization", abs(result["contact_cross_section_recoil_integral_ratio"]-1) < 1e-12)
    for p in result["profiles"]:
        q = p["Q"]
        source = PARENT / "data" / f"profile_Q{q}.csv"
        check(f"profile {q}: frozen hash", hashlib.sha256(source.read_bytes()).hexdigest() == p["source_sha256"])
        check(f"profile {q}: depletion and zero-momentum normalization", abs(p["normalization_relative_residual"]) < 1e-9 and p["F_at_qhat"][0] == 1. and max(abs(x) for x in p["F_at_qhat"]) <= 1.)
        check(f"profile {q}: moment expansion", abs(p["small_q_moment_residual"]) < 1e-6)
        check(f"profile {q}: decimation", p["decimation_max_absolute_difference"] < 1e-7)
        check(f"profile {q}: alternative interpolation/quadrature", abs(p["qhat1_Pchip_quad_difference"]) < 1e-6)
        check(f"profile {q}: radial-source EOM", abs(p["radial_source_EOM_integral_relative"]) < 1e-9)
        check(f"profile {q}: radial-source correction identity", abs(p["radial_source_identity_relative"]) < 1e-9 and p["radial_source_extra_over_E_over_v"] > 0)
    for r in result["examples"]:
        name = f"Q{r['Q']} v{r['v_phi_GeV']:g}"
        check(name+": two Higgs-source routes", abs(r["source_integral_vs_mass_relative"]) < 1e-9 and abs(r["uniform_Higgs_derivative_vs_source_relative"]) < 1e-12)
        check(name+": sufficient Born control at declared coupling", r["Born_kernel_sufficient_upper_bound"] < .001)
        check(name+": positive rate and normalized suppression", r["ideal_events_tonne_year_5_40keV"] > 0 and 0 < r["extended_to_point_window_ratio"] <= 1.0000000001)
        check(name+": weak mixing and bounded quartic", abs(r["small_mixing_theta"]) < 1e-6 and 0 < r["quartic_stability_ratio"] < 1)
        check(name+": formal amplitude scaling", math.isclose(r["formal_kappa_for_one_ideal_event_per_tonne_year"]**2*r["ideal_events_tonne_year_5_40keV"], 1e-20, rel_tol=1e-12))
    for scale in (1., 10., 1000., 2.435e13):
        normalized = []
        for r in result["examples"]:
            if r["v_phi_GeV"] == scale:
                mu = .939*r["mass_GeV"]/(.939+r["mass_GeV"])
                normalized.append(r["sigma_n_zero_cm2"]/(mu*r["mass_GeV"])**2)
        check(f"v{scale:g}: common mass-normalized cross section", max(normalized)/min(normalized)-1 < 1e-12)
    raw = np.loadtxt(ROOT / "data/recoil_Q39944_v1.csv", delimiter=",", skiprows=1)
    chosen = next(x for x in result["examples"] if x["Q"] == 39944 and x["v_phi_GeV"] == 1.)
    check("raw spectrum integral", math.isclose(simpson(raw[:, 5], x=raw[:, 0]), chosen["ideal_events_tonne_year_5_40keV"], rel_tol=1e-12))
    check("raw spectrum point/object factorization", np.allclose(raw[:, 5], raw[:, 4]*raw[:, 2]**2, rtol=1e-13, atol=0))
    check("recoil-grid decimation", abs(simpson(raw[::2, 5], x=raw[::2, 0])/simpson(raw[:, 5], x=raw[:, 0])-1) < 1e-7)
    for script in (ROOT/"code").glob("*.py"):
        ast.parse(script.read_text())
    check("Python syntax", True)
    manifest = ROOT / "MANIFEST.json"
    if manifest.exists():
        entries = json.loads(manifest.read_text())["files"]
        actual = sorted(str(p.relative_to(ROOT)) for p in ROOT.rglob("*") if p.is_file() and p.name != "MANIFEST.json" and "__pycache__" not in p.parts)
        check("manifest exact file set", actual == sorted(entries))
        check("manifest hashes", all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == meta["sha256"] and (ROOT/p).stat().st_size == meta["bytes"] for p, meta in entries.items()))
    if replay:
        with tempfile.TemporaryDirectory(prefix="grut-portal-replay-") as tmp:
            here = Path(tmp)/ROOT.name
            there = Path(tmp)/PARENT.name
            (here/"code").mkdir(parents=True)
            (there/"data").mkdir(parents=True)
            shutil.copy2(ROOT/"code/calculate.py", here/"code/calculate.py")
            shutil.copy2(PARENT/"RESULTS.json", there/"RESULTS.json")
            for q in (19972, 39944, 79988):
                shutil.copy2(PARENT/"data"/f"profile_Q{q}.csv", there/"data"/f"profile_Q{q}.csv")
            subprocess.run([sys.executable, str(here/"code/calculate.py")], check=True, stdout=subprocess.DEVNULL)
            expected, fresh = leaves(result), leaves(json.loads((here/"RESULTS.json").read_text()))
            check("fresh replay schema", expected.keys() == fresh.keys())
            numeric = [p for p in expected if isinstance(expected[p], (int, float))]
            check(f"fresh replay {len(numeric)} numerical leaves", all(math.isclose(expected[p], fresh[p], rel_tol=1e-12, abs_tol=0) for p in numeric))
            check("fresh replay nonnumerical leaves", all(expected[p] == fresh[p] for p in expected if p not in numeric))
            check("fresh raw spectrum replay", (here/"data/recoil_Q39944_v1.csv").read_bytes() == (ROOT/"data/recoil_Q39944_v1.csv").read_bytes())
    failed = [name for name, okay in checks if not okay]
    print(json.dumps({"conditions": len(checks), "failed": failed,
                      "scope": "Same-author internal verification; not an observational constraint or outside replication."}))
    return bool(failed)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--replay", action="store_true")
    raise SystemExit(main(parser.parse_args().replay))
