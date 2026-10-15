"""Reproduce this packet without overwriting its saved scientific outputs.

Default: validate saved-result assertions, parent reproduction, sources and links.
--rerun: also execute every new numerical/exact calculation in a temporary mirror
and compare science outputs. --baseline additionally repeats the earlier 1801afe
Node/Python suite, which was already executed without caches during this sprint.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse
import hashlib
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np
import scipy

P = Path(__file__).resolve().parent
ENV = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
ENV.pop("CP12_FORMATION_CACHE", None)


def execute(base, script):
    start = time.monotonic()
    result = subprocess.run([sys.executable, str(base / script)], cwd=base,
                            env=ENV, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"{script}: {result.stderr[-6000:]}\n{result.stdout[-2000:]}")
    return dict(script=script, exit_code=0, seconds=time.monotonic()-start,
                stdout_sha256=hashlib.sha256(result.stdout.encode()).hexdigest(),
                stderr=result.stderr)


def compare(a, b, path="", differences=None):
    if differences is None:
        differences = []
    if isinstance(a, dict) and isinstance(b, dict):
        assert a.keys() == b.keys(), (path, "keys")
        for key in a:
            if key not in {"seconds", "numeric_runtime_seconds", "platform", "python"}:
                compare(a[key], b[key], path+"."+key, differences)
    elif isinstance(a, list) and isinstance(b, list):
        assert len(a) == len(b), (path, "length")
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, path+f"[{i}]", differences)
    elif isinstance(a, float) or isinstance(b, float):
        assert math.isclose(a, b, rel_tol=2e-8, abs_tol=2e-12), (path, a, b)
        if a != b:
            differences.append(dict(path=path, saved=a, fresh=b))
    else:
        assert a == b, (path, a, b)
    return differences


def check_saved():
    charged = json.loads((P/"charged_access/RESULTS.json").read_text())
    assert charged["status"] == "EXACT_GRAPH_CHOICE_ACCESS_PASS"
    assert len(charged["independent_DP_checks"]) == 16
    assert {x["N"] for x in charged["hosts"]} == {5, 6, 7}
    gauge = json.loads((P/"induced_gauge/RESULTS.json").read_text())
    assert "PASS" in gauge["status"]
    for relative in ["protected_geometry/RESULTS.json", "protected_geometry/EXACT_ONSET.json",
                     "protected_geometry/NEUTRAL_RECORD.json", "independent_review/ALGEBRA_CHECK.json",
                     "theory_context/LOCALITY_BRIDGE_CHECK.json", "PREPARATION.json"]:
        data = json.loads((P/relative).read_text())
        assert "PASS" in data["status"], (relative, data["status"])
    expected = {
        "cp12_ground_connectivity/frozen_model.py": "e896fce8a3e61948ff41891576b4c3ac028ad26785cafabb30bddef93c5adb9a",
        "cp12_ground_connectivity/exact_destination.py": "1a1ccb4f06602179f14c2834dea3566fada022ca6a23f58370abe0d4f1e5b35b",
        "cp12_distance_history/evolution.py": "672233f2d9a934af047cfb3a9c5a24a2eaf35469584a4c761e43febddc3bc578",
    }
    for relative, digest in expected.items():
        assert hashlib.sha256((P.parent/relative).read_bytes()).hexdigest() == digest
    broken = []
    for path in P.rglob("*.md"):
        if "unpublished_local_backup" in path.parts:
            continue
        for target in re.findall(r"\]\(([^)\n]+)\)", path.read_text()):
            if not target.startswith(("https:", "http:", "#", "mailto:")):
                local = target.split("#")[0]
                if local and not (path.parent/local).exists():
                    broken.append(str(path.relative_to(P))+": "+target)
    # VERIFICATION.json is this run's output and can be absent on first execution.
    broken = [x for x in broken if not x.endswith(": VERIFICATION.json")]
    assert not broken, broken
    return expected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rerun", action="store_true")
    parser.add_argument("--baseline", action="store_true")
    args = parser.parse_args()
    started = time.monotonic()
    sources = check_saved()
    saved_checks = [execute(P, "protected_geometry/verify.py"),
                    execute(P, "reproduction/compare_results.py")]
    fresh, comparisons = [], {}
    if args.rerun:
        with tempfile.TemporaryDirectory(prefix="grut-sprint-reproduction-") as temporary:
            rf = Path(temporary)/"reality_framework"
            rf.mkdir()
            for sibling in P.parent.iterdir():
                if sibling.name != P.name:
                    (rf/sibling.name).symlink_to(sibling, target_is_directory=sibling.is_dir())
            mirror = rf/P.name
            shutil.copytree(P, mirror, ignore=shutil.ignore_patterns("__pycache__", "unpublished_local_backup"))
            scripts = ["charged_access/run.py", "induced_gauge/run.py",
                       "protected_geometry/run.py", "protected_geometry/exact_onset.py",
                       "protected_geometry/neutral_record.py", "theory_context/locality_bridge_check.py",
                       "independent_review/algebra_check.py", "preparation.py"]
            # Keep the electric evolution and exact combinatorics sequential:
            # their temporary workspaces can otherwise overlap substantially.
            for script in scripts:
                print("Reproducing "+script, flush=True)
                fresh.append(execute(mirror, script))
            for relative in ["charged_access/RESULTS.json", "induced_gauge/RESULTS.json",
                             "protected_geometry/RESULTS.json", "protected_geometry/EXACT_ONSET.json",
                             "protected_geometry/NEUTRAL_RECORD.json",
                             "theory_context/LOCALITY_BRIDGE_CHECK.json",
                             "independent_review/ALGEBRA_CHECK.json", "PREPARATION.json"]:
                differences = compare(json.loads((P/relative).read_text()),
                                      json.loads((mirror/relative).read_text()))
                comparisons[relative] = dict(status="PASS", nonbitwise_float_entries=len(differences),
                                              comparison_rtol=2e-8, comparison_atol=2e-12)
            fresh.append(execute(mirror, "protected_geometry/verify.py"))
            if args.baseline:
                # This earlier suite writes only the mirror's reproduction logs.
                fresh.append(execute(mirror, "reproduction/run_reproduction.py"))
                fresh.append(execute(mirror, "reproduction/compare_results.py"))
    result = dict(status="PASS", baseline="a2f447dc37ea0a71f62c6e990a4630d727ebd85e",
        parent="1801afee27096c03c2a3b64cc811ac0be0d2a4a7", law_changed=False,
        new_candidate_hamiltonian=False, cache_used=False,
        fresh_new_calculation_reproduction=args.rerun, fresh_parent_reproduction_this_invocation=args.baseline,
        parent_reproduction_completed_during_sprint=True,
        saved_assertions=saved_checks, fresh_executions=fresh, scientific_output_comparisons=comparisons,
        frozen_source_sha256=sources,
        environment=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform()),
        seconds=time.monotonic()-started,
        limitations="Ordinary floating-point evaluation of analytical bounds; no formal outward-rounding or external certification. See individual derivations for approximation scopes.")
    (P/"VERIFICATION.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(dict(status="PASS", fresh_scripts=len(fresh), comparisons=len(comparisons), seconds=result["seconds"])))


if __name__ == "__main__":
    main()
