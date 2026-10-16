#!/usr/bin/env python3
"""Verify the published packet or rerun it in a disposable copy."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import scipy

import campaign
import controls
import probability
import verify

BASE = Path(__file__).resolve().parent
CHARTERS = {
    "CHARTER.md": "58b94aa9d0c4b797aa961de7ac94f0317c2a3c66f4d2a2167736b131468e6fa6",
    "PROBABILITY_CHARTER.md": "499851a26aebebb17eaa6f1d14ea8023e7860fd0bed50c338564435ee0dbee42"
}


def save(path, value):
    path = Path(path)
    temporary = path.with_name(path.name+".partial")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")
    temporary.replace(path)


def check_charters():
    for name, digest in CHARTERS.items():
        assert verify.sha(BASE/name) == digest, "Frozen charter changed: "+name


def same_evidence(left, right):
    if isinstance(left, dict):
        return left.keys() == right.keys() and all(same_evidence(left[k],right[k]) for k in left)
    if isinstance(left, list):
        return len(left) == len(right) and all(same_evidence(a,b) for a,b in zip(left,right))
    if isinstance(left, (int,float)) and not isinstance(left,bool):
        return math.isclose(left,right,rel_tol=2e-12,abs_tol=2e-12)
    return left == right


def extra_checks():
    extra = dict(sign_symmetry=controls.sign_symmetry_checks(),
                 physical_receiver_reference=controls.physical_reference_check())
    maximum, count = 0.0, 0
    for a in range(101):
        p = a/100
        angle = math.acos(math.sqrt(p))
        for k in range(-50,51):
            rotation = k/100
            actual = math.cos(angle+rotation)**2
            error = 2*abs(math.sin(rotation/2))
            maximum = max(maximum, probability.lower_bound(p,error)-actual)
            count += 1
    assert maximum < 2e-12
    extra["probability_bound"] = dict(analytic_two_dimensional_controls=count, maximum_violation=maximum)
    return extra


def worker():
    check_charters()
    result = campaign.run()
    result["extra_checks"] = extra_checks()
    save(BASE/"RESULTS.json", result)
    extension = probability.run()
    audit = verify.audit()
    return dict(S5=result, P1=extension, audited=audit, verification=verify.summary(audit),
                runtime=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                             platform=platform.platform(), at_utc=datetime.now(timezone.utc).isoformat()),
                independent_validation=False, outward_rounded=False)


def packet_manifest():
    files = sorted(p for p in BASE.iterdir() if p.is_file() and p.suffix in (".md", ".py", ".json")
                   and p.name not in ("manifest.json", "VERIFICATION.json", "RERUN.json"))
    return dict(base_commit="c6ac0c68480d4164ad620bbf2a4b3d8c65ac85b1",
                S5_charter_commit="64378573e8251224dd563928667e961a0effb3da",
                P1_charter_commit="364781b841121e3169041b121e661048fd793524",
                files={p.name:dict(sha256=verify.sha(p),bytes=p.stat().st_size) for p in files})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--recompute", action="store_true", help="rerun all bounded controls in a temporary copy")
    parser.add_argument("--output", type=Path, help="optional JSON evidence destination for a recomputation")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    check_charters()
    if args.worker:
        save(BASE/"RERUN.json", worker())
        return
    if args.recompute:
        before = packet_manifest()
        with tempfile.TemporaryDirectory(prefix="grut_s5_rerun_") as directory:
            parent = Path(directory)/"reality_framework"
            parent.mkdir()
            target = parent/BASE.name
            ignore = shutil.ignore_patterns("__pycache__", "*.partial")
            shutil.copytree(BASE,target,ignore=ignore)
            prior = BASE.parent/"wm2_prepared_disorder_20261016"
            shutil.copytree(prior,parent/prior.name,ignore=ignore)
            environment = os.environ.copy()
            environment["OPENBLAS_NUM_THREADS"] = "1"
            subprocess.run([sys.executable,str(target/"reproduce.py"),"--worker"], check=True,env=environment)
            evidence = json.loads((target/"RERUN.json").read_text())
        assert packet_manifest() == before, "Archived S5 files changed during recomputation"
        if args.output:
            args.output.parent.mkdir(parents=True,exist_ok=True)
            save(args.output,evidence)
        print(json.dumps(evidence["verification"],indent=2))
        return
    audit = verify.audit()
    assert same_evidence(audit,json.loads((BASE/"AUDITED_RESULTS.json").read_text()))
    assert packet_manifest() == json.loads((BASE/"manifest.json").read_text())
    verification = verify.summary(audit)
    verification["manifest_sha256"] = verify.sha(BASE/"manifest.json")
    assert same_evidence(verification,json.loads((BASE/"VERIFICATION.json").read_text()))
    print(json.dumps(verification,indent=2))


if __name__ == "__main__":
    main()
