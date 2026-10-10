#!/usr/bin/env python3
"""Two-axis engineering/scientific reporting for issue #2. No acceptance operations."""
import argparse
import ast
from datetime import datetime, timezone
from importlib.metadata import version
import json
import hashlib
from collections import Counter
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET

from freeze_manifest import PINNED_SOURCE
from integrity import (adjudication_problems, baseline_deltas, classifier_agrees,
    json_read, locked_manifest, outcome_returncode_agrees, profile_baseline,
    profiles_env, read_junit, source_identity, verify_assets)

ROOT = Path(__file__).resolve().parents[1]


def numerical_environment(requirements):
    observed = {}
    for line in Path(requirements).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if not re.fullmatch(r'[A-Za-z0-9_-]+==[0-9]+(?:\.[0-9]+)+', line):
            raise ValueError('Expensive-profile dependency must have an exact version pin')
        name, expected = line.split('==')
        installed = version(name)
        if installed != expected:
            raise ValueError(f'Expensive-profile dependency mismatch: {name} {installed} != {expected}')
        observed[name] = installed
    if not observed:
        raise ValueError('Missing numerical dependency contract')
    return observed


def command_status(name, returncode, output):
    if returncode != 0:
        return "FAIL"
    if name == "bank_gate":
        verdicts = re.findall(r"^OVERALL: (\S+)\s*$", output, re.MULTILINE)
        return "PASS" if verdicts == ["CLEAN"] else "REVIEW_REQUIRED"
    if name == "validate":
        return "PASS" if len(re.findall(r"^PASS:", output, re.MULTILINE)) == 1 and not re.search(
            r"^FAIL:", output, re.MULTILINE) else "FAIL"
    if name == "validate_scoped":
        scopes = re.findall(r"^SCOPE: (\S+)", output, re.MULTILINE)
        passes = re.findall(r"^\s*PASS:", output, re.MULTILINE)
        return "PASS" if scopes == ["grut", "vacuum-cluster"] and len(passes) == 2 and not re.search(
            r"^\s*FAIL", output, re.MULTILINE) else "FAIL"
    return "PASS"


def run_check(name, argv, cwd, outdir, env, timeout=900, fresh_reports=()):
    for path in fresh_reports:
        Path(path).unlink(missing_ok=True)
    started = time.monotonic()
    try:
        proc = subprocess.run(argv, cwd=cwd, env=env, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout)
        output, rc = proc.stdout, proc.returncode
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ""
        if isinstance(output, bytes):
            output = output.decode("utf-8", errors="replace")
        output, rc = output + "\nCHECK TIMEOUT\n", 124
    except OSError as exc:
        output, rc = str(exc) + "\n", 127
    (outdir / (name + ".log")).write_text(output, encoding="utf-8")
    result = {"name": name, "command": argv, "returncode": rc,
              "status": command_status(name, rc, output),
              "seconds": round(time.monotonic() - started, 3), "log": name + ".log"}
    print(f"{name}: {result['status']} (exit {rc})", flush=True)
    return result


def syntax_check(outdir):
    paths = sorted((ROOT / "provenance").glob("*.py")) + sorted((ROOT / "development").glob("*.py"))
    errors = []
    for path in paths:
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except (SyntaxError, UnicodeError) as exc:
            errors.append(f"{path.relative_to(ROOT)}: {exc}")
    output = "\n".join(errors) if errors else f"Parsed {len(paths)} Python files.\n"
    (outdir / "syntax.log").write_text(output, encoding="utf-8")
    result = {"name": "syntax", "status": "FAIL" if errors else "PASS",
              "returncode": 1 if errors else 0, "files": len(paths), "log": "syntax.log"}
    print(f"syntax: {result['status']}", flush=True)
    return result


def junit_counts(path):
    """Legacy diagnostic API. Decisions use the stricter, case-validating read_junit."""
    if not Path(path).is_file():
        return None
    root = ET.parse(path).getroot()
    return {key: sum(int(s.get(key, 0)) for s in root.iter("testsuite"))
            for key in ("tests", "failures", "errors", "skipped")}


def choose_axes(infrastructure_errors, deltas, problems, classifier_agreement, state, cases, profile):
    if infrastructure_errors or deltas or (profile == "default" and not classifier_agreement):
        engineering = "FAIL"
    elif profile != "default":
        engineering = "NOT_ESTABLISHED"
    elif problems:
        engineering = "BLOCKED"
    else:
        engineering = "PASS"
    bank = state.get("bank_inventory", {}).get("report", {}).get("overall", "UNKNOWN")
    open_passes = [pid for pid, item in state.get("open_passes", {}).items() if item["status"] == "OPEN"]
    scientific = "REVIEW_REQUIRED" if bank != "CLEAN" or open_passes or any(
        x != "PASS" for x in cases.values()) else "NO_RECORDED_OPEN_FLAGS"
    return {"engineering_integrity": engineering, "scientific_status": scientific,
            "profile_execution_integrity": "FAIL" if infrastructure_errors or deltas else "PASS",
            "bank_gate": bank, "open_passes": sorted(open_passes),
            "scientific_approval": False, "external_review": "NOT_PERFORMED_BY_THIS_RUN"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "development/results")
    parser.add_argument("--profile", choices=("default", "full-mutation", "slow-falsifiers"), default="default")
    parser.add_argument("--reconciliation", type=Path,
                        help="Explicit locked owner transition; original manifest is never refreshed")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    outdir = Path(tempfile.mkdtemp(prefix=args.profile + "-", dir=args.output.resolve()))
    env = profiles_env(os.environ, args.profile)
    checks, errors, deltas, problems, cases, counts, state = [], [], [], [], {}, None, {}
    agreement = False
    numerical_packages = {}
    identity = source_identity(ROOT)
    manifest_path = ROOT / "development/expected_red_manifest.json"
    manifest = None
    try:
        manifest = locked_manifest(manifest_path)
        if manifest["source_commit"] != PINNED_SOURCE:
            raise ValueError("Manifest source differs from issue #2's frozen source")
        if args.reconciliation:
            from owner_transition import reconciled_manifest
            manifest = reconciled_manifest(manifest, args.reconciliation, ROOT)
        if args.profile != 'default':
            numerical_packages = numerical_environment(ROOT / 'development/requirements-numerical.txt')
        checks.append(syntax_check(outdir))
        capture = ROOT / "development/capture_state.py"
        checks.append(run_check("state_before", [sys.executable, str(capture), "--root", str(ROOT)], ROOT, outdir, env))
        before = json_read(outdir / "state_before.log")
        if before["source_commit"] != identity["source_commit"]:
            raise ValueError("Captured source SHA does not match checkout")
        for name, script in (("validate", "validate.py"), ("validate_scoped", "validate_scoped.py"),
                             ("bank_gate", "bankgate.py")):
            checks.append(run_check(name, [sys.executable, "provenance/" + script], ROOT, outdir, env))
        checks.append(run_check("runner_tests", [sys.executable, "-m", "unittest", "discover",
                                "-s", "development", "-p", "test_*.py"], ROOT, outdir, env))
        junit = outdir / "provenance.xml"
        argv = [sys.executable, "-m", "pytest", "-q", "--junitxml=" + str(junit)]
        if args.profile == "full-mutation":
            argv.append("test_mutation_battery.py")
        raw = run_check("provenance_suite", argv, ROOT / "provenance", outdir, env,
                        timeout=900 if args.profile == "default" else 10800, fresh_reports=[junit])
        checks.append(raw)
        cases, counts = read_junit(junit)
        if not outcome_returncode_agrees(cases, raw["returncode"]):
            raise ValueError("Raw pytest exit code disagrees with case outcomes")
        if args.profile == "default":
            classifier = run_check("expected_red_classifier", [sys.executable, "provenance/expected_red.py"], ROOT, outdir, env)
            checks.append(classifier)
        checks.append(run_check("state_after", [sys.executable, str(capture), "--root", str(ROOT)], ROOT, outdir, env))
        state = json_read(outdir / "state_after.log")
        if state != before or source_identity(ROOT) != identity:
            raise ValueError("Inputs/source identity changed during execution")
        expected = dict(manifest)
        expected["pytest_cases"] = profile_baseline(manifest, args.profile)
        deltas = baseline_deltas(expected, cases, state)
        audit_cases = cases if args.profile == "default" else manifest["pytest_cases"]
        problems = adjudication_problems(audit_cases, state)
        if args.profile == "default":
            agreement = classifier_agrees(cases, problems,
                (outdir / "expected_red_classifier.log").read_text(), classifier["returncode"])
        verdicts = re.findall(r"^OVERALL: (\S+)\s*$", (outdir / "bank_gate.log").read_text(), re.MULTILINE)
        if verdicts != [state["bank_inventory"]["report"]["overall"]]:
            raise ValueError("Structured bank inventory and CLI verdict disagree")
        required = [c["log"] for c in checks] + ["provenance.xml"]
        assets = verify_assets(outdir, required)
        (outdir / "assets.json").write_text(json.dumps(assets, indent=2) + "\n")
    except Exception as exc:
        errors.append(f"{type(exc).__name__}: {exc}")
    critical = {"syntax", "state_before", "state_after", "validate", "validate_scoped", "runner_tests"}
    errors.extend(c["name"] + " did not pass" for c in checks if c["name"] in critical and c["status"] != "PASS")
    axes = choose_axes(errors, deltas, problems, agreement, state, cases, args.profile)
    result = {"schema_version": 2, "repository": "ryangrvr/TestingGRUT", **identity,
              "utc": datetime.now(timezone.utc).isoformat(), "profile": args.profile,
              "owner_reconciliation": str(args.reconciliation) if args.reconciliation else None,
              "baseline_source": manifest["source_commit"] if manifest else None,
              "numerical_packages": numerical_packages,
              "checks": checks, "provenance_tests": counts, "axes": axes,
              "baseline_deltas": deltas, "adjudication_problems": problems,
              "classifier_agreement_with_raw_and_independent_audit": agreement if args.profile == "default" else "NOT_RUN_SUBSET_OR_SLOW_PROFILE",
              "infrastructure_errors": errors, "case_outcomes": cases}
    (outdir / "summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    (outdir / "summary.md").write_text(
        f"Engineering integrity: **{axes['engineering_integrity']}**\n\n"
        f"Scientific status: **{axes['scientific_status']}**\n\n"
        f"Profile: `{args.profile}`; source `{identity['source_commit']}`; dirty tree: {identity['working_tree_dirty']}\n\n"
        f"Baseline deltas: {len(deltas)}; adjudication blockers: {len(problems)}; infrastructure errors: {len(errors)}.\n\n"
        "Observed failures are not new expected-red authorizations. No scientific approval is supplied.\n")
    try:
        verify_assets(outdir, ["summary.json", "summary.md", "assets.json"])
    except ValueError as exc:
        print(f"REPORT FAILURE: {exc}", flush=True)
        return 1
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a", encoding="utf-8") as target:
            target.write((outdir / "summary.md").read_text())
    github_output = os.environ.get("GITHUB_OUTPUT")
    if github_output:
        with open(github_output, "a", encoding="utf-8") as target:
            target.write(f"engineering_integrity={axes['engineering_integrity']}\n")
            target.write(f"scientific_status={axes['scientific_status']}\n")
    compact = {"source_commit": identity["source_commit"], "dirty": identity["working_tree_dirty"],
        "profile": args.profile, "axes": axes, "counts": counts,
        "baseline_delta_count": len(deltas), "classifier_agreement": agreement,
        "infrastructure_errors": errors, "blocker_counts": dict(Counter(p["kind"] for p in problems)),
        "case_digest": hashlib.sha256(json.dumps(cases, sort_keys=True).encode()).hexdigest(),
        "expected_case_digest": hashlib.sha256(json.dumps(profile_baseline(manifest, args.profile), sort_keys=True).encode()).hexdigest() if manifest else None}
    print("GRUT_RESULT_SUMMARY=" + json.dumps(compact, sort_keys=True), flush=True)
    print(f"Engineering: {axes['engineering_integrity']}; scientific: {axes['scientific_status']}; report: {outdir / 'summary.json'}", flush=True)
    if args.profile != "default":
        return 0 if axes["profile_execution_integrity"] == "PASS" else 1
    return 0 if axes["engineering_integrity"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
