#!/usr/bin/env python3
"""Read-only development checks. A successful command is not necessarily a clean gate."""
import argparse
import ast
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]


def command_status(name, returncode, output):
    if returncode != 0:
        return "FAIL"
    if name == "bank_gate":
        verdicts = re.findall(r"^OVERALL: (\S+)\s*$", output, re.MULTILINE)
        return "PASS" if verdicts == ["CLEAN"] else "REVIEW_REQUIRED"
    if name == "validate" and not re.search(r"^PASS:", output, re.MULTILINE):
        return "FAIL"
    if name == "validate_scoped" and ("PASS:" not in output or "FAIL" in output):
        return "FAIL"
    return "PASS"


def run_check(name, argv, cwd, outdir, env, timeout=900):
    started = time.monotonic()
    try:
        proc = subprocess.run(argv, cwd=cwd, env=env, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              timeout=timeout)
        output, rc = proc.stdout, proc.returncode
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ""
        if isinstance(output, bytes):
            output = output.decode("utf-8", errors="replace")
        output += "\nCHECK TIMEOUT\n"
        rc = 124
    except OSError as exc:
        output, rc = str(exc) + "\n", 127
    (outdir / (name + ".log")).write_text(output, encoding="utf-8")
    result = {"name": name, "command": argv, "returncode": rc,
              "status": command_status(name, rc, output),
              "seconds": round(time.monotonic() - started, 3),
              "log": name + ".log"}
    print(f"{name}: {result['status']} (exit {rc})", flush=True)
    return result


def syntax_check(outdir):
    paths = sorted((ROOT / "provenance").glob("*.py")) + sorted(
        (ROOT / "development").glob("*.py"))
    errors = []
    for path in paths:
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except (SyntaxError, UnicodeError) as exc:
            errors.append(f"{path.relative_to(ROOT)}: {exc}")
    output = "\n".join(errors) if errors else f"Parsed {len(paths)} Python files.\n"
    (outdir / "syntax.log").write_text(output, encoding="utf-8")
    result = {"name": "syntax", "status": "FAIL" if errors else "PASS",
              "returncode": 1 if errors else 0, "files": len(paths),
              "log": "syntax.log"}
    print(f"syntax: {result['status']}", flush=True)
    return result


def junit_counts(path):
    if not path.is_file():
        return None
    root = ET.parse(path).getroot()
    suites = list(root.iter("testsuite"))
    return {key: sum(int(s.get(key, 0)) for s in suites)
            for key in ("tests", "failures", "errors", "skipped")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "development/results")
    args = parser.parse_args()
    outdir = args.output.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    junit = outdir / "provenance.xml"
    junit.unlink(missing_ok=True)  # A failed new run must not inherit an old report.
    env = os.environ.copy()
    env.update({"PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1", "PYTEST_ADDOPTS": ""})
    for key in ("GRUT_FULL_MUTATION", "GRUT_RUN_SLOW"):
        env.pop(key, None)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True, check=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                           capture_output=True, text=True, check=True).stdout
    checks = [syntax_check(outdir)]
    for name, script in (("validate", "validate.py"),
                         ("validate_scoped", "validate_scoped.py"),
                         ("bank_gate", "bankgate.py")):
        checks.append(run_check(name, [sys.executable, "provenance/" + script],
                                ROOT, outdir, env))
    checks.append(run_check("runner_tests", [sys.executable, "-m", "unittest",
                            "discover", "-s", "development", "-p", "test_*.py"],
                            ROOT, outdir, env))
    checks.append(run_check("provenance_suite", [sys.executable, "-m", "pytest", "-q",
                            "--junitxml=" + str(junit)], ROOT / "provenance", outdir, env))
    try:
        counts = junit_counts(junit)
    except (ET.ParseError, ValueError):
        counts = None
    if counts is None or counts["tests"] == 0:
        checks[-1]["status"] = "FAIL"
        checks[-1]["report_error"] = "No valid nonempty JUnit report was produced."
    result = {"schema_version": 1, "repository": "ryangrvr/TestingGRUT",
              "source_commit": head, "working_tree_dirty": bool(dirty),
              "python": sys.version, "utc": datetime.now(timezone.utc).isoformat(),
              "profile": "default: slow falsifiers and slow mutants excluded",
              "checks": checks, "provenance_tests": counts,
              "status": "PASS" if all(c["status"] == "PASS" for c in checks) else "FAIL",
              "scientific_approval": False}
    (outdir / "summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Overall: {result['status']}; report: {outdir / 'summary.json'}", flush=True)
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
