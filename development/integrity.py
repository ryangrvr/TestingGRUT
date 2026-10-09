"""Development-only observation and comparison; no adjudication authority."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET

ENGINEERING_AREAS = ("development/", ".github/workflows/")


def json_read(path):
    def unique(pairs):
        out = {}
        for key, value in pairs:
            if key in out:
                raise ValueError(f"Duplicate JSON key: {key}")
            out[key] = value
        return out
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=unique)


def locked_manifest(path):
    path = Path(path)
    lock = path.with_suffix(path.suffix + ".sha256").read_text().strip()
    if not re.fullmatch("[0-9a-f]{64}", lock) or hashlib.sha256(path.read_bytes()).hexdigest() != lock:
        raise ValueError("Frozen observational manifest failed its hash lock")
    manifest = json_read(path)
    if manifest.get("authorization") != "OBSERVATIONS_ONLY_NO_ADJUDICATION":
        raise ValueError("Manifest lacks the observational-only scope")
    return manifest


def verify_assets(directory, required):
    root = Path(directory).resolve()
    hashes = {}
    for name in required:
        path = (root / name).resolve()
        if path.parent != root or not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing/nonlocal/empty report asset: {name}")
        hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return hashes


def source_identity(root):
    root = Path(root)
    def git(*args):
        return subprocess.check_output(["git", *args], cwd=root).decode().strip()
    paths = git("ls-files", "-z").split("\0")
    hashes = {p: hashlib.sha256((root / p).read_bytes()).hexdigest()
              for p in sorted(paths) if p and p != "AGENT_COORDINATION.md"
              and not p.startswith(ENGINEERING_AREAS)}
    untracked = [p for p in git("ls-files", "--others", "--exclude-standard", "-z").split("\0") if p]
    unknown = [p for p in untracked if p != "AGENT_COORDINATION.md" and not p.startswith(ENGINEERING_AREAS)]
    engineering = {p: hashlib.sha256((root / p).read_bytes()).hexdigest()
                   for p in sorted(set(paths + untracked)) if p and p.startswith(ENGINEERING_AREAS)}
    return {"source_commit": git("rev-parse", "HEAD"),
            "working_tree_dirty": bool(git("status", "--porcelain")),
            "protected_inputs": {"algorithm": "sha256-of-sorted-path-sha256-map-v1",
                "digest": hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest(),
                "tracked_file_count": len(hashes)},
            "engineering_inputs": {"digest": hashlib.sha256(json.dumps(engineering, sort_keys=True).encode()).hexdigest(),
                                   "file_count": len(engineering)},
            "untracked_protected_paths": sorted(unknown)}


def read_junit(path):
    root = ET.parse(path).getroot()
    if root.tag not in ("testsuite", "testsuites"):
        raise ValueError("Invalid JUnit root")
    suites = [root] if root.tag == "testsuite" else list(root)
    if not suites or any(s.tag != "testsuite" for s in suites):
        raise ValueError("No valid test suites")
    cases = {}
    for suite in suites:
        if list(suite.iter("testsuite")) != [suite]:
            raise ValueError("Nested suites cannot be counted unambiguously")
        for item in suite.findall("testcase"):
            classname, name = item.get("classname", ""), item.get("name", "")
            if not re.fullmatch(r"test_[A-Za-z0-9_]+(?:\.[A-Za-z0-9_]+)*", classname) or not name:
                raise ValueError("Missing or unsupported exact pytest node identity")
            module, *classes = classname.split(".")
            node = module + ".py::" + "::".join(classes + [name])
            if node in cases:
                raise ValueError(f"Duplicate pytest node: {node}")
            marks = [tag for tag in ("failure", "error", "skipped") if item.find(tag) is not None]
            if len(marks) > 1:
                raise ValueError(f"Conflicting outcomes for {node}")
            cases[node] = {"failure": "FAIL", "error": "ERROR", "skipped": "SKIP"}.get(
                marks[0] if marks else "", "PASS")
    counts = Counter(cases.values())
    want = {"tests": len(cases), "failures": counts["FAIL"],
            "errors": counts["ERROR"], "skipped": counts["SKIP"]}
    got = {key: sum(int(s.attrib[key]) for s in suites) for key in want}
    if not cases or got != want:
        raise ValueError(f"JUnit counts do not certify the case list: {got} vs {want}")
    return dict(sorted(cases.items())), want


def outcome_returncode_agrees(cases, returncode):
    expected = 1 if any(x in ("FAIL", "ERROR") for x in cases.values()) else 0
    return returncode == expected


def baseline_deltas(manifest, cases, state):
    deltas = []
    old = manifest["pytest_cases"]
    for node in sorted(set(old) | set(cases)):
        if old.get(node) != cases.get(node):
            deltas.append({"kind": "PYTEST_CASE", "node": node,
                           "baseline": old.get(node), "observed": cases.get(node)})
    for key in ("declarations", "open_passes", "bank_inventory", "protected_inputs", "tool_versions"):
        if manifest[key] != state[key]:
            deltas.append({"kind": key.upper(), "baseline": manifest[key], "observed": state[key]})
    if state.get("untracked_protected_paths"):
        deltas.append({"kind": "UNTRACKED_PROTECTED", "paths": state["untracked_protected_paths"]})
    return deltas


def adjudication_problems(cases, state):
    """Independently mirror the existing classifier's obligations; never broaden them."""
    declared = state["declarations"]
    failing = {n for n, value in cases.items() if value in ("FAIL", "ERROR")}
    problems = []
    for node in sorted(failing - set(declared)):
        problems.append({"kind": "UNDECLARED_FAILING_TEST", "node": node})
    for node in sorted(set(declared) - failing):
        problems.append({"kind": "STALE_DECLARATION", "node": node})
    cited = set()
    for node, entry in sorted(declared.items()):
        live, want = set(entry["live_cases"]), set(entry["cases"])
        if node in failing:
            if not live:
                problems.append({"kind": "UNMODELLED_FAILURE", "node": node})
            for case in sorted(live - want):
                problems.append({"kind": "NEW_DECLARED_SET_MEMBER", "node": node, "case": case})
            for case in sorted(want - live):
                problems.append({"kind": "STALE_DECLARED_CASE", "node": node, "case": case})
        for case, pid in sorted(entry["cases"].items()):
            cited.add(pid)
            if state["open_passes"].get(pid, {}).get("status") != "OPEN":
                problems.append({"kind": "NON_OPEN_OR_UNKNOWN_ADJUDICATION", "pass": pid, "case": case})
    for pid, entry in sorted(state["open_passes"].items()):
        if entry["status"] == "OPEN" and not entry["symptomless"] and pid not in cited:
            problems.append({"kind": "ORPHANED_OPEN_PASS", "pass": pid})
    return problems


def classifier_observed_failures(output):
    return set(re.findall(r"^\s*EXPECTED-RED\s+(\S+)\s*$", output, re.MULTILINE)) | set(
        re.findall(r"^\s*\*\*\* NEW RED: (\S+)\s*$", output, re.MULTILINE))


def classifier_agrees(cases, problems, output, returncode):
    failures = {node for node, status in cases.items() if status in ("FAIL", "ERROR")}
    verdict = "FAIL: the failing set is not the declared set." if problems else "No new red."
    diagnostics = {
        "UNDECLARED_FAILING_TEST": {(n,) for n in re.findall(r"^\s*\*\*\* NEW RED: (\S+)", output, re.MULTILINE)},
        "STALE_DECLARATION": {(n,) for n in re.findall(r"^\s*\*\*\* STALE DECLARATION: (\S+)", output, re.MULTILINE)},
        "ORPHANED_OPEN_PASS": {(p,) for p in re.findall(r"^\s*\*\*\* ORPHANED OPEN PASS '([^']+)'", output, re.MULTILINE)},
        "NEW_DECLARED_SET_MEMBER": set(re.findall(r"^\s*\*\*\* NEW RED \(case\): (\S+)\n[ \t]+([^\n]+)", output, re.MULTILINE)),
        "STALE_DECLARED_CASE": set(re.findall(r"^\s*\*\*\* STALE CASE: (\S+)\n[ \t]+([^\n]+)", output, re.MULTILINE)),
        "UNMODELLED_FAILURE": {(n,) for n in re.findall(r"^\s*\*\*\* UNMODELLED FAILURE: (\S+)", output, re.MULTILINE)},
    }
    for kind, observed in diagnostics.items():
        expected = {(p["node"], p["case"]) if "case" in p else (p.get("node", p.get("pass")),)
                    for p in problems if p["kind"] == kind}
        if observed != expected:
            return False
    invalid_passes = {p["pass"] for p in problems if p["kind"] == "NON_OPEN_OR_UNKNOWN_ADJUDICATION"}
    if set(re.findall(r"^\s*\*\*\* (?:UNKNOWN|NON-OPEN) PASS '([^']+)'", output, re.MULTILINE)) != invalid_passes:
        return False
    return (classifier_observed_failures(output) == failures
            and returncode == (1 if problems else 0) and verdict in output)


def profiles_env(base, profile):
    env = dict(base)
    env.update({"PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1", "PYTEST_ADDOPTS": ""})
    env.pop("GRUT_FULL_MUTATION", None)
    env.pop("GRUT_RUN_SLOW", None)
    if profile == "full-mutation":
        env["GRUT_FULL_MUTATION"] = "1"
    elif profile == "slow-falsifiers":
        env["GRUT_RUN_SLOW"] = "1"
    elif profile != "default":
        raise ValueError("Unknown execution profile")
    return env


def profile_baseline(manifest, profile):
    out = dict(manifest["pytest_cases"])
    if profile == "full-mutation":
        return {n: v for n, v in out.items() if n.startswith("test_mutation_battery.py::")}
    if profile == "slow-falsifiers":
        node = "test_layer5_overturning.py::TestLayer5::test_every_cited_falsifier_exits_zero"
        if out.get(node) != "SKIP":
            raise ValueError("Frozen slow-falsifier skip is absent")
        out[node] = "PASS"  # Explicit required outcome when the documented guard is enabled.
    return out
