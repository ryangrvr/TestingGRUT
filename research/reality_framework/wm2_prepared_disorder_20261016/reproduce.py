#!/usr/bin/env python3
"""Verify saved S4 evidence, or rerun with --recompute. No prior file is edited."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import platform
from pathlib import Path
import sys
import time

import numpy as np
import scipy

import controls as C

BASE = Path(__file__).resolve().parent
S3 = "80f5d6ab2e5d74307a9445cbbe9182ba2d4c4131"
CHARTER_COMMIT = "4c4d7a8d590b1f8587a571ece32bf0709ff06d9b"


def save(name, value):
    (BASE / name).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def run():
    result = dict(base_commit=S3, charter_commit=CHARTER_COMMIT,
                  executed_at_utc=datetime.now(timezone.utc).isoformat(),
                  runtime=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                               platform=platform.platform()),
                  arithmetic="IEEE-754 ordinary double precision; analytic cut errors exclude roundoff",
                  external_independent_validation=False)
    tasks = [("S3_Node", lambda: C.node_reproduction(BASE)),
             ("independent_matrix_construction", lambda: C.cross_runtime_matrices(BASE)),
             ("reduction", C.reduction_checks), ("first_error", C.first_error_checks),
             ("physical_receiver_formation", C.primary_formation),
             ("direct_receiver_cut_comparison", C.receiver_cut_comparison),
             ("small_mode_controls", C.small_mode_formation),
             ("independent_detuning_Q_packets", C.pole_packet_checks)]
    for key, function in tasks:
        start = time.monotonic()
        print("Running", key, flush=True)
        result[key] = function()
        print("Completed", key, "in", round(time.monotonic() - start, 3), "seconds", flush=True)
        save("CHECKPOINT.json", result)
    probes = []
    for n, delta, word in [(3, [0, -.25, .25], (0, 1, 2)),
                            (4, [0, -.25, .25, -.25], (0, 1, 2, 3)),
                            (5, [0, -.25, .25, -.25, .25], (0, 1, 2, 3, 4))]:
        for placement in ["interior", "end"]:
            print("Running conditional probe", n, placement, flush=True)
            probes.append(C.conditional_probe(n, 2, delta, word, placement))
    result["conditional_probes"] = probes
    save("RESULTS.json", result)
    (BASE / "CHECKPOINT.json").unlink(missing_ok=True)
    return result


def verify(r):
    assert r["base_commit"] == S3 and r["charter_commit"] == CHARTER_COMMIT
    freeze = json.loads((BASE / "CHARTER_FREEZE.json").read_text())
    body = freeze["raw_git_commit_object"].encode()
    assert hashlib.sha1(b"commit " + str(len(body)).encode() + b"\0" + body).hexdigest() == CHARTER_COMMIT
    charter = (BASE / "CHARTER.md").read_bytes()
    assert hashlib.sha1(b"blob " + str(len(charter)).encode() + b"\0" + charter).hexdigest() == freeze["charter_git_blob_sha"]
    assert r["S3_Node"]["complete_recomputation"] == "PASS"
    assert r["S3_Node"]["archived_packet_unchanged"]
    first = r["first_error"]
    assert sum(x["complete_consistent_sources_checked"] for x in first) == 195
    reduction = r["reduction"]["exact_type_word_checks"]
    assert len(reduction) == 12
    assert max(x["isometry_wavefunction_error"] for x in reduction) < 1e-10
    physical = r["physical_receiver_formation"]["cases"]
    observations = [x for a in physical for x in a["rows"]]
    assert len(physical) == 30 and len(observations) == 90
    assert max(x["unitary_norm_error"] for x in observations) < 2e-9
    assert min(x["macro_probability"] for x in observations) > .5
    assert max(x["physical_receiver_probability_error_upper"] for x in observations) < 1e-15
    assert len(r["physical_receiver_formation"]["all_binary_corner_root_averages"]) == 72
    assert len(r["small_mode_controls"]) == 20
    assert len(r["independent_matrix_construction"]) == 3
    assert len(r["direct_receiver_cut_comparison"]["rows"]) == 2
    assert len(r["independent_detuning_Q_packets"]) == 4
    probes = r["conditional_probes"]
    assert len(probes) == 6 and sum(len(x["rows"]) for x in probes) == 24
    pulse = max(x["pulse_identity_error"] for a in probes for x in a["rows"])
    assert pulse < 1e-10
    return dict(execution_consistency="PASS",
                physical_status="Conditional bulk dynamics proved; independent-disorder all-size formation OPEN",
                first_storage_error_sources_checked=195,
                exact_type_word_wavefunction_checks=12,
                physical_receiver_parameter_sets=len(physical),
                physical_receiver_observations=len(observations),
                binary_corner_root_averages=72,
                small_receiver_pattern_sets=20,
                conditional_probe_controls=6, phase_pulse_checks=24,
                independent_detuning_Q_packet_checks=4,
                largest_formation_matrix=max(x["exact_type_word_dimension"] for x in physical),
                largest_direct_physical_receiver_matrix=r["direct_receiver_cut_comparison"]["physical_dimension"],
                largest_full_one_probe_matrix=max(x["full_dimension"] for x in probes),
                independent_complete_matrix_entry_checks=sum(x["entries_compared"] for x in r["independent_matrix_construction"]),
                maximum_direct_receiver_cut_norm_discrepancy=max(x["directly_observed_state_difference"] for x in r["direct_receiver_cut_comparison"]["rows"]),
                maximum_type_word_wavefunction_error=max(x["isometry_wavefunction_error"] for x in reduction),
                maximum_formation_norm_error=max(x["unitary_norm_error"] for x in observations),
                maximum_receiver_cut_probability_error_upper=max(x["physical_receiver_probability_error_upper"] for x in observations),
                minimum_observed_primary_formation=min(x["macro_probability"] for x in observations),
                maximum_pulse_identity_error=pulse,
                independent_external_validation=False,
                all_N_nonshrinking_disorder_formation_proved=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args()
    r = run() if args.recompute else json.loads((BASE / "RESULTS.json").read_text())
    v = verify(r)
    files = sorted(p for p in BASE.iterdir() if p.suffix in (".py", ".md", ".json")
                   and p.name not in ("manifest.json", "VERIFICATION.json", "CHECKPOINT.json"))
    entries = {p.name: dict(sha256=hashlib.sha256(p.read_bytes()).hexdigest(), bytes=p.stat().st_size)
               for p in files}
    manifest = dict(research_base_commit=S3, frozen_charter_commit=CHARTER_COMMIT, files=entries)
    if args.recompute:
        save("manifest.json", manifest)
    else:
        assert manifest == json.loads((BASE / "manifest.json").read_text()), "Content changed"
    v["manifest_sha256"] = hashlib.sha256((BASE / "manifest.json").read_bytes()).hexdigest()
    if args.recompute:
        save("VERIFICATION.json", v)
    else:
        assert v == json.loads((BASE / "VERIFICATION.json").read_text())
    print(json.dumps(v, indent=2))


if __name__ == "__main__":
    main()
