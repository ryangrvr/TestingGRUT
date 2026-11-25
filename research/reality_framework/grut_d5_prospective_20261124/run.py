"""Evaluate the frozen D5 comparison arms and exactly grouped pulse scores."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from model import (
    BLANK, MODES, archived_flip, prefixes, pulse_result, route, seed, target,
    window_certificate,
)

ROOT = Path(__file__).resolve().parent
FREEZE_COMMIT = "9767724c9ee8302df9c29c8545744728989db136"


def display(word: bytes) -> str:
    return "".join("01e"[symbol] for symbol in word)


def occupation(word: bytes) -> int:
    return len(word) - word.count(BLANK)


def snapshots(states, edges, size, length):
    by_link = {edge.link: state for edge, state in zip(edges, states)}
    last_rebuild = max(index for index, edge in enumerate(edges) if edge.operation == "rebuild")
    first, second = by_link[length // 4], by_link[length // 2]
    assert occupation(first.a) == size // 4
    assert occupation(second.a) == size // 2
    assert states[0].memory == BLANK and first.memory == second.memory == states[0].h
    assert all(occupation(state.a) < size for state in states[:last_rebuild + 1])
    assert occupation(states[last_rebuild + 1].a) == size
    return {
        "first_early_tape": display(first.a),
        "second_early_tape": display(second.a),
        "filled_before_first_full_write": occupation(states[last_rebuild].a),
        "first_full_tape": display(states[last_rebuild + 1].a),
        "first_full_clock_position": edges[last_rebuild].link + 1,
        "memory_before_assay": states[0].memory,
        "memory_at_first_early_stage": first.memory,
        "memory_at_second_early_stage": second.memory,
    }


def run() -> dict:
    protocol_bytes = (ROOT / "PROTOCOL.json").read_bytes()
    protocol = json.loads(protocol_bytes)
    gates = protocol["thresholds"]
    histories = list(product((0, 1), repeat=3))
    summary = []
    for size in protocol["sizes"]:
        length = protocol["head_edges_per_cell"] * size
        edges = route(size, length)
        central_index = next(index for index, edge in enumerate(edges) if edge.operation == "control")
        exact_arms, controls, pulses = [], [], []
        for h, g, coin in histories:
            word = target(size, h, g)
            for fault in range(4):
                initial = seed(size, h, g, coin, fault)
                states = prefixes(initial, edges)
                endpoint = states[-1]
                witness = fault == 0 or archived_flip(endpoint, word, fault)
                assert endpoint.a == word and witness
                assert states == prefixes(initial, edges, null=True)
                exact_arms.append({
                    "h": h, "gene": g, "coin": coin, "fault": fault,
                    "score_head_and_target": 1,
                    "damage_archive_matches": witness if fault else None,
                    "displayed_word": display(endpoint.a),
                    "old_prefix_archive": display(endpoint.b),
                    "maintenance_archive": display(endpoint.c),
                    "observationally_identical_null": True,
                    "path_snapshots": snapshots(states, edges, size, length),
                })

            damaged = prefixes(seed(size, h, g, coin, 2), edges)[central_index + 1]
            for mode in MODES:
                initial = seed(size, h, g, coin, 2, mode)
                states = prefixes(initial, edges)
                checkpoint = states[central_index + 1]
                assert (checkpoint.a, checkpoint.b, checkpoint.c, checkpoint.phase) == (
                    damaged.a, damaged.b, damaged.c, damaged.phase)
                endpoint = states[-1]
                null_end = prefixes(initial, edges, null=True)[-1]
                controls.append({
                    "h": h, "gene": g, "coin": coin, "mode": mode,
                    "score_head_and_original_target": int(endpoint.a == word),
                    "null_score_head_and_original_target": int(null_end.a == word),
                    "damaged_tape_at_control": display(checkpoint.a),
                    "target_after_control": display(endpoint.a),
                    "null_target_after_control": display(null_end.a),
                    "gene_after_control": checkpoint.gene,
                    "memory_after_control": checkpoint.memory,
                    "gene_backup": checkpoint.gene_backup,
                    "memory_backup": checkpoint.memory_backup,
                    "same_damaged_substrate_at_control": True,
                    "clock_energy_distribution": "Binomial(L,1/2) over omega*(j-L/2)",
                })

            for stage, p in enumerate(protocol["external_trials"]["probabilities_p"], start=1):
                result = pulse_result(seed(size, h, g, coin), length, stage, Fraction(p))
                assert result["joint_completion_and_archive_probability"] >= gates["minimum_external_joint_success"]
                assert result["correct_stage_probability"] >= gates["minimum_external_correct_stage_probability"]
                pulses.append({"h": h, "gene": g, "coin": coin, "stage": stage, "p": p, **result})

        mean_scores = {}
        for mode in MODES:
            rows = [row for row in controls if row["mode"] == mode]
            mean = sum(row["score_head_and_original_target"] for row in rows) / len(rows)
            null_mean = sum(row["null_score_head_and_original_target"] for row in rows) / len(rows)
            expected = 1 if mode == "intact" else (0.5 if mode.endswith("scrambled") else 0)
            assert mean == expected
            assert null_mean == (1 if mode.startswith("M ") else expected)
            mean_scores[mode] = {"read_record": mean, "pilot_null": null_mean}
        summary.append({
            "tape_sites": size,
            "clock_edges": length,
            "declared_update_links": len(edges),
            "certificate": window_certificate(size, length),
            "internal_intact_rows": exact_arms,
            "internal_control_rows": controls,
            "external_single_pulse_rows": pulses,
            "aggregate_control_scores": mean_scores,
            "worst_external_log_failure": max(
                row["log_joint_completion_and_archive_failure"] for row in pulses
            ),
        })

    return {
        "protocol_id": protocol["id"],
        "protocol_sha256": sha256(protocol_bytes).hexdigest(),
        "freeze_commit": FREEZE_COMMIT,
        "sizes": summary,
        "negative_controls": {
            "inverted_inputs": "No original whole-word success.",
            "disabled_inputs": "No original whole-word success.",
            "scrambled_inputs": "Half original whole-word success.",
            "null": "Identical intact trajectory and acquired record, yet memory interventions leave the null's output unchanged.",
            "recurrence": "All one-head states recur at 2T for the supplied finite clock.",
            "unspecified_future": "For identical accessible conditional states and four equal-prior whole words, success is at most 1/4."
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "RESULTS.json")
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n")
    print(json.dumps([
        {"N": row["tape_sites"], "controls": row["aggregate_control_scores"],
         "worst_external_log_failure": row["worst_external_log_failure"]}
        for row in result["sizes"]
    ], indent=2))


if __name__ == "__main__":
    main()
