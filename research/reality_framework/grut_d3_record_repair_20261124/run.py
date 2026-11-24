"""Run the published finite GRUT-D3 design and write its numerical record."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model import (
    outcome_probabilities,
    polynomial_propagation,
    repair_adjacency,
    state_after_fault,
    window_certificate,
)

ROOT = Path(__file__).resolve().parent
FREEZE_COMMIT = "53fd7c9e13185e79f0983932f32807f543591e2a"


def write_diagnostics(state: np.ndarray, adjacency) -> dict:
    blank = float(abs(state[0]) ** 2)
    applied = adjacency @ state
    hopping_mean = float(np.vdot(state, applied).real)
    variance = float(np.vdot(applied, applied).real - hopping_mean**2)
    return {
        "written_record_probability": 1 - blank,
        "blank_probability": blank,
        "intact_pointer_information_bits": 1 - blank,
        "scrambled_pointer_information_bits": 0.0,
        "before_fault_energy_mean": 4.0,
        "before_fault_energy_variance": 1.0,
        "after_fault_energy_mean": 8 + hopping_mean,
        "after_fault_energy_variance": variance,
        "fault_mean_work": 4 + hopping_mean,
        "norm_error": float(abs(np.vdot(state, state) - 1)),
    }


def short_chain_control() -> dict:
    write_time = float(np.pi / 2)
    repair_time = float(np.pi / np.sqrt(2))
    vector = state_after_fault(1, 1, write_time)
    adjacency = repair_adjacency(1, 1)
    final = expm_multiply(-1j * repair_time * adjacency, vector, traceA=0.0)
    probabilities = outcome_probabilities(final, 1, 1)
    result = {
        "writer_length": 1,
        "repair_length": 1,
        "write_time": write_time,
        "repair_time": repair_time,
        "ordinary_configuration_probability": float(probabilities[:, 1].sum()),
        "restore_and_retain_probability": float(probabilities[1, 1]),
        "M_blank_probability": float(probabilities[0, :].sum()),
        "joint_gate": "fails",
    }
    assert abs(result["ordinary_configuration_probability"] - 1) < 1e-10
    assert result["restore_and_retain_probability"] < 1e-10
    return result


def run() -> dict:
    protocol_bytes = (ROOT / "PROTOCOL.json").read_bytes()
    protocol = json.loads(protocol_bytes)
    proof = window_certificate(protocol)
    assert all(proof["strict_gates_pass"].values())

    write_times = protocol["writer"]["sample_gw"]
    repair_times = protocol["repair"]["sample_gu"]
    writer_length = protocol["writer"]["L"]
    repair_length = protocol["repair"]["L"]
    states = np.column_stack([
        state_after_fault(writer_length, repair_length, write_time) for write_time in write_times
    ])
    adjacency = repair_adjacency(writer_length, repair_length)
    propagated = expm_multiply(
        -1j * adjacency, states, start=repair_times[0], stop=repair_times[-1],
        num=len(repair_times), traceA=0.0,
    )
    polynomial_times = [8, 12, 16]
    polynomial_vectors, polynomial_tail = polynomial_propagation(adjacency, states, polynomial_times)
    maximum_vector_discrepancy = max(
        float(np.linalg.norm(propagated[repair_times.index(time)] - predicted, axis=0).max())
        for time, predicted in zip(polynomial_times, polynomial_vectors)
    )
    tolerance = protocol["verification"]["absolute_match_tolerance"]
    assert maximum_vector_discrepancy < tolerance

    results = []
    writer_samples = []
    certificate_lower = proof["probability_lower_bounds"]["restore_and_retain"]
    for column, write_time in enumerate(write_times):
        diagnostics = write_diagnostics(states[:, column], adjacency)
        assert abs(diagnostics["fault_mean_work"] - 4) < tolerance
        assert diagnostics["written_record_probability"] > protocol["analytic_gates"]["window_write_accuracy_strict_lower"]
        writer_samples.append({"write_time": write_time, **diagnostics})
        for index, repair_time in enumerate(repair_times):
            probabilities = outcome_probabilities(propagated[index, :, column], writer_length, repair_length)
            intact_joint = float(probabilities[1, 1])
            restored = float(probabilities[:, 1].sum())
            norm_error = float(abs(probabilities.sum() - 1))
            assert norm_error < tolerance
            assert intact_joint > certificate_lower
            results.append({
                "write_time": write_time,
                "repair_time": repair_time,
                "intact_restore_probability": restored,
                "intact_restore_and_retain_probability": intact_joint,
                "scrambled_restore_probability": restored / 2,
                "scrambled_restore_and_retain_probability": intact_joint / 2,
                "wrong_restore_probability": 0.0,
                "wrong_restore_and_retain_probability": 0.0,
                "ordinary_configuration_probability_all_arms": restored,
                "outcomes_by_record_availability_and_repair": probabilities.tolist(),
                "norm_error": norm_error,
            })
    return {
        "protocol_id": protocol["id"],
        "protocol_sha256": hashlib.sha256(protocol_bytes).hexdigest(),
        "freeze_commit": FREEZE_COMMIT,
        "grid_dimension_for_one_selected_label": int(adjacency.shape[0]),
        "window_certificate": proof,
        "writer_samples": writer_samples,
        "time_rows": results,
        "minimum_sampled_joint_success": min(row["intact_restore_and_retain_probability"] for row in results),
        "maximum_sampled_joint_success": max(row["intact_restore_and_retain_probability"] for row in results),
        "polynomial_check": {
            "times": polynomial_times,
            "degree": 192,
            "matched_vectors": len(polynomial_times) * len(write_times),
            "maximum_L2_difference": maximum_vector_discrepancy,
            "analytic_series_tail_bound": float(polynomial_tail),
            "scope": "Tail controls exact-arithmetic series truncation. Cross-check measures floating-point agreement.",
        },
        "short_chain_control": short_chain_control(),
        "unflagged_color_flip": "Exact charge obstruction to repair with B vacuum.",
        "verdicts": {
            "historical_record_selects_repair": "derived_for_this_finite_model",
            "record_retained_on_scored_window": "derived_for_this_finite_model",
            "target_class": "supplied",
            "external_time_and_fault_controller": "supplied",
            "repeated_repair": "unestablished",
            "new_physical_law": "unestablished",
            "cosmological_selection": "unestablished",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "RESULTS.json")
    args = parser.parse_args()
    results = run()
    args.output.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "samples": len(results["time_rows"]),
        "minimum_sampled_joint_success": results["minimum_sampled_joint_success"],
        "polynomial_check": results["polynomial_check"],
        "uniform_joint_success_lower": results["window_certificate"]["probability_lower_bounds"]["restore_and_retain"],
    }, indent=2))


if __name__ == "__main__":
    main()
