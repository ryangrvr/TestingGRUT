"""Run the frozen finite comparisons and produce the D4 evidence ledger."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
import json
from math import log2, pi, sqrt
from pathlib import Path

import numpy as np

from model import DELTA, KAPPA, MASKS, ready_statistics, summarize, uniform_certificate

ROOT = Path(__file__).resolve().parent
FREEZE_COMMIT = "5c59a2dce5858752a785f563dd779fa1034a9371"


def expected_counts(sequence: tuple[str, ...]) -> tuple[int, int]:
    return tuple(sum(MASKS[site][i] for site in sequence) for i in range(2))


def scaling_results() -> dict:
    rows = []
    for labels in [2, 3, 4, 8, 16, 64]:
        syndrome_outputs = 1 + 3 * (labels - 1)
        rows.append({
            "historical_configurations": labels,
            "uniform_one_bit_side_information_success_ceiling": float(min(Fraction(1), Fraction(2, labels))),
            "bits_for_perfect_uniform_history": log2(labels),
            "orthogonal_outputs_needed_for_no_or_one_replica_error": syndrome_outputs,
            "minimum_pure_environment_capacity_bits_for_one_exact_repair": log2(syndrome_outputs),
            "three_replica_at_most_two_errors_minimax_ceiling_on_worst_tuple": 1 / min(labels, 3),
        })
    enumeration = []
    for labels in [2, 3, 4, 8]:
        for copies, max_faults in [(2, 1), (3, 1), (3, 2)]:
            worst = max(
                sum(word.count(symbol) >= copies - max_faults for symbol in range(labels))
                for word in product(range(labels), repeat=copies)
            )
            forecast = min(labels, copies // (copies - max_faults))
            assert worst == forecast
            enumeration.append({
                "labels": labels,
                "copies": copies,
                "allowed_faults": max_faults,
                "max_compatible_histories": worst,
                "minimax_success_ceiling_at_worst_tuple": 1 / worst,
            })
    return {"information_bound_rows": rows, "finite_tuple_checks": enumeration}


def short_chain_controls() -> dict:
    # Full writer-on invariant sector with one M1 defect and L1=L2=1:
    # [defect, aligned with excitation in B1, blank with excitation in B1].
    hopping = np.array([[0, 1, 0], [1, 0, KAPPA], [0, KAPPA, 0]], dtype=float)
    omega = sqrt(1 + KAPPA**2)
    values, vectors = np.linalg.eigh(hopping)

    def state(source: int, time: float) -> np.ndarray:
        return vectors @ (np.exp(-1j * values * time) * vectors[source].conj())

    recurrence_time = 2 * pi / omega
    recovery_time = pi / (2 * omega)
    recurrence = float(abs(state(0, recurrence_time)[0]) ** 2)
    recovery = float(abs(state(0, recovery_time)[1]) ** 2)
    loaded = float(abs(state(1, recovery_time)[0]) ** 2)
    assert abs(recurrence - 1) < 1e-10
    assert abs(recovery - 1 / (1 + KAPPA**2)) < 1e-10
    assert abs(loaded - recovery) < 1e-10
    return {
        "omega_over_g": omega,
        "recurrence_time": recurrence_time,
        "defect_probability_at_recurrence": recurrence,
        "early_joint_repair": recovery,
        "loaded_interface_discrepancy_without_replica_flip": loaded,
        "loaded_interface_scope": "Diagnostic starts with an excited bath site. The vacuum no-fault control remains dark.",
    }


def run() -> dict:
    protocol_bytes = (ROOT / "PROTOCOL.json").read_bytes()
    protocol = json.loads(protocol_bytes)
    certificate = uniform_certificate(protocol)
    assert certificate["strict_write_gate_pass"]
    assert all(row["strict_codeword_gate_pass"] and row["strict_counts_gate_pass"]
               for row in certificate["after_faults"])
    lengths = tuple(protocol["correction_chain_spins"])
    sites = protocol["fault_sites"]
    tolerance = protocol["verification"]["tolerance"]
    cycles = protocol["repeated"]["max_faults"]
    last_count_bound = certificate["after_faults"][-1]["historical_codeword_and_counts_lower_bound"]
    longest_time = cycles * protocol["repeated"]["gaps_and_terminal_wait_range"][1]

    repeated_rows = []
    max_mode_error = 0.0
    max_normalization_error = 0.0
    min_counts_comparator = 1.0
    for gap in protocol["repeated"]["sample_equal_gaps"]:
        for sequence in product(sites, repeat=cycles):
            pulses = tuple((gap * j, (site,)) for j, site in enumerate(sequence))
            result = ready_statistics(lengths, pulses, gap * cycles)
            joint = result["joint"]
            measures = summarize(joint)
            expected = expected_counts(sequence)
            count_success = float(joint[0, 0, 0, expected[0], expected[1]])
            min_counts_comparator = min(min_counts_comparator, count_success)
            assert measures["normalization_error"] < tolerance
            max_normalization_error = max(max_normalization_error, measures["normalization_error"])
            max_mode_error = max(max_mode_error, *[
                item["largest_mode_gram_error"] for item in result["mode_diagnostics"]
            ])
            repeated_rows.append({
                "fault_sequence": list(sequence),
                "wait_after_each_fault": gap,
                "expected_bath_counts": list(expected),
                "historical_joint_and_exact_counts": count_success,
                "readout": measures,
            })

    single = {}
    for site in sites:
        pulse = ((0.0, (site,)),)
        joint = ready_statistics(lengths, pulse, 136.0)["joint"]
        immediate = ready_statistics(lengths, pulse, 0.0)["joint"]
        mask = MASKS[site]
        assert abs(immediate[0, mask[0], mask[1], 0, 0] - 1) < tolerance
        single[site] = {
            "internal_parity_syndrome": list(mask),
            "exact_fault_work_at_nominal_write": DELTA * sum(mask),
            "readout": summarize(joint),
        }

    simultaneous = []
    for size in (2, 3):
        for at_once in combinations(sites, size):
            pulses = ((0.0, at_once),)
            immediate = summarize(ready_statistics(lengths, pulses, 0.0)["joint"])
            output = summarize(ready_statistics(lengths, pulses, 136.0)["joint"])
            assert output["historical_joint"] < tolerance
            simultaneous.append({
                "sites_at_once": list(at_once),
                "parity_probabilities_immediately": immediate["parity_probabilities"],
                "readout_in_ready_comparison": output,
                "exact_writer_on_historical_codeword_at_nominal_write": 0,
                "verdict": "mismatch but opposite majority" if size == 2 else "opposite codeword undetected",
            })

    single_state = ready_statistics(lengths, ((0.0, ("S",)),), 136.0)["joint"]
    opposite_history_state = ready_statistics(lengths, ((0.0, ("M1", "M2")),), 136.0)["joint"][::-1]
    collision_error = float(np.max(abs(single_state - opposite_history_state)))
    assert collision_error < tolerance
    assert max_mode_error < tolerance
    maximum_count = 2 * cycles
    return {
        "protocol_id": protocol["id"],
        "protocol_sha256": sha256(protocol_bytes).hexdigest(),
        "freeze_commit": FREEZE_COMMIT,
        "certificate_for_physical_writer_on_hamiltonian": certificate,
        "comparison_scope": "Ready-sector calculation uses exact physical fault maps with the writer exchange omitted. The uniform proof bounds the full writer-on system.",
        "samples_in_ready_comparison": {
            "four_fault_rows": repeated_rows,
            "minimum_joint": min(row["readout"]["historical_joint"] for row in repeated_rows),
            "minimum_joint_with_exact_bath_counts": min_counts_comparator,
            "maximum_mode_gram_error": max_mode_error,
            "maximum_normalization_error": max_normalization_error,
            "no_reservoir_reset": True,
            "no_intermediate_measurement": True,
        },
        "single_faults": single,
        "simultaneous_faults": simultaneous,
        "ambiguous_pair": {
            "history_zero_one_system_error": [1, 0, 0],
            "history_one_two_memory_errors": [1, 0, 0],
            "full_state_collision_exact_after_nominal_write": True,
            "maximum_post_evolution_probability_discrepancy": collision_error,
            "equal_prior_historical_success_ceiling": 0.5,
        },
        "no_fault": {
            "duration_after_nominal_write": longest_time,
            "false_repair_bath_activation": 0,
            "historical_system": 1,
            "historical_joint": float(np.cos(KAPPA * longest_time)**2),
            "writer_blank_due_to_reversion": float(np.sin(KAPPA * longest_time)**2),
            "writer_reversion_time_after_nominal_write": pi / (2 * KAPPA),
        },
        "short_reservoir": short_chain_controls(),
        "maximum_load_four_system_faults": {
            "expected_bath_counts": [cycles, cycles],
            "physical_mean_total_absorbed_quanta_interval": [last_count_bound * maximum_count, maximum_count],
            "physical_mean_net_pulse_work_interval": [
                (DELTA - 2) * last_count_bound * maximum_count - KAPPA,
                (DELTA + 2) * maximum_count + KAPPA,
            ],
            "work_scope": "A conservative bound from charge conservation, occupancy, and hopping spectral width. Work is supplied by the fault pulse source.",
        },
        "scaling": scaling_results(),
        "verdicts": {
            "historical_restoration": "pass for separated single faults",
            "internal_fault_detection": "pass for ready-state mismatch; culprit identification conditional on one fault",
            "repeated_self_maintenance": "pass for four separated faults with the same prepared baths",
            "autonomous_timing": "correction dynamics always on; writing window and diagnostic faults timed externally",
            "physical_novelty": "standard quantum mechanics and a supplied repetition-code repair rule",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "RESULTS.json")
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    samples = result["samples_in_ready_comparison"]
    print(json.dumps({
        "repeated_timing_rows": len(samples["four_fault_rows"]),
        "minimum_comparison_joint": samples["minimum_joint"],
        "minimum_comparison_joint_and_counts": samples["minimum_joint_with_exact_bath_counts"],
        "physical_four_fault_certificate": result["certificate_for_physical_writer_on_hamiltonian"]["after_faults"][-1],
    }, indent=2))


if __name__ == "__main__":
    main()
