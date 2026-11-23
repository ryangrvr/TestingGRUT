"""Run only the published GRUT-D2-01 family. No fitted parameters."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import numpy as np
import scipy
from model import (WallPath, size_card, bounds, count_information, count_joint,
                   order_probability, initial_resources, kicked_expected_walls,
                   bath_echo_block, echo_energy_variance, subsystem_entropies)

ROOT = Path(__file__).resolve().parent


def run():
    protocol_bytes = (ROOT / "PROTOCOL.json").read_bytes()
    protocol = json.loads(protocol_bytes)
    output = {"protocol_id": protocol["id"],
              "freeze_commit": "8e3a0773dbb4d8d54cb9261f12e76e2129fae403",
              "protocol_sha256": hashlib.sha256(protocol_bytes).hexdigest(),
              "runtime": {"python": platform.python_version(), "numpy": np.__version__,
                          "scipy": scipy.__version__, "platform": platform.system()},
              "cases": [], "negative_results": {
                  "strict_target_free_emergence": "not_established",
                  "record_dependent_regulation": "not_established",
                  "derived_internal_time": "not_established",
                  "new_GRUT_physics": "not_established",
                  "cosmological_selection": "not_established",
                  "consciousness": "not_established"}}
    for n in protocol["system_spin_sizes"]:
        card = size_card(n)
        m, L = card["m"], card["L"]
        path = WallPath(m, L)
        ts = np.linspace(card["window_start"], card["window_end"],
                         protocol["family"]["time_samples"])
        rows = []
        for t in ts:
            a = path.block(t)
            joint = count_joint(a)
            information = count_information(joint)
            order = order_probability(a)
            rows.append({"t": float(t), "order_probability": order,
                         "expected_walls": 0.5 * float(np.linalg.norm(a)**2),
                         "count_record_success": information["success"],
                         "count_mutual_information_bits": information["mutual_information_bits"],
                         "joint_normalization_error": information["normalization_error"],
                         "order_count_identity_error": abs(order - information["success"]),
                         "joint_K_Nbath": joint.tolist()})
        tau = card["window_start"]
        a = path.block(tau)
        echo = order_probability(bath_echo_block(a))
        recovery = kicked_expected_walls(path, tau, tau)
        kick_wall = 0.5 * float(np.linalg.norm(a[-1])**2)
        item = {**card, "resources": initial_resources(m, L),
                "entropy_at_window_start": subsystem_entropies(a),
                "certificate": bounds(m, L, card["window_start"], card["window_end"]),
                "sampled": {"minimum_order_probability": min(r["order_probability"] for r in rows),
                            "minimum_count_mutual_information_bits": min(r["count_mutual_information_bits"] for r in rows),
                            "maximum_expected_walls": max(r["expected_walls"] for r in rows),
                            "maximum_normalization_error": max(r["joint_normalization_error"] for r in rows),
                            "maximum_order_count_identity_error": max(r["order_count_identity_error"] for r in rows)},
                "disturbance": {"time": tau, "recovery_time": 2*tau,
                                "kick_mean_work": 4 * (1 - 2*kick_wall),
                                "mean_work_formula": "Delta*(1-2*<w_last>) ; real hopping expectations vanish",
                                "expected_walls_after_recovery": recovery,
                                "order_probability_union_lower": 1 - recovery,
                                "shifted_count_record_union_lower": 1 - recovery - kick_wall,
                                "bath_reset": False},
                "bath_phase_echo": {"kick_time": tau, "read_time": 2*tau,
                                    "current_system_marginal_change": 0,
                                    "current_count_readout_change": 0,
                                    "baseline_final_order_probability": rows[-1]["order_probability"],
                                    "echo_final_order_probability": echo,
                                    "echo_operation": "Z on bath sites whose full path index is odd",
                                    "mean_work": 0,
                                    "energy_variance_after_kick": echo_energy_variance(path, tau),
                                    "energy_variance_before_kick": initial_resources(m, L)["energy_variance"],
                                    "work_scope": "zero mean work for this bipartite real preparation; coherent controller and timing supplied",
                                    "causal_scope": "quantum bath phase affects future, no classical-count feedback demonstrated"},
                "time_rows": rows}
        output["cases"].append(item)
        print(f"n={n}, L={L}: sampled order >= {item['sampled']['minimum_order_probability']:.12f}; "
              f"certified > {item['certificate']['order_lower_display']:.8f}; "
              f"repair union >= {1-recovery:.12f}", flush=True)
    output["short_bath_recurrence"] = {
        "n": 2, "L": 1,
        "times": [0, float(np.pi/2), float(np.pi)],
        "exact_order_probabilities": [0.5, 1.0, 0.5],
        "interpretation": "Rabi transfer returns to initial disorder; no attractor"}
    output["maximally_mixed_active_control"] = {
        "n": 2, "L": 1, "initial_reference": "K=initial wall",
        "active_state": "I/4 on wall+bath, orientation independent fair bit",
        "full_reference_active_state": "correlated, not maximally mixed",
        "test_time_gt": float(np.pi/2),
        "initial_order_probability": 0.5, "final_order_probability": 0.5,
        "initial_count_record_success": 0.5, "final_count_record_success": 1.0,
        "frozen_no_new_record_forecast": "failed",
        "correction": "Stationarity of active marginal does not prevent transfer of initial reference correlations.",
        "source": "CONTROL_CORRECTION_01.md; exact exchange checked independently"}
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "RESULTS.json")
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
