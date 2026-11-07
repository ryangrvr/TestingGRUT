"""Run the frozen calibration or open the withheld disturbance program."""

import argparse
import hashlib
import json
import platform
from math import comb
from pathlib import Path

import numpy as np
import scipy

import model as m
import reduced
from reference import FullEvolution, clock_binomial, interaction_supports
from reference import path_history, predictive_rank, primitive_components


MID = np.pi / (2 * m.OMEGA)
END = np.pi / m.OMEGA
TIME_GRID = np.array(m.CARD["clock_sample_times_over_pi"]) * END
CLOCK_FLOOR = m.CARD["condition_probability_minimum"]
CLOCK_LEVELS = np.arange(m.LINKS + 1)
ENERGY_WEIGHTS = np.array([comb(m.LINKS, k) / 2 ** m.LINKS for k in CLOCK_LEVELS])
READOUT = np.indices((2, 2, 2))
RECORD_AGREEMENT = READOUT[0] == READOUT[1]
ACTION_AGREEMENT = READOUT[0] == READOUT[2]
FIXED_TARGET = READOUT[2] == 0


def distribution_at_all_ticks(state):
    return np.array([m.probability(m.conditioned_work(state, k)) for k in CLOCK_LEVELS])


def max_tv(first, second):
    return float(np.max(np.sum(abs(first - second), axis=(-3, -2, -1)) / 2))


def distance(first, second):
    return float(np.sum(abs(np.linalg.eigvalsh(first - second))) / 2)


def entropy(density):
    eigenvalues = np.linalg.eigvalsh(density)
    occupied = eigenvalues[eigenvalues > 1e-14]
    return float(-np.sum(occupied * np.log2(occupied)))


def classical_information(distribution):
    joint = distribution.sum(axis=2)
    independent = joint.sum(axis=1)[:, None] * joint.sum(axis=0)
    positive = joint > 1e-14
    return float(np.sum(joint[positive] * np.log2(joint[positive] / independent[positive])))


def readout_series(distributions):
    return [{
        "tick": k,
        "joint_probabilities_s_parity_a": p.tolist(),
        "record_agreement": float(p[RECORD_AGREEMENT].sum()),
        "action_agreement": float(p[ACTION_AGREEMENT].sum()),
        "signal_parity_mutual_information_bits": classical_information(p),
    } for k, p in enumerate(distributions)]


def agrees(values, target):
    return bool(np.all(abs(np.asarray(values) - target) < m.TOL))


def clock_diagnostics(evolution, links, initial, work):
    history = path_history(links, work)
    timeseries = []
    for time in TIME_GRID:
        state = evolution.evolve(initial, time)
        clock = m.clock_distribution(state)
        condition_errors = []
        for k in np.flatnonzero(clock >= CLOCK_FLOOR):
            conditioned = m.conditioned_work(state, k)
            ideal = history[k]
            # Orthogonal component gives trace distance without subtracting near-unit squares.
            condition_errors.append(float(np.linalg.norm(conditioned - ideal * np.vdot(ideal, conditioned))))
        flat = state.reshape(m.LINKS + 1, -1)
        clock_state = flat @ flat.conj().T
        mean = float(clock @ CLOCK_LEVELS)
        timeseries.append({
            "time_over_pi": float(time / END),
            "tick_probabilities": clock.tolist(),
            "max_binomial_probability_error": float(np.max(abs(clock - clock_binomial(time)))),
            "max_conditioned_work_trace_distance": max(condition_errors),
            "mean_tick": mean,
            "tick_variance": float(clock @ (CLOCK_LEVELS - mean) ** 2),
            "clock_work_entanglement_entropy_bits": entropy(clock_state),
            "clock_purity": float(np.trace(clock_state @ clock_state).real),
        })

    lo, hi = m.CARD["clock_rms_interval_over_pi"]
    rms_times = END * np.linspace(lo, hi, m.CARD["clock_rms_grid_points"])
    estimates = 2 * np.arcsin(np.sqrt(CLOCK_LEVELS / m.LINKS)) / m.OMEGA
    rms = np.array([np.sqrt(clock_binomial(time) @ (estimates - time) ** 2) for time in rms_times])
    mid_state = evolution.evolve(initial, MID)
    read_state = np.concatenate([m.at_tick(k, mid_state.reshape(m.LINKS + 1, m.DIM, -1)[k])
                                 for k in CLOCK_LEVELS], axis=1)
    after_read = m.clock_distribution(evolution.evolve(read_state, MID))
    unread = m.clock_distribution(evolution.evolve(initial, END))

    eigenhistory = np.concatenate([np.sqrt(p) * state for p, state in zip(ENERGY_WEIGHTS, history)])
    eigen_residual = np.linalg.norm(evolution.h @ eigenhistory - (m.OMEGA * m.LINKS / 2) * eigenhistory)
    eigen_final = evolution.evolve(eigenhistory, MID)
    nominal = distribution_at_all_ticks(mid_state)
    return {
        "time_samples": timeseries,
        "single_read_rms_at_midpoint": float(np.sqrt(clock_binomial(MID) @ (estimates - MID) ** 2)),
        "maximum_single_read_rms_on_grid": float(rms.max()),
        "single_read_rms_grid": np.column_stack((rms_times / END, rms)).tolist(),
        "final_tick_probability_unread": float(unread[-1]),
        "final_tick_probability_after_ignored_clock_read": float(after_read[-1]),
        "expected_final_tick_probability_after_read": comb(2 * m.LINKS, m.LINKS) / 4 ** m.LINKS,
        "stationary_history_eigen_residual": float(eigen_residual),
        "stationary_vs_localized_conditional_tv": max_tv(distribution_at_all_ticks(eigen_final), nominal),
    }


def algebra_diagnostics(evolution):
    groups = primitive_components()
    commutant = ((m.GAUGE_X, m.GAUGE_Y, m.GAUGE_Z), (m.TAU_X, m.TAU_Y, m.P))
    pauli_residual = max(np.linalg.norm(x @ y - 1j * z) for x, y, z in commutant)
    cross_residual = max(np.linalg.norm(a @ b - b @ a) for a in commutant[0] for b in commutant[1])
    symmetries = {}
    for register in ("S", "M1", "D1", "D2"):
        diagonal = np.tile(np.diag(m.work_operator(**{register: m.Z})), m.LINKS + 1)
        symmetries[register] = float(np.linalg.norm(evolution.h * (diagonal[None, :] - diagonal[:, None])))
    parity = np.repeat((-1.0) ** CLOCK_LEVELS, m.DIM)
    chiral = np.linalg.norm(evolution.h * (1 + parity[:, None] * parity[None, :]))
    return {
        "primitive_noncommutation_components": groups,
        "interaction_supports": interaction_supports(evolution.h),
        "memory_pauli_multiplication_residual": float(pauli_residual),
        "memory_algebra_cross_commutator_residual": float(cross_residual),
        "work_symmetry_commutator_norms": symmetries,
        "clock_parity_anticommutator_norm": float(chiral),
        "hermiticity_residual": float(np.linalg.norm(evolution.h - evolution.h.conj().T)),
    }


def program_test(evolution, links, program):
    work = m.preparation(program)
    initial = m.at_tick(0, work)
    state = evolution.evolve(initial, MID)
    exact = distribution_at_all_ticks(state)
    compressed = reduced.solve(work, program)
    comparator = reduced.populations(work, program)
    mixed = m.mixed_preparation(program)
    mixed_exact = distribution_at_all_ticks(evolution.evolve(m.at_tick(0, mixed), MID))
    return state, {
        "program": list(program),
        "energy_probabilities": evolution.energy_distribution(initial).tolist(),
        "exact": readout_series(exact),
        "reduced": readout_series(compressed),
        "population_only": readout_series(comparator),
        "pure_reduced_max_tv": max_tv(exact, compressed),
        "mixed_reduced_max_tv": max_tv(mixed_exact, reduced.solve(mixed, program)),
        "population_only_max_tv": max_tv(exact, comparator),
        "backward_span": predictive_rank(links, program),
    }


def record_and_action(evolution, state, work):
    clean = distribution_at_all_ticks(state)
    one = m.conditioned_work(state, 1)
    two = m.conditioned_work(state, 2)
    prior = m.conditioned_work(state, 5)
    local_memory_residual = 0.0
    for k in range(1, m.LINKS + 1):
        conditioned = m.conditioned_work(state, k)
        for source in (0, 1):
            branch = conditioned * (m.S_BIT == source)[:, None]
            branch /= np.linalg.norm(branch)
            for carrier in ("M1", "M2"):
                local_memory_residual = max(local_memory_residual, distance(m.marginal(branch, (carrier,)), m.I / 2))
    projected = m.unobserved_read(one, m.M1_BIT)
    resumed = evolution.evolve(m.at_tick(1, projected), MID)
    resumed_clock = m.clock_distribution(resumed)
    after_read = []
    for k in range(2, m.LINKS + 1):
        if resumed_clock[k] >= CLOCK_FLOOR:
            row = m.probability(m.conditioned_work(resumed, k))
            after_read.append({"tick": k, "probability_of_tick": float(resumed_clock[k]),
                               "record_agreement": float(row[RECORD_AGREEMENT].sum())})

    variants = (
        ("scrambled", m.preparation((0, 0), scrambled=True), evolution, {}),
        ("writer_removed", work, FullEvolution(m.gates(no_writer=True)), {"no_writer": True}),
        ("fixed_action", work, FullEvolution(m.gates(fixed_action=True)), {"fixed_action": True}),
    )
    controls = {}
    for name, input_work, solver, flags in variants:
        initial = m.at_tick(0, input_work)
        control_state = solver.evolve(initial, MID)
        exact = distribution_at_all_ticks(control_state)
        controls[name] = {
            "exact": readout_series(exact),
            "reduced_max_tv": max_tv(exact, reduced.solve(input_work, (0, 0), **flags)),
            "energy_probabilities": solver.energy_distribution(initial).tolist(),
            "pure_state_purity": float(np.linalg.norm(control_state) ** 4),
        }
        if name == "scrambled":
            scrambled_prior = m.conditioned_work(control_state, 5)
            matched_marginals = {
                "+".join(registers): distance(m.marginal(prior, registers),
                                             m.marginal(scrambled_prior, registers))
                for registers in (("S",), ("M1", "M2"), ("A",), ("D1", "D2"))
            }
        if name == "fixed_action":
            controls[name]["max_target_bit_one_probability"] = float(np.max(1 - exact[:, FIXED_TARGET].sum(axis=1)))
    before_action = clean[5]
    # Use separated S, full dimer, and target marginals at the handoff.
    product = (before_action.sum(axis=(1, 2))[:, None, None]
               * before_action.sum(axis=(0, 2))[None, :, None]
               * before_action.sum(axis=(0, 1))[None, None, :])
    product_success = reduced.act(product, False)[ACTION_AGREEMENT].sum()
    clean_success = clean[7][ACTION_AGREEMENT].sum()
    scrambled_success = controls["scrambled"]["exact"][7]["action_agreement"]
    return {
        "resource_entropies_bits": [{
            "tick": k,
            "memory": entropy(m.marginal(w, ("M1", "M2"))),
            "signal": entropy(m.marginal(w, ("S",))),
            "target": entropy(m.marginal(w, ("A",))),
            "signal_memory_quantum_mutual_information": (
                entropy(m.marginal(w, ("M1", "M2"))) + entropy(m.marginal(w, ("S",)))
                - entropy(m.marginal(w, ("S", "M1", "M2")))
            ),
        } for k, w in ((0, work), (1, one), (7, m.conditioned_work(state, 7)))],
        "gauge_x_at_tick_1": m.expectation(one, m.GAUGE_X),
        "gauge_x_at_tick_2": m.expectation(two, m.GAUGE_X),
        "gauge_y_at_tick_2": m.expectation(two, m.GAUGE_Y),
        "gauge_x_after_memory_read": m.expectation(projected, m.GAUGE_X),
        "max_source_conditioned_single_carrier_distance_to_mixed": local_memory_residual,
        "parity_preserving_read_future_ticks": after_read,
        "primary_purity": float(np.linalg.norm(state) ** 4),
        "scrambled_marginal_trace_distances_at_tick_5": matched_marginals,
        "product_of_marginals_action_prediction": float(product_success),
        "marginals_only_worst_case_lower_bound": float(abs(clean_success - scrambled_success) / 2),
        "controls": controls,
    }


def calibration():
    links = m.gates()
    evolution = FullEvolution(links)
    nominal_program = tuple(m.CARD["calibration_programs"][0])
    work = m.preparation(nominal_program)
    initial = m.at_tick(0, work)
    algebra = algebra_diagnostics(evolution)
    clock = clock_diagnostics(evolution, links, initial, work)
    tests = {}
    for program in m.CARD["calibration_programs"]:
        state, report = program_test(evolution, links, tuple(program))
        tests["".join(map(str, program))] = report
        if tuple(program) == nominal_program:
            nominal = state
    record = record_and_action(evolution, nominal, work)
    clean = tests["00"]["exact"]
    controls = record["controls"]
    interaction = algebra["interaction_supports"]
    expected_edges = {"+".join(sorted(("C", name))) for name in m.REGISTERS}
    expected_edges.update(("M2+S", "D1+M2", "D2+M2", "M1+M2", "A+M1", "A+M2"))
    energy_rows = [row["energy_probabilities"] for row in [*tests.values(), *controls.values()]]
    components = {frozenset(g) for g in algebra["primitive_noncommutation_components"]}
    expected_components = {frozenset((name + ":X", name + ":Z")) for name in m.REGISTERS}
    expected_components.add(frozenset([f"Q{k}" for k in CLOCK_LEVELS]
                                     + [f"B{k}" for k in range(m.LINKS)]))
    checks = {
        "microscopic_access_components": components == expected_components,
        "interaction_graph_and_parseval": set(interaction["pair_graph"]) == expected_edges
            and agrees(interaction["sum_weights_squared"], interaction["direct_normalized_trace_h_squared"]),
        "symmetries_and_algebras": agrees([
            algebra["memory_pauli_multiplication_residual"],
            algebra["memory_algebra_cross_commutator_residual"],
            algebra["clock_parity_anticommutator_norm"],
            algebra["hermiticity_residual"],
            *algebra["work_symmetry_commutator_norms"].values(),
        ], 0),
        "binomial_clock_and_exact_conditioning": all(max(
            row["max_binomial_probability_error"], row["max_conditioned_work_trace_distance"]
        ) < m.TOL for row in clock["time_samples"]),
        "stationary_history": max(clock["stationary_history_eigen_residual"],
                                 clock["stationary_vs_localized_conditional_tv"]) < m.TOL,
        "clock_read_backaction": agrees(clock["final_tick_probability_unread"], 1)
            and agrees(clock["final_tick_probability_after_ignored_clock_read"],
                       clock["expected_final_tick_probability_after_read"]),
        "matched_full_energy_distributions": all(agrees(row, ENERGY_WEIGHTS) for row in energy_rows),
        "record_written_and_persistent": agrees(clean[0]["record_agreement"], 0.5)
            and agrees([row["record_agreement"] for row in clean[1:]], 1)
            and agrees([row["signal_parity_mutual_information_bits"] for row in clean[1:]], 1)
            and agrees(record["max_source_conditioned_single_carrier_distance_to_mixed"], 0),
        "persistent_record_with_changed_dimer_and_read_disturbance": agrees(record["gauge_x_at_tick_1"], 1)
            and agrees(record["gauge_x_at_tick_2"], 0) and agrees(record["gauge_y_at_tick_2"], 1)
            and agrees(record["gauge_x_after_memory_read"], 0)
            and len(record["parity_preserving_read_future_ticks"]) == m.LINKS - 1
            and agrees([row["record_agreement"] for row in record["parity_preserving_read_future_ticks"]], 1),
        "writer_removed": agrees([row["record_agreement"] for row in controls["writer_removed"]["exact"]], 0.5),
        "informative_action": agrees([row["action_agreement"] for row in clean[6:]], 1),
        "scrambled_with_matched_marginals": agrees(list(record["scrambled_marginal_trace_distances_at_tick_5"].values()), 0)
            and agrees([record["primary_purity"], controls["scrambled"]["pure_state_purity"]], 1)
            and agrees(controls["scrambled"]["exact"][5]["signal_parity_mutual_information_bits"], 0)
            and agrees([row["action_agreement"] for row in controls["scrambled"]["exact"][6:]], 0.5),
        "fixed_action": agrees(controls["fixed_action"]["max_target_bit_one_probability"], 0)
            and agrees([row["action_agreement"] for row in controls["fixed_action"]["exact"][6:]], 0.5),
        "marginals_only_error_bound": agrees(record["product_of_marginals_action_prediction"], 0.5)
            and agrees(record["marginals_only_worst_case_lower_bound"], 0.25),
        "backward_ranks": [tests[k]["backward_span"]["rank_including_identity"] for k in ("00", "10", "01")] == [8, 12, 12],
        "reduced_matches_reference": all(max(x["pure_reduced_max_tv"], x["mixed_reduced_max_tv"]) < m.TOL
                                         for x in tests.values())
            and all(x["reduced_max_tv"] < m.TOL for x in controls.values()),
        "population_only_matches_calibration": all(x["population_only_max_tv"] < m.TOL for x in tests.values())
            and agrees([tests[k]["exact"][7]["action_agreement"] for k in ("10", "01")], 0.75),
    }
    return {"checks": checks, "algebra": algebra, "clock": clock, "programs": tests, "record": record,
            "energy_levels": (m.OMEGA * (CLOCK_LEVELS - m.LINKS / 2)).tolist(),
            "expected_energy_probabilities": ENERGY_WEIGHTS.tolist()}


def holdout():
    links = m.gates()
    evolution = FullEvolution(links)
    program = tuple(m.CARD["withheld_program"])
    state, result = program_test(evolution, links, program)
    middle = m.conditioned_work(state, 3)
    pinched = m.unobserved_read(middle, m.PARITY_BIT)
    readout_change = max_tv(m.probability(middle), m.probability(pinched))
    start = m.at_tick(3, middle)
    pinched_start = m.at_tick(3, pinched)
    continued = evolution.evolve(pinched_start, MID)
    endpoint = m.conditioned_work(continued, 7)
    final = m.probability(endpoint)
    reference_energy_change = float(np.max(abs(
        evolution.energy_distribution(start) - evolution.energy_distribution(pinched_start)
    )))
    exact_success = result["exact"][7]["action_agreement"]
    comparator_success = result["population_only"][7]["action_agreement"]
    pinched_success = float(final[ACTION_AGREEMENT].sum())
    separation = exact_success - comparator_success
    checks = {
        "frozen_coherent_prediction": agrees(exact_success, m.CARD["withheld_success_probability"]),
        "frozen_population_prediction": agrees(comparator_success, m.CARD["withheld_population_only_success_probability"]),
        "coherence_decision_margin": separation >= m.CARD["withheld_minimum_separation"],
        "pinched_full_state_matches_population_only": agrees(pinched_success, comparator_success)
            and max_tv(final, np.array(result["population_only"][7]["joint_probabilities_s_parity_a"])) < m.TOL,
        "pinching_matches_immediate_readout_and_energy_at_tick_3": max(readout_change, reference_energy_change) < m.TOL,
        "backward_rank": result["backward_span"]["rank_including_identity"] == 16,
        "reduced_matches_reference": max(result["pure_reduced_max_tv"], result["mixed_reduced_max_tv"]) < m.TOL,
        "withheld_total_energy_distribution": agrees(result["energy_probabilities"], ENERGY_WEIGHTS),
    }
    return {
        "checks": checks,
        "program": result,
        "coherent_minus_population_agreement": separation,
        "pinched_full_state_final_agreement": pinched_success,
        "pinched_final_distribution": final.tolist(),
        "pinched_tick_7_probability": float(m.clock_distribution(continued)[7]),
        "pinching_immediate_readout_tv": readout_change,
        "pinching_energy_distribution_max_difference": float(reference_energy_change),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("calibration", "holdout"))
    parser.add_argument("--calibration", type=Path, help="Passed calibration JSON required for holdout")
    args = parser.parse_args()
    calibration_sha = None
    if args.stage == "holdout":
        if args.calibration is None:
            parser.error("Run and supply the calibration before opening the withheld program.")
        calibration_bytes = args.calibration.read_bytes()
        prior = json.loads(calibration_bytes)
        if (prior.get("stage") != "calibration" or prior.get("status") != "passed"
                or prior.get("parameter_card_sha256") != m.CARD_SHA
                or prior.get("specification_sha256") != m.SPEC_SHA
                or not prior.get("checks", {}).get("population_only_matches_calibration")):
            parser.error("The passed calibration must use the frozen model and population comparator.")
        calibration_sha = hashlib.sha256(calibration_bytes).hexdigest()
    result = calibration() if args.stage == "calibration" else holdout()
    passed = all(result["checks"].values())
    print(json.dumps({
        "stage": args.stage,
        "status": "passed" if passed else "failed",
        "specification_sha256": m.SPEC_SHA,
        "parameter_card_sha256": m.CARD_SHA,
        "calibration_sha256": calibration_sha,
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        **result,
    }, indent=2, allow_nan=False))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
