"""Same-author mathematical verification, with independent matrix constructions."""
import math
import numpy as np
from common import ROOT, error, save
import json


def main():
    load = lambda name: json.loads((ROOT/name).read_text())
    model = np.load(ROOT / "HISTORY_MODEL.npz")
    result = load("WHOLE_PARTS_RESULTS.json")
    checks = []

    def group(name, predicate, diagnostics):
        checks.append({"name": name, "pass": bool(predicate), "diagnostics": diagnostics})

    # Build the clock graph independently and conjugate by its block-diagonal W.
    clock = np.roll(np.eye(4), 1, axis=0)
    w = model["unitaries"]
    block = np.zeros((64, 64), complex)
    for t in range(4):
        block[16*t:16*(t+1), 16*t:16*(t+1)] = w[t]
    expected_c = block @ np.kron(2*np.eye(4)-clock-clock.T, np.eye(16)) @ block.conj().T
    c, g, j = model["constraint"], model["global_neutrality"], model["logical_isometry"]
    spectrum = np.linalg.eigvalsh(c)
    target = np.repeat([0., 2., 2., 4.], 16)
    matrix_error = error(c-expected_c)
    spectrum_error = error(np.sort(spectrum)-np.sort(target))
    neutral_kernel = int(np.sum(np.linalg.eigvalsh(c+g@g) < 1e-10))
    group("constraints_and_three_logical_factors", matrix_error < 1e-13 and spectrum_error < 1e-12
          and neutral_kernel == 8 and error(j.conj().T@j-np.eye(8)) < 1e-13
          and error(j@j.conj().T-model["physical_projector"]) < 1e-13,
          {"independent_clock_graph_error": matrix_error, "spectrum_error": spectrum_error,
           "neutral_physical_dimension": neutral_kernel})

    fields = ("constraint_residual", "neutrality_residual", "conditional_dictionary_error",
              "physical_effect_commutator_error", "effect_completeness_error",
              "physical_effect_positivity_error", "neutral_probe_pauli_algebra_error",
              "logical_probe_dictionary_error")
    worst = max(result[f] for f in fields)
    # Direct physical probe effect at tick zero versus the group-average dictionary.
    local = np.kron(np.eye(8), np.array([[.5, .5], [.5, .5]]))
    direct = block @ np.kron(np.eye(4), local) @ block.conj().T
    state = model["state"]
    probe_error = abs(np.vdot(state, direct@state).real-result["probe_evolution"][0]["plus_probability"])
    group("physical_observation_dictionary", worst < 1e-13 and probe_error < 1e-13,
          {"maximum_stored_residual": worst, "direct_logical_effect_error": probe_error})

    # Explicit amplitude support: at final two ticks only X=R is supported.
    record_support_error = 0.
    for t in (2, 3):
        amplitudes = (w[t] @ (model["embedding"].conj().T @ state)).reshape(2,2,2,2)
        record_support_error = max(record_support_error, error(amplitudes[0,:,1,:]),
                                   error(amplitudes[1,:,0,:]))
    r = result["records"]
    record_error = max(abs(r[t]["X_charge_conditioned_on_record"][b]-[6,-6][b])
                       for t in (2,3) for b in (0,1))
    group("internal_records_and_global_coherence", record_support_error < 1e-13 and record_error < 1e-13
          and abs(result["whole_coherence_pure"]-1) < 1e-13
          and abs(result["whole_coherence_mixture"]) < 1e-13,
          {"wrong_record_support": record_support_error, "conditional_charge_error": record_error,
           "pure_whole_coherence": result["whole_coherence_pure"],
           "dephased_whole_coherence": result["whole_coherence_mixture"]})

    group("internal_asymmetry_and_clock_resolution", abs(result["internal_asymmetric_X_charge"]-6) < 1e-13
          and result["internal_asymmetric_neutrality_error"] < 1e-13
          and result["internal_asymmetric_constraint_error"] < 1e-13
          and result["clock_formula_error"] < 1e-13
          and all(abs(row["probe_purity"]-row["predicted_purity"]) < 1e-13
                  for row in result["clock_readouts"]),
          {"asymmetric_local_charge": result["internal_asymmetric_X_charge"],
           "clock_formula_error": result["clock_formula_error"],
           "ideal_tick_probe_purity": result["clock_readouts"][0]["probe_purity"]})

    group("dimension_one_and_finite_CCR_lemmas", result["dimension_one_pauli_anticommutator_min_residual"] == 2
          and abs(result["finite_ladder_commutator_trace"]) < 1e-13
          and error(np.array(result["finite_ladder_commutator_diagonal"])-[1,1,1,-3]) < 1e-13,
          {"one_dimensional_Pauli_residual": 2, "finite_ladder_trace": result["finite_ladder_commutator_trace"]})

    vacuum = load("VACUUM_COUNTEREXAMPLE.json")
    phases = vacuum["phases"]
    difference = phases[0]["gravitating_density_epsilon_minus_mu_q"]-phases[1]["gravitating_density_epsilon_minus_mu_q"]
    group("ground_reference_does_not_identify_present_phase", abs(difference-vacuum["delta"]) < 1e-14
          and all(p["matter_second_derivative"] > 0 and p["matter_derivative"] == 0 for p in phases)
          and abs(phases[1]["gravitating_density_epsilon_minus_mu_q"]) < 1e-14,
          {"positive_false_vacuum_gap": difference, "true_reference_density": phases[1]["gravitating_density_epsilon_minus_mu_q"]})

    portal = load("PORTAL_JOINT_GATE.json")
    # Eliminate N*lambda^2 between the width and beta; independent of multiplicity.
    allowed = portal["allowed_extra_width_GeV"]
    closed = allowed*portal["mH_GeV"]*portal["log_span"]/(math.pi*portal["vEW_GeV"]**2*portal["phase_space_factor"])
    delta_error = max(abs(row["frozen_leading_log_shift_at_max"]/closed-1) for row in portal["rows"])
    group("portal_width_beta_joint_identity", delta_error < 1e-13 and portal["hessian_trace_max_relative_error"] < 1e-13,
          {"multiplicity_independent_shift": closed, "relative_error": delta_error,
           "hessian_trace_relative_error": portal["hessian_trace_max_relative_error"]})

    tcc, uv = load("TCC_JOINT_GATE.json"), load("UV_COUNTEREXAMPLE.json")
    group("conditional_TCC_and_UV_free_direction", tcc["bisection_formula_relative_error"] < 1e-13
          and tcc["efold_gap"] > 40 and not tcc["finite_Hilbert_space_implies_TCC"]
          and max(abs(row["flow"][-1]["y"]-2) for row in uv["flows"]) < 2e-9
          and len(set(row["flow"][0]["y"] for row in uv["flows"])) == 3,
          {"TCC_analytic_vs_bisection_error": tcc["bisection_formula_relative_error"],
           "heavy_point_efold_gap": tcc["efold_gap"], "UV_regular_family_members": 3})
    save("VERIFICATION.json", {"authorship": "same-author; no external replication claimed",
                               "groups": checks, "all_pass": all(x["pass"] for x in checks)})
    print(json.dumps({"groups": len(checks), "all_pass": all(x["pass"] for x in checks),
                      "failed": [x["name"] for x in checks if not x["pass"]]}))
    if not all(x["pass"] for x in checks):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
