#!/usr/bin/env python3
"""Tests measurement bookkeeping; no additional oscillator-physics claim."""
import copy
import hashlib
import json
import tempfile
from pathlib import Path

from joint_contract import ContractError, DESIGN, build

ROOT = Path(__file__).resolve().parent


def run():
    cal = json.loads((DESIGN/"software_control_calibration.json").read_text())
    protocol = json.loads((DESIGN/"software_control_protocol.json").read_text())
    measurement = json.loads((ROOT/"software_control_measurement.json").read_text())
    checks = []

    def check(name, value):
        if not value:
            raise AssertionError(name)
        checks.append({"name": name, "pass": True})

    result = build(cal, protocol, measurement)
    check("fixture_stays_software_control", result["status"] == "SOFTWARE_CONTROL_ONLY")
    check("spin_and_energy_budget_added", abs(result["joint_false_rejection_bound_conditional"]-0.02) < 1e-15)
    check("power_miss_budget_split_over_both_contrasts", result["miss_probability_allocated_per_contrast"] == 0.025)
    check("declared_fixture_contrasts_have_conditional_power", result["target_contrasts_power_status"] ==
          "ALL_TARGET_CONTRASTS_RESOLVABLE_CONDITIONAL")
    check("no_science_or_experiment_self_certificate", result["physical_certificates_assessed_by_tool"] is False
          and result["experimental_validation_performed"] is False)

    asymmetric = copy.deepcopy(measurement)
    case = asymmetric["energy_tests"][1]
    case["baseline"].update(unresolved_tail_mean_max_phonons=0.01)
    case["driven"].update(unresolved_tail_mean_max_phonons=0.04)
    case["energy_model_bias_abs_phonons"] = 0.003
    a = build(cal, protocol, asymmetric)["energy_tests"][1]
    d_lo, d_hi = a["physical_uptake_prediction_interval_phonons"]
    lo, hi = a["reported_uptake_mean_interval_phonons"]
    check("driven_tail_lowers_capped_uptake_baseline_tail_raises_it",
          abs(lo-(d_lo-0.047)) < 1e-13 and abs(hi-(d_hi+0.017)) < 1e-13)

    wide = copy.deepcopy(measurement)
    for case in wide["energy_tests"]:
        for role in ("baseline", "driven"):
            case[role].update(estimator_min=-20, estimator_max=80)
    w = build(cal, protocol, wide)
    check("ill_conditioned_bounded_decoder_loses_certified_power",
          w["energy_contrast_power"][0]["status"] == "INSUFFICIENT_CERTIFIED_POWER")
    check("signed_estimator_range_is_explicitly_supported", w["energy_tests"][0]["groups"][0]["range"] == [-20, 80])

    tiny = copy.deepcopy(measurement)
    for case in tiny["energy_tests"]:
        case["baseline"]["shots"] = case["driven"]["shots"] = 1
    t = build(cal, protocol, tiny)
    check("energy_test_accepting_entire_support_is_flagged",
          all(x["accepts_entire_measurement_support"] for x in t["energy_tests"]))
    spin_wide = copy.deepcopy(cal)
    spin_wide["readout"]["absolute_probability_model_error"] = 1
    s = build(spin_wide, protocol, measurement)
    check("spin_test_accepting_every_possible_fraction_is_flagged",
          len(s["spin_settings_accepting_all_bright_fractions"]) == len(protocol["settings"]))
    check("vacuous_spin_calibration_cannot_pass_power", s["target_contrasts_power_status"] ==
          "INSUFFICIENT_CERTIFIED_POWER")

    def rejects(name, mutate):
        c, p, m = copy.deepcopy(cal), copy.deepcopy(protocol), copy.deepcopy(measurement)
        mutate(c, p, m)
        try:
            build(c, p, m)
        except ContractError:
            check(name, True)
        else:
            check(name, False)

    rejects("missing_energy_contract_blocks_joint_test", lambda c, p, m: m.update(energy_tests=[]))
    rejects("missing_energy_tail_moment_blocks_joint_test", lambda c, p, m:
            m["energy_tests"][0]["driven"].pop("unresolved_tail_mean_max_phonons"))
    rejects("unbounded_energy_estimator_rejected", lambda c, p, m:
            m["energy_tests"][0]["driven"].update(estimator_max=float('inf')))
    rejects("negative_tail_moment_rejected", lambda c, p, m:
            m["energy_tests"][0]["driven"].update(unresolved_tail_mean_max_phonons=-0.1))
    rejects("silent_baseline_reuse_rejected", lambda c, p, m:
            m["energy_tests"][1]["baseline"].update(id=m["energy_tests"][0]["baseline"]["id"]))
    rejects("holdout_in_prospective_contract_rejected", lambda c, p, m: m.update(outcomes=[1, 2, 3]))
    rejects("postselection_rejected", lambda c, p, m: m.update(no_postselection_assumed=False))
    rejects("dependent_trials_rejected", lambda c, p, m: m.update(independent_shots_assumed=False))
    rejects("zero_sampling_budget_rejected", lambda c, p, m: m.update(energy_sampling_failure_probability=0))
    rejects("total_probability_budget_over_one_rejected", lambda c, p, m:
            m.update(energy_sampling_failure_probability=0.999))
    rejects("zero_energy_cap_rejected", lambda c, p, m: m["energy_tests"][0].update(occupation_cap=0))
    rejects("unknown_mode_rejected", lambda c, p, m: m["energy_tests"][0].update(mode_id="missing"))
    rejects("identical_settings_cannot_be_contrast", lambda c, p, m:
            m["spin_contrasts"][0].update(setting_b=m["spin_contrasts"][0]["setting_a"]))
    rejects("purpose_cannot_upgrade_software_values", lambda c, p, m:
            m.update(purpose="INDEPENDENT_CALIBRATION_FORECAST"))

    # Only artificial byte fixtures: successful hashing does not certify the
    # estimator, tail bound, apparatus calibration or experiment.
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        blob = b"SOFTWARE CONTRACT TEST ONLY; NO ENERGY CALIBRATION\n"
        (root/"fixture.txt").write_bytes(blob)
        m = copy.deepcopy(measurement)
        m["evidence"] = [{"id": "fixture", "role": "bounded_energy_estimator", "path": "fixture.txt",
                          "sha256": hashlib.sha256(blob).hexdigest()}]
        valid = build(cal, protocol, m, measurement_root=root)
        check("energy_evidence_hash_verified_without_science_certificate",
              bool(valid["energy_evidence_byte_identity_verified"]) and
              valid["physical_certificates_assessed_by_tool"] is False)
        (root/"fixture.txt").write_bytes(blob+b"modified")
        try:
            build(cal, protocol, m, measurement_root=root)
        except ContractError:
            check("energy_evidence_hash_change_rejected", True)
        else:
            check("energy_evidence_hash_change_rejected", False)

    status = json.loads((ROOT/"platform_readiness.json").read_text())
    check("actual_platform_remains_blocked_despite_fixture_power",
          status["status"].startswith("BLOCKED") and
          status["prospective_forecast_for_actual_platform_generated"] is False)
    return {"status": "JOINT_MEASUREMENT_CONTRACT_SOFTWARE_CHECKS_PASS",
            "checks_passed": len(checks), "checks": checks,
            "new_grut_law": False, "card_2_opened": False,
            "actual_platform_status": status["status"],
            "actual_experimental_validation": "NOT_PERFORMED",
            "independent_review": "PENDING",
            "oscillator_equations_or_design_packet_changed": False}


if __name__ == "__main__":
    r = run()
    (ROOT/"results.json").write_text(json.dumps(r, indent=2)+"\n")
    print(r["status"])
    print(f"{r['checks_passed']} contract checks; actual experiment BLOCKED")
