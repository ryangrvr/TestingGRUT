#!/usr/bin/env python3
"""Joint spin/energy measurement and power contract; no observation or fit input."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

DESIGN = Path(__file__).resolve().parent.parent / "constructive_decision"
FORECAST_PATH = DESIGN / "forecast.py"
spec = importlib.util.spec_from_file_location("fixed_approach_b_forecast", FORECAST_PATH)
forecast = importlib.util.module_from_spec(spec)
spec.loader.exec_module(forecast)
ContractError, number, identifier = forecast.ContractError, forecast.number, forecast.identifier


def sha(data):
    return hashlib.sha256(data).hexdigest()


def count(value, name):
    if type(value) is not int or value <= 0:
        raise ContractError(f"{name}: require a positive integer")
    return value


def allowance(width, trials, failure, groups):
    return width * math.sqrt(math.log(2 * groups / failure) / (2 * trials))


def power(interval_a, interval_b, groups, null_allowance, beta):
    gap = max(0.0, interval_a[0]-interval_b[1], interval_b[0]-interval_a[1])
    alt = sum(allowance(w, n, beta, len(groups)) for w, n in groups)
    supported = gap > null_allowance + alt
    return {"operational_equality_null": True, "minimum_admitted_mean_separation": gap,
            "null_sampling_threshold": null_allowance,
            "alternative_sampling_allowance": alt,
            "requested_conditional_power": 1-beta,
            "status": "TARGET_CONTRAST_RESOLVABLE_CONDITIONAL" if supported else
                      "INSUFFICIENT_CERTIFIED_POWER",
            "physical_certificates_assessed_by_tool": False}


def build(calibration, protocol, measurement, calibration_root=None, measurement_root=None):
    # Frozen design bytes are checked; its equations and code are not revised.
    expected = json.loads((DESIGN / "MANIFEST.json").read_text())["sha256"]["forecast.py"]
    if sha(FORECAST_PATH.read_bytes()) != expected:
        raise ContractError("fixed forecast differs from its design manifest")
    spin = forecast.predict(calibration, protocol, calibration_root)
    if type(measurement.get("schema_version")) is not int or measurement["schema_version"] != 1:
        raise ContractError("measurement schema_version must be integer 1")
    if measurement.get("purpose") != calibration["purpose"]:
        raise ContractError("measurement and calibration purposes must match")
    if any(k in measurement for k in ("outcomes", "holdout", "observations", "fit")):
        raise ContractError("a prospective contract cannot contain holdout observations")
    if measurement.get("independent_shots_assumed") is not True or measurement.get("no_postselection_assumed") is not True:
        raise ContractError("independent trials and frozen complete trial retention required")
    mid = identifier(measurement.get("id"), "measurement.id")
    ae = number(measurement.get("energy_sampling_failure_probability"), "energy sampling failure", 0, 1)
    ad = number(measurement.get("claimed_simultaneous_energy_certificate_failure_probability"), "energy certificate failure", 0, 1)
    beta = number(measurement.get("power_miss_probability"), "power miss probability", 0, 1)
    if not 0 < ae < 1 or not 0 < beta < 1:
        raise ContractError("energy sampling and power probabilities must lie in (0,1)")
    total = spin["conditional_total_false_rejection_bound"] + ae + ad
    if total >= 1:
        raise ContractError("joint false-rejection probability budget must be less than one")

    evidence = {}
    for entry in measurement.get("evidence", []):
        eid = identifier(entry.get("id"), "energy evidence id")
        role = identifier(entry.get("role"), "energy evidence role")
        digest = identifier(entry.get("sha256"), "energy evidence sha256")
        path = identifier(entry.get("path"), "energy evidence path")
        if eid in evidence or len(digest) != 64 or any(x not in "0123456789abcdef" for x in digest):
            raise ContractError("duplicate or malformed energy evidence")
        if measurement_root is None or not (measurement_root / path).is_file() or sha((measurement_root / path).read_bytes()) != digest:
            raise ContractError("energy evidence missing or hash mismatched")
        evidence[eid] = {"role": role, "sha256": digest, "path": path}

    def source(obj, key, role):
        if calibration["purpose"] == "SOFTWARE_CONTROL_ONLY":
            return
        ref = obj.get(key)
        if ref not in evidence or evidence[ref]["role"] != role:
            raise ContractError(f"{key}: require {role} evidence")

    settings = {x["id"]: x for x in spin["settings"]}
    cases = measurement.get("energy_tests", [])
    if not cases:
        raise ContractError("joint test blocked: energy measurement contract absent")
    energy, used_groups = {}, set()
    for case in cases:
        cid = identifier(case.get("id"), "energy test id")
        sid, mode_id = case.get("setting_id"), case.get("mode_id")
        if cid in energy or sid not in settings:
            raise ContractError("duplicate energy test or unknown spin setting")
        modes = {x["id"]: x for x in settings[sid]["modes"]}
        if mode_id not in modes:
            raise ContractError("energy test references unknown mode")
        count(case.get("occupation_cap"), "occupation cap")
        source(case, "tail_model_evidence_id", "energy_tail_and_model")
        bm = number(case.get("energy_model_bias_abs_phonons"), "energy model discrepancy", 0)
        groups = []
        for key in ("baseline", "driven"):
            group = case.get(key, {})
            gid = identifier(group.get("id"), "energy group id")
            if gid in used_groups:
                raise ContractError("energy groups must be distinct; silent baseline reuse rejected")
            used_groups.add(gid)
            trials = count(group.get("shots"), "energy shots")
            lo = number(group.get("estimator_min"), "estimator lower bound")
            hi = number(group.get("estimator_max"), "estimator upper bound")
            if hi <= lo or not math.isfinite(hi-lo):
                raise ContractError("estimator requires a finite positive range width")
            bias = number(group.get("mean_readout_bias_abs_phonons"), "mean readout bias", 0)
            tail = number(group.get("unresolved_tail_mean_max_phonons"), "unresolved tail moment", 0)
            source(group, "readout_evidence_id", "bounded_energy_estimator")
            groups.append({"id": gid, "shots": trials, "range": [lo, hi], "width": hi-lo,
                           "mean_bias": bias, "tail_mean_max": tail,
                           "sampling_allowance": allowance(hi-lo, trials, ae, 2*len(cases))})
        b, a = groups
        mode = modes[mode_id]
        alpha, radius = abs(complex(*mode["alpha"])), mode["alpha_disk_radius"]
        d_lo, d_hi = max(0.0, alpha-radius)**2, (alpha+radius)**2
        bias_sum = a["mean_bias"]+b["mean_bias"]+bm
        interval = [d_lo-a["tail_mean_max"]-bias_sum, d_hi+b["tail_mean_max"]+bias_sum]
        eps = a["sampling_allowance"]+b["sampling_allowance"]
        support = [a["range"][0]-b["range"][1], a["range"][1]-b["range"][0]]
        accepted = [max(support[0], interval[0]-eps), min(support[1], interval[1]+eps)]
        if accepted[0] > accepted[1]:
            raise ContractError("energy certificates conflict with their estimator support")
        energy[cid] = {"id": cid, "setting_id": sid, "mode_id": mode_id,
                       "physical_uptake_prediction_phonons": alpha**2,
                       "physical_uptake_prediction_interval_phonons": [d_lo, d_hi],
                       "reported_uptake_mean_interval_phonons": interval,
                       "accepted_observed_estimator_uptake_phonons": accepted,
                       "measurement_support": support,
                       "accepts_entire_measurement_support": accepted == support,
                       "groups": groups}

    spin_pairs, energy_pairs = measurement.get("spin_contrasts", []), measurement.get("energy_contrasts", [])
    if not spin_pairs or not energy_pairs:
        raise ContractError("freeze at least one spin and one energy contrast for this comparison")
    beta_each = beta/(len(spin_pairs)+len(energy_pairs))
    spin_power, energy_power = [], []
    for pair in spin_pairs:
        aa, bb = pair.get("setting_a"), pair.get("setting_b")
        if aa not in settings or bb not in settings or aa == bb:
            raise ContractError("spin contrast needs two distinct known settings")
        a, b = settings[aa], settings[bb]
        p = power(a["bright_probability_interval"], b["bright_probability_interval"],
                  [(1, a["shots"]), (1, b["shots"])],
                  a["sampling_fraction_allowance"]+b["sampling_fraction_allowance"], beta_each)
        p.update(setting_a=aa, setting_b=bb)
        spin_power.append(p)
    for pair in energy_pairs:
        aa, bb = pair.get("test_a"), pair.get("test_b")
        if aa not in energy or bb not in energy or aa == bb:
            raise ContractError("energy contrast needs two distinct known tests")
        a, b = energy[aa], energy[bb]
        groups = a["groups"]+b["groups"]
        p = power(a["reported_uptake_mean_interval_phonons"], b["reported_uptake_mean_interval_phonons"],
                  [(g["width"], g["shots"]) for g in groups],
                  sum(g["sampling_allowance"] for g in groups), beta_each)
        p.update(test_a=aa, test_b=bb)
        energy_power.append(p)
    if not spin_power or not energy_power:
        raise ContractError("freeze at least one spin and one energy contrast for this comparison")
    informative = all(p["status"] == "TARGET_CONTRAST_RESOLVABLE_CONDITIONAL" for p in spin_power+energy_power)
    return {"schema_version": 1, "measurement_id": mid,
            "status": "SOFTWARE_CONTROL_ONLY" if calibration["purpose"] == "SOFTWARE_CONTROL_ONLY" else
                      "CONDITIONAL_JOINT_CONTRACT_NOT_INDEPENDENTLY_CERTIFIED",
            "new_grut_law": False, "card_2_opened": False, "experimental_validation_performed": False,
            "design_forecast_code_sha256": expected,
            "joint_contract_code_sha256": sha(Path(__file__).read_bytes()),
            "certified_interval_arithmetic": False,
            "joint_false_rejection_bound_conditional": total,
            "probability_allocations": {"spin_sampling": spin["sampling_failure_probability"],
                                        "energy_sampling": ae,
                                        "force_state_readout_calibration": spin["claimed_calibration_failure_probability"],
                                        "energy_certificate": ad},
            "energy_evidence_byte_identity_verified": evidence,
            "physical_certificates_assessed_by_tool": False,
            "spin_settings_accepting_all_bright_fractions": [s["id"] for s in spin["settings"]
                if s["accepted_observed_bright_fraction"] == [0, 1]],
            "energy_tests": list(energy.values()), "spin_contrast_power": spin_power,
            "energy_contrast_power": energy_power,
            "requested_joint_conditional_power": 1-beta,
            "miss_probability_allocated_per_contrast": beta_each,
            "target_contrasts_power_status": "ALL_TARGET_CONTRASTS_RESOLVABLE_CONDITIONAL" if informative else
                                              "INSUFFICIENT_CERTIFIED_POWER"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("calibration", "protocol", "measurement", "output"):
        parser.add_argument("--"+name, required=True, type=Path)
    args = parser.parse_args()
    try:
        inputs = {k: getattr(args, k).read_bytes() for k in ("calibration", "protocol", "measurement")}
        result = build(*(json.loads(inputs[k]) for k in ("calibration", "protocol", "measurement")),
                       args.calibration.parent, args.measurement.parent)
        result["input_sha256"] = {k: sha(v) for k, v in inputs.items()}
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
        print(result["status"])
        print(result["target_contrasts_power_status"])
    except (ContractError, KeyError, TypeError, ValueError, OverflowError, OSError) as e:
        parser.exit(2, f"JOINT_CONTRACT_REJECTED: {e}\n")


if __name__ == "__main__":
    main()
