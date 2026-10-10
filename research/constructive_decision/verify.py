#!/usr/bin/env python3
"""Author software checks of a known model, not an independent science review."""
import copy
import hashlib
import json
import math
import platform
import tempfile
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import expm

from forecast import ContractError, predict

ROOT = Path(__file__).resolve().parent


def run() -> dict:
    cal = json.loads((ROOT / "software_control_calibration.json").read_text())
    protocol = json.loads((ROOT / "software_control_protocol.json").read_text())
    checks = []

    def checked(name, condition, discrepancy=None):
        if not condition:
            raise AssertionError(name)
        checks.append({"name": name, "pass": True, "discrepancy": discrepancy})

    output = predict(cal, protocol)
    for inp, out in zip(protocol["settings"], output["settings"]):
        mode = cal["modes"][0]
        om, n = mode["detuning_rad_s"], mode["occupation"]
        g = mode["drives"][0]["real_rad_s"]
        t = sum(s["duration_s"] for s in inp["segments"])
        expected = (4*(g/om)**2*(2*n+1)*(1-math.cos(om*t)) if
                    len(inp["segments"]) == 1 else
                    32*(g/om)**2*(2*n+1)*math.sin(om*t/4)**4)
        checked(f"closed_form_{inp['id']}", abs(out["gamma"]-expected) < 1e-12,
                abs(out["gamma"]-expected))

    # Independent calculation: exponentiate the conditional oscillator
    # Hamiltonians Omega*N +/- f*(a+a†), not the displacement formula.
    discrepancies = []
    for inp, out in zip(protocol["settings"], output["settings"]):
        dimension_results = []
        for dim in (24, 40):
            a = np.diag(np.sqrt(np.arange(1, dim)), 1)
            num = a.T @ a
            n = cal["modes"][0]["occupation"]
            weights = (n/(n+1))**np.arange(dim)/(n+1)
            rho = np.diag(weights/weights.sum())
            branches = []
            for z in (1, -1):
                u = np.eye(dim, dtype=complex)
                for seg in inp["segments"]:
                    sign = seg["drives"][0]["sign"]
                    force = sign * cal["modes"][0]["drives"][0]["real_rad_s"]
                    h = cal["modes"][0]["detuning_rad_s"]*num + z*force*(a.T+a)
                    u = expm(-1j*h*seg["duration_s"]) @ u
                branches.append(u)
            coherence = np.trace(branches[0] @ rho @ branches[1].conj().T)
            uptake = np.trace(num @ (branches[0] @ rho @ branches[0].conj().T-rho)).real
            dimension_results.append((coherence, uptake))
        high, low = dimension_results[1], dimension_results[0]
        discrepancy = max(abs(high[0]-out["visibility"]),
                          abs(high[1]-abs(complex(*out["modes"][0]["alpha"]))**2))
        discrepancies.append(float(discrepancy))
        checked(f"independent_fock_evolution_{inp['id']}", discrepancy < 2e-10,
                float(discrepancy))
        checked(f"fock_cutoff_comparison_{inp['id']}",
                max(abs(high[0]-low[0]), abs(high[1]-low[1])) < 2e-10)

    # Resonance must be handled by the stable integral rather than g/Omega.
    resonant = copy.deepcopy(cal)
    resonant["modes"][0]["detuning_rad_s"] = 0.0
    rz = predict(resonant, protocol)
    for inp, out in zip(protocol["settings"], rz["settings"]):
        signed_duration = sum(seg["duration_s"]*seg["drives"][0]["sign"] for seg in inp["segments"])
        alpha_sq = (cal["modes"][0]["drives"][0]["real_rad_s"]*signed_duration)**2
        expected = 2*(2*cal["modes"][0]["occupation"]+1)*alpha_sq
        checked(f"zero_detuning_{inp['id']}", abs(out["gamma"]-expected) < 1e-12)

    # Stress-check the analytic whole-domain error formula, not a convergence
    # certificate or evidence that the physical uncertainty claims are true.
    uncertain = copy.deepcopy(cal)
    m = uncertain["modes"][0]
    m["detuning_abs_error_rad_s"] = 0.04
    m["occupation_abs_error"] = 0.08
    m["physical_frequency_abs_error_rad_s"] = 0.01*m["physical_frequency_rad_s"]
    m["drives"][0]["absolute_error_rad_s"] = 0.01
    up = copy.deepcopy(protocol)
    for s in up["settings"]:
        for seg in s["segments"]:
            seg["duration_abs_error_s"] = 0.001
    bounded = predict(uncertain, up)
    rng = np.random.default_rng(20261010)
    all_inside = True
    for _ in range(64):
        sample_c, sample_p = copy.deepcopy(uncertain), copy.deepcopy(up)
        mode = sample_c["modes"][0]
        mode["detuning_rad_s"] += rng.uniform(-0.04, 0.04)
        mode["occupation"] += rng.uniform(-0.08, 0.08)
        mode["physical_frequency_rad_s"] *= 1+rng.uniform(-0.01, 0.01)
        r, phase = 0.01*math.sqrt(rng.uniform()), rng.uniform(0, math.tau)
        mode["drives"][0]["real_rad_s"] += r*math.cos(phase)
        mode["drives"][0]["imag_rad_s"] += r*math.sin(phase)
        for s in sample_p["settings"]:
            for seg in s["segments"]:
                seg["duration_s"] += rng.uniform(-0.001, 0.001)
        actual = predict(sample_c, sample_p)
        for observed, enclosure in zip(actual["settings"], bounded["settings"]):
            for key, interval in (("gamma", "gamma_interval"),
                                  ("energy_uptake_joule", "energy_uptake_interval_joule")):
                all_inside &= enclosure[interval][0] <= observed[key] <= enclosure[interval][1]
    checked("64_uncertainty_stress_draws_inside_analytic_enclosure", all_inside)

    # Two modes with different detunings and complex pulse phases test the
    # implemented domain, including coherent addition across segments.
    multi = copy.deepcopy(cal)
    first = multi["modes"][0]
    first["drives"][0]["imag_rad_s"] = 0.07
    first["drives"].append({"id": "force_2", "real_rad_s": 0.13,
                            "imag_rad_s": -0.08, "absolute_error_rad_s": 0.0})
    second = copy.deepcopy(first)
    second.update(id="mode_2", detuning_rad_s=-1.3, occupation=0.05)
    multi["modes"].append(second)
    mp = copy.deepcopy(protocol)
    mp["settings"] = [{"id": "complex_two_mode_control", "shots": 3000, "segments": [
        {"duration_s": dt, "duration_abs_error_s": 0.0, "drives": [
            {"mode_id": "mode_1", "drive_id": did, "sign": sign},
            {"mode_id": "mode_2", "drive_id": did, "sign": -sign}]}
        for dt, did, sign in ((0.4, "force_1", 1), (0.2, "force_2", 1), (0.7, "force_1", -1))]}]
    pred = predict(multi, mp)["settings"][0]
    full_coherence = 1.0+0j
    dimension = 40
    a = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    num = a.T @ a
    max_energy_error = 0.0
    for mode, mode_pred in zip(multi["modes"], pred["modes"]):
        n = mode["occupation"]
        weights = (n/(n+1))**np.arange(dimension)/(n+1)
        rho = np.diag(weights/weights.sum())
        branches = []
        for z in (1, -1):
            u = np.eye(dimension, dtype=complex)
            for seg in mp["settings"][0]["segments"]:
                command = next(x for x in seg["drives"] if x["mode_id"] == mode["id"])
                drive = next(x for x in mode["drives"] if x["id"] == command["drive_id"])
                force = command["sign"]*complex(drive["real_rad_s"], drive["imag_rad_s"])
                h = mode["detuning_rad_s"]*num + z*(force*a.T+force.conjugate()*a)
                u = expm(-1j*h*seg["duration_s"]) @ u
            branches.append(u)
        full_coherence *= np.trace(branches[0] @ rho @ branches[1].conj().T)
        energy_n = np.trace(num @ (branches[0] @ rho @ branches[0].conj().T-rho)).real
        max_energy_error = max(max_energy_error,
                               abs(energy_n-abs(complex(*mode_pred["alpha"]))**2))
    error = float(max(abs(full_coherence-pred["visibility"]), max_energy_error))
    discrepancies.append(error)
    checked("complex_two_mode_independent_fock_evolution", error < 2e-10, error)

    def rejects(name, mutate):
        c, p = copy.deepcopy(cal), copy.deepcopy(protocol)
        mutate(c, p)
        try:
            predict(c, p)
        except ContractError:
            checked(name, True)
        else:
            checked(name, False)

    rejects("nonfinite_force_rejected", lambda c, p: c["modes"][0]["drives"][0].update(real_rad_s=float('nan')))
    rejects("negative_occupation_rejected", lambda c, p: c["modes"][0].update(occupation=-0.1))
    rejects("nonpositive_physical_frequency_rejected", lambda c, p: c["modes"][0].update(physical_frequency_rad_s=0))
    rejects("duplicate_mode_rejected", lambda c, p: c["modes"].append(copy.deepcopy(c["modes"][0])))
    rejects("unknown_drive_rejected", lambda c, p: p["settings"][0]["segments"][0]["drives"][0].update(drive_id="unmeasured"))
    rejects("unbounded_timing_interval_rejected", lambda c, p: p["settings"][0]["segments"][0].update(duration_abs_error_s=99))
    rejects("holdout_in_forecast_input_rejected", lambda c, p: p.update(observations=[0.9]))
    rejects("independent_forecast_without_evidence_rejected", lambda c, p: c.update(purpose="INDEPENDENT_CALIBRATION_FORECAST"))
    rejects("readout_outside_probability_range_rejected", lambda c, p: c["readout"].update(bright_given_plus_x=1.1))
    rejects("dependent_shots_rejected", lambda c, p: p.update(independent_shots_assumed=False))

    # Provenance bytes must match. This is manufactured evidence for a software
    # contract test only; it is never a scientific calibration certificate.
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        blob = b"SOFTWARE CONTRACT FIXTURE: NOT LABORATORY CALIBRATION\n"
        (root/"fixture.txt").write_bytes(blob)
        cc = copy.deepcopy(cal)
        cc["evidence"] = [{"id": "fixture", "role": "force", "path": "fixture.txt",
                           "sha256": hashlib.sha256(blob).hexdigest()}]
        checked("evidence_hash_validated", bool(predict(cc, protocol, root)["evidence_byte_identity_verified"]))
        forecast_fixture = copy.deepcopy(cal)
        forecast_fixture["purpose"] = "INDEPENDENT_CALIBRATION_FORECAST"
        roles = ("frequency", "occupation_and_state", "force", "readout", "model_error", "clock_and_waveform")
        forecast_fixture["evidence"] = [{"id": role, "role": role, "path": "fixture.txt",
                                        "sha256": hashlib.sha256(blob).hexdigest()} for role in roles]
        forecast_fixture["clock_evidence_id"] = "clock_and_waveform"
        forecast_fixture["modes"][0].update(frequency_evidence_id="frequency",
                                           occupation_evidence_id="occupation_and_state")
        forecast_fixture["modes"][0]["drives"][0]["force_evidence_id"] = "force"
        forecast_fixture["readout"].update(readout_evidence_id="readout", model_error_evidence_id="model_error")
        result = predict(forecast_fixture, protocol, root)
        checked("conditional_forecast_contract_complete_fixture", result["status"] ==
                "CONDITIONAL_FORECAST_NOT_EXPERIMENTALLY_CERTIFIED")
        checked("byte_verification_does_not_self_certify_science", result["evidence_scientifically_certified_by_this_tool"] is False)
        (root/"fixture.txt").write_bytes(blob+b"changed")
        try:
            predict(cc, protocol, root)
        except ContractError:
            checked("evidence_hash_mismatch_rejected", True)
        else:
            checked("evidence_hash_mismatch_rejected", False)

    return {"status": "CONVENTIONAL_EFFECTIVE_FORECAST_SOFTWARE_CONTROLS_PASS",
            "new_grut_law": False, "independent_scientific_review": "PENDING",
            "experimental_calibration": "ABSENT", "card_2_opened": False,
            "checks_passed": len(checks), "checks": checks,
            "independent_fock_max_discrepancy": max(discrepancies),
            "fock_dimensions": [24, 40], "uncertainty_stress_draws": 64,
            "stress_draws_are_not_the_error_certificate": True,
            "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__}


if __name__ == "__main__":
    result = run()
    (ROOT/"results.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(result["status"])
    print(f"{result['checks_passed']} checks; Fock discrepancy {result['independent_fock_max_discrepancy']:.3g}")
