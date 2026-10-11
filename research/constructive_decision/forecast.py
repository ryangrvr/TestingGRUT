#!/usr/bin/env python3
"""Prospective conventional force/thermal-state forecast; never fits holdouts."""
from __future__ import annotations

import argparse
import cmath
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any

HBAR = 6.62607015e-34 / math.tau
ROUNDING_REL_ALLOWANCE = 1e-12  # disclosed allowance, not certified interval arithmetic


class ContractError(ValueError):
    pass


def number(value: Any, name: str, lower: float | None = None,
           upper: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError(f"{name}: expected a finite number")
    value = float(value)
    if not math.isfinite(value) or (lower is not None and value < lower) or (
            upper is not None and value > upper):
        raise ContractError(f"{name}: outside its finite admissible interval")
    return value


def identifier(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ContractError(f"{name}: expected a nonempty identifier")
    return value


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def integral(omega: float, a: float, b: float) -> complex:
    """Exact integral exp(i omega t), evaluated stably at omega=0."""
    d = b - a
    x = omega * d / 2
    sinc = 1.0 if x == 0 else math.sin(x) / x
    return d * sinc * cmath.exp(1j * omega * (a + b) / 2)


def predict(calibration: dict, protocol: dict,
            evidence_root: Path | None = None) -> dict:
    if any(type(obj.get("schema_version")) is not int or obj["schema_version"] != 1
           for obj in (calibration, protocol)):
        raise ContractError("both inputs require schema_version=1")
    purpose = calibration.get("purpose")
    if purpose not in ("SOFTWARE_CONTROL_ONLY", "INDEPENDENT_CALIBRATION_FORECAST"):
        raise ContractError("declare the calibration purpose explicitly")
    if any(key in obj for obj in (calibration, protocol)
           for key in ("observations", "holdout", "fit", "outcomes")):
        raise ContractError("forecast inputs cannot contain holdout observations or fits")
    calibration_id = identifier(calibration.get("id"), "calibration.id")
    protocol_id = identifier(protocol.get("id"), "protocol.id")

    # Local byte identity verifies evidence provenance, not its scientific adequacy.
    provenance = {}
    for entry in calibration.get("evidence", []):
        eid = identifier(entry.get("id"), "evidence.id")
        if eid in provenance:
            raise ContractError(f"duplicate evidence id: {eid}")
        role = identifier(entry.get("role"), "evidence.role")
        path_text = identifier(entry.get("path"), "evidence.path")
        digest = identifier(entry.get("sha256"), "evidence.sha256")
        if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ContractError("evidence sha256 must be 64 lowercase hexadecimal digits")
        if evidence_root is None:
            raise ContractError("evidence_root required for provenance verification")
        path = (evidence_root / path_text).resolve()
        if not path.is_file() or sha256_bytes(path.read_bytes()) != digest:
            raise ContractError(f"missing or hash-mismatched evidence: {eid}")
        provenance[eid] = {"role": role, "sha256": digest, "path": path_text}

    def source(obj: dict, key: str, role: str) -> None:
        if purpose == "SOFTWARE_CONTROL_ONLY":
            return
        sid = obj.get(key)
        if sid not in provenance or provenance[sid]["role"] != role:
            raise ContractError(f"{key}: require independently supplied {role} evidence")

    modes = {}
    for mode in calibration.get("modes", []):
        mid = identifier(mode.get("id"), "mode.id")
        if mid in modes:
            raise ContractError(f"duplicate mode id: {mid}")
        w = number(mode.get("physical_frequency_rad_s"), "physical frequency", 0)
        ew = number(mode.get("physical_frequency_abs_error_rad_s"), "frequency error", 0)
        if w <= ew:
            raise ContractError("physical mode frequency interval must be strictly positive")
        om = number(mode.get("detuning_rad_s"), "detuning")
        eo = number(mode.get("detuning_abs_error_rad_s"), "detuning error", 0)
        n = number(mode.get("occupation"), "occupation", 0)
        en = number(mode.get("occupation_abs_error"), "occupation error", 0)
        source(mode, "frequency_evidence_id", "frequency")
        source(mode, "occupation_evidence_id", "occupation_and_state")
        drives = {}
        for drive in mode.get("drives", []):
            did = identifier(drive.get("id"), "drive.id")
            if did in drives:
                raise ContractError(f"duplicate drive id in mode {mid}: {did}")
            f = complex(number(drive.get("real_rad_s"), "force real"),
                        number(drive.get("imag_rad_s"), "force imaginary"))
            ef = number(drive.get("absolute_error_rad_s"), "force error", 0)
            source(drive, "force_evidence_id", "force")
            drives[did] = (f, ef)
        if not drives:
            raise ContractError(f"mode {mid} has no independently calibrated drive")
        modes[mid] = (w, ew, om, eo, n, en, drives)
    if not modes:
        raise ContractError("a nonempty finite frozen mode inventory is required")

    if calibration.get("state_model") != "CENTERED_PRODUCT_ISOTROPIC_GAUSSIAN":
        raise ContractError("the declared Gaussian product state model is required")
    if calibration.get("initial_system_environment_product") is not True:
        raise ContractError("initial system/environment independence must be explicit")
    readout = calibration.get("readout", {})
    pp = number(readout.get("bright_given_plus_x"), "plus readout", 0, 1)
    ep = number(readout.get("plus_abs_error"), "plus readout error", 0, 1)
    pm = number(readout.get("bright_given_minus_x"), "minus readout", 0, 1)
    em = number(readout.get("minus_abs_error"), "minus readout error", 0, 1)
    model_error = number(readout.get("absolute_probability_model_error"), "model error", 0, 1)
    source(readout, "readout_evidence_id", "readout")
    source(readout, "model_error_evidence_id", "model_error")
    source(calibration, "clock_evidence_id", "clock_and_waveform")
    ac = number(calibration.get("claimed_simultaneous_calibration_failure_probability"),
                "calibration failure probability", 0, 1)
    ass = number(protocol.get("sampling_failure_probability"), "sampling failure probability", 0, 1)
    if ass <= 0 or ass + ac >= 1:
        raise ContractError("require 0 < sampling failure and total failure < 1")
    if protocol.get("independent_shots_assumed") is not True:
        raise ContractError("independent shot assumption must be explicit")
    settings = protocol.get("settings", [])
    if not settings:
        raise ContractError("freeze a nonempty finite setting chart")
    result_settings, seen = [], set()
    for setting in settings:
        sid = identifier(setting.get("id"), "setting.id")
        if sid in seen:
            raise ContractError(f"duplicate setting id: {sid}")
        seen.add(sid)
        shots = setting.get("shots")
        if isinstance(shots, bool) or not isinstance(shots, int) or shots <= 0:
            raise ContractError("shots must be a positive integer")
        segments = setting.get("segments", [])
        if not segments:
            raise ContractError("each setting requires at least one measured segment")
        alphas = {mid: 0j for mid in modes}
        radii = {mid: 0.0 for mid in modes}
        absolute_integrals = {mid: 0.0 for mid in modes}
        a, endpoint_error = 0.0, 0.0
        for segment in segments:
            dt = number(segment.get("duration_s"), "duration", 0)
            edt = number(segment.get("duration_abs_error_s"), "duration error", 0)
            if dt <= edt:
                raise ContractError("duration interval must be strictly positive")
            b, eb = a + dt, endpoint_error + edt
            if not math.isfinite(b) or not math.isfinite(eb):
                raise ContractError("timing overflow")
            used = set()
            for command in segment.get("drives", []):
                mid, did = command.get("mode_id"), command.get("drive_id")
                if mid not in modes or did not in modes[mid][6] or mid in used:
                    raise ContractError("unknown drive/mode or duplicate command for a mode")
                used.add(mid)
                sign = command.get("sign")
                if isinstance(sign, bool) or sign not in (-1, 1):
                    raise ContractError("drive sign must be +1 or -1")
                _, _, om, eo, _, _, drives = modes[mid]
                f, ef = drives[did]
                f *= sign
                alphas[mid] += -1j * f * integral(om, a, b)
                radii[mid] += (ef * (dt + edt) + abs(f) * (
                    endpoint_error + eb + min(2 * dt, eo * dt * (a + b) / 2)))
                absolute_integrals[mid] += abs(f) * dt
            a, endpoint_error = b, eb
        gamma, glo, ghi, energy, elo, ehi = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        mode_predictions = []
        for mid, (w, ew, _, _, n, en, _) in modes.items():
            alpha = alphas[mid]
            # Accounts for cancellation sensitivity; not a formal libm certificate.
            radius = radii[mid] + ROUNDING_REL_ALLOWANCE * max(1.0, absolute_integrals[mid])
            amin, amax = max(0.0, abs(alpha) - radius), abs(alpha) + radius
            lo, hi = 2 * (2 * max(0.0, n - en) + 1) * amin**2, 2 * (2 * (n + en) + 1) * amax**2
            gj = 2 * (2 * n + 1) * abs(alpha)**2
            ej = HBAR * w * abs(alpha)**2
            ejlo, ejhi = HBAR * (w - ew) * amin**2, HBAR * (w + ew) * amax**2
            gamma += gj
            glo += lo
            ghi += hi
            energy += ej
            elo += ejlo
            ehi += ejhi
            mode_predictions.append({"id": mid, "alpha": [alpha.real, alpha.imag],
                                     "alpha_disk_radius": radius,
                                     "energy_uptake_joule": ej,
                                     "energy_uptake_interval_joule": [ejlo, ejhi]})
        if not all(math.isfinite(x) for x in (gamma, glo, ghi, energy, elo, ehi)):
            raise ContractError("forecast overflow: reduce the supported numeric range")
        v, vlo, vhi = math.exp(-gamma), math.exp(-ghi), math.exp(-glo)
        corners = [0.5 * (x * (1 + vv) + y * (1 - vv)) for x, y, vv in itertools.product(
            (max(0, pp - ep), min(1, pp + ep)),
            (max(0, pm - em), min(1, pm + em)), (vlo, vhi))]
        p = 0.5 * (pp * (1 + v) + pm * (1 - v))
        plo, phi = max(0.0, min(corners) - model_error), min(1.0, max(corners) + model_error)
        eps = math.sqrt(math.log(2 * len(settings) / ass) / (2 * shots))
        result_settings.append({"id": sid, "duration_s": a, "shots": shots,
                                "modes": mode_predictions, "gamma": gamma,
                                "gamma_interval": [glo, ghi], "visibility": v,
                                "visibility_interval": [vlo, vhi],
                                "energy_uptake_joule": energy,
                                "energy_uptake_interval_joule": [elo, ehi],
                                "bright_probability": p,
                                "bright_probability_interval": [plo, phi],
                                "sampling_fraction_allowance": eps,
                                "accepted_observed_bright_fraction": [max(0, plo-eps), min(1, phi+eps)]})
    return {"schema_version": 1, "calibration_id": calibration_id, "protocol_id": protocol_id,
            "status": ("SOFTWARE_CONTROL_ONLY" if purpose == "SOFTWARE_CONTROL_ONLY" else
                       "CONDITIONAL_FORECAST_NOT_EXPERIMENTALLY_CERTIFIED"),
            "new_grut_law": False, "card_2_opened": False,
            "sampling_failure_probability": ass,
            "claimed_calibration_failure_probability": ac,
            "conditional_total_false_rejection_bound": ass + ac,
            "rounding_relative_allowance": ROUNDING_REL_ALLOWANCE,
            "certified_interval_arithmetic": False,
            "evidence_byte_identity_verified": provenance,
            "evidence_scientifically_certified_by_this_tool": False,
            "settings": result_settings}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calibration", required=True, type=Path)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        cb, pb = args.calibration.read_bytes(), args.protocol.read_bytes()
        prediction = predict(json.loads(cb), json.loads(pb), args.calibration.parent)
        prediction["input_sha256"] = {"calibration": sha256_bytes(cb), "protocol": sha256_bytes(pb)}
        prediction["forecast_code_sha256"] = sha256_bytes(Path(__file__).read_bytes())
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(prediction, indent=2, allow_nan=False) + "\n")
        print(prediction["status"])
    except (ContractError, KeyError, TypeError, OverflowError, OSError, ValueError) as exc:
        parser.exit(2, f"FORECAST_REJECTED: {exc}\n")


if __name__ == "__main__":
    main()
