#!/usr/bin/env python3
"""Bounded Newton stationarity check for the existing independent FV results."""
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.sparse import bmat, diags
from scipy.sparse.linalg import splu

from independent_minimize import RadialEnergy


HERE = Path(__file__).resolve().parent


def hessian_parts(model, x):
    f, s = model.unpack(x)
    n = model.n
    inertia = float(np.dot(x[:n], x[:n]))
    omega = model.charge / inertia
    diag = np.zeros(n)
    diag[:-1] += model.edge_coefficient / model.volume[:-1]
    diag[1:] += model.edge_coefficient / model.volume[1:]
    diag[-1] += model.outer_coefficient / model.volume[-1]
    off = -model.edge_coefficient / (model.sqrt_volume[:-1] * model.sqrt_volume[1:])
    lap = diags([off, diag, off], [-1, 0, 1], format="csc")
    hff = lap + diags(s * s + 3 * model.b * f * f - omega**2)
    hss = lap + diags(model.a * (3 * s * s - 1) + f * f)
    cross = diags(2 * f * s)
    local = bmat([[hff, cross], [cross, hss]], format="csc")
    z = np.concatenate((x[:n], np.zeros(n)))
    alpha = 4 * omega**2 / inertia
    return local, z, alpha


def preflight(a, b):
    model = RadialEnergy(32, 8.0, 13.0, a, b)
    x = model.pack(0.4 + 0.2 * np.exp(-model.r**2 / 9), 0.65 + 0.1 * np.exp(-model.r**2 / 16))
    direction = np.random.default_rng(20261019).normal(size=x.shape)
    direction /= np.linalg.norm(direction)
    local, z, alpha = hessian_parts(model, x)
    exact = local @ direction + alpha * z * np.dot(z, direction)
    checks = []
    for step in (1e-5, 3e-6):
        _, gp = model.evaluate(x + step * direction)
        _, gm = model.evaluate(x - step * direction)
        numerical = (gp - gm) / (2 * step)
        relative = float(np.linalg.norm(numerical - exact) / max(1.0, np.linalg.norm(exact)))
        checks.append({"step": step, "relative_Hessian_vector_discrepancy": relative})
    if min(x["relative_Hessian_vector_discrepancy"] for x in checks) > 1e-5:
        raise AssertionError(f"Hessian preflight failed: {checks}")
    return checks


def project(model, x):
    y = x.copy()
    y[:model.n] = np.maximum(y[:model.n], 0.0)
    y[model.n:] = np.clip(y[model.n:], -model.sqrt_volume, 0.0)
    return y


def diagnostics(model, x):
    e, g = model.evaluate(x)
    _, parts = model.evaluate(x, gradient=False)
    f, s = model.unpack(x)
    pg = model.projected_gradient(x, g)
    rf = g[:model.n] / model.sqrt_volume
    rs = g[model.n:] / model.sqrt_volume
    maskf = f > 1e-4
    masks = (s > 1e-4) & (s < 1 - 1e-4)
    return dict(parts, projected_gradient_inf=float(np.max(np.abs(pg))),
                full_gradient_inf=float(np.max(np.abs(g))),
                F_core_EL_residual_inf=float(np.max(np.abs(rf[maskf]))) if np.any(maskf) else None,
                S_core_EL_residual_inf=float(np.max(np.abs(rs[masks]))) if np.any(masks) else None), g


def main():
    protocol_path = HERE / "independent_refinement_protocol.json"
    protocol_bytes = protocol_path.read_bytes()
    protocol = json.loads(protocol_bytes)
    output_path = HERE / "independent_refinement_results.json"
    log_path = HERE / "independent_refinement_execution.jsonl"
    profiles_path = HERE / "independent_refinement_profiles.npz"
    if output_path.exists() or log_path.exists() or profiles_path.exists():
        raise SystemExit("Refinement outputs already exist; refusing to overwrite.")
    original_path = HERE / "independent_results.json"
    original = json.loads(original_path.read_text())
    profiles = np.load(HERE / "independent_profiles.npz")
    a, b = original["protocol"]["a"], original["protocol"]["b"]
    checks = preflight(a, b)
    result = {"protocol": protocol,
              "protocol_sha256": hashlib.sha256(protocol_bytes).hexdigest(),
              "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "original_results_sha256": hashlib.sha256(original_path.read_bytes()).hexdigest(),
              "Hessian_preflight": checks, "runs": []}
    saved_profiles = {}
    started = time.perf_counter()
    for fraction in protocol["charges"]:
        label = f"q{fraction:.2f}_N{protocol['N']}_independent_gaussian"
        model = RadialEnergy(protocol["N"], protocol["L"], fraction * original["protocol"]["target_Qbar"], a, b)
        x = model.pack(profiles[label + "_F"], profiles[label + "_S"])
        row = {"label": label, "fraction": fraction, "Qbar": model.charge, "steps": []}
        current, g = diagnostics(model, x)
        row["initial"] = current
        status = "step_limit"
        for step_index in range(protocol["max_newton_steps"]):
            if current["projected_gradient_inf"] <= protocol["requested_projected_gradient_inf"]:
                status = "requested_gradient_reached"
                break
            if time.perf_counter() - started > protocol["wall_time_limit_total_seconds"]:
                status = "total_wall_limit"
                break
            local, z, alpha = hessian_parts(model, x)
            try:
                lu = splu(local)
                solved = lu.solve(np.column_stack((-g, z)))
                v, w = solved[:, 0], solved[:, 1]
                denominator = 1 + alpha * np.dot(z, w)
                if abs(denominator) < 1e-12:
                    raise ArithmeticError("Near-singular Sherman-Morrison denominator")
                delta = v - (alpha * np.dot(z, v) / denominator) * w
            except (RuntimeError, ArithmeticError) as exc:
                status = "linear_solve_failure"
                row["failure"] = str(exc)
                break
            linear_residual = local @ delta + alpha * z * np.dot(z, delta) + g
            step = {"step": step_index + 1, "Sherman_Morrison_denominator": float(denominator),
                    "linear_residual_inf": float(np.max(np.abs(linear_residual))), "trials": []}
            accepted = False
            for k in range(11):
                scale = 2.0**(-k)
                raw = x + scale * delta
                trial = project(model, raw)
                metrics, trial_g = diagnostics(model, trial)
                step["trials"].append({"scale": scale, "projection_inf": float(np.max(np.abs(trial - raw))),
                                        **metrics})
                if metrics["projected_gradient_inf"] < current["projected_gradient_inf"] or metrics["projected_gradient_inf"] <= protocol["requested_projected_gradient_inf"]:
                    x, current, g = trial, metrics, trial_g
                    step["accepted_scale"] = scale
                    accepted = True
                    break
            row["steps"].append(step)
            if not accepted:
                status = "backtracking_failed_to_reduce_gradient"
                break
        if current["projected_gradient_inf"] <= protocol["requested_projected_gradient_inf"]:
            status = "requested_gradient_reached"
        row["status"] = status
        row["final"] = current
        f, s = model.unpack(x)
        saved_profiles[label + "_r"] = model.r
        saved_profiles[label + "_F"] = f
        saved_profiles[label + "_S"] = s
        result["runs"].append(row)
        result["total_wall_seconds"] = time.perf_counter() - started
        output_path.write_text(json.dumps(result, indent=2) + "\n")
        np.savez_compressed(profiles_path, **saved_profiles)
        with log_path.open("a") as handle:
            handle.write(json.dumps(row) + "\n")
        print(json.dumps({"label": label, "status": status, **current}), flush=True)


if __name__ == "__main__":
    main()
