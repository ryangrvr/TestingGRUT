#!/usr/bin/env python3
"""Local fixed-Q finite-volume energy minimization; no stationary ODE solver.

Read independent_protocol.json before interpreting outputs. This implementation
does not import the existing BVP or branch-scanning implementation.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import minimize


HERE = Path(__file__).resolve().parent


class RuntimeLimit(Exception):
    pass


class RadialEnergy:
    """Cell-centered shell-volume functional and its exact discrete gradient."""

    def __init__(self, n: int, length: float, charge: float, a: float, b: float):
        self.n = n
        self.length = length
        self.h = length / n
        self.charge = charge
        self.a = a
        self.b = b
        edges = np.arange(n + 1, dtype=float) * self.h
        self.r = (np.arange(n, dtype=float) + 0.5) * self.h
        self.volume = (4.0 * np.pi / 3.0) * np.diff(edges**3)
        self.sqrt_volume = np.sqrt(self.volume)
        self.edge_coefficient = 4.0 * np.pi * edges[1:-1] ** 2 / self.h
        self.outer_coefficient = 4.0 * np.pi * length**2 / (0.5 * self.h)
        self.calls = 0

    def pack(self, f: np.ndarray, s: np.ndarray) -> np.ndarray:
        return np.concatenate((self.sqrt_volume * f, self.sqrt_volume * (s - 1.0)))

    def unpack(self, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return x[: self.n] / self.sqrt_volume, 1.0 + x[self.n :] / self.sqrt_volume

    def bounds(self) -> list[tuple[float | None, float | None]]:
        return [(0.0, None)] * self.n + [(-float(z), 0.0) for z in self.sqrt_volume]

    def evaluate(self, x: np.ndarray, gradient: bool = True):
        self.calls += 1
        f, s = self.unpack(x)
        inertia = float(np.dot(x[: self.n], x[: self.n]))
        if inertia <= 1e-200:
            raise FloatingPointError("Vanishing inertia reached; fixed-Q energy is singular.")
        omega = self.charge / inertia
        kinetic = self.charge**2 / (2.0 * inertia)

        df = np.diff(f)
        ds = np.diff(s)
        grad_energy = 0.5 * float(np.dot(self.edge_coefficient, df * df + ds * ds))
        grad_energy += 0.5 * self.outer_coefficient * (f[-1] ** 2 + (s[-1] - 1.0) ** 2)
        potential_density = (
            0.25 * self.a * (s * s - 1.0) ** 2
            + 0.5 * s * s * f * f
            + 0.25 * self.b * f**4
        )
        potential = float(np.dot(self.volume, potential_density))
        energy = kinetic + grad_energy + potential
        if not gradient:
            return energy, {
                "Ibar": inertia,
                "Omega": omega,
                "Ebar": energy,
                "Etime_bar": kinetic,
                "Egrad_bar": grad_energy,
                "EU_bar": potential,
                "Ebar_over_Qbar": energy / self.charge,
            }

        # Differentiate the face energies with respect to the unweighted fields.
        gf = np.zeros(self.n)
        gs = np.zeros(self.n)
        edge_f = self.edge_coefficient * df
        edge_s = self.edge_coefficient * ds
        gf[:-1] -= edge_f
        gf[1:] += edge_f
        gs[:-1] -= edge_s
        gs[1:] += edge_s
        gf[-1] += self.outer_coefficient * f[-1]
        gs[-1] += self.outer_coefficient * (s[-1] - 1.0)

        # Potential and charge terms, transformed to u=sqrt(V)F,w=sqrt(V)(S-1).
        gu = gf / self.sqrt_volume
        gu += self.sqrt_volume * (s * s * f + self.b * f**3)
        gu -= omega**2 * x[: self.n]
        gw = gs / self.sqrt_volume
        gw += self.sqrt_volume * (self.a * (s * s - 1.0) * s + s * f * f)
        return energy, np.concatenate((gu, gw))

    def projected_gradient(self, x: np.ndarray, g: np.ndarray) -> np.ndarray:
        pg = g.copy()
        at_lower_u = x[: self.n] <= 1e-12
        pg[: self.n][at_lower_u & (g[: self.n] > 0.0)] = 0.0
        w = x[self.n :]
        at_lower_w = w <= -self.sqrt_volume + 1e-12 * np.maximum(1.0, self.sqrt_volume)
        at_upper_w = w >= -1e-12
        pg[self.n :][at_lower_w & (g[self.n :] > 0.0)] = 0.0
        pg[self.n :][at_upper_w & (g[self.n :] < 0.0)] = 0.0
        return pg


def gradient_preflight(a: float, b: float) -> list[dict]:
    model = RadialEnergy(32, 8.0, 13.0, a, b)
    f = 0.4 + 0.2 * np.exp(-model.r**2 / 9.0)
    s = 0.65 + 0.1 * np.exp(-model.r**2 / 16.0)
    x = model.pack(f, s)
    _, g = model.evaluate(x)
    rng = np.random.default_rng(20261018)
    rows = []
    for direction_index in range(4):
        direction = rng.normal(size=x.shape)
        direction /= np.linalg.norm(direction)
        exact = float(np.dot(g, direction))
        discrepancies = []
        for step in (1e-5, 3e-6):
            ep, _ = model.evaluate(x + step * direction, gradient=False)
            em, _ = model.evaluate(x - step * direction, gradient=False)
            numerical = (ep - em) / (2.0 * step)
            relative = abs(numerical - exact) / max(1.0, abs(exact))
            discrepancies.append(relative)
            rows.append({
                "direction": direction_index,
                "step": step,
                "analytic": exact,
                "finite_difference": numerical,
                "relative_discrepancy": relative,
            })
        if min(discrepancies) > 1e-5:
            raise AssertionError(f"Gradient preflight failed in direction {direction_index}: {rows[-2:]}")
    return rows


def initial_fields(model: RadialEnergy, name: str, fraction: float, source) -> tuple[np.ndarray, np.ndarray]:
    scale = fraction ** (1.0 / 3.0)
    if name == "bvp_shape":
        r = model.r / scale
        f = np.interp(r, source["rho"], source["y"][0], left=float(source["y"][0, 0]), right=0.0)
        s = np.interp(r, source["rho"], source["y"][2], left=float(source["y"][2, 0]), right=1.0)
    elif name == "independent_gaussian":
        sigma = 4.0 * scale
        f = 1.25 * np.exp(-model.r**2 / (2.0 * sigma**2))
        s = 1.0 - 0.8 * np.exp(-model.r**2 / (2.0 * (1.1 * sigma) ** 2))
    else:
        raise ValueError(name)
    return np.maximum(f, 0.0), np.clip(s, 0.0, 1.0)


def run_minimization(model: RadialEnergy, x0: np.ndarray, options: dict, limit: float) -> tuple[dict, np.ndarray]:
    started = time.perf_counter()
    last_x = x0.copy()
    last_iteration = 0

    def callback(x):
        nonlocal last_x, last_iteration
        last_x = x.copy()
        last_iteration += 1
        if time.perf_counter() - started > limit:
            raise RuntimeLimit(f"Per-run wall limit {limit}s exceeded.")

    initial_energy, initial_parts = model.evaluate(x0, gradient=False)
    try:
        result = minimize(
            model.evaluate, x0, method="L-BFGS-B", jac=True,
            bounds=model.bounds(), callback=callback, options=options,
        )
        x = result.x.copy()
        outcome = {
            "optimizer_success": bool(result.success),
            "optimizer_status": int(result.status),
            "message": str(result.message),
            "iterations": int(result.nit),
            "function_evaluations_reported": int(result.nfev),
        }
    except (RuntimeLimit, FloatingPointError) as exc:
        x = last_x.copy()
        outcome = {
            "optimizer_success": False,
            "optimizer_status": "runtime_limit" if isinstance(exc, RuntimeLimit) else "floating_point_failure",
            "message": str(exc),
            "iterations": last_iteration,
        }

    energy, gradient = model.evaluate(x)
    _, parts = model.evaluate(x, gradient=False)
    f, s = model.unpack(x)
    pg = model.projected_gradient(x, gradient)
    projected_inf = float(np.max(np.abs(pg)))
    # An unweighted Euler residual is informative on populated interior cells;
    # tiny tail cells and explicitly active bounds need a separate KKT reading.
    rf = gradient[: model.n] / model.sqrt_volume
    rs = gradient[model.n :] / model.sqrt_volume
    core_f = f > 1e-4
    core_s = (s > 1e-4) & (s < 1.0 - 1e-4)
    outcome.update(parts)
    outcome.update({
        "initial_Ebar": initial_energy,
        "initial_parts": initial_parts,
        "energy_change": energy - initial_energy,
        "projected_gradient_inf": projected_inf,
        "projected_gradient_meets_requested_gtol": projected_inf <= options["gtol"],
        "full_gradient_inf": float(np.max(np.abs(gradient))),
        "F_core_EL_residual_inf": float(np.max(np.abs(rf[core_f]))) if np.any(core_f) else None,
        "S_core_EL_residual_inf": float(np.max(np.abs(rs[core_s]))) if np.any(core_s) else None,
        "F_first_cell": float(f[0]),
        "S_first_cell": float(s[0]),
        "F_last_cell": float(f[-1]),
        "S_last_cell_minus_one": float(s[-1] - 1.0),
        "function_evaluations_actual": model.calls,
        "wall_seconds": time.perf_counter() - started,
    })
    return outcome, x


def main() -> None:
    protocol_path = HERE / "independent_protocol.json"
    protocol_bytes = protocol_path.read_bytes()
    protocol = json.loads(protocol_bytes)
    results_path = HERE / "independent_results.json"
    log_path = HERE / "independent_execution.jsonl"
    profiles_path = HERE / "independent_profiles.npz"
    if results_path.exists() or log_path.exists() or profiles_path.exists():
        raise SystemExit("Independent outputs already exist; refusing to overwrite the recorded experiment.")
    a, b = protocol["a"], protocol["b"]
    source_path = (HERE / protocol["starts"][0]["source"]).resolve()
    source = np.load(source_path)
    preflight = gradient_preflight(a, b)
    experiment = {
        "protocol_sha256": hashlib.sha256(protocol_bytes).hexdigest(),
        "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "profile_source": str(source_path),
        "profile_source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "gradient_preflight": preflight,
        "protocol": protocol,
        "runs": [],
    }
    results_path.write_text(json.dumps(experiment, indent=2) + "\n")
    all_started = time.perf_counter()
    profiles = {}

    for fraction in protocol["charge_fractions"]:
        for n in protocol["cell_counts"]:
            for start_spec in protocol["starts"]:
                name = start_spec["name"]
                label = f"q{fraction:.2f}_N{n}_{name}"
                remaining = protocol["wall_time_limit_total_seconds"] - (time.perf_counter() - all_started)
                if remaining <= 0.0:
                    row = {"label": label, "fraction": fraction, "N": n, "start": name, "skipped": True,
                           "message": "Total frozen wall-time budget exhausted before this run."}
                else:
                    charge = fraction * protocol["target_Qbar"]
                    model = RadialEnergy(n, protocol["radial_domain_L"], charge, a, b)
                    f0, s0 = initial_fields(model, name, fraction, source)
                    row, x = run_minimization(
                        model, model.pack(f0, s0), protocol["options"],
                        min(protocol["wall_time_limit_per_run_seconds"], remaining),
                    )
                    row.update({"label": label, "fraction": fraction, "Qbar": charge,
                                "N": n, "L": model.length, "h": model.h, "start": name})
                    f, s = model.unpack(x)
                    profiles[label + "_r"] = model.r
                    profiles[label + "_F"] = f
                    profiles[label + "_S"] = s
                experiment["runs"].append(row)
                with log_path.open("a") as handle:
                    handle.write(json.dumps(row) + "\n")
                print(json.dumps({k: row.get(k) for k in (
                    "label", "optimizer_success", "message", "iterations", "Ebar", "Omega",
                    "projected_gradient_inf", "wall_seconds", "skipped")}), flush=True)
                experiment["total_wall_seconds"] = time.perf_counter() - all_started
                results_path.write_text(json.dumps(experiment, indent=2) + "\n")
                np.savez_compressed(profiles_path, **profiles)

    fieldnames = ["label", "fraction", "Qbar", "N", "h", "start", "optimizer_success", "message",
                  "iterations", "Ebar", "Omega", "Ebar_over_Qbar", "projected_gradient_inf",
                  "F_core_EL_residual_inf", "S_core_EL_residual_inf", "wall_seconds", "skipped"]
    with (HERE / "independent_results.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(experiment["runs"])


if __name__ == "__main__":
    main()
