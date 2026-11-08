"""M1-01 massive lattice comparator; reviewed ED Hamiltonian reused explicitly."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import time

import numpy as np
from scipy.sparse.linalg import eigsh, expm_multiply

from continuum import key

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "grut_qed_detector_histories_20261107/ed/schwinger.py"
spec = importlib.util.spec_from_file_location("archived_schwinger", SOURCE)
schwinger = importlib.util.module_from_spec(spec)
spec.loader.exec_module(schwinger)


def run(card, n, bare_zero=False, refined=False):
    start = time.monotonic()
    a = card["length"] / (n + 1)
    x = a ** -2
    mass = card["mass_R_over_g"] - a / 8
    w, electric, states = schwinger.hamiltonian(n, x, mass)
    h = w / (2 * np.sqrt(x))
    tol = card["eigensolver_refined_tol"] if refined else card["eigensolver_tol"]
    v0 = np.random.default_rng(card["eigensolver_seed_base"] + n).standard_normal(h.shape[0])
    energies, vectors = eigsh(h, k=2, which="SA", tol=tol, v0=v0, maxiter=10000)
    order = np.argsort(energies)
    energies, vectors = energies[order], vectors[:, order]
    residuals = [float(np.linalg.norm(h @ vectors[:, j] - energies[j] * vectors[:, j]))
                 for j in range(2)]
    vacuum = vectors[:, 0].astype(complex)
    z = (np.arange(n - 1) + 1.5) * a
    weights = np.exp(-0.5 * ((z - card["length"] / 2) / card["sigma"]) ** 2)
    weights /= weights.sum()
    f = electric @ weights
    mean = float(np.vdot(vacuum, f * vacuum).real)
    delta = (f - mean) * vacuum
    trace_h = float(h.diagonal().sum())
    substeps = card["refined_substeps"] if refined else 1
    def evolve(state, dt):
        for _ in range(substeps):
            if dt:
                step = dt / substeps
                state = expm_multiply((-1j * step) * h, state,
                                      traceA=-1j * step * trace_h)
        return state
    histories = {}
    max_free_norm_error = 0.0
    for times in card["protocols"]:
        branches = vacuum[:, None]
        labels = [""]
        previous = 0.0
        phase_w = np.exp(1j * (card["phase"] + card["coupling"] * f))[:, None]
        for moment in times:
            before = np.linalg.norm(branches) ** 2
            branches = evolve(branches, moment - previous)
            max_free_norm_error = max(max_free_norm_error,
                                      abs(np.linalg.norm(branches) ** 2 - before))
            shifted = phase_w * branches
            branches = np.concatenate(((branches + shifted) / 2,
                                       (branches - shifted) / 2), axis=1)
            labels = [label + "+" for label in labels] + [label + "-" for label in labels]
            previous = moment
        histories[key(times)] = {label: float(np.vdot(branches[:, j], branches[:, j]).real)
                                 for j, label in enumerate(labels)}
    covariance = {}
    for moment in (0.0, 1.0, 2.0):
        value = np.vdot(delta, evolve(delta, moment)) * np.exp(1j * energies[0] * moment)
        covariance[f"{moment:g}"] = [float(value.real), float(value.imag)]
    one = histories["0"]
    disturbance = {}
    for label in ("0,1", "0,2"):
        p = histories[label]
        disturbance[label] = 0.5 * sum(abs(sum(v for k, v in p.items() if k[-1] == s) - one[s])
                                          for s in ("+", "-"))
    return {"N": n, "a": a, "x": x, "mass_lat_over_g": mass,
            "bare_zero_control": bare_zero, "solver_refinement": refined,
            "dimension": int(h.shape[0]), "ground_energy": float(energies[0]),
            "gap": float(energies[1] - energies[0]), "mean": mean,
            "covariance": covariance, "histories": histories, "disturbance": disturbance,
            "eigen_residuals": residuals, "free_norm_error": float(max_free_norm_error),
            "probability_norm_error": max(abs(sum(p.values()) - 1) for p in histories.values()),
            "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            "seconds": time.monotonic() - start}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, required=True)
    parser.add_argument("--bare-zero", action="store_true")
    parser.add_argument("--refined", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    card_path = HERE / "protocol_m1_01.json"
    card = json.loads(card_path.read_text())
    data = run(card, args.N, args.bare_zero, args.refined)
    data["protocol_sha256"] = hashlib.sha256(card_path.read_bytes()).hexdigest()
    with args.output.open("x") as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: data[k] for k in ("N", "mass_lat_over_g", "mean", "gap", "seconds")}))
