"""One frozen finite experiment. No external service, fitting, or file writes."""
import importlib.util
import json
import math
from pathlib import Path
import resource
import sys
import time

import numpy as np
from scipy.linalg import expm
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def load_file(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def maximum(matrix):
    return float(np.max(np.abs(matrix.data), initial=0))


def evaluate():
    started = time.monotonic()
    protocol = json.loads((HERE / "PROTOCOL.json").read_text())
    ledger = json.loads((HERE / "claims.json").read_text())
    auditor = load_file("discovery_auditor", ROOT / "provenance/auditor.py")
    audit = auditor.audit(ledger["claims"], ledger["sources"])
    if not audit.ok:
        raise ValueError(audit.blocking)
    p1 = load_file("discovery_p1", HERE.parent / "grut_p1_preparation_20261017/model.py")
    model = p1.Model(3, modes=2)
    assert model.dim <= protocol["budgets"]["max_matrix_dimension"]
    h = model.H0 + model.C + model.V
    tol = protocol["numerical_tolerance"]
    hermitian = maximum(h - h.getH())
    conserved = model.fcount + model.number
    commutator = maximum(h @ diags(conserved) - diags(conserved) @ h)
    # P1's kernel has equal species counts; check that directly, not by naming.
    assert all(len(s[3]) == len(s[4]) for s in model.states)
    assert hermitian < tol and commutator < tol
    seed_vacuum = model.initial()
    seed = np.zeros(model.dim, complex)
    for i in np.flatnonzero(abs(seed_vacuum) > 0):
        word, q, receiver, plus, minus, fuel = model.states[i]
        seed[model.index[(word, q, receiver, plus, minus, fuel | 1)]] = seed_vacuum[i]
    initial_fuel = float(np.vdot(seed, model.fuel * seed).real)
    assert abs(initial_fuel - 3) < tol
    times = np.linspace(0, protocol["time_max"], protocol["sample_count"])
    trajectories = expm_multiply(-1j * h, seed, start=0, stop=times[-1], num=len(times))
    probabilities = abs(trajectories) ** 2
    pairs = probabilities @ model.number
    fuels = probabilities @ model.fcount
    norms = probabilities.sum(axis=1)
    energies = np.array([np.vdot(v, h @ v).real for v in trajectories])
    norm_error = float(max(abs(norms - 1)))
    energy_error = float(max(abs(energies - energies[0])))
    conservation_error = float(max(abs(fuels + pairs - 1)))
    assert max(norm_error, energy_error, conservation_error) < tol
    # One independently implemented dense spectral evolution at two times.
    eigenvalues, eigenvectors = np.linalg.eigh(h.toarray())
    coefficients = eigenvectors.conj().T @ seed
    dense_errors = []
    for index in (16, 32):
        dense = eigenvectors @ (np.exp(-1j * times[index] * eigenvalues) * coefficients)
        dense_errors.append(float(np.linalg.norm(dense - trajectories[index])))
    assert max(dense_errors) < tol
    vacuum = expm_multiply(-1j * h * math.pi, seed_vacuum)
    baseline = expm_multiply(-1j * (model.H0 + model.V) * math.pi, seed)
    vacuum_pairs = float(np.vdot(vacuum, model.number * vacuum).real)
    baseline_pairs = float(np.vdot(baseline, model.number * baseline).real)
    assert max(abs(vacuum_pairs), abs(baseline_pairs)) < tol
    assert model.C.nnz and model.T.nnz and model.G.nnz and model.B.nnz
    assert np.linalg.norm(model.C @ seed) > 0
    # Exercise geometry and hopping on actual states, not nnz counts alone.
    geometry_action = float(np.linalg.norm((model.G + model.B) @ seed))
    hopping_action = float(np.linalg.norm(model.T @ trajectories[16]))
    assert geometry_action > 1e-8 and hopping_action > 1e-8
    controls = []
    for detuning in (-.5, 0., .5):
        matrix = np.array([[3+detuning, 1], [1, 3]], complex)
        omega = math.sqrt(1 + detuning**2 / 4)
        for t in (math.pi / 2, math.pi):
            observed = abs((expm(-1j * t * matrix) @ [1, 0])[1])**2
            analytic = math.sin(t * omega)**2 / omega**2
            assert abs(observed - analytic) < tol
            controls.append({"detuning": detuning, "time": t,
                             "probability": float(observed), "formula": analytic})
    low, high = protocol["persistence_window"]
    window = np.flatnonzero((times >= low - 1e-14) & (times <= high + 1e-14))
    worst = int(window[np.argmin(pairs[window])])
    refuted = bool(pairs[worst] < protocol["pair_threshold"] - tol)
    return {
        "candidate": protocol["candidate"],
        "decision": "refuted finite persistence target" if refuted else "sampled support only; continuous target unresolved",
        "dimension": model.dim,
        "scope": protocol["scope_limits"],
        "audit": {"passed": audit.ok, "warnings": audit.warnings, "net_supplied_inputs": audit.net,
                  "certifies_physical_truth": False},
        "structural_checks": {"hermiticity_max": hermitian, "conserved_commutator_max": commutator,
                              "norm_max_error": norm_error, "energy_max_error": energy_error,
                              "fuel_plus_pair_max_error": conservation_error,
                              "dense_sparse_state_errors": dense_errors,
                              "vacuum_pairs_at_pi": vacuum_pairs, "no_conversion_pairs_at_pi": baseline_pairs,
                              "initial_supplied_fuel_energy": initial_fuel,
                              "geometry_action_norm": geometry_action, "hopping_action_norm": hopping_action},
        "persistence_witness": {"time": float(times[worst]), "expected_pairs": float(pairs[worst]),
                                "threshold": protocol["pair_threshold"], "refuted": refuted},
        "peak_sampled_expected_pairs": float(max(pairs)),
        "isolated_rabi_controls": controls,
        "series": [{"time": float(t), "expected_pairs": float(p), "expected_fuel": float(f)}
                   for t, p, f in zip(times, pairs, fuels)],
        "resources": {"wall_seconds": time.monotonic() - started,
                      "peak_rss_MiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024},
        "independent_validation": False, "experimental_observations": 0,
        "next_decision": "stop; seek a separately frozen mechanism with a persistence bound, not another tuned pulse"
    }


if __name__ == "__main__":
    json.dump(evaluate(), sys.stdout, indent=2, allow_nan=False)
    print()
