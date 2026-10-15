"""Independent finite physical controls of the full-H/fixed-graph bridge.

The theorem is stated without an electric cutoff in synthesis.md.  This script
checks its algebra against the unchanged raw frozen CP12 builder, including the
rewiring escape diagonal.  L is a finite regulator for this check only.
"""
from __future__ import annotations

import importlib.util
import hashlib
import itertools
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad
from scipy.sparse.linalg import expm_multiply

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / "cp12_ground_connectivity" / "frozen_model.py"
spec = importlib.util.spec_from_file_location("cp12_frozen_bridge_control", SOURCE)
assert spec is not None and spec.loader is not None
frozen = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = frozen
spec.loader.exec_module(frozen)


def run() -> dict:
    rng = np.random.default_rng(20261015)
    n, q, cutoff = 4, 4, 1
    edges = list(itertools.combinations(range(n), 2))
    m = len(edges)
    matrix, states, _, _ = frozen.build(n, q, cutoff)
    graph_masks = sorted({s[0] for s in states})
    c = frozen.GAMMA / m * np.sqrt((m - q) * q)
    rows = []
    worst_gram = 0.0
    worst_state = 0.0
    worst_response = 0.0
    max_escape_diagonal_error = 0.0
    for graph in graph_masks:
        inside = np.array([i for i, s in enumerate(states) if s[0] == graph])
        outside = np.array([i for i, s in enumerate(states) if s[0] != graph])
        compression = matrix[inside][:, inside].toarray()
        outgoing = matrix[outside][:, inside]
        removable = np.array([
            sum(
                1 for e in range(m)
                if (graph >> e) & 1 and (states[i][1] >> (n + e)) & 1
                and states[i][2][2 * e] == states[i][2][2 * e + 1] == 0
            )
            for i in inside
        ])
        predicted = np.diag((frozen.GAMMA / m) ** 2 * (m - q) * removable)
        gram = (outgoing.T @ outgoing).toarray()
        gram_error = float(np.max(np.abs(gram - predicted)))
        worst_gram = max(worst_gram, gram_error)
        assert gram_error < 1e-14
        assert np.sqrt(np.linalg.eigvalsh(gram)[-1]) <= c + 1e-13
        # Independently rebuild the hopping+electric/mass diagonal and isolate
        # the compression's mandatory escape potential.
        nonrewiring, _, _, _ = frozen.build(n, q, cutoff, masks=[graph], rewire=False)
        escape = compression - nonrewiring.toarray()
        escape_error = float(np.max(np.abs(
            escape - np.diag(frozen.GAMMA / m * (m - q) * removable)
        )))
        max_escape_diagonal_error = max(max_escape_diagonal_error, escape_error)
        assert escape_error < 1e-13
        if graph != graph_masks[0]:
            continue
        psi_local = rng.normal(size=len(inside)) + 1j * rng.normal(size=len(inside))
        psi_local /= np.linalg.norm(psi_local)
        psi = np.zeros(len(states), complex)
        psi[inside] = psi_local
        source_diag = np.array([(s[1] & 1) for s in states], float)
        detector_diag = np.array([((s[1] >> 2) & 1) for s in states], float)
        bpsi = source_diag * psi
        for time in (0.1, 0.5, 1.0, 2.0):
            v_local = expm(-1j * compression * time) @ psi_local
            vb_local = expm(-1j * compression * time) @ bpsi[inside]
            v = np.zeros(len(states), complex)
            vb = np.zeros(len(states), complex)
            v[inside] = v_local
            vb[inside] = vb_local
            pair = expm_multiply(-1j * matrix * time, np.column_stack((psi, bpsi)))
            upsi, ubpsi = pair[:, 0], pair[:, 1]
            state_error = float(np.linalg.norm(upsi - v))
            integral = quad(lambda s: float(np.linalg.norm(
                outgoing @ (expm(-1j * compression * s) @ psi_local)
            )), 0, time, epsabs=1e-12)[0]
            integral_b = quad(lambda s: float(np.linalg.norm(
                outgoing @ (expm(-1j * compression * s) @ bpsi[inside])
            )), 0, time, epsabs=1e-12)[0]
            response = float(2 * np.imag(np.vdot(upsi, detector_diag * ubpsi)))
            local_response = float(2 * np.imag(np.vdot(v, detector_diag * vb)))
            theta = 1e-5
            plus = upsi + np.expm1(-1j * theta) * ubpsi
            minus = upsi + np.expm1(1j * theta) * ubpsi
            finite_difference = float((
                np.vdot(plus, detector_diag * plus).real
                - np.vdot(minus, detector_diag * minus).real
            ) / (2 * theta))
            assert abs(finite_difference - response) < 1e-9
            response_error = abs(response - local_response)
            response_bound = 2 * (np.linalg.norm(bpsi) * integral + integral_b)
            worst_state = max(worst_state, state_error - integral)
            worst_response = max(worst_response, response_error - response_bound)
            assert state_error <= integral + 1e-12
            assert integral <= c * time + 1e-12
            assert response_error <= response_bound + 1e-12
            rows.append({
                "time": time,
                "state_difference_norm": state_error,
                "integrated_residual_bound": integral,
                "uniform_residual_bound": float(c * time),
                "full_kubo_response": response,
                "source_pulse_convention": "exp(-i theta n_vertex_0)",
                "direct_central_source_derivative": finite_difference,
                "compressed_kubo_response": local_response,
                "response_difference": response_error,
                "response_residual_bound": float(response_bound),
            })
    return {
        "status": "PASS",
        "frozen_builder_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "scope": "Regulated raw physical finite-H control; theorem itself is unregulated analytical.",
        "N": n,
        "Q": q,
        "electric_cutoff": cutoff,
        "dimension": len(states),
        "fixed_graph_blocks_checked": len(graph_masks),
        "max_gram_identity_error": worst_gram,
        "maximum_positive_state_bound_violation": max(0.0, worst_state),
        "maximum_positive_response_bound_violation": max(0.0, worst_response),
        "maximum_escape_diagonal_reconstruction_error": max_escape_diagonal_error,
        "fixed_graph_coefficient_C_N": float(c),
        "time_controls": rows,
    }


if __name__ == "__main__":
    result = run()
    (HERE / "LOCALITY_BRIDGE_CHECK.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
