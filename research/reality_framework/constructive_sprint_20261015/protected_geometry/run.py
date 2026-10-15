"""Full-law finite-time control on an explicitly supplied flux-carrying cycle.

The electric box is a numerical retained space, never a changed physical law.
Every omitted outgoing transition enters a Duhamel error certificate.  Geometry,
orientation and the localized initial charge are supplied diagnostics.
"""
from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from math import factorial, lcm, sqrt
from pathlib import Path
import argparse
import json
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix

HERE = Path(__file__).resolve().parent
RF = HERE.parents[1]
FLOATING_ALLOWANCE = 5e-12  # reported numerical tolerance, not an interval proof
sys.path.insert(0, str(RF / "cp12_ground_connectivity"))
sys.path.insert(0, str(RF / "cp12_distance_history"))
from frozen_model import build, gauss, hop
from exact_destination import transitions
from evolution import Evolution


def local_norm_checks():
    """Independent local CAR calculation, with arbitrary unitary rotor phases."""
    rows = []
    for left, right in [(0., 0.), (.3, -.8), (1.7, 2.4)]:
        h = np.zeros((8, 8), dtype=complex)
        for b in range(8):
            for v, phase, orient in [(0, left, 1), (1, right, -1)]:
                for dst, src, z in [(v, 2, np.exp(1j * phase)),
                                    (2, v, np.exp(-1j * phase))]:
                    result = hop(b, dst, src)
                    if result:
                        nb, sign = result
                        h[nb, b] += orient * sign * z
        vals = np.linalg.eigvalsh(h)
        assert np.max(abs(h - h.conj().T)) < 1e-14
        assert abs(np.max(abs(vals)) - sqrt(2)) < 1e-14
        rows.append(dict(rotor_angles=[left, right], operator_norm=float(np.max(abs(vals)))))
    return rows


def initial(n, L):
    es = list(combinations(range(n), 2))
    directed = [(0, 1), (1, 2), (2, 3), (3, 0)]
    g = 0
    E = [0] * (2 * len(es))
    for a, b in directed:
        e = es.index(tuple(sorted((a, b))))
        g |= 1 << e
        E[2 * e] = L if a < b else -L
        E[2 * e + 1] = -E[2 * e]
    bits = sum(1 << (n + e) for e in range(len(es)))
    e = es.index((0, 1))
    bits, sign = hop(bits, 0, n + e)
    E[2 * e] += 1
    state = g, bits, tuple(E)
    assert not np.any(gauss(n, state))
    return state, sign


def full_boundary(n, H, states, lookup):
    """Compare every retained column independently and retain every exit row."""
    M = n * (n - 1) // 2
    scale = lcm(8, 10 * M)
    outside = {}
    br, bc, bv = [], [], []
    rr, cc, vv = [], [], []
    for col, state in enumerate(states):
        for dest, amp in transitions(n, state, scale):
            a = amp / scale
            if dest in lookup:
                rr.append(lookup[dest]); cc.append(col); vv.append(a)
            else:
                row = outside.setdefault(dest, len(outside))
                br.append(row); bc.append(col); bv.append(a)
    other = coo_matrix((vv, (rr, cc)), shape=H.shape).tocsr()
    delta = H - other
    discrepancy = float(np.max(abs(delta.data))) if delta.nnz else 0.
    assert discrepancy < 1e-12
    B = coo_matrix((bv, (br, bc)), shape=(len(outside), len(states))).tocsr()
    return B, discrepancy


def interval(p, error):
    e = error + FLOATING_ALLOWANCE
    return [min(1., max(0., sqrt(max(0., p)) - e)**2),
            min(1., (sqrt(max(0., p)) + e)**2)]


def observables(states, host):
    host_value = np.array([g == host for g, b, E in states], float)
    protection = []
    circulation = []
    charge = [[] for _ in range(4)]
    for g, b, E in states:
        present = [e for e in range(6) if g >> e & 1]
        protection.append(all(E[2 * e] or E[2 * e + 1] for e in present))
        # On this supplied directed square, require its original four through-flows.
        es = list(combinations(range(4), 2))
        original_flow = g == host
        for a, c in [(0, 1), (1, 2), (2, 3), (3, 0)]:
            e = es.index(tuple(sorted((a, c))))
            leg = 0 if a < c else 1
            original_flow = original_flow and E[2 * e + leg] > 0 and E[2 * e + 1 - leg] < 0
        circulation.append(original_flow)
        for j in range(4):
            charge[j].append((b >> j) & 1)
    return dict(host=host_value, protection=np.array(protection, float),
                original_circulation=np.array(circulation, float),
                density=np.array(charge, float))


def run_box(cutoff, fluxes, waits_lags, steps):
    start = time.monotonic()
    H, states, lookup, _ = build(4, 4, cutoff)
    B, discrepancy = full_boundary(4, H, states, lookup)
    host = initial(4, 1)[0][0]
    O = observables(states, host)
    rows = []
    for L in fluxes:
        key, sign = initial(4, L)
        assert key in lookup
        v = np.zeros(len(states), dtype=complex)
        v[lookup[key]] = sign
        energy = float(np.vdot(v, H @ v).real)
        energy_variance = float(np.linalg.norm(H @ v - energy * v)**2)
        E = Evolution(H, B, energy, steps=steps)
        for wait, lag in waits_lags:
            w, _, ew = E(v, wait)
            pair = np.column_stack([w, O["density"][0] * w])
            z, lag_errors, _ = E(pair, lag)
            x, y = z[:, 0], z[:, 1]
            ex, ey = float(ew + lag_errors[0]), float(ew + lag_errors[1])
            error = 2 * (ex + ey + ex * ey)
            probability = abs(x)**2
            values = {}
            for name in ["host", "protection", "original_circulation"]:
                p = float(O[name] @ probability)
                values[name] = dict(value=p, interval=interval(p, ex))
            density = [float(d @ probability) for d in O["density"]]
            response = [float(2 * np.vdot(x, d * y).imag) for d in O["density"]]
            # The source perturbation is exp(-i theta n_0), a gauge-invariant unitary.
            # The final detector remains distance-independent local vertex density.
            source = []
            for j, chi in enumerate(response):
                distance = 0 if j == 0 else (2 if j == 2 else 1)
                joint = float(2 * np.vdot(x, O["host"] * O["density"][j] * y).imag)
                # Exact phase-shift identity for the projector source n_0.
                plus = x + (-1j-1) * y
                minus = x + (1j-1) * y
                phase_identity = float((np.vdot(plus,O["density"][j]*plus).real
                                        -np.vdot(minus,O["density"][j]*minus).real)/2)
                assert abs(phase_identity-chi) < 1e-12
                source.append(dict(receiver=j, initial_host_distance=distance,
                    response=chi, error_bound=error,
                    floating_point_allowance=FLOATING_ALLOWANCE,
                    interval=[chi-error-FLOATING_ALLOWANCE, chi+error+FLOATING_ALLOWANCE],
                    zero_excluded=abs(chi) > error+FLOATING_ALLOWANCE, same_host_joint_response=joint,
                    escaped_host_response=chi-joint,
                    exact_phase_shift_formula_difference=abs(phase_identity-chi)))
            T = wait + lag
            edge_amplitude_bound = (sqrt(2) * T)**(2 * L) / factorial(2 * L)
            row = dict(N=4, Q=4, flux=L, cutoff=cutoff, wait=wait, lag=lag,
                time=T, initial_energy=energy, initial_energy_variance=energy_variance,
                initial_electric_square=8 * L * L + 2 * L + 1,
                forward_state_error=ex, perturbed_path_error=ey,
                graph_change_probability_upper=min(1., 4 * edge_amplitude_bound**2),
                final_density=density, detector_responses=source, **values)
            rows.append(row)
            print(json.dumps(dict(cutoff=cutoff, flux=L, time=T, host=values["host"]["value"],
                response_d2=response[2], response_error=error, state_error=ex)), flush=True)
    return dict(cutoff=cutoff, dimension=len(states), nonzeros=H.nnz,
        outgoing_dimension=B.shape[0], outgoing_nonzeros=B.nnz,
        independent_H_column_maximum_error=discrepancy,
        seconds=time.monotonic()-start, results=rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cutoffs", default="6,10")
    parser.add_argument("--fluxes", default="1,2")
    parser.add_argument("--steps", type=int, default=64)
    parser.add_argument("--out", default="RESULTS.json")
    args = parser.parse_args()
    fluxes = [int(x) for x in args.fluxes.split(",")]
    cutoffs = [int(x) for x in args.cutoffs.split(",")]
    checks = local_norm_checks()
    boxes = [run_box(c, fluxes, [(.25,.25),(.5,.5)], args.steps) for c in cutoffs]
    convergence = []
    for small, large in zip(boxes, boxes[1:]):
        for a, b in zip(small["results"], large["results"]):
            assert (a["flux"],a["time"]) == (b["flux"],b["time"])
            diff = max(abs(x["response"]-y["response"])
                       for x,y in zip(a["detector_responses"],b["detector_responses"]))
            allowance = (a["detector_responses"][0]["error_bound"]
                         + b["detector_responses"][0]["error_bound"])
            assert diff <= allowance + 5e-12
            convergence.append(dict(flux=a["flux"],time=a["time"],
                maximum_detector_difference=diff,analytic_allowance=allowance,
                floating_point_tolerance=5e-12))
    result = dict(status="PASS", local_hopping_norm_checks=checks,
        truncation_error_method="Full outgoing residual; fifth-order midpoint upper sums for leakage, unitary Duhamel norm error; ordinary floating-point evaluation, not formal interval arithmetic.",
        reported_floating_point_allowance=FLOATING_ALLOWANCE,
        reported_floating_point_allowance_scope="A conservative numerical reporting tolerance added to probability amplitudes and response intervals. It exceeds the observed cutoff discrepancies; it is not a mathematically certified roundoff bound. error_bound fields separately report the analytic electric-truncation estimate.",
        source="Gauge-invariant phase kick exp(-i theta n_0) after unforced waiting time.",
        supplied_structure="Directed four-vertex square with integer flux L, plus one initially localized Gauss-admissible vertex fermion and matching edge hole. This is neither formation from the original neutral ray nor an emergent-dimension calculation.",
        same_state="Host persistence, protected edges, circulation and receiver density responses are measured on the same evolved preparation at each L and time. Responses are unconditional, with no loop filter.",
        hopping_off_control="Exactly zero response: all vertex densities commute with H_diagonal+H_rewire. This control is an operator identity, not a numerical subtraction.",
        boxes=boxes, cutoff_comparisons=convergence)
    (HERE / args.out).write_text(json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    main()
