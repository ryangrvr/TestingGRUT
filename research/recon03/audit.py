"""Disclosed standard-quantum controls; no GRUT candidate is searched or tested."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import expm, schur, svdvals

ROOT = Path(__file__).resolve().parent
SEED = 20261010
TOL = 1e-10  # Floating diagnostic tolerance, not an exact interval certificate.
I2 = np.eye(2, dtype=complex)
Z = np.diag([1, -1]).astype(complex)
HAD = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
CNOT = np.array([[1, 0, 0, 0], [0, 1, 0, 0],
                 [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex)
RCNOT = np.array([[1, 0, 0, 0], [0, 0, 0, 1],
                  [0, 0, 1, 0], [0, 1, 0, 0]], dtype=complex)
BELL = CNOT @ np.kron(HAD @ np.diag([1, 1j]), I2)


def enc(a):
    a = np.asarray(a)
    return {"real": a.real.tolist(), "imag": a.imag.tolist()}


def log_h(u):
    """Hermitian principal log, u=exp(-i*h); bounded by pi in pulse units."""
    t, q = schur(u, output="complex")
    h = q @ np.diag(-np.angle(np.diag(t))) @ q.conj().T
    return (h + h.conj().T) / 2


def gate(i, u, h=None):
    return {"edge": i, "u": u, "h": log_h(u) if h is None else h}


def apply_pair(state, u, i, n):
    perm = [i, i + 1] + [j for j in range(n) if j not in (i, i + 1)]
    a = state.reshape([2] * n).transpose(perm).reshape(4, -1)
    return (u @ a).reshape([2] * n).transpose(np.argsort(perm)).reshape(-1)


def evolve(state, layers, n):
    out = state.copy()
    for layer in layers:
        for g in layer:
            out = apply_pair(out, g["u"], g["edge"], n)
    return out


def reduced(state, sites, n):
    sites = list(sites)
    other = [i for i in range(n) if i not in sites]
    a = state.reshape([2] * n).transpose(sites + other)
    a = a.reshape(2 ** len(sites), -1)
    return a @ a.conj().T


def entropy_from_amplitudes(a):
    p = svdvals(a) ** 2
    p = p[p > 1e-15]
    return float(-np.sum(p * np.log2(p)))


def information(states, sites, n):
    """Pointer Holevo information and actual I(S:A); neither replaces trace D."""
    sites = list(sites)
    other = [i for i in range(n) if i not in sites]
    a0, a1 = [s.reshape([2]*n).transpose(sites+other).reshape(2**len(sites), -1)
              for s in states]
    sa = entropy_from_amplitudes(np.concatenate((a0, a1), axis=1) / np.sqrt(2))
    sb = entropy_from_amplitudes(np.concatenate((a0.T, a1.T), axis=1) / np.sqrt(2))
    mean_branch_entropy = (entropy_from_amplitudes(a0) + entropy_from_amplitudes(a1)) / 2
    return sa - mean_branch_entropy, 1 + sa - sb  # S entropy exactly 1.


def td(rho, sigma):
    return float(np.abs(np.linalg.eigvalsh(rho - sigma)).sum() / 2)


def packing(distances, n, threshold=1.0):
    dp = [0] * (n + 1)
    for end in range(n):
        dp[end + 1] = dp[end]
        for start in range(end + 1):
            if distances[start, end] >= threshold - TOL:
                dp[end + 1] = max(dp[end + 1], dp[start] + 1)
    return dp[n]


def cone_layers(layers, center):
    support = {center}
    selected = []
    for layer in layers:
        chosen = []
        old = support.copy()
        for g in layer:
            edge = {g["edge"], g["edge"] + 1}
            if old & edge:
                chosen.append(g)
                support |= edge
        selected.append(chosen)
    return selected, sorted(support)


def inverse_decode(states, selected, center, n):
    decoded = [s.copy() for s in states]
    for layer in reversed(selected):
        for g in reversed(layer):
            decoded = [apply_pair(s, g["u"].conj().T, g["edge"], n)
                       for s in decoded]
    out = []
    for state in decoded:
        rho = reduced(state, [center], n)
        out.append(float(np.real(np.trace(Z @ rho))))
    return out


def random_layers(n, d, rng):
    layers = []
    for depth in range(d):
        layer = []
        for i in range(depth % 2, n - 1, 2):
            a = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
            h = (a + a.conj().T) / 2
            h *= np.pi * rng.uniform(0.15, 1) / np.linalg.norm(h, 2)
            layer.append(gate(i, expm(-1j * h), h))
        layers.append(layer)
    return layers


def compression(n, d):
    b = d + 1
    layers = [[] for _ in range(d)]
    retained = []
    for start in range(0, n, b):
        stop = min(start + b, n)
        retained.append(start)
        for step, i in enumerate(range(stop - 1, start, -1)):
            layers[step].append(gate(i - 1, CNOT))
    return layers, retained


def balanced_compression(n, d):
    """Inward CNOT sweeps in blocks of at most 2*d; last shared-root gates serial."""
    assert d >= 1
    layers = [[] for _ in range(d)]
    retained = []
    for start in range(0, n, 2*d):
        stop = min(start + 2*d, n)
        root = (start + stop - 1) // 2
        retained.append(root)
        for step, i in enumerate(range(start, root)):
            layers[step].append(gate(i, RCNOT))
        for step, i in enumerate(range(stop - 1, root, -1)):
            g = gate(i - 1, CNOT)
            # An odd block's two final pulses share the root; delay the right one.
            occupied = {j for h in layers[step] for j in (h["edge"], h["edge"]+1)}
            if {i-1, i} & occupied:
                step += 1
            assert step < d
            layers[step].append(g)
    return layers, retained


def main():
    rng = np.random.default_rng(SEED)
    specs = []
    for n in (4, 6, 8):
        for d in (0, 1, 2, 3):
            for rep in range(2):
                specs.append((f"random_N{n}_d{d}_{rep}", n, d,
                              random_layers(n, d, rng), "random", None))
    for d in (1, 2, 3):
        specs.append((f"identity_N8_d{d}", 8, d, [[] for _ in range(d)],
                      "identity", None))
    for n in (4, 6, 8):
        specs.append((f"bell_N{n}_d1", n, 1,
                      [[gate(i, BELL) for i in range(0, n, 2)]], "bell", None))
        for d in (0, 1, 2, 3):
            layers, kept = compression(n, d)
            specs.append((f"compress_N{n}_d{d}", n, d, layers, "compression", kept))
        layers, kept = compression(n, n - 1)
        specs.append((f"full_compress_N{n}_d{n-1}", n, n - 1, layers,
                      "longer_budget_full_compression", kept))
        for d in (1, 2, 3):
            layers, kept = balanced_compression(n, d)
            specs.append((f"balanced_compress_N{n}_d{d}", n, d, layers,
                          "balanced_compression", kept))
        layers, kept = balanced_compression(n, n//2)
        specs.append((f"balanced_full_N{n}_d{n//2}", n, n//2, layers,
                      "longer_budget_balanced_full", kept))
    specs.append(("raw_detector_failure_N8_d1", 8, 1,
                  [[gate(i, np.kron(HAD, HAD)) for i in range(0, 8, 2)]],
                  "raw_detector_failure", None))

    fixtures = []
    max_error = 0.0
    counts = {"intervals": 0, "predeclared_blocks": 0, "decoded_blocks": 0,
              "channel_matrix_units": 0, "local_pulses": 0,
              "preparation_error_controls": 0}
    for name, n, d, layers, kind, kept in specs:
        assert len(layers) <= d
        pulse_norms, pulse_errors = [], []
        for layer in layers:
            used = set()
            for g in layer:
                i = g["edge"]
                assert 0 <= i < n - 1 and not ({i, i + 1} & used)
                used |= {i, i + 1}
                norm = float(np.linalg.norm(g["h"], 2))
                assert norm <= np.pi + TOL
                error = float(np.max(np.abs(expm(-1j * g["h"]) - g["u"])))
                assert error < TOL
                pulse_norms.append(norm)
                pulse_errors.append(error)
                counts["local_pulses"] += 1
        e0, e1 = np.zeros(2 ** n, complex), np.zeros(2 ** n, complex)
        e0[0], e1[-1] = 1, 1
        states = [evolve(e, layers, n) for e in (e0, e1)]
        orth_error = float(abs(np.vdot(states[1], states[0])))
        norm_error = max(float(abs(np.vdot(s, s) - 1)) for s in states)
        assert orth_error < TOL and norm_error < TOL

        # Full reduced S channel: trace each conditional matrix-unit dilation.
        channel_errors = []
        for b in (0, 1):
            for c in (0, 1):
                actual = np.trace(np.outer(states[b], states[c].conj()))
                channel_errors.append(float(abs(actual - (1 if b == c else 0))))
                counts["channel_matrix_units"] += 1
        assert max(channel_errors) < TOL

        distances, interval_rows = {}, []
        rho_pairs = {}
        for start in range(n):
            for end in range(start, n):
                rho0, rho1 = [reduced(s, range(start, end + 1), n) for s in states]
                value = td(rho0, rho1)
                assert -TOL <= value <= 1 + TOL
                distances[start, end] = value
                rho_pairs[start, end] = (rho0, rho1)
                chi, mutual = information(states, range(start, end+1), n)
                assert -TOL <= chi <= 1+TOL and -TOL <= mutual <= 2+TOL
                if value >= 1-TOL:
                    assert abs(chi-1) < TOL
                interval_rows.append({"start": start, "end": end, "D": value,
                                      "pointer_Holevo_bits": chi, "I_S_A_bits": mutual})
                counts["intervals"] += 1
        block_rows = []
        for m in range(1, n + 1):
            values = [distances[start, start + m - 1] for start in range(0, n-m+1, m)]
            block_rows.append({"size": m, "count": len(values), "D": values,
                               "perfect_count": sum(v >= 1-TOL for v in values)})
            counts["predeclared_blocks"] += len(values)

        # Conventional averaged Holevo redundancy, with an expressly fixed
        # distribution over contiguous intervals (not over arbitrary subsets).
        averaged_chart = []
        for m in range(1, n+1):
            rows = [r for r in interval_rows if r["end"]-r["start"]+1 == m]
            averaged_chart.append({"size": m,
                                   "mean_pointer_Holevo_bits": float(np.mean([r["pointer_Holevo_bits"] for r in rows])),
                                   "mean_I_S_A_bits": float(np.mean([r["I_S_A_bits"] for r in rows]))})
        m90 = next(r["size"] for r in averaged_chart if r["mean_pointer_Holevo_bits"] >= 0.9-TOL)

        length = 2 * d + 1
        certified = []
        for start in range(0, n - length + 1, length):
            center, end = start + d, start + length - 1
            selected, support = cone_layers(layers, center)
            assert all(start <= i <= end for i in support)
            decoded = inverse_decode(states, selected, center, n)
            error = max(abs(decoded[0] - 1), abs(decoded[1] + 1),
                        abs(distances[start, end] - 1))
            assert error < TOL
            certified.append({"start": start, "end": end, "center": center,
                              "cone_support": support, "decoded_Z": decoded,
                              "D": distances[start, end]})
            counts["decoded_blocks"] += 1
            for p in (0.01, 0.2):
                r0, r1 = rho_pairs[start, end]
                measured = td((1-p)*r0 + p*r1, (1-p)*r1 + p*r0)
                assert abs(measured - (1 - 2*p)) < TOL
                counts["preparation_error_controls"] += 1

        capacity = packing(distances, n)
        lower = max(1, n // length, 2 if 2*d < n-1 else 1)
        assert capacity >= lower
        endpoint_decoders = []
        for center in (0, n-1):
            selected, support = cone_layers(layers, center)
            assert all(abs(i-center) <= d for i in support)
            decoded = inverse_decode(states, selected, center, n)
            assert max(abs(decoded[0]-1), abs(decoded[1]+1)) < TOL
            endpoint_decoders.append({"center": center, "support": support,
                                      "decoded_Z": decoded})
        raw_single = []
        for i in range(n):
            r0, r1 = rho_pairs[i, i]
            raw_single.append(float(abs(np.real(r0[0, 0] - r1[0, 0]))))
        single = [distances[i, i] for i in range(n)]
        if kind == "identity":
            assert capacity == n and max(abs(v-1) for v in single) < TOL
        elif kind == "bell":
            assert capacity == n // 2 and max(abs(v) for v in single) < TOL
            assert max(abs(distances[i, i+1]-1) for i in range(0, n, 2)) < TOL
        elif kind in ("compression", "longer_budget_full_compression",
                      "balanced_compression", "longer_budget_balanced_full"):
            expected_idx = sum(1 << (n-1-i) for i in kept)
            assert abs(abs(states[0][0])-1) < TOL
            assert abs(abs(states[1][expected_idx])-1) < TOL
            assert capacity == len(kept)
            assert max(abs(distances[a,b] - int(any(a <= i <= b for i in kept)))
                       for a,b in distances) < TOL
        elif kind == "raw_detector_failure":
            assert capacity == n and max(abs(v-1) for v in single) < TOL
            assert max(abs(v) for v in raw_single) < TOL
        max_error = max(max_error, norm_error, orth_error, max(channel_errors),
                        max(pulse_errors, default=0),
                        max((max(abs(c["decoded_Z"][0]-1),
                                 abs(c["decoded_Z"][1]+1), abs(c["D"]-1))
                             for c in certified), default=0))
        fixtures.append({"name": name, "kind": kind, "N": n, "depth_ceiling": d,
                         "used_layers": len(layers), "packing_capacity": capacity,
                         "universal_packing_lower_bound": lower,
                         "single_D": single, "raw_single_Z_D": raw_single,
                         "fixed_fragment_charts": block_rows,
                         "fixed_contiguous_interval_averaging": averaged_chart,
                         "averaged_Holevo_R_delta_0p1": n/m90,
                         "certified_blocks": certified, "all_intervals": interval_rows,
                         "endpoint_decoders": endpoint_decoders,
                         "channel_error": max(channel_errors),
                         "pulse_norms_hbar_over_tau": pulse_norms,
                         "layers": [[{"edge": g["edge"], "u": enc(g["u"]),
                                      "h": enc(g["h"])} for g in layer] for layer in layers]})

    # Explicit outside-baseline reset: each local channel has Kraus K0,K1.
    k0, k1 = np.array([[1, 0], [0, 0]]), np.array([[0, 1], [0, 0]])
    assert np.allclose(k0.T @ k0 + k1.T @ k1, I2)
    reset_outputs = [sum(k @ r @ k.T for k in (k0, k1))
                     for r in (np.diag([1, 0]), np.diag([0, 1]))]
    assert td(*reset_outputs) == 0
    result = {"status": "STANDARD_LOCAL_RECORD_CAPACITY_CONTROLS_PASS_NO_GRUT_LAW",
              "date": "2026-10-10", "seed": SEED, "numpy": np.__version__,
              "scipy": scipy.__version__, "floating_tolerance": TOL,
              "not_an_interval_or_experimental_certificate": True,
              "fixture_count": len(fixtures), "counts": counts,
              "maximum_identity_or_decoder_discrepancy": max_error,
              "reset_control": {"inside_baseline": False,
                                "resource_change": "one fresh blank ancilla per E site, local SWAP then trace",
                                "E_record_D": 0, "S_channel": "unchanged complete dephasing"},
              "budget": {"cards_consumed": 1, "cards_remaining": 2, "card2": "UNOPENED"},
              "fixtures": fixtures}
    # One compact line per fixture avoids tens of thousands of matrix lines
    # while retaining every numerical cell in ordinary machine-readable JSON.
    header = json.dumps({k: v for k, v in result.items() if k != "fixtures"}, indent=2)
    body = ",\n".join("    " + json.dumps(f, separators=(",", ":")) for f in fixtures)
    (ROOT / "results.json").write_text(header[:-2] + ',\n  "fixtures": [\n' + body + '\n  ]\n}\n')
    lookup = {f["name"]: f for f in fixtures}
    lines = ["# Generated finite-control results", "",
             "Generated by audit.py from results.json; synthetic diagnostics, not a device certificate.", "",
             "| Control, N=8 | Depth ceiling | Packing capacity | Perfect singleton records | Averaged Holevo R_0.1 |",
             "|---|---:|---:|---:|---:|"]
    names = [("identity_N8_d1", "Identity"), ("bell_N8_d1", "Bell pair encoding"),
             ("balanced_compress_N8_d1", "Balanced block compression"),
             ("balanced_compress_N8_d2", "Balanced block compression"),
             ("balanced_full_N8_d4", "Full balanced compression; larger budget"),
             ("raw_detector_failure_N8_d1", "Hadamards; raw Z detector loses all singleton records")]
    for name, label in names:
        f = lookup[name]
        perfect = sum(v >= 1-TOL for v in f["single_D"])
        lines.append(f"| {label} | {f['depth_ceiling']} | {f['packing_capacity']} | {perfect} | {f['averaged_Holevo_R_delta_0p1']:.8g} |")
    lines += ["", f"Disclosed circuit fixtures: **{len(fixtures)}**.", ""]
    for key, val in counts.items():
        lines.append(f"- {key}: {val}")
    lines += ["", f"Maximum identity/decoder discrepancy: `{max_error:.16g}`.",
              f"Diagnostic tolerance: `{TOL}`; no interval certification is claimed.", "",
              "At N=8, d=1 the identity, Bell and balanced-block controls use the same ceiling.",
              "d=2 and d=4 rows are separately labelled larger-budget comparisons.",
              "Hadamards preserve all optimal singleton D values at 1, but raw Z distinguishability is 0.",
              "Averaged Holevo R_0.1 uses uniform contiguous intervals at each size; packing does not.", ""]
    (ROOT / "TABLES.md").write_text("\n".join(lines))
    print(result["status"])
    print(json.dumps({"fixtures": len(fixtures), "counts": counts,
                      "maximum_discrepancy": max_error}, indent=2))


if __name__ == "__main__":
    main()
