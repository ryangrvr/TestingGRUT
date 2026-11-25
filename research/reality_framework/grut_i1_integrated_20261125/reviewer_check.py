#!/usr/bin/env python3
"""Separate density-matrix cross-check of the frozen I1 coarsest grid.

Registry order is B,M,Q,R. Each source preparation is retained as a
conditional density matrix. Bath Kraus operators come from the explicitly
exponentiated two-site exchange or XX coupling.
"""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix, identity, kron
from scipy.sparse.linalg import expm_multiply


HERE = Path(__file__).resolve().parent
ID = np.eye(2)
X = np.array([[0, 1], [1, 0]])
Z = np.diag([1, -1])
LOWER = np.array([[0, 1], [0, 0]])
ONE = np.diag([0, 1])
BELL = np.array([1, 0, 0, 1]) / math.sqrt(2)
EFFECT = np.outer(BELL, BELL)


def product(matrices):
    result = np.array([[1]])
    for matrix in matrices:
        result = np.kron(result, matrix)
    return result


def on_site(matrix, index):
    return product([matrix if k == index else ID for k in range(4)])


SOURCE_SIGN = on_site(Z, 0)
MEMORY_SIGN = on_site(Z, 1)
Q_EXCITED = on_site(ONE, 2)
FIELD = product([Z, ID, Z, ID])
RESONANT_EXCHANGE = np.kron(LOWER.T, LOWER) + np.kron(LOWER, LOWER.T)


def initial_states(contrast):
    theta = math.asin(contrast) / 2
    decode = expm(-1j * math.pi * X / 4)
    plus = np.ones(2) / math.sqrt(2)
    states = []
    for bit in range(2):
        probe = decode @ expm(-1j * (-1) ** bit * theta * Z) @ plus
        source = np.diag([1 - bit, bit])
        memory = np.diag(abs(probe) ** 2)
        states.append(product([source, memory, EFFECT]))
    return np.array(states)


def bath_kraus(generator, angle, excited, target):
    unitary = expm(-1j * angle * generator).reshape(2, 2, 2, 2)
    return np.array([
        [on_site(unitary[:, out, :, before] * math.sqrt(weight), target)
         for before, weight in enumerate((1 - excited, excited))]
        for out in range(2)
    ])


def apply_collision(rhos, kraus):
    result = np.zeros_like(rhos)
    bath = np.zeros((2, 2), dtype=complex)
    mean = rhos.mean(axis=0)
    for out, initial in itertools.product(range(2), repeat=2):
        left = kraus[out, initial]
        result += left @ rhos @ left.conj().T
        for other in range(2):
            bath[out, other] += np.trace(left @ mean @ kraus[other, initial].conj().T)
    return result, bath


def entropy(rho):
    values = np.linalg.eigvalsh(rho)
    return -math.fsum(p * math.log(p) for p in values if p > 1e-15)


def expectation(effect, state):
    return float(np.trace(effect @ state).real)


def exactly_evolve(initial, severity, steps, duration, regime, hot, cold):
    dt = duration / steps
    field = expm(-1j * severity * dt / duration * FIELD)
    flips = []
    for key, target in (("source_rate", 0), ("memory_rate", 1)):
        angle = math.asin(math.sqrt(regime[key] * dt))
        flips.append(bath_kraus(np.kron(X, X), angle, 0.5, target))
    thermal = []
    for key, excited in (("hot_rate", hot), ("cold_rate", cold)):
        angle = math.asin(math.sqrt(regime[key] * dt))
        thermal.append((bath_kraus(RESONANT_EXCHANGE, angle, excited, 2), excited))
    rhos = initial.copy()
    q_energy = [0.0, 0.0]
    sigma_min = 0.0
    energy_error = 0.0
    entropy_identity_error = 0.0
    cold_count = round(steps * regime["cold_fraction"])
    for k in range(steps):
        rhos = field @ rhos @ field.conj().T
        for kraus in flips:
            rhos, _ = apply_collision(rhos, kraus)
        for stream, (kraus, excited) in enumerate(thermal):
            if stream == 1 and k >= cold_count:
                continue
            before = rhos.mean(axis=0)
            rhos, bath = apply_collision(rhos, kraus)
            after = rhos.mean(axis=0)
            delta = float(bath[1, 1].real) - excited
            q_energy[stream] += delta
            bare_balance = delta + expectation(Q_EXCITED, after - before)
            energy_error = max(energy_error, abs(bare_balance))
            beta = math.log((1 - excited) / excited)
            sigma = entropy(after) - entropy(before) + beta * delta
            sigma_min = min(sigma_min, sigma)
            relative_entropy = -entropy(bath) - float(np.dot(
                bath.diagonal().real, np.log([1 - excited, excited])))
            mutual_information = (entropy(after) + entropy(bath)
                                  - entropy(before) - entropy(np.diag([1 - excited, excited])))
            entropy_identity_error = max(entropy_identity_error,
                                         abs(sigma - relative_entropy - mutual_information))
    return rhos, dict(hot_heat=q_energy[0], cold_heat=q_energy[1],
                      energy_balance_error=energy_error, entropy_min=sigma_min,
                      entropy_identity_error=entropy_identity_error)


def scores(rhos, jitters):
    state = rhos.mean(axis=0)
    blocks = state.reshape((2, 2, 4, 2, 2, 4))
    branches = [sum(blocks[b, m, :, b, m, :] for b in range(2)) for m in range(2)]
    unflagged = sum(branches)
    # For Bell inputs, only the unnormalised |00>, |11> coherence
    # blocks and those two populations enter a diagonal recovery.
    # The review proves this recovery optimal when both labels share
    # the same physical population map with nonnegative contraction.
    blind = expectation(EFFECT, unflagged)
    conditional = []
    for rho_qr in branches:
        c = rho_qr[0, 3]
        angle = float(np.angle(c)) / 2
        branch_idle = expectation(EFFECT, rho_qr)
        branch_clean = 0.5 * (float((rho_qr[0, 0] + rho_qr[3, 3]).real) + 2 * abs(c))
        branch_jitter = []
        for eps in jitters:
            values = []
            for noise in (-eps, eps):
                unitary = np.kron(expm(-1j * (angle + noise) * Z), ID)
                values.append(expectation(EFFECT, unitary @ rho_qr @ unitary.conj().T))
            branch_jitter.append(max(branch_idle, math.fsum(values) / 2))
        conditional.append((branch_clean, branch_jitter))
    assisted = math.fsum(clean for clean, _ in conditional)
    cur = expectation(SOURCE_SIGN @ MEMORY_SIGN, state)
    retained = (expectation(MEMORY_SIGN, rhos[0]) - expectation(MEMORY_SIGN, rhos[1])) / 2
    return dict(assisted=assisted, blind=blind, ideal_gain=assisted - blind,
                excited_probability=expectation(Q_EXCITED, state),
                memory_initial_contrast=retained, memory_current_contrast=cur,
                raw_offdiag=[float(unflagged[0, 3].real), float(unflagged[0, 3].imag)],
                raw_conditional_offdiag=[[float(r[0, 3].real), float(r[0, 3].imag)] for r in branches],
                jitter=[dict(epsilon=eps, assisted=math.fsum(s[i] for _, s in conditional),
                             blind=blind) for i, eps in enumerate(jitters)])


def dissipator(operator):
    jump = csr_matrix(operator)
    rate = jump.conj().T @ jump
    unit = identity(jump.shape[0], dtype=complex)
    return kron(jump.conj(), jump) - (kron(unit, rate) + kron(rate.T, unit)) / 2


def continuous(initial, severity, duration, regime, hot, cold):
    field = csr_matrix(severity / duration * FIELD)
    unit = identity(16, dtype=complex)
    generator = -1j * (kron(unit, field) - kron(field.T, unit))
    for key, target in (("source_rate", 0), ("memory_rate", 1)):
        generator += regime[key] * dissipator(on_site(X, target))
    down = dissipator(on_site(LOWER, 2))
    up = dissipator(on_site(LOWER.T, 2))
    generator += regime["hot_rate"] * ((1 - hot) * down + hot * up)
    with_cold = generator + regime["cold_rate"] * ((1 - cold) * down + cold * up)
    vector = np.column_stack([r.reshape(-1, order="F") for r in initial])
    cold_time = duration * regime["cold_fraction"]
    if cold_time:
        vector = expm_multiply(cold_time * with_cold, vector)
    if duration > cold_time:
        vector = expm_multiply((duration - cold_time) * generator, vector)
    return np.array([column.reshape(16, 16, order="F") for column in vector.T])


def run(spec):
    duration = spec["episode_time"]
    n = spec["direct_density_slices"]
    hot, cold = [spec[f"{stream}_excited_fraction"] for stream in ("hot", "cold")]
    jitters = spec["phase_error_radians"]
    rows = []
    for (regime_name, regime), (assay, contrast), (phase, fraction_of_pi) in itertools.product(
        spec["regimes"].items(), spec["assays"].items(), spec["phase_over_pi"].items()
    ):
        angle = math.pi * fraction_of_pi
        initial = initial_states(contrast)
        final, ledger = exactly_evolve(initial, angle, n, duration, regime, hot, cold)
        continuum = continuous(initial, angle, duration, regime, hot, cold)
        row = dict(regime=regime_name, assay=assay, phase=phase, slices=n,
                   **scores(final, jitters), **ledger)
        row["continuum"] = scores(continuum, jitters)
        row["max_trace_error"] = float(np.max(abs(np.trace(final, axis1=1, axis2=2) - 1)))
        row["min_density_eigenvalue"] = min(float(np.linalg.eigvalsh(x)[0]) for x in final)
        row["continuum_trace_error"] = float(np.max(abs(np.trace(continuum, axis1=1, axis2=2) - 1)))
        rows.append(row)
    f = {r["phase"]: r for r in rows if (r["regime"], r["assay"]) == ("anchor", "weak")}
    anchor_error = max(abs(f["shallow"]["assisted"] - (1 + math.sqrt(5 / 8)) / 2),
                       abs(f["deep"]["assisted"] - 0.75),
                       abs(f["shallow"]["blind"] - (1 + 1 / math.sqrt(2)) / 2),
                       abs(f["deep"]["blind"] - 0.5))
    return dict(rows=rows, count=len(rows), anchor_error=anchor_error,
                max_energy_error=max(r["energy_balance_error"] for r in rows),
                max_entropy_identity_error=max(r["entropy_identity_error"] for r in rows),
                min_entropy_production=min(r["entropy_min"] for r in rows),
                min_eigenvalue=min(r["min_density_eigenvalue"] for r in rows),
                max_trace_error=max(r["max_trace_error"] for r in rows),
                max_continuum_trace_error=max(r["continuum_trace_error"] for r in rows))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    opts = parser.add_mutually_exclusive_group()
    opts.add_argument("--output", type=Path)
    opts.add_argument("--verify", type=Path)
    args = parser.parse_args()
    spec = json.loads((HERE / "PARAMETERS.json").read_text())
    result = run(spec)
    serialized = json.dumps(result, indent=2, allow_nan=False) + "\n"
    tol = spec["thresholds"]["density_crosscheck_abs"]
    if max(result["anchor_error"], result["max_energy_error"], result["max_entropy_identity_error"],
           result["max_trace_error"], result["max_continuum_trace_error"], -result["min_eigenvalue"],
           -result["min_entropy_production"]) > tol:
        raise SystemExit("Independent density check failed")
    if args.output:
        args.output.write_text(serialized)
    if args.verify and json.loads(args.verify.read_text()) != result:
        raise SystemExit("Independent density replay differs")
    print(json.dumps({k: v for k, v in result.items() if k != "rows"}, indent=2))


if __name__ == "__main__":
    main()
