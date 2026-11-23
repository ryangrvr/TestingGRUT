"""Independent physical-Pauli oracle for GRUT-D2; prints its review record.

This reviewer starts with MODEL.md's physical-spin terms and explicit mixtures
of computational basis states. It uses full finite matrices on small systems,
not the builder's wall reduction, to assess the reduced formulas. No external
RAI service or human reviewer is represented by this script.
"""
from __future__ import annotations

import json
import math
import numpy as np
from scipy.linalg import eigh, expm

import model


I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
N = (I - Z) / 2


def tensor_term(qubits, operators):
    answer = np.ones((1, 1), dtype=complex)
    for site in range(qubits):
        answer = np.kron(answer, operators.get(site, I))
    return answer


def physical_hamiltonian(n, L, delta=4.0, connected=True):
    q = n + L
    identity = np.eye(2**q, dtype=complex)
    h = np.zeros_like(identity)
    for i in range(n - 1):
        h += delta * (identity - tensor_term(q, {i: Z, i + 1: Z})) / 2
    for j in range(L):
        h += delta * tensor_term(q, {n + j: N})
    for i in range(n - 2):
        h += (tensor_term(q, {i + 1: X})
              - tensor_term(q, {i: Z, i + 1: X, i + 2: Z})) / 2
    if connected:
        h += (tensor_term(q, {n - 1: X, n: X})
              + tensor_term(q, {n - 2: Z, n - 1: Y, n: Y})) / 2
    for j in range(L - 1):
        h += (tensor_term(q, {n + j: X, n + j + 1: X})
              + tensor_term(q, {n + j: Y, n + j + 1: Y})) / 2
    return h


def binary_rows(qubits):
    return ((np.arange(2**qubits)[:, None]
             >> np.arange(qubits - 1, -1, -1)) & 1).astype(int)


def basis_dictionary(n, L):
    bits = binary_rows(n + L)
    walls = bits[:, :n - 1] ^ bits[:, 1:n]
    encoded = np.column_stack([bits[:, 0], walls, bits[:, n:]])
    indices = encoded @ (2**np.arange(n + L - 1, -1, -1))
    return np.argsort(indices), walls.sum(1), bits[:, n:].sum(1), bits


def wall_chain_hamiltonian(m, L, delta=4.0):
    q = m + L
    h = sum((delta * tensor_term(q, {j: N}) for j in range(q)),
            start=np.zeros((2**q, 2**q), complex))
    for j in range(q - 1):
        h += (tensor_term(q, {j: X, j + 1: X})
              + tensor_term(q, {j: Y, j + 1: Y})) / 2
    return np.kron(I, h)


def statistics(columns, walls, bath, source_k, weights):
    probability = np.abs(columns)**2
    final = probability @ weights
    joint = np.zeros((max(source_k) + 1, max(max(source_k), max(bath)) + 1))
    for col, k in enumerate(source_k):
        for nb in np.unique(bath):
            joint[k, nb] += weights[col] * probability[bath == nb, col].sum()
    return float(final[walls == 0].sum()), float(final @ walls), joint


def energy_moments(h, columns, weights):
    hc = h @ columns
    mean = float(np.sum(weights * np.sum(columns.conj() * hc, axis=0)).real)
    second = float(np.sum(weights * np.sum(abs(hc)**2, axis=0)).real)
    return mean, second - mean**2


def reduced_physical(columns, n, L, weights):
    shaped = columns.reshape(2**n, 2**L, -1)
    return np.einsum("ibc,jbc,c->ij", shaped, shaped.conj(), weights)


def matrix_entropy(density):
    eigenvalues = np.linalg.eigvalsh(density)
    positive = eigenvalues[eigenvalues > 1e-14]
    return float(-positive @ np.log2(positive))


def main():
    worst = {
        "physical_to_wall_matrix": 0.0,
        "global_flip_commutator": 0.0,
        "anchor_commutator": 0.0,
        "organization_probability": 0.0,
        "initial_count_bath_joint": 0.0,
        "kick_expected_walls": 0.0,
        "echo_organization": 0.0,
        "echo_present_system_state": 0.0,
        "echo_present_count_distribution": 0.0,
        "echo_mean_work": 0.0,
        "echo_energy_variance": 0.0,
        "subsystem_entropies": 0.0,
        "initial_energy_moments": 0.0,
        "disconnected_energy_moments": 0.0,
        "source_block_vs_matrix_exponential": 0.0,
    }
    numerical_cases = []
    determinant_counterexample = None
    echo_variance_examples = []

    def update(name, error):
        worst[name] = max(worst[name], float(error))

    for n in [2, 3, 4]:
        m = n - 1
        for L in [1, 2, 3]:
            h = physical_hamiltonian(n, L)
            inverse_map, walls, bath, bits = basis_dictionary(n, L)
            update("physical_to_wall_matrix", np.max(abs(
                h[np.ix_(inverse_map, inverse_map)] - wall_chain_hamiltonian(m, L))))
            flip = tensor_term(n + L, {j: X for j in range(n)})
            anchor = tensor_term(n + L, {0: Z})
            update("global_flip_commutator", np.max(abs(h @ flip - flip @ h)))
            update("anchor_commutator", np.max(abs(h @ anchor - anchor @ h)))
            eigenvalues, vectors = eigh(h)

            def propagator(t):
                return (vectors * np.exp(-1j * eigenvalues * t)) @ vectors.conj().T

            initial_indices = np.arange(2**n) * 2**L
            initial = np.eye(2**(n + L), dtype=complex)[:, initial_indices]
            source_k = walls[initial_indices]
            weights = np.full(2**n, 1 / 2**n)
            resources = model.initial_resources(m, L)
            mean, variance = energy_moments(h, initial, weights)
            update("initial_energy_moments", max(abs(mean - resources["mean_energy"]),
                abs(variance - resources["energy_variance"])))
            hd = physical_hamiltonian(n, L, connected=False)
            mean_d, variance_d = energy_moments(hd, initial, weights)
            update("disconnected_energy_moments", max(
                abs(mean_d - resources["disconnected_mean_energy"]),
                abs(variance_d - resources["disconnected_energy_variance"])))

            path = model.WallPath(m, L)
            hopping = np.diag(np.ones(m + L - 1), 1) + np.diag(np.ones(m + L - 1), -1)
            for t in [0.0, 0.413, 0.731, 1.618]:
                columns = propagator(t) @ initial
                p_order, expected, joint = statistics(columns, walls, bath, source_k, weights)
                a = path.block(t)
                update("organization_probability", abs(p_order - model.order_probability(a)))
                update("initial_count_bath_joint", np.max(abs(joint[:, :m + 1]
                    - model.count_joint(a))))
                update("source_block_vs_matrix_exponential", np.max(abs(
                    a - expm(-1j * t * hopping)[:m, :m])))

            tau, s = 0.731, 0.413
            before = propagator(tau) @ initial
            kick = tensor_term(n + L, {n - 1: X})
            after = propagator(s) @ kick @ before
            p_kick, e_kick, _ = statistics(after, walls, bath, source_k, weights)
            predicted = model.kicked_expected_walls(path, tau, s)
            update("kick_expected_walls", abs(e_kick - predicted))
            if p_kick + 1e-11 < 1 - predicted:
                raise AssertionError("Expected-wall union bound failed")

            # Deliberately test the invalid number-conserving determinant
            # after X. The normal covariance remains useful for N_wall but
            # generally loses anomalous information needed for P(N_wall=0).
            u1, u2 = expm(-1j * tau * hopping), expm(-1j * s * hopping)
            covariance = 0.5 * u1[:, :m] @ u1[:, :m].conj().T
            signs = np.ones(m + L)
            signs[m - 1], signs[m:] = 0, -1
            changed = signs[:, None] * covariance * signs[None, :]
            changed[m - 1, m - 1] = 1 - covariance[m - 1, m - 1]
            evolved = u2 @ changed @ u2.conj().T
            incorrect_determinant = float(np.linalg.det(np.eye(m) - evolved[:m, :m]).real)
            if m > 1 and abs(incorrect_determinant - p_kick) > 1e-7:
                determinant_counterexample = {
                    "n": n, "L": L, "tau": tau, "recovery_time": s,
                    "exact_organization": p_kick,
                    "invalid_normal_covariance_determinant": incorrect_determinant,
                    "absolute_mismatch": abs(incorrect_determinant - p_kick),
                }

            # Apply physical bath-only staggered Z phases; no source operation.
            selected = [n + j for j in range(L) if (m + j) % 2]
            phase = tensor_term(n + L, {j: Z for j in selected})
            phased = phase @ before
            before_order, _, before_joint = statistics(before, walls, bath, source_k, weights)
            phase_order, _, phase_joint = statistics(phased, walls, bath, source_k, weights)
            update("echo_present_system_state", np.max(abs(
                reduced_physical(phased, n, L, weights)
                - reduced_physical(before, n, L, weights))))
            update("echo_present_count_distribution", np.max(abs(before_joint - phase_joint)))
            e_before, var_before = energy_moments(h, before, weights)
            e_after, var_after = energy_moments(h, phased, weights)
            update("echo_mean_work", abs(e_after - e_before))
            update("echo_energy_variance", abs(var_after - model.echo_energy_variance(path, tau)))
            shaped = before.reshape(2**n, 2**L, -1)
            bath_density = np.einsum("ibc,idc,c->bd", shaped, shaped.conj(), weights)
            direct_system_entropy = matrix_entropy(reduced_physical(before, n, L, weights))
            direct_bath_entropy = matrix_entropy(bath_density)
            predicted_entropies = model.subsystem_entropies(path.block(tau))
            update("subsystem_entropies", max(
                abs(direct_system_entropy - predicted_entropies["physical_system_entropy_bits"]),
                abs(direct_bath_entropy - predicted_entropies["bath_entropy_bits"]),
                abs(direct_system_entropy + direct_bath_entropy - n
                    - predicted_entropies["system_bath_mutual_information_bits"])))
            echo_columns = propagator(tau) @ phased
            echo_order, _, echo_joint = statistics(echo_columns, walls, bath, source_k, weights)
            echo_a = model.bath_echo_block(path.block(tau))
            update("echo_organization", abs(echo_order - model.order_probability(echo_a)))
            update("initial_count_bath_joint", np.max(abs(echo_joint[:, :m + 1]
                - model.count_joint(echo_a))))
            echo_variance_examples.append({"n": n, "L": L,
                "energy_variance_change": var_after - var_before})
            numerical_cases.append({"n": n, "L": L, "dimension": len(h)})

    if determinant_counterexample is None:
        raise AssertionError("Failed to expose invalid post-kick determinant")
    if max(worst.values()) >= 1e-11:
        raise AssertionError(f"Full physical-oracle discrepancy: {worst}")

    # Falsify the preregistered 'no new count record' statement for a hot bath.
    # The active maximally mixed state is stationary, but its inert reference
    # tagging the initial wall is correlated, so those correlations can move.
    n, L = 2, 1
    h = physical_hamiltonian(n, L)
    _, walls, bath, _ = basis_dictionary(n, L)
    initial = np.eye(2**(n + L), dtype=complex)
    weights = np.full(len(h), 1 / len(h))
    _, _, hot_before = statistics(initial, walls, bath, walls, weights)
    evolved = expm(-1j * math.pi / 2 * h)
    hot_order, _, hot_after = statistics(evolved, walls, bath, walls, weights)
    hot_mi_before = model.count_information(hot_before)["mutual_information_bits"]
    hot_mi_after = model.count_information(hot_after)["mutual_information_bits"]
    assert abs(hot_order - 0.5) < 1e-12
    assert abs(hot_mi_before) < 1e-12 and abs(hot_mi_after - 1) < 1e-12
    assert np.max(abs(evolved @ evolved.conj().T / len(h) - np.eye(len(h)) / len(h))) < 1e-12

    record = {
        "reviewer": "separately tasked internal AI, not external RAI or human peer review",
        "oracle": "full physical Pauli Hamiltonian, dense eigensolver/exponential, explicit initial basis mixtures",
        "case_count": len(numerical_cases), "cases": numerical_cases,
        "worst_absolute_errors": worst,
        "all_oracle_checks_pass": True,
        "invalid_post_kick_determinant_counterexample": determinant_counterexample,
        "hot_bath_preregistered_no_record_claim": {
            "verdict": "falsified; ordering stays unchanged but reference information can transfer",
            "n": n, "L": L, "time": math.pi / 2,
            "active_state_stationary": True, "organization_probability": hot_order,
            "count_information_before_bits": hot_mi_before,
            "count_information_after_bits": hot_mi_after,
            "joint_K_NB_before": hot_before.tolist(),
            "joint_K_NB_after": hot_after.tolist(),
        },
        "echo_variance_examples": echo_variance_examples,
    }
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
