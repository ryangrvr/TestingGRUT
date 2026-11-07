"""Reproduce finite access, response, and compression checks with NumPy."""

import json
from itertools import combinations
from math import factorial, sqrt

import numpy as np


TOLERANCE = 5e-11
SEED = 20261106
IDENTITY = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
PAULIS = (X, Y, Z)


def embed(site_count, factors):
    operator = np.ones((1, 1), dtype=complex)
    for site in range(site_count):
        operator = np.kron(operator, factors.get(site, IDENTITY))
    return operator


def commutator(a, b):
    return a @ b - b @ a


def inner(a, b):
    return np.vdot(a, b) / a.shape[0]


def norm(a):
    return float(np.linalg.norm(a) / sqrt(a.shape[0]))


def recover_factors(controls):
    """Recover qubit factors from an unlabeled list of candidate involutions."""
    dimension = controls[0].shape[0]
    identity = np.eye(dimension)
    for a in controls:
        if norm(a - a.conj().T) > TOLERANCE or norm(a @ a - identity) > TOLERANCE:
            raise ValueError("A candidate control is not a Hermitian involution.")

    neighbors = [set() for _ in controls]
    for i, j in combinations(range(len(controls)), 2):
        if norm(commutator(controls[i], controls[j])) > TOLERANCE:
            neighbors[i].add(j)
            neighbors[j].add(i)

    groups = []
    visited = set()
    for i, adjacent in enumerate(neighbors):
        if i in visited:
            continue
        if len(adjacent) != 1:
            raise ValueError("The noncommutation graph is not a disjoint union of pairs.")
        j = next(iter(adjacent))
        if neighbors[j] != {i}:
            raise ValueError("The noncommutation graph is not a disjoint union of pairs.")
        a, b = controls[i], controls[j]
        if norm(a @ b + b @ a) > TOLERANCE:
            raise ValueError("A recovered pair does not anticommute.")
        groups.append(((i, j), (a, 1j * a @ b, b)))
        visited.update((i, j))

    if 2 ** len(groups) != dimension:
        raise ValueError("The pairs do not account for the Hilbert-space dimension.")
    return groups


def interaction_weight(hamiltonian, bases):
    nested = [hamiltonian]
    for basis in bases:
        nested = [commutator(term, a) for term in nested for a in basis]
    return sqrt(sum(norm(term) ** 2 for term in nested) / 8 ** len(bases))


def exact_supports(hamiltonian, bases):
    """Subtract inclusive supersets, in descending support size."""
    energies = {}
    sites = range(len(bases))
    for size in range(len(bases), 0, -1):
        for sites_in_support in combinations(sites, size):
            support = frozenset(sites_in_support)
            inclusive = interaction_weight(hamiltonian, [bases[i] for i in sites_in_support])
            energies[support] = inclusive ** 2 - sum(
                energy for superset, energy in energies.items() if support < superset
            )
    return energies


def response_derivatives(hamiltonian, sender, receiver, max_order):
    terms = list(receiver)
    strengths = []
    for order in range(max_order + 1):
        strengths.append(sqrt(sum(norm(commutator(a, b)) ** 2 for a in terms for b in sender)))
        if order < max_order:
            terms = [1j * commutator(hamiltonian, a) for a in terms]
    return strengths


def finite_time_response(hamiltonian, sender, receiver, times):
    energies, vectors = np.linalg.eigh(hamiltonian)
    values = []
    for time in times:
        unitary = (vectors * np.exp(-1j * time * energies)) @ vectors.conj().T
        evolved = [unitary.conj().T @ a @ unitary for a in receiver]
        values.append(sqrt(sum(norm(commutator(a, b)) ** 2 for a in evolved for b in sender)))
    return np.array(values)


def order_check(hamiltonian, sender, receiver, expected_order):
    derivatives = response_derivatives(hamiltonian, sender, receiver, 4)
    order = next(n for n, value in enumerate(derivatives) if value > TOLERANCE)
    assert order == expected_order
    times = np.array([0.01, 0.02, 0.04, 0.08])
    values = finite_time_response(hamiltonian, sender, receiver, times)
    slope = float(np.polyfit(np.log(times), np.log(values), deg=1)[0])
    leading = derivatives[order] / factorial(order)
    relative_error = abs(values[0] / (leading * times[0] ** order) - 1)
    assert abs(slope - order) < 0.01
    assert relative_error < 0.001
    return {
        "first_nonzero_derivative_order": order,
        "derivative_norms_orders_0_to_4": derivatives,
        "leading_response_coefficient": leading,
        "times": times.tolist(),
        "responses": values.tolist(),
        "finite_time_log_slope": slope,
        "leading_term_relative_error_at_first_time": float(relative_error),
    }


def access_check():
    sites = 5
    dimension = 2 ** sites
    rng = np.random.default_rng(SEED)
    random_matrix = rng.normal(size=(dimension, dimension)) + 1j * rng.normal(
        size=(dimension, dimension)
    )
    rotation, _ = np.linalg.qr(random_matrix)
    chain_bonds = (0.31, 0.43, 0.37, 0.29)
    fields = (0.13, -0.19, 0.07, 0.11, -0.17)
    hidden_hamiltonian = sum(
        coupling * embed(sites, {i: axis, i + 1: axis})
        for i, coupling in enumerate(chain_bonds)
        for axis in (X, Y)
    ) + sum(field * embed(sites, {i: Z}) for i, field in enumerate(fields))
    primitives = [embed(sites, {i: axis}) for i in range(sites) for axis in (X, Z)]
    permutation = rng.permutation(len(primitives))

    def rotated(operator):
        return rotation @ operator @ rotation.conj().T

    def shuffled(operators):
        return [rotated(operators[i]) for i in permutation]

    def validate_and_label(groups):
        labeled = {}
        for (i, j), basis in groups:
            # The planted labels are only consulted after the recovery.
            label = int(permutation[i] // 2)
            assert label == permutation[j] // 2
            labeled[label] = basis
        assert len(labeled) == sites
        return [labeled[i] for i in range(sites)]

    hamiltonian = rotated(hidden_hamiltonian)
    groups = recover_factors(shuffled(primitives))
    bases = validate_and_label(groups)
    weights = {
        (i, j): interaction_weight(hamiltonian, [bases[i], bases[j]])
        for i, j in combinations(range(sites), 2)
    }
    pair_residual = max(
        abs(weight - (sqrt(2) * chain_bonds[i] if j == i + 1 else 0.0))
        for (i, j), weight in weights.items()
    )
    assert pair_residual < TOLERANCE

    h_direct = hamiltonian + 0.11 * rotated(embed(sites, {0: Z, 4: Z}))
    h_triple = hamiltonian + 0.09 * rotated(embed(sites, {0: Z, 2: X, 4: Y}))
    end_sites = (bases[0], bases[4])
    direct_weight = interaction_weight(h_direct, end_sites)
    triple_pair_weight = interaction_weight(h_triple, end_sites)
    assert abs(direct_weight - 0.11) < TOLERANCE
    assert abs(triple_pair_weight - 0.09) < TOLERANCE
    exact_triple_supports = exact_supports(h_triple, bases)
    triple_exact_pair = exact_triple_supports[frozenset({0, 4})]
    triple_exact_triple = exact_triple_supports[frozenset({0, 2, 4})]
    assert abs(triple_exact_pair) < TOLERANCE
    assert abs(triple_exact_triple - 0.09 ** 2) < TOLERANCE

    # Change the supplied access algebras with H and its Gibbs state held fixed.
    labels = np.arange(dimension)
    endpoint_cnot = np.eye(dimension)[labels ^ (labels >> (sites - 1))]
    alternative = [endpoint_cnot @ a @ endpoint_cnot.T for a in primitives]
    alternate_bases = validate_and_label(recover_factors(shuffled(alternative)))
    alternate_end_sites = (alternate_bases[0], alternate_bases[4])
    alternate_weight = interaction_weight(hamiltonian, alternate_end_sites)
    expected_alternate = sqrt(2 * chain_bonds[0] ** 2 + chain_bonds[-1] ** 2 + fields[-1] ** 2)
    assert abs(alternate_weight - expected_alternate) < TOLERANCE

    energies, vectors = np.linalg.eigh(hamiltonian)
    populations = np.exp(-energies)
    state = (vectors * (populations / populations.sum())) @ vectors.conj().T
    phases = np.exp(1j * rng.uniform(0, 2 * np.pi, size=dimension))
    stabilizer = (vectors * phases) @ vectors.conj().T
    preserved_h = norm(stabilizer @ hamiltonian @ stabilizer.conj().T - hamiltonian)
    preserved_state = norm(stabilizer @ state @ stabilizer.conj().T - state)
    assert max(preserved_h, preserved_state) < TOLERANCE
    probe = stabilizer @ bases[0][0] @ stabilizer.conj().T
    basis_with_identity = (np.eye(dimension), *bases[0])
    projection = sum(inner(a, probe) * a for a in basis_with_identity)
    leakage = norm(probe - projection)
    assert leakage > 0.1

    readout = rotated(embed(sites, {0: X}))
    kick = rotated(embed(sites, {4: X}))
    identity = np.eye(dimension)
    probe_state = (
        (identity + rotated(embed(sites, {0: Y})))
        @ (identity + rotated(embed(sites, {4: Y})))
        / dimension
    )

    def kick_slope(candidate):
        readout_derivative = 1j * commutator(candidate, readout)
        return float(np.trace(probe_state @ (1j * commutator(kick, readout_derivative))).real)

    chain_slope, direct_slope = (kick_slope(candidate) for candidate in (hamiltonian, h_direct))
    assert abs(chain_slope) < TOLERANCE
    assert abs(direct_slope + 0.44) < TOLERANCE

    broken = list(primitives)
    broken[3] = embed(sites, {0: Z, 1: Z})
    try:
        recover_factors(shuffled(broken))
    except ValueError as error:
        rejection = str(error)
    else:
        raise AssertionError("A joint control was accepted as a set of local primitives.")

    return {
        "dimension": dimension,
        "recovered_groups_input_indices": [list(pair) for pair, _ in groups],
        "chain_pair_weights": {f"{i}-{j}": weight for (i, j), weight in weights.items()},
        "max_chain_pair_weight_residual": pair_residual,
        "direct_end_weight": direct_weight,
        "three_site_control_inclusive_end_weight": triple_pair_weight,
        "three_site_control_exact_end_pair_weight_squared": triple_exact_pair,
        "three_site_control_exact_triple_weight_squared": triple_exact_triple,
        "alternate_access_same_hamiltonian_end_weight": alternate_weight,
        "stabilizer_hamiltonian_residual": preserved_h,
        "stabilizer_stationary_state_residual": preserved_state,
        "stabilizer_rotated_local_operator_leakage": leakage,
        "broken_control_rejection": rejection,
        "specific_kick_initial_time_slopes": {"chain": chain_slope, "direct_end_term": direct_slope},
        "chain_endpoint_response": order_check(hamiltonian, *end_sites, expected_order=4),
        "direct_endpoint_response": order_check(h_direct, *end_sites, expected_order=1),
        "three_site_endpoint_response": order_check(h_triple, *end_sites, expected_order=1),
        "alternate_access_endpoint_response": order_check(
            hamiltonian, *alternate_end_sites, expected_order=1
        ),
    }


def compression_check():
    hamiltonian = np.kron(X, X)
    local = [np.kron(axis, IDENTITY) for axis in PAULIS]
    local += [np.kron(IDENTITY, axis) for axis in PAULIS]
    basis = [np.eye(4), *local]
    index = 0
    while index < len(basis):
        generated = 1j * commutator(hamiltonian, basis[index])
        for a in basis:
            generated -= inner(a, generated) * a
        length = norm(generated)
        if length > TOLERANCE:
            basis.append(generated / length)
        index += 1
    assert len(basis) == 11

    correlation = np.kron(Y, X)
    states = [(np.eye(4) + sign * correlation) / 4 for sign in (1, -1)]
    assert all(np.min(np.linalg.eigvalsh(state)) >= -TOLERANCE for state in states)
    difference = states[0] - states[1]
    assert max(abs(np.trace(a @ difference)) for a in local) < TOLERANCE

    distances = []
    for time in (np.pi / 8, np.pi / 4):
        unitary = np.cos(time) * np.eye(4) - 1j * np.sin(time) * hamiltonian
        evolved_difference = unitary @ difference @ unitary.conj().T
        # Trace distance on a qubit is half the Euclidean Bloch-vector separation.
        bloch_difference = [np.trace(a @ evolved_difference).real for a in local[:3]]
        distance = float(np.linalg.norm(bloch_difference) / 2)
        assert abs(distance - abs(np.sin(2 * time))) < TOLERANCE
        distances.append(distance)

    return {
        "local_expectations_excluding_trace": len(local),
        "minimal_predictive_expectations_excluding_trace": len(basis) - 1,
        "initial_both_marginals_match": True,
        "receiver_trace_distances_at_pi_over_8_and_pi_over_4": distances,
        "common_prediction_worst_case_lower_bounds": [distance / 2 for distance in distances],
    }


if __name__ == "__main__":
    print(json.dumps({
        "seed": SEED,
        "matrix_tolerance": TOLERANCE,
        "access": access_check(),
        "compression": compression_check(),
    }, indent=2, allow_nan=False))
