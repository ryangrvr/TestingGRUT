"""Full-matrix spectral evolution and diagnostics of the supplied access."""

from itertools import combinations
from math import comb

import numpy as np
from scipy.sparse import csr_matrix, eye, kron

import model as m


class FullEvolution:
    def __init__(self, links):
        self.h = m.hamiltonian(links)
        self.energies, self.vectors = np.linalg.eigh(self.h)

    def evolve(self, state, time):
        phase = np.exp(-1j * self.energies * time)
        return self.vectors @ (phase[:, None] * (self.vectors.conj().T @ state))

    def energy_distribution(self, state):
        amplitudes = self.vectors.conj().T @ state
        levels = np.rint(self.energies / m.OMEGA + m.LINKS / 2).astype(int)
        return np.bincount(levels, np.sum(abs(amplitudes) ** 2, axis=1), minlength=m.LINKS + 1)


def path_history(links, work):
    """Conditional-state identity, kept separate from full evolution."""
    history = [work]
    for link in links:
        history.append(link @ history[-1])
    return history


def clock_binomial(time):
    p = np.sin(m.OMEGA * time / 2) ** 2
    return np.array([comb(m.LINKS, k) * p ** k * (1 - p) ** (m.LINKS - k)
                     for k in range(m.LINKS + 1)])


def primitive_components():
    clock_dim = m.LINKS + 1
    probes = []
    for k in range(clock_dim):
        project = csr_matrix(([1], ([k], [k])), shape=(clock_dim, clock_dim))
        probes.append((f"Q{k}", kron(project, eye(m.DIM), format="csr")))
    for k in range(m.LINKS):
        link = csr_matrix(([1, 1], ([k, k + 1], [k + 1, k])), shape=(clock_dim, clock_dim))
        probes.append((f"B{k}", kron(link, eye(m.DIM), format="csr")))
    for name in m.REGISTERS:
        for axis, pauli in (("X", m.X), ("Z", m.Z)):
            probes.append((name + ":" + axis, kron(
                eye(clock_dim), csr_matrix(m.work_operator(**{name: pauli})), format="csr"
            )))
    edges = [set() for _ in probes]
    for i, j in combinations(range(len(probes)), 2):
        a, b = probes[i][1], probes[j][1]
        if np.linalg.norm((a @ b - b @ a).data) > m.TOL:
            edges[i].add(j)
            edges[j].add(i)
    components = []
    visited = set()
    for start in range(len(probes)):
        if start in visited:
            continue
        stack = [start]
        visited.add(start)
        group = []
        while stack:
            i = stack.pop()
            group.append(probes[i][0])
            for j in edges[i] - visited:
                visited.add(j)
                stack.append(j)
        components.append(sorted(group))
    return components


def pauli_coefficients(matrix):
    count = len(m.REGISTERS)
    # Interleave each qubit's row and column so each transform is a local trace.
    axes = [axis for i in range(count) for axis in (i, count + i)]
    coefficients = matrix.reshape((2,) * (2 * count)).transpose(axes).reshape((4,) * count)
    transform = np.stack((m.I, m.X, m.Y, m.Z)).conj().reshape(4, 4) / 2
    for axis in range(count):
        coefficients = np.moveaxis(np.tensordot(transform, coefficients, axes=(1, axis)), 0, axis)
    return coefficients


def interaction_supports(h):
    squared_coefficients = np.zeros((4,) * len(m.REGISTERS))
    for k in range(m.LINKS):
        clock_block = h[(k + 1) * m.DIM:(k + 2) * m.DIM, k * m.DIM:(k + 1) * m.DIM]
        squared_coefficients += 2 * abs(pauli_coefficients(clock_block)) ** 2 / (m.LINKS + 1)
    weights_squared = {}
    for axes in np.ndindex(squared_coefficients.shape):
        support = ("C",) + tuple(site for site, axis in zip(m.REGISTERS, axes) if axis)
        weights_squared[support] = weights_squared.get(support, 0.0) + squared_coefficients[axes]
    supports = {"+".join(s): float(np.sqrt(w)) for s, w in weights_squared.items()
                if w > m.CARD["support_weight_minimum"] ** 2}
    edges = sorted({"+".join(sorted(pair)) for s in supports for pair in combinations(s.split("+"), 2)})
    return {
        "weights": supports,
        "pair_graph": edges,
        "sum_weights_squared": float(sum(weights_squared.values())),
        "direct_normalized_trace_h_squared": float(np.vdot(h, h).real / h.shape[0]),
    }


def predictive_rank(links, program):
    sector = m.program_mask(program)
    bits = np.stack(m.READOUT_BITS, axis=1)[sector]
    effects = [np.diag(np.all(bits == signs, axis=1)) for signs in np.ndindex(2, 2, 2)]
    history = path_history([u[np.ix_(sector, sector)] for u in links], np.eye(len(bits)))
    vectors = np.array([(v.conj().T @ effect @ v).ravel() for v in history for effect in effects])
    singular_values = np.linalg.svd(vectors, compute_uv=False)
    rank = int(np.count_nonzero(singular_values > m.CARD["svd_rank_tolerance"]))
    return {
        "rank_including_identity": rank,
        "smallest_retained_singular_value": float(singular_values[rank - 1]),
        "largest_discarded_singular_value": float(singular_values[rank]),
    }
