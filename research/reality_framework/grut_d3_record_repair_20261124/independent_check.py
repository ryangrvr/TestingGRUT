"""Independent physical tensor-Hamiltonian audit of GRUT-D3.

The physical reference uses two qutrits and explicit reservoir spins. Inert
D and C are accounted for by their classical conditional histories, without
building redundant identity tensor factors. The tested implementation is
model.py; its reduced sector is never used to generate the physical oracle.
"""
from __future__ import annotations

import json
import math

import numpy as np
from scipy.linalg import eigh, expm

import model


DELTA = 4.0
F = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]], complex)
b = np.array([[0, 1], [0, 0]], complex)
number = b.conj().T @ b


def ketbra(row, col, dimension=3):
    result = np.zeros((dimension, dimension), complex)
    result[row, col] = 1
    return result


def tensor_term(dimensions, local):
    result = np.ones((1, 1), complex)
    for site, dimension in enumerate(dimensions):
        result = np.kron(result, local.get(site, np.eye(dimension)))
    return result


def physical_hamiltonian(Lw, Lb):
    dimensions = [3, 3] + [2] * (Lw + Lb)
    d = math.prod(dimensions)
    qw = tensor_term(dimensions, {1: ketbra(2, 2)})
    qb = tensor_term(dimensions, {0: ketbra(2, 2)})
    for j in range(Lw):
        qw += tensor_term(dimensions, {2 + j: number})
    for j in range(Lb):
        qb += tensor_term(dimensions, {2 + Lw + j: number})
    h = DELTA * (qw + qb)
    for start, length in [(2, Lw), (2 + Lw, Lb)]:
        for j in range(length - 1):
            term = tensor_term(dimensions, {start + j: b.conj().T, start + j + 1: b})
            h += term + term.conj().T
    for r in [0, 1]:
        write = tensor_term(dimensions, {0: ketbra(r, r), 1: ketbra(r, 2),
                                        2: b.conj().T})
        repair = tensor_term(dimensions, {0: ketbra(r, 2), 1: ketbra(r, r),
                                         2 + Lw: b.conj().T})
        h += write + write.conj().T + repair + repair.conj().T
    fm = tensor_term(dimensions, {1: F})
    fs = tensor_term(dimensions, {0: F})
    assert h.shape == (d, d)
    return h, qw, qb, fm, fs


def basis_index(S, M, Lw, Lb, x=0, y=0):
    spin_bits = 0
    if x:
        spin_bits += 1 << (Lw + Lb - x)
    if y:
        spin_bits += 1 << (Lb - y)
    return (3 * S + M) * 2**(Lw + Lb) + spin_bits


def grid_embedding(Lw, Lb, selected):
    dimension = 9 * 2**(Lw + Lb)
    embedding = np.zeros((dimension, (Lw + 1) * (Lb + 1)), complex)
    for x in range(Lw + 1):
        for y in range(Lb + 1):
            S = selected if y else 2
            M = selected if x else 2
            embedding[basis_index(S, M, Lw, Lb, x, y), x * (Lb + 1) + y] = 1
    return embedding


def writer_embedding(Lw, Lb, selected):
    dimension = 9 * 2**(Lw + Lb)
    embedding = np.zeros((dimension, Lw + 1), complex)
    for x in range(Lw + 1):
        M = selected if x else 2
        embedding[basis_index(selected, M, Lw, Lb, x), x] = 1
    return embedding


def independent_grid(Lw, Lb):
    side = Lb + 1
    dimension = (Lw + 1) * side
    result = np.zeros((dimension, dimension))
    for x in range(Lw + 1):
        for y in range(Lb + 1):
            origin = x * side + y
            if x < Lw and not (x == 0 and y == 0):
                target = (x + 1) * side + y
                result[origin, target] = result[target, origin] = 1
            if y < Lb and not (x == 0 and y == 0):
                target = x * side + y + 1
                result[origin, target] = result[target, origin] = 1
    return result


def physical_scores(vector, historical):
    probabilities = abs(vector.reshape(3, 3, -1))**2
    sm = probabilities.sum(axis=2)
    flags = np.zeros((2, 2))
    for S in range(3):
        for M in range(3):
            flags[int(M != 2), int(S != 2)] += sm[S, M]
    return {
        "joint": float(sm[historical, historical]),
        "restore": float(sm[historical].sum()),
        "ordinary": float(sm[:2].sum()),
        "flags": flags,
    }


def nonmemory_density(vector):
    shaped = vector.reshape(3, 3, -1)
    result = np.einsum("smb,tmd->sbtd", shaped, shaped.conj())
    return result.reshape(3 * shaped.shape[2], 3 * shaped.shape[2])


def energy_groups(eigenvalues, tolerance=1e-9):
    groups = []
    for index, value in enumerate(eigenvalues):
        if not groups or abs(value - eigenvalues[groups[-1][0]]) > tolerance:
            groups.append([index])
        else:
            groups[-1].append(index)
    return groups


def spectral_law(vector, vectors, groups):
    probabilities = abs(vectors.conj().T @ vector)**2
    return np.array([probabilities[group].sum() for group in groups])


def expectation(h, vector):
    return float(np.vdot(vector, h @ vector).real)


def main():
    error = {name: 0.0 for name in [
        "charge_commutators", "color_symmetry_commutator", "grid_matrix",
        "grid_invariance", "common_corner_intersection", "writer_matrix",
        "writer_invariance", "writer_state", "write_accuracy", "path_state",
        "post_fault_state", "conditional_energy_spectral_distribution",
        "control_equals_symmetry_on_erased_support", "unconditional_active_plus_coin_state",
        "nonmemory_joint_marginal", "common_fault_mean_work", "terminal_flags",
        "terminal_intact_scores", "wrong_historical_scores", "coin_half_score_law",
        "ordinary_score_equality", "probability_normalization",
    ]}

    def update(name, value):
        error[name] = max(error[name], float(value))

    cases = []
    point_count = 0
    writing_times = [0.37, 1.11, 2.37]
    repair_times = [0.29, 0.83, 1.41]
    for Lw, Lb in [(1, 1), (1, 2), (2, 1), (2, 2), (3, 2), (2, 3)]:
        h, qw, qb, fm, fs = physical_hamiltonian(Lw, Lb)
        global_flip = fs @ fm
        for q in [qw, qb]:
            update("charge_commutators", np.max(abs(h @ q - q @ h)))
        update("color_symmetry_commutator", np.max(abs(h @ global_flip - global_flip @ h)))
        grid = independent_grid(Lw, Lb)
        supplied_grid = model.repair_adjacency(Lw, Lb)
        supplied_grid = supplied_grid.toarray() if hasattr(supplied_grid, "toarray") else np.asarray(supplied_grid)
        update("grid_matrix", np.max(abs(grid - supplied_grid)))
        embeddings = [grid_embedding(Lw, Lb, r) for r in [0, 1]]
        reduced_h = grid + 2 * DELTA * np.eye(len(grid))
        for embedding in embeddings:
            update("grid_matrix", np.max(abs(embedding.conj().T @ h @ embedding - reduced_h)))
            update("grid_invariance", np.max(abs(h @ embedding - embedding @ reduced_h)))
        intersection = embeddings[0].conj().T @ embeddings[1]
        correct_intersection = np.zeros_like(intersection)
        correct_intersection[0, 0] = 1
        update("common_corner_intersection", np.max(abs(intersection - correct_intersection)))
        path_h = np.diag(np.ones(Lw), 1) + np.diag(np.ones(Lw), -1)
        for r in [0, 1]:
            writer = writer_embedding(Lw, Lb, r)
            writer_h = path_h + DELTA * np.eye(Lw + 1)
            update("writer_matrix", np.max(abs(writer.conj().T @ h @ writer - writer_h)))
            update("writer_invariance", np.max(abs(h @ writer - writer @ writer_h)))

        eigenvalues, eigenvectors = eigh(h)
        groups = energy_groups(eigenvalues)

        def evolve(vector, t):
            return eigenvectors @ (np.exp(-1j * eigenvalues * t) * (eigenvectors.conj().T @ vector))

        for w in writing_times:
            source_path = expm(-1j * w * path_h)[:, 0]
            update("path_state", np.max(abs(source_path - model.path_state(Lw, w))))
            grid_initial = np.zeros(len(grid), complex)
            grid_initial[::Lb + 1] = source_path
            update("post_fault_state", np.max(abs(grid_initial - model.state_after_fault(Lw, Lb, w))))
            post = []
            for r in [0, 1]:
                initial = writer_embedding(Lw, Lb, r)[:, 0]
                written = evolve(initial, w)
                predicted = np.exp(-1j * DELTA * w) * writer_embedding(Lw, Lb, r) @ source_path
                update("writer_state", np.max(abs(written - predicted)))
                sm = (abs(written.reshape(3, 3, -1))**2).sum(2)
                update("write_accuracy", abs(sm[:, r].sum() - (1 - abs(source_path[0])**2)))
                # This is the exact conditional representation of SWAP(S,D)
                # with D initially e and S known to be r in this history.
                fault_map = tensor_term([3, 3] + [2] * (Lw + Lb), {0: ketbra(2, r)})
                disturbed = fault_map @ written
                update("post_fault_state", np.max(abs(disturbed
                    - np.exp(-1j * DELTA * w) * embeddings[r] @ grid_initial)))
                update("common_fault_mean_work", abs(expectation(h, disturbed)
                    - expectation(h, written) - DELTA))
                mirrored = fm @ disturbed
                update("control_equals_symmetry_on_erased_support", np.max(abs(
                    mirrored - global_flip @ disturbed)))
                update("conditional_energy_spectral_distribution", np.max(abs(
                    spectral_law(disturbed, eigenvectors, groups)
                    - spectral_law(mirrored, eigenvectors, groups))))
                update("nonmemory_joint_marginal", np.max(abs(nonmemory_density(disturbed)
                    - nonmemory_density(mirrored))))
                post.append(disturbed)

            baseline_density = sum(np.outer(v, v.conj()) for v in post) / 2
            flipped_density = fm @ baseline_density @ fm.conj().T
            update("unconditional_active_plus_coin_state", np.max(abs(
                baseline_density - flipped_density)))
            # For a fair classical C, the controlled arm is block diagonal:
            # each C block is one half of one of these equal active densities.
            for coin in [0, 1]:
                controlled_density = baseline_density if coin == 0 else flipped_density
                update("unconditional_active_plus_coin_state", np.max(abs(
                    0.5 * controlled_density - 0.5 * baseline_density)))

            for u in repair_times:
                reference_grid = expm(-1j * u * grid) @ grid_initial
                expected_flags = np.zeros((2, 2))
                for x in range(Lw + 1):
                    for y in range(Lb + 1):
                        expected_flags[int(x > 0), int(y > 0)] += abs(
                            reference_grid[x * (Lb + 1) + y])**2
                supplied_flags = model.outcome_probabilities(reference_grid, Lw, Lb)
                update("terminal_flags", np.max(abs(expected_flags - supplied_flags)))
                for r in [0, 1]:
                    intact = physical_scores(evolve(post[r], u), r)
                    wrong = physical_scores(evolve(fm @ post[r], u), r)
                    for scores in [intact, wrong]:
                        update("terminal_flags", np.max(abs(scores["flags"] - expected_flags)))
                        update("probability_normalization", abs(scores["flags"].sum() - 1))
                    update("terminal_intact_scores", max(
                        abs(intact["joint"] - expected_flags[1, 1]),
                        abs(intact["restore"] - expected_flags[:, 1].sum())))
                    update("wrong_historical_scores", max(abs(wrong["joint"]), abs(wrong["restore"])))
                    scrambled = {key: (intact[key] + wrong[key]) / 2
                                 for key in ["joint", "restore", "ordinary"]}
                    update("coin_half_score_law", max(
                        abs(scrambled["joint"] - intact["joint"] / 2),
                        abs(scrambled["restore"] - intact["restore"] / 2)))
                    update("ordinary_score_equality", max(
                        abs(intact["ordinary"] - wrong["ordinary"]),
                        abs(intact["ordinary"] - scrambled["ordinary"])))
                point_count += 1
        cases.append({"writer_length": Lw, "repair_length": Lb,
                      "physical_active_dimension": len(h),
                      "reduced_conditional_dimension": len(grid),
                      "energy_eigenspaces": len(groups)})

    # Exact forecasted negative control: a three-vertex path transfers the
    # repair excitation while returning the writer to its blank endpoint.
    h, _, _, _, fs = physical_hamiltonian(1, 1)
    initial = writer_embedding(1, 1, 0)[:, 0]
    w, u = math.pi / 2, math.pi / math.sqrt(2)
    written = expm(-1j * w * h) @ initial
    fault_map = tensor_term([3, 3, 2, 2], {0: ketbra(2, 0)})
    repaired = expm(-1j * u * h) @ fault_map @ written
    negative = physical_scores(repaired, 0)
    assert abs(negative["restore"] - 1) < 1e-10 and abs(negative["joint"]) < 1e-10
    assert abs(negative["flags"][0, 1] - 1) < 1e-10

    # Additional adversarial diagnostic, outside the scored erasure protocol:
    # an unflagged 0<->1 fault leaves S ordinary and is not corrected. Here the
    # original preparation label r, rather than a SWAP-populated D, is scored.
    unflagged = expm(-1j * 1.234 * h) @ fs @ written
    unflagged_score = physical_scores(unflagged, 0)
    assert abs(unflagged_score["restore"]) < 1e-10
    assert abs(unflagged_score["ordinary"] - 1) < 1e-10

    if max(error.values()) >= 1e-10:
        raise AssertionError(f"Physical oracle disagreement: {error}")
    result = {
        "reviewer": "separately tasked internal AI; no external RAI execution or human peer review",
        "oracle": "independent physical tensor Hamiltonians; dense eigensolver/exponentials; explicit two-label and coin histories",
        "cases": cases, "terminal_points": point_count,
        "nonspecial_write_times": writing_times, "nonspecial_repair_times": repair_times,
        "worst_absolute_errors": error, "all_small_oracle_checks_pass": True,
        "short_chain_joint_gate_negative_control": {
            "writer_length": 1, "repair_length": 1, "writing_time": w,
            "repair_duration": u, "restore": negative["restore"],
            "restore_and_retain": negative["joint"],
            "memory_blank_and_system_restored": float(negative["flags"][0, 1]),
        },
        "additional_unflagged_fault_diagnostic": {
            "scope": "alternate fault, scored against known initial preparation label; not the common SWAP-erasure experiment",
            "restore": unflagged_score["restore"], "ordinary": unflagged_score["ordinary"],
            "verdict": "no correction of an unflagged color flip",
        },
        "shared_corner": "color-sector embeddings intersect in exactly the common isolated corner; they are not an orthogonal direct sum",
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
