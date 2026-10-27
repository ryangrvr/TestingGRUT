"""Exact finite constrained history: one whole, distinguishable internal systems.

The charge constraint is only GLOBAL neutrality. This is not lattice gauge theory.
All logical gates and the initial state are supplied. See MATHEMATICS.md.
"""
import numpy as np
from common import ROOT, density, entropy, error, expectation, kron, purity, reduced, save

I2 = np.eye(2, dtype=complex)
X2 = np.array([[0, 1], [1, 0]], complex)
Y2 = np.array([[0, -1j], [1j, 0]], complex)
Z2 = np.diag([1.0, -1.0]).astype(complex)
P0 = np.diag([1.0, 0.0]).astype(complex)
P1 = I2 - P0
I16, I64 = np.eye(16, dtype=complex), np.eye(64, dtype=complex)


def build():
    # Ordering within the nonclock system: X, compensator, record, neutral probe.
    qx = 6 * kron(Z2, I2, I2, I2)
    qc = 6 * kron(I2, Z2, I2, I2)
    g = qx + qc
    copy = kron(P0, I2, I2, I2) + kron(P1, I2, X2, I2)
    phase = np.diag([1.0, 1j])
    w = np.array([
        (copy if t >= 2 else I16) @ kron(I2, I2, I2, np.linalg.matrix_power(phase, t))
        for t in range(4)
    ])
    shift = np.zeros((64, 64), complex)
    embedding = np.zeros((64, 16), complex)
    for t in range(4):
        nxt = (t + 1) % 4
        shift[nxt*16:(nxt+1)*16, t*16:(t+1)*16] = w[nxt] @ w[t].conj().T
        embedding[t*16:(t+1)*16] = w[t] / 2
    constraint = 2 * I64 - shift - shift.conj().T
    gauss = kron(np.eye(4), g)
    neutral = np.diag((np.abs(np.diag(g)) < 1e-14).astype(float))
    physical = embedding @ neutral @ embedding.conj().T
    return w, shift, embedding, constraint, gauss, physical, qx, qc


def relational(operator, tick, shift):
    """Group-averaged physical effect; sums all cyclic gauge translates."""
    seed = np.zeros((64, 64), complex)
    seed[tick*16:(tick+1)*16, tick*16:(tick+1)*16] = operator
    result, power = np.zeros((64, 64), complex), I64.copy()
    for _ in range(4):
        result += power @ seed @ power.conj().T
        power = shift @ power
    return result


def main():
    w, shift, embedding, constraint, gauss, physical, qx, qc = build()
    e = np.eye(2, dtype=complex)
    plus = (e[0] + e[1]) / np.sqrt(2)
    base_a = kron(e[0], e[1], e[0], plus)
    base_b = kron(e[1], e[0], e[0], plus)
    base = (base_a + base_b) / np.sqrt(2)
    psi = embedding @ base
    rho = density(psi)
    asymmetric = embedding @ base_a
    mixture = (density(embedding @ base_a) + density(embedding @ base_b)) / 2
    conditional = np.array([density(w[t] @ base) for t in range(4)])
    record = [kron(I2, I2, r, I2) for r in (P0, P1)]
    probe_plus = kron(I2, I2, I2, density(plus))
    probe_ops = [kron(I2, I2, I2, a) for a in (X2, Y2, Z2)]
    transfer = np.zeros((4, 4), complex)
    transfer[1, 2] = transfer[2, 1] = 1
    coherence = relational(kron(transfer, I2, I2), 0, shift)

    conditional_error = completeness_error = commutator_error = positivity_error = 0.0
    records, probe = [], []
    for t in range(4):
        row = {"tick": t,
               "X_entropy_bits": entropy(reduced(conditional[t], [2]*4, [0])),
               "record_entropy_bits": entropy(reduced(conditional[t], [2]*4, [2])),
               "conditional_remainder_purity": purity(conditional[t]),
               "X_purity": purity(reduced(conditional[t], [2]*4, [0]))}
        record_effects = [relational(r, t, shift) for r in record]
        row["record_probabilities"] = [expectation(rho, z) for z in record_effects]
        row["X_charge_conditioned_on_record"] = [
            float(np.trace(conditional[t] @ r @ qx).real / np.trace(conditional[t] @ r).real)
            if np.trace(conditional[t] @ r).real > 1e-14 else None for r in record
        ]
        records.append(row)
        p_effect = relational(probe_plus, t, shift)
        probe.append({"tick": t,
                      "bloch": [expectation(conditional[t], o) for o in probe_ops],
                      "plus_probability": expectation(rho, p_effect),
                      "purity": purity(reduced(conditional[t], [2]*4, [3]))})
        pairs = [(probe_plus, p_effect)] + list(zip(record, record_effects))
        for local, effect in pairs:
            conditional_error = max(conditional_error, abs(
                expectation(rho, effect) - expectation(conditional[t], local)))
            commutator_error = max(commutator_error, error(constraint @ effect - effect @ constraint),
                                   error(gauss @ effect - effect @ gauss))
            positivity_error = max(positivity_error, max(0., -float(np.linalg.eigvalsh(effect).min())))
        completeness_error = max(completeness_error, error(sum(record_effects) - I64))

    # Change only the clock readout, leaving the same global pure state and dynamics.
    clock_results = []
    clock_error = 0.0
    for eta in (0., .1, .25, .5, 1.):
        noise = (1-eta) * np.eye(4)[0] + eta/4
        post = sum(noise[t] * conditional[t] / 4 for t in range(4))
        post /= np.trace(post).real
        probe_rho = reduced(post, [2]*4, [3])
        predicted = (1-eta) * reduced(conditional[0], [2]*4, [3]) + eta*I2/2
        clock_error = max(clock_error, error(probe_rho - predicted))
        clock_results.append({"eta": eta, "probe_purity": purity(probe_rho),
                              "predicted_purity": 1-eta+eta**2/2,
                              "relative_coherence": float(abs(probe_rho[0, 1])*2)})

    # Genuine charge-flip symmetry of THIS finite model, not relativistic CP or T.
    flip_base = kron(X2, X2, I2, I2)
    flip = sum(np.kron(np.diag(np.eye(4)[t]), w[t] @ flip_base @ w[t].conj().T)
               for t in range(4))
    # Neutral probe has a full physical Pauli algebra after relational dressing.
    dx, dy, dz = [relational(o, 0, shift) for o in probe_ops]
    probe_algebra = max(error(dx @ dy - 1j*dz), error(dx @ dx - I64),
                        error(dy @ dy - I64), error(dz @ dz - I64))

    # Finite oscillator: exact canonical CCR fails in the highest state.
    annihilator = np.diag(np.sqrt(np.arange(1, 4)), 1)
    ladder_commutator = annihilator @ annihilator.T - annihilator.T @ annihilator
    eigen = np.linalg.eigvalsh(constraint + gauss @ gauss)
    # The physical space contains three logical factors: neutral XC pair, R, P.
    logical = embedding[:, 4:12]
    logical_probe_error = error(logical.conj().T @ dx @ logical - kron(I2, I2, X2))
    results = {
        "status": "exact finite example; supplied gates and state, not new empirical physics",
        "kinematic_dimension": 64, "physical_dimension": int(np.sum(np.abs(eigen) < 1e-10)),
        "global_purity": purity(rho), "constraint_residual": error(constraint @ psi),
        "neutrality_residual": error(gauss @ psi),
        "shift_unitarity_error": error(shift @ shift.conj().T - I64),
        "constraint_neutrality_commutator": error(constraint @ gauss - gauss @ constraint),
        "physical_projector_error": error(physical @ physical - physical),
        "embedding_isometry_error": error(embedding.conj().T @ embedding - I16),
        "physical_effect_commutator_error": commutator_error,
        "physical_effect_positivity_error": positivity_error,
        "effect_completeness_error": completeness_error,
        "conditional_dictionary_error": conditional_error,
        "neutral_probe_pauli_algebra_error": probe_algebra,
        "physical_logical_dimensions": [2, 2, 2],
        "logical_isometry_error": error(logical.conj().T @ logical - np.eye(8)),
        "logical_projector_error": error(logical @ logical.conj().T - physical),
        "logical_probe_dictionary_error": logical_probe_error,
        "records": records, "probe_evolution": probe,
        "whole_coherence_pure": expectation(rho, coherence),
        "whole_coherence_mixture": expectation(mixture, coherence),
        "whole_coherence_constraint_commutator": error(constraint @ coherence - coherence @ constraint),
        "whole_coherence_neutrality_commutator": error(gauss @ coherence - coherence @ gauss),
        "unconditional_X_charge": expectation(rho, kron(np.eye(4), qx)),
        "internal_asymmetric_X_charge": expectation(density(asymmetric), kron(np.eye(4), qx)),
        "internal_asymmetric_constraint_error": error(constraint @ asymmetric),
        "internal_asymmetric_neutrality_error": error(gauss @ asymmetric),
        "charge_flip_constraint_error": error(constraint @ flip - flip @ constraint),
        "symmetric_state_flip_error": error(flip @ psi - psi),
        "clock_readouts": clock_results, "clock_formula_error": clock_error,
        "dimension_one_pauli_anticommutator_min_residual": min(abs(2*x*z) for x in (-1, 1) for z in (-1, 1)),
        "finite_ladder_commutator_diagonal": np.diag(ladder_commutator).tolist(),
        "finite_ladder_commutator_trace": float(np.trace(ladder_commutator)),
        "physical_nbar": None, "vacuum_density_prediction": None,
        "full_GR_or_local_gauge_model": False,
    }
    save("WHOLE_PARTS_RESULTS.json", results)
    np.savez_compressed(ROOT / "HISTORY_MODEL.npz", shift=shift, constraint=constraint,
                        global_neutrality=gauss, physical_projector=physical, embedding=embedding,
                        state=psi, conditional_states=conditional, unitaries=w, logical_isometry=logical)
    print({key: results[key] for key in ("kinematic_dimension", "physical_dimension", "constraint_residual",
                                      "conditional_dictionary_error", "whole_coherence_pure")})


if __name__ == "__main__":
    main()
