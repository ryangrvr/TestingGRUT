#!/usr/bin/env python3
"""Finite checks for PHOTON_BRIDGE_REVIEW.md; no continuum numerical claim."""

import argparse
import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import expm
from scipy.special import eval_genlaguerre, gammaln


def displacement(z, cutoff):
    """Exact block <m|exp(z b^dagger-z* b)|n>, NOT exp of a cut generator."""
    out = np.empty((cutoff + 1, cutoff + 1), dtype=complex)
    mu = abs(z) ** 2
    for m in range(cutoff + 1):
        for n in range(cutoff + 1):
            lo, hi = min(m, n), max(m, n)
            scale = np.exp(-mu / 2 + (gammaln(lo + 1) - gammaln(hi + 1)) / 2)
            power = z ** (m - n) if m >= n else (-z.conjugate()) ** (n - m)
            out[m, n] = scale * power * eval_genlaguerre(lo, hi - lo, mu)
    return out


def instrument(w, phase):
    ident = np.eye(w.shape[0], dtype=complex)
    kp = (ident + np.exp(1j * phase) * w) / 2
    km = (ident - np.exp(1j * phase) * w) / 2
    deficit = ident - kp.conj().T @ kp - km.conj().T @ km
    deficit = (deficit + deficit.conj().T) / 2
    values, vectors = np.linalg.eigh(deficit)
    assert values.min() > -3e-13, values.min()
    ko = (vectors * np.sqrt(np.maximum(values, 0))) @ vectors.conj().T
    closure = kp.conj().T @ kp + km.conj().T @ km + ko.conj().T @ ko
    assert np.max(np.abs(closure - ident)) < 3e-13
    return {"+": kp, "-": km, "o": ko}


def histories(instruments, initial, propagators=None):
    if propagators is None:
        propagators = [np.eye(len(initial)) for _ in instruments]
    branches = {"": initial}
    for inst, u in zip(instruments, propagators):
        branches = {
            word + outcome: op @ (u @ vec)
            for word, vec in branches.items()
            for outcome, op in inst.items()
        }
    probs = {word: float(np.vdot(vec, vec).real) for word, vec in branches.items()}
    assert abs(sum(probs.values()) - 1) < 2e-12
    return probs, branches


def coherent_overlap(beta, gamma):
    return np.exp(-(abs(beta) ** 2 + abs(gamma) ** 2) / 2 + beta.conjugate() * gamma)


def exact_branches(zs, phases):
    """Uncut reference from coherent overlaps and Weyl phases, without matrices."""
    branches = {"": [(1 + 0j, 0j)]}
    for z, phi in zip(zs, phases):
        nxt = {}
        for word, terms in branches.items():
            for label, sign in [("+", 1), ("-", -1)]:
                out = []
                for coeff, beta in terms:
                    out.append((coeff / 2, beta))
                    weyl_phase = np.exp((z * beta.conjugate() - z.conjugate() * beta) / 2)
                    out.append((sign * np.exp(1j * phi) * coeff * weyl_phase / 2, beta + z))
                nxt[word + label] = out
        branches = nxt
    probs = {}
    for word, terms in branches.items():
        prob = sum(
            b.conjugate() * c * coherent_overlap(beta, gamma)
            for b, beta in terms
            for c, gamma in terms
        )
        assert abs(prob.imag) < 1e-12 and prob.real > -1e-13
        probs[word] = float(prob.real)
    assert abs(sum(probs.values()) - 1) < 2e-12
    return probs, branches


def coherent_block(terms, cutoff):
    vec = np.zeros(cutoff + 1, dtype=complex)
    for coeff, beta in terms:
        component = np.exp(-abs(beta) ** 2 / 2)
        vec[0] += coeff * component
        for n in range(1, cutoff + 1):
            component *= beta / math.sqrt(n)
            vec[n] += coeff * component
    return vec


def check_displacement():
    dimension, cutoff, z = 64, 12, 0.8 + 0.35j
    lowering = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    large = expm(z * lowering.conj().T - z.conjugate() * lowering)
    exact = displacement(z, cutoff)
    err = float(np.max(np.abs(large[: cutoff + 1, : cutoff + 1] - exact)))
    assert err < 1e-12
    small = lowering[: cutoff + 1, : cutoff + 1]
    wrong = expm(z * small.conj().T - z.conjugate() * small)
    wrong_err = float(np.linalg.norm(wrong - exact))
    assert wrong_err > 0.1
    return {"large_exponential_entry_error": err, "exponentiate_compressed_generator_error": wrong_err}


def check_axial():
    spacing = 2 * np.pi / 9
    u = (2 * np.sin(spacing / 2) / spacing) ** 2
    omega_ax = np.sqrt(u)
    dts = [0.3, 0.7, 0.2]
    zs = [0.75 + 0.20j, -0.45 + 0.60j, 0.1 - 0.8j]
    phases = [0.4, -0.8, np.pi / 2]
    rows = []
    for cutoff in [1, 2, 4]:
        d = cutoff + 1
        for spectator in ["vacuum", "nonstationary"]:
            distributions = []
            for eta in [0, 1]:
                omega_rest = np.sqrt(2 * u * (1 + eta * spacing ** 2 * u))
                energies = (omega_ax * np.arange(d)[:, None]
                            + omega_rest * np.arange(d)[None, :]).ravel()
                initial_ax = np.eye(d, dtype=complex)[:, 0]
                initial_rest = initial_ax.copy()
                if spectator == "nonstationary":
                    initial_rest = np.array([(0.3 + 0.1j) ** n for n in range(d)])
                    initial_rest /= np.linalg.norm(initial_rest)
                initial = np.kron(initial_ax, initial_rest)
                # Independently complete the FULL finite instrument, including the spectator.
                instruments = [
                    instrument(np.kron(displacement(z, cutoff), np.eye(d)), phi)
                    for z, phi in zip(zs, phases)
                ]
                propagators = [np.diag(np.exp(-1j * energies * dt)) for dt in dts]
                probs, _ = histories(instruments, initial, propagators)
                distributions.append(probs)
            max_diff = max(abs(distributions[0][w] - distributions[1][w]) for w in distributions[0])
            overflow = sum(p for w, p in distributions[0].items() if "o" in w)
            assert max_diff < 2e-12 and overflow > 1e-6
            row = {"cutoff": cutoff, "spectator": spectator, "max_history_difference": max_diff,
                   "any_overflow_probability": overflow}
            if cutoff == 2 and spectator == "vacuum":
                row["all_history_probabilities_eta_0"] = distributions[0]
                row["all_history_probabilities_eta_1"] = distributions[1]
            rows.append(row)
    return rows


def check_basis_control():
    z, cutoff = 1 + 0j, 1
    initial = np.eye(4, dtype=complex)[:, 0]
    w_original = np.kron(displacement(z, cutoff), np.eye(2))
    w_rotated = np.kron(displacement(z / np.sqrt(2), cutoff),
                        displacement(z / np.sqrt(2), cutoff))
    probs = []
    for w in [w_original, w_rotated]:
        op = instrument(w, 0.2)["o"]
        probs.append(float(np.linalg.norm(op @ initial) ** 2))
    expected = [0.5 * (1 - np.exp(-1) * 2), 0.5 * (1 - np.exp(-1) * 2.25)]
    assert np.max(np.abs(np.array(probs) - expected)) < 2e-14
    assert abs(probs[0] - probs[1]) > 0.04
    return {"overflow_original": probs[0], "overflow_rotated": probs[1],
            "analytic_values": expected, "difference": probs[0] - probs[1]}


def check_occupation_convergence():
    times = [0, 0.6, 1.1]
    zs = [z * np.exp(1j * t) for z, t in zip(
        [1.1 + 0.4j, -0.4 + 0.9j, 0.8 - 0.7j], times)]
    phases = [0.2, np.pi / 2, -0.3]
    ref_probs, ref_terms = exact_branches(zs, phases)
    b = float(sum(np.cumsum(np.abs(zs))))
    constant = 2 ** (len(zs) + 1) * b
    rows = []
    for cutoff in [0, 1, 2, 4, 8, 16, 24]:
        initial = np.eye(cutoff + 1, dtype=complex)[:, 0]
        instruments = [instrument(displacement(z, cutoff), p) for z, p in zip(zs, phases)]
        probs, branches = histories(instruments, initial)
        tv = 0.5 * sum(abs(p - ref_probs.get(w, 0)) for w, p in probs.items())
        max_vec_error = 0
        for word, terms in ref_terms.items():
            inner = np.vdot(branches[word], coherent_block(terms, cutoff))
            squared = probs[word] + ref_probs[word] - 2 * inner.real
            assert squared > -2e-12
            max_vec_error = max(max_vec_error, np.sqrt(max(0, squared)))
        vector_bound = b / np.sqrt(cutoff + 1)
        assert max_vec_error <= vector_bound + 3e-8
        assert tv <= min(1, constant / np.sqrt(cutoff + 1)) + 2e-12
        rows.append({"cutoff": cutoff, "tv_to_uncut": float(tv),
                     "any_overflow_probability": sum(p for w, p in probs.items() if "o" in w),
                     "max_ordinary_branch_vector_error": float(max_vec_error),
                     "vector_bound": float(vector_bound),
                     "tv_bound_unclipped": float(constant / np.sqrt(cutoff + 1))})
    assert rows[-1]["tv_to_uncut"] < 2e-9
    return {"reference_method": "exact coherent-state overlaps, no occupation cutoff",
            "reference_probabilities": ref_probs, "C_occ": constant, "rows": rows}


def check_commutators():
    s, t = 1.1, 0.9
    rows = []
    for cutoff in [1, 2, 4, 8, 12]:
        d = cutoff + 1
        # W=exp(i s q_1), V=exp(i t p_2), uncut [q_1,p_2]=0.
        w = np.kron(displacement(0.5j * s, cutoff), displacement(0.5j * s, cutoff))
        v = np.kron(displacement(-0.5 * t + 0j, cutoff),
                    displacement(0.5 * t + 0j, cutoff))
        comm = w @ v - v @ w
        vacuum_norm = float(np.linalg.norm(comm[:, 0]))
        op_norm = float(np.linalg.norm(comm, 2))
        bound = (s + t) / np.sqrt(2 * (cutoff + 1))
        assert vacuum_norm <= bound + 1e-13
        lower = np.diag(np.sqrt(np.arange(1, d)), 1).astype(complex)
        q = (lower + lower.conj().T) / np.sqrt(2)
        p = -1j * (lower - lower.conj().T) / np.sqrt(2)
        ident = np.eye(d)
        q1 = (np.kron(q, ident) + np.kron(ident, q)) / np.sqrt(2)
        p2 = (np.kron(p, ident) - np.kron(ident, p)) / np.sqrt(2)
        boundary = np.diag([0.] * cutoff + [1.])
        predicted = 0.5j * (cutoff + 1) * (np.kron(ident, boundary) - np.kron(boundary, ident))
        residual = float(np.max(np.abs(q1 @ p2 - p2 @ q1 - predicted)))
        assert residual < 2e-13
        rows.append({"cutoff": cutoff, "compressed_weyl_commutator_on_vacuum": vacuum_norm,
                     "compressed_weyl_commutator_operator_norm": op_norm,
                     "vacuum_bound": float(bound),
                     "field_commutator_identity_residual": residual,
                     "field_commutator_operator_norm": (cutoff + 1) / 2})
    assert rows[0]["compressed_weyl_commutator_on_vacuum"] > 1e-3
    assert rows[-1]["compressed_weyl_commutator_on_vacuum"] < 1e-7
    return {"uncut_weyl_symplectic_pairing": 0.0, "rows": rows}


def stiffness(k, a, eta=0, zeta=0):
    kh = 2 * np.sin(a * np.array(k) / 2) / a
    x, y, z = kh
    curl = 1j * np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    weight = np.eye(3) + eta * a * a * np.diag(kh ** 2)
    isotropic_extra = zeta * a ** 4 * sum(kh[i] ** 2 * kh[j] ** 2
                                         for i in range(3) for j in range(i + 1, 3))
    weight += isotropic_extra * np.eye(3)
    return curl.conj().T @ weight @ curl


def check_spectra():
    a = 2 * np.pi / 9
    u = (2 * np.sin(a / 2) / a) ** 2
    rows = []
    for k in [(1, 0, 0), (1, 1, 0)]:
        row = {"wavevector": k}
        for label, eta, zeta in [("plaquette", 0, 0), ("eta_1", 1, 0), ("zeta_1", 0, 1)]:
            values = np.linalg.eigvalsh(stiffness(k, a, eta, zeta))
            assert abs(values[0]) < 1e-12 and values[1] > 0
            if k == (1, 0, 0):
                expected = np.array([u, u])
            elif eta:
                expected = np.array([2 * u, 2 * u * (1 + a * a * u)])
            else:
                expected = np.repeat(2 * u * (1 + zeta * a ** 4 * u ** 2), 2)
            assert np.max(np.abs(values[1:] - expected)) < 2e-12
            row[label] = values[1:].tolist()
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = {
        "scope": "same-author finite verification and negative controls; no continuum simulation",
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "displacement": check_displacement(),
        "axial_complete_histories": check_axial(),
        "basis_mismatch_control": check_basis_control(),
        "occupation_convergence": check_occupation_convergence(),
        "compression_and_commutation": check_commutators(),
        "nonselection_spectra": check_spectra(),
        "all_assertions_passed": True,
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "all_assertions_passed": True,
        "axial_max_probability_difference": max(x["max_history_difference"] for x in result["axial_complete_histories"]),
        "basis_control_difference": result["basis_mismatch_control"]["difference"],
        "K24_history_tv_error": result["occupation_convergence"]["rows"][-1]["tv_to_uncut"],
        "output": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
