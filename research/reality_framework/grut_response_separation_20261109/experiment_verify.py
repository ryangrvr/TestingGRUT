"""EX1-01 ideal verifier and prospective raw-count scorer; Python standard library.

No hardware API, simulated shots, fitted parameters, or empirical pass is produced.
Run: python experiment_verify.py --output experiment_ideal.json
Optional future data: --counts counts.json (format documented in the protocol).
"""
import argparse
import cmath
import hashlib
import json
import math
from pathlib import Path

NQ = 4
DIM = 1 << NQ
S, M1, M2, A = range(NQ)
ALPHA_VALIDATION = 0.005
ALPHA_PER_OBJECT = ALPHA_VALIDATION / 4
COUNTS_REQUIRED = {"C11": 20000, "P11": 20000, "I_C": 50000, "I_P": 50000}
IDEAL_ACTION = {"00": 1.0, "10": 0.75, "01": 0.75, "11": 13 / 16}


def bit(i, q):
    return (i >> (NQ - 1 - q)) & 1


def single(v, q, m):
    out = list(v)
    mask = 1 << (NQ - 1 - q)
    for i in range(DIM):
        if not i & mask:
            j = i | mask
            out[i] = m[0][0] * v[i] + m[0][1] * v[j]
            out[j] = m[1][0] * v[i] + m[1][1] * v[j]
    return out


def cnot(v, control, target):
    out = [0j] * DIM
    mask = 1 << (NQ - 1 - target)
    for i, amplitude in enumerate(v):
        out[i ^ mask if bit(i, control) else i] = amplitude
    return out


def rx(v, q, angle):
    c, s = math.cos(angle / 2), -1j * math.sin(angle / 2)
    return single(v, q, ((c, s), (s, c)))


def rz(v, q, angle):
    return [x * cmath.exp(-0.5j * angle * (-1) ** bit(i, q)) for i, x in enumerate(v)]


def rzz(v, angle):
    return [x * cmath.exp(-0.5j * angle * (-1) ** (bit(i, M1) ^ bit(i, M2)))
            for i, x in enumerate(v)]


def parity(v):
    return [x * (-1) ** (bit(i, M1) ^ bit(i, M2)) for i, x in enumerate(v)]


def prep():
    v = [0j] * DIM
    v[0] = 1
    h = tuple(tuple(x / math.sqrt(2) for x in row) for row in ((1, 1), (1, -1)))
    return cnot(single(single(v, S, h), M1, h), M1, M2)


def before_dephaser(program, no_writer=False):
    v = prep()
    if not no_writer:
        v = cnot(v, S, M2)
    v = rz(v, M1, math.pi / 2)
    return rx(v, M2, math.pi / 3) if program[0] == "1" else v


def run(program, dephase=False, beta=2 * math.pi / 3, no_writer=False,
        fixed_action=False, immediate=False, wrong_dephaser=False):
    v = before_dephaser(program, no_writer)
    states = [(1.0, v)]
    if dephase:
        pv = ([x * (-1) ** bit(i, M1) for i, x in enumerate(v)]
              if wrong_dephaser else parity(v))
        states = [(0.5, v), (0.5, pv)]
    if immediate:
        return states
    out = []
    for weight, state in states:
        state = rzz(state, beta)
        if program[1] == "1":
            state = rx(state, M2, math.pi / 3)
        if not fixed_action:
            state = cnot(cnot(state, M1, A), M2, A)
        out.append((weight, state))
    return out


def raw_distribution(states):
    return [sum(weight * abs(state[i]) ** 2 for weight, state in states) for i in range(DIM)]


def registered_distribution(raw):
    out = [0.0] * 9
    for i, p in enumerate(raw[:DIM]):
        j = 4 * bit(i, S) + 2 * (bit(i, M1) ^ bit(i, M2)) + bit(i, A)
        out[j] += p
    if len(raw) == DIM + 1:
        out[8] = raw[DIM]
    return out


def agreement(raw):
    return sum(p for i, p in enumerate(raw[:DIM]) if bit(i, S) == bit(i, A))


def tv(p, q):
    return sum(abs(a - b) for a, b in zip(p, q)) / 2


def density(states):
    return [[sum(w * v[i] * v[j].conjugate() for w, v in states)
             for j in range(DIM)] for i in range(DIM)]


def max_matrix_difference(a, b):
    return max(abs(a[i][j] - b[i][j]) for i in range(DIM) for j in range(DIM))


def scalar_radius(n):
    return math.sqrt(math.log(2 / ALPHA_PER_OBJECT) / (2 * n))


def tv_radius(n, k=9):
    # TV = max_A |p(A)-q(A)|; union-bound Hoeffding over nontrivial subsets.
    return math.sqrt(math.log(2 * (2 ** k - 2) / ALPHA_PER_OBJECT) / (2 * n))


def stats_profile():
    return {
        "validation_family_failure_probability": ALPHA_VALIDATION,
        "calibration_family_failure_probability_reserved": 0.005,
        "alpha_each_of_four_confidence_objects": ALPHA_PER_OBJECT,
        "validation_shots": COUNTS_REQUIRED,
        "validation_total_shots": sum(COUNTS_REQUIRED.values()),
        "calibration_total_cap": 1000000,
        "primary_scalar_radius": scalar_radius(20000),
        "contrast_radius": 2 * scalar_radius(20000),
        "each_immediate_tv_radius": tv_radius(50000),
        "immediate_comparison_tv_radius": 2 * tv_radius(50000),
        "ideal_primary_equivalence_tolerance": 0.06,
        "contrast_lower_bound_required": 0.10,
        "immediate_tv_upper_bound_required": 0.05,
        "scope": "independent bounded outcomes; time-averaged probabilities if nonidentical",
    }


def ideal_verify():
    checks, programs = {}, {}
    for program, target in IDEAL_ACTION.items():
        raw = raw_distribution(run(program))
        result = agreement(raw)
        programs[program] = {"action_agreement": result, "target": target,
                             "raw_bits_S_M1_M2_A": raw,
                             "joint_S_parity_A_plus_erasure": registered_distribution(raw)}
        checks["program_" + program] = abs(result - target) < 1e-12
        checks["normalization_" + program] = abs(sum(raw) - 1) < 1e-12
        theta1 = math.pi / 3 if program[0] == "1" else 0
        theta2 = math.pi / 3 if program[1] == "1" else 0
        analytic = (1 + math.cos(theta1) * math.cos(theta2)
                    - math.sin(theta1) * math.sin(theta2) * math.cos(2 * math.pi / 3)) / 2
        checks["independent_logical_rotation_" + program] = abs(result - analytic) < 1e-12
    pinched_raw = raw_distribution(run("11", dephase=True))
    pinched = agreement(pinched_raw)
    checks["pinched_prediction"] = abs(pinched - 5 / 8) < 1e-12
    immediate_c = registered_distribution(raw_distribution(run("11", immediate=True)))
    immediate_p = registered_distribution(raw_distribution(run("11", immediate=True, dephase=True)))
    checks["immediate_distribution_unchanged"] = tv(immediate_c, immediate_p) < 1e-12
    v = before_dephaser("11")
    projected = [[x if (bit(i, M1) ^ bit(i, M2)) == b else 0j
                  for i, x in enumerate(v)] for b in (0, 1)]
    luders = density([(1, projected[0]), (1, projected[1])])
    randomized = density(run("11", immediate=True, dephase=True))
    checks["random_parity_channel_equals_luders"] = max_matrix_difference(luders, randomized) < 1e-12
    blank = prep()
    blank_rho = density([(1, blank)])
    parity_blank = density([(0.5, blank), (0.5, parity(blank))])
    local_z_blank = [x * (-1) ** bit(i, M1) for i, x in enumerate(blank)]
    locally_pinched_blank = density([(0.5, blank), (0.5, local_z_blank)])
    checks["parity_channel_preserves_within_sector_Bell_coherence"] = max_matrix_difference(blank_rho, parity_blank) < 1e-12
    checks["individual_Z_measurement_is_stronger_than_parity_channel"] = max_matrix_difference(blank_rho, locally_pinched_blank) > 0.1
    phase_error = 0.0
    for i in range(DIM):
        basis = [0j] * DIM
        basis[i] = 1
        direct = rzz(basis, 2 * math.pi / 3)
        decomposed = cnot(rz(cnot(basis, M1, M2), M2, 2 * math.pi / 3), M1, M2)
        phase_error = max(phase_error, max(abs(a - b) for a, b in zip(direct, decomposed)))
    checks["RZZ_decomposition"] = phase_error < 1e-12
    checks["writer_removed_control"] = abs(agreement(raw_distribution(run("00", no_writer=True))) - 0.5) < 1e-12
    checks["actuator_removed_control"] = abs(agreement(raw_distribution(run("00", fixed_action=True))) - 0.5) < 1e-12
    wrong_phase = agreement(raw_distribution(run("11", beta=0)))
    wrong_pinch = agreement(raw_distribution(run("11", dephase=True, wrong_dephaser=True)))
    checks["negative_control_missing_phase_rejected"] = abs(wrong_phase - 13 / 16) > 0.06
    checks["negative_control_wrong_dephaser_rejected"] = abs(wrong_pinch - 5 / 8) > 0.06
    return {"status": "ideal_reference_checks_passed" if all(checks.values()) else "failed",
            "empirical_status": "not_run", "checks": checks, "programs": programs,
            "pinched_action_agreement": pinched,
            "coherent_minus_pinched": programs["11"]["action_agreement"] - pinched,
            "immediate_registered_tv": tv(immediate_c, immediate_p),
            "dephasing_channel_identity_residual": max_matrix_difference(luders, randomized),
            "RZZ_decomposition_residual": phase_error,
            "negative_controls": {"missing_beta_phase_action": wrong_phase,
                                  "wrong_local_M1_dephaser_action": wrong_pinch},
            "prospective_statistics": stats_profile()}


def score_counts(path):
    """Only scores statistical gates. This deliberately cannot certify the apparatus."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if set(data) != {"conditions"} or set(data["conditions"]) != set(COUNTS_REQUIRED):
        raise ValueError("Expected only 'conditions', containing exactly C11, P11, I_C, I_P")
    keys = [format(i, "04b") for i in range(DIM)] + ["erasure"]
    records = {}
    for condition, n in COUNTS_REQUIRED.items():
        counts = data["conditions"][condition]
        if set(counts) != set(keys):
            raise ValueError(f"{condition}: require all 16 S,M1,M2,A bitstrings and erasure")
        if any(type(x) is not int or x < 0 for x in counts.values()) or sum(counts.values()) != n:
            raise ValueError(f"{condition}: nonnegative integer counts must total {n}")
        raw = [counts[key] / n for key in keys]
        records[condition] = {"shots": n, "action_agreement": agreement(raw),
                              "joint_S_parity_A_plus_erasure": registered_distribution(raw),
                              "erasure_rate": raw[-1]}
    p_c, p_p = (records[x]["action_agreement"] for x in ("C11", "P11"))
    r = scalar_radius(20000)
    immediate_tv = tv(records["I_C"]["joint_S_parity_A_plus_erasure"],
                      records["I_P"]["joint_S_parity_A_plus_erasure"])
    metrics = {"coherent_ideal_error_upper": abs(p_c - 13 / 16) + r,
               "pinched_ideal_error_upper": abs(p_p - 5 / 8) + r,
               "contrast_lower": p_c - p_p - 2 * r,
               "immediate_tv_upper": min(1.0, immediate_tv + 2 * tv_radius(50000))}
    gates = {"coherent_equivalence": metrics["coherent_ideal_error_upper"] <= 0.06,
             "pinched_equivalence": metrics["pinched_ideal_error_upper"] <= 0.06,
             "coherence_contrast": metrics["contrast_lower"] >= 0.10,
             "immediate_registered_equivalence": metrics["immediate_tv_upper"] <= 0.05}
    return {"status": "statistical_gates_pass" if all(gates.values()) else "statistical_gates_fail",
            "overall_empirical_claim": "not_evaluated: apparatus_calibration_and_model_qualification_require_separate_review",
            "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "records": records, "metrics": metrics, "gates": gates,
            "sampling_assumptions_verified_by_script": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--counts", type=Path)
    args = parser.parse_args()
    result = ideal_verify()
    result["verifier_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.counts:
        result["future_real_count_scoring"] = score_counts(args.counts)
    with args.output.open("x", encoding="utf-8") as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write("\n")
    print(json.dumps({"output": str(args.output), "status": result["status"],
                      "empirical_status": result["empirical_status"],
                      "checks": len(result["checks"])}))
    return 0 if all(result["checks"].values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
