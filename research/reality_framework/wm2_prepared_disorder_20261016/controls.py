"""Bounded S4 experiments. Every finite receiver is explicitly identified."""
from collections import deque
from itertools import product
from math import ceil, exp, factorial, sqrt
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
from scipy.linalg import eigh
from scipy.special import gammainc

import model as M


def dyson_tail(r, x):
    if r <= 0:
        return 1.0
    if x == 0:
        return 0.0
    return min(1.0, float(exp(x) * gammainc(r, x)))


def macro_population(a, x):
    threshold = max(2, ceil(a["n"] / 2)) - 1
    return float(sum(abs(x[j]) ** 2 for j, (_, o) in enumerate(a["states"])
                     if len(o) >= threshold))


def node_reproduction(base):
    archived = base.parent / "wm2_detuning_receiver_20261016"
    before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in archived.iterdir() if p.is_file()}
    node = shutil.which("node")
    assert node is not None, "Node is needed for the archived reproduction"
    with tempfile.TemporaryDirectory(prefix="grut_s3_node_") as directory:
        temp = Path(directory)
        for p in archived.iterdir():
            if p.is_file():
                shutil.copy2(p, temp / p.name)
        subprocess.run([node, str(temp / "reproduce.cjs")], check=True,
                       capture_output=True, text=True)
        completed = subprocess.run([node, str(temp / "reproduce.cjs"), "--recompute"],
                                   check=True, capture_output=True, text=True)
        verification = json.loads(completed.stdout)
        original = json.loads((archived / "RESULTS.json").read_text())
        recomputed = json.loads((temp / "RESULTS.json").read_text())
        differences = []
        numeric_leaves = 0

        def compare(a, b, path=""):
            nonlocal numeric_leaves
            if path == "execution":
                return
            if isinstance(a, dict):
                assert a.keys() == b.keys()
                for k in a:
                    compare(a[k], b[k], path + "." + k if path else k)
            elif isinstance(a, list):
                assert len(a) == len(b)
                for k, (x, y) in enumerate(zip(a, b)):
                    compare(x, y, f"{path}[{k}]")
            elif isinstance(a, (int, float)) and not isinstance(a, bool):
                numeric_leaves += 1
                error = abs(a - b)
                assert error <= 2e-10 * max(1, abs(a)), (path, a, b)
                if error:
                    differences.append(dict(path=path, absolute_difference=error))
            else:
                assert a == b, (path, a, b)

        compare(original, recomputed)
    after = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
             for p in archived.iterdir() if p.is_file()}
    assert before == after
    return dict(node_version=subprocess.check_output([node, "--version"], text=True).strip(),
                saved_evidence_driver="PASS", complete_recomputation="PASS",
                archived_packet_unchanged=True, numeric_leaves_compared=numeric_leaves,
                largest_numeric_difference=max((x["absolute_difference"] for x in differences), default=0),
                changed_numeric_leaves=len(differences), recomputed_verification=verification,
                independent_external_validation=False)


def reduction_checks():
    rows = []
    for n, modes in [(3, 3), (4, 3), (5, 2), (6, 2)]:
        delta = [0.0] + [-0.25 if v % 2 else 0.25 for v in range(1, n)]
        full = M.vacuum(n, modes, delta, classes=False)
        reduced = M.vacuum(n, modes, delta, classes=True)
        target = reduced["index"]
        mapping, factors = [], []
        for w, o in full["states"]:
            a = tuple(reduced["values"].index(delta[v]) for v in w[1:])
            multiplicity = 1
            for k, count in enumerate(reduced["counts"]):
                multiplicity *= factorial(count) // factorial(count - a.count(k))
            mapping.append(target[a, o])
            factors.append(1 / sqrt(multiplicity))
        mapping, factors = np.asarray(mapping), np.asarray(factors)
        for t in [0.25, 1.0, 4.0]:
            x = M.evolve(full["H"], M.unit(len(full["states"]), full["seed"]), t)
            z = M.evolve(reduced["H"], M.unit(len(reduced["states"]), reduced["seed"]), t)
            error = float(np.linalg.norm(x - z[mapping] * factors))
            assert error < 1e-10
            rows.append(dict(N=n, modes=modes, t=t, full_dimension=len(x),
                             exact_type_word_dimension=len(z), isometry_wavefunction_error=error))

    n, modes, delta = 4, 2, [0, -0.25, 0.25, -0.25]
    a = M.vacuum(n, modes, delta)
    x = np.zeros(len(a["states"]), dtype=complex)
    for v in range(1, n):
        x[a["index"][(0, v), (0,)]] = 1 / sqrt(n - 1)
    y = a["H"] @ x
    radial = y.copy()
    for k in range(n):
        for o in set(o for _, o in a["states"] if len(o) == k):
            indices = [j for j, (w, oo) in enumerate(a["states"]) if oo == o]
            radial[indices] = np.mean(y[indices])
    variance = float(np.std(delta[1:]))
    leakage = float(np.linalg.norm(y - radial))
    assert abs(leakage - variance) < 1e-12
    seed = M.unit(len(a["states"]), a["seed"])
    cubic = a["H"] @ (a["H"] @ (a["H"] @ seed))
    difference = cubic[a["index"][(0, 1, 2), (0, 0)]] - cubic[a["index"][(0, 2, 1), (0, 0)]]
    predicted = sqrt(2) / n ** 2 * (delta[1] - delta[2])
    assert abs(difference - predicted) < 1e-12
    return dict(exact_type_word_checks=rows, radial_invariance_counterexample=dict(
        N=n, delta=delta, leakage_norm=leakage, detuning_standard_deviation=variance,
        ordered_word_H3_difference=float(difference.real), exact_predicted_difference=predicted))


def cross_runtime_matrices(base):
    script = """
const S=require(process.argv[1]);
const h=S.make(Number(process.argv[2]),2,{delta:.5,g:.9});
console.log(JSON.stringify({states:h.states,cols:h.cols}));
"""
    rows = []
    source = base.parent / "wm2_detuning_receiver_20261016/model.cjs"
    for n in [2, 3, 4]:
        raw = json.loads(subprocess.check_output(["node", "-e", script, str(source), str(n)], text=True))
        python = M.full_probe(n, 2, [.5] * n, .9)
        keys = [(tuple(s["w"]), s["q"], tuple(r for r, k in enumerate(s["o"]) for _ in range(k)), s["p"])
                for s in raw["states"]]
        permutation = [python["index"][s] for s in keys]
        columns = python["H"].tocsc()
        maximum = 0
        entries = 0
        for j, column in enumerate(raw["cols"]):
            expected = {}
            for entry in column:
                i, value = entry[:2]
                expected[permutation[i]] = expected.get(permutation[i], 0) + value
            k = permutation[j]
            actual = {int(i): float(v) for i, v in zip(columns.indices[columns.indptr[k]:columns.indptr[k+1]],
                                                     columns.data[columns.indptr[k]:columns.indptr[k+1]])}
            for i in expected.keys() | actual.keys():
                maximum = max(maximum, abs(expected.get(i, 0) - actual.get(i, 0)))
                entries += 1
        assert maximum < 1e-13
        rows.append(dict(N=n, modes=2, complete_probe_dimension=len(keys), entries_compared=entries,
                         maximum_matrix_entry_difference=maximum))
    return rows


def first_error_checks():
    rows = []
    for n, modes in [(2, 2), (3, 2), (4, 2)]:
        a = M.full_probe(n, modes, [0] + [0.25] * (n - 1))
        adjacency = [[] for _ in a["states"]]
        for j, i, cost, kind in a["support"]:
            adjacency[j].append((i, cost))
        qset = set(a["Q"])
        checked, minimum_margin = 0, 1000
        for start in a["P"]:
            w, _, _, p = a["states"][start]
            if p not in w:
                continue
            expected = max(1, len(w) - 1 - w.index(p))
            distance = {start: 0}
            queue = deque([start])
            best = 10000
            while queue:
                j = queue.popleft()
                d = distance[j]
                if d >= best:
                    continue
                for i, cost in adjacency[j]:
                    dd = d + cost
                    if i in qset:
                        best = min(best, dd)
                    elif dd < distance.get(i, 10000):
                        distance[i] = dd
                        (queue.append if cost else queue.appendleft)(i)
            assert best == expected
            minimum_margin = min(minimum_margin, best - expected)
            checked += 1
        rows.append(dict(N=n, modes=modes, full_dimension=len(a["states"]),
                         complete_consistent_sources_checked=checked,
                         minimum_measured_minus_claimed=minimum_margin,
                         receiver_moves_cost=0))
    return rows


def primary_formation():
    """All binary corners for N=2,3, averaged over the original roots.

    Equal root/non-root detunings produce only a conserved root phase.
    There are at most two non-root detuning classes. A controlled mode
    cut is used, with a norm error against the physical L_N receiver.
    """
    cases, lookup, ensemble_rows = [], {}, []
    for n, modes in [(2, 128), (3, 224)]:
        for bar, g in product([-0.5, 0.0, 0.5], [0.9, 1.1]):
            needed = sorted(set(tuple(sorted(x)) for x in product([bar - .25, bar + .25], repeat=n - 1)))
            for nonroot in needed:
                a = M.vacuum(n, modes, [bar] + list(nonroot), g, classes=True)
                x = M.unit(len(a["states"]), a["seed"])
                rows = []
                previous = 0
                for scale in [4, 6, 8]:
                    t = scale * n * n
                    x = M.evolve(a["H"], x, t - previous)
                    previous = t
                    p = macro_population(a, x)
                    norm_error = abs(M.norm2(x) - 1)
                    assert norm_error < 2e-9
                    error = M.receiver_cut_error(n, modes, t)
                    assert error < 1e-15
                    populations = [float(sum(abs(x[j]) ** 2 for j, (_, o) in enumerate(a["states"]) if len(o) == k))
                                   for k in range(n)]
                    hx = a["H"] @ x
                    energy = float(np.vdot(x, hx).real)
                    variance = M.norm2(hx) - energy ** 2
                    expected_variance = (g / n) ** 2 * (n - 1)
                    assert abs(energy) < 1e-10
                    assert abs(variance - expected_variance) < 1e-10
                    work = -float(np.vdot(x, a["H_emission"] @ x).real)
                    word_prob = {}
                    for j, (w, o) in enumerate(a["states"]):
                        word_prob[w] = word_prob.get(w, 0) + abs(x[j]) ** 2
                    entropy = 0
                    for w, probability in word_prob.items():
                        multiplicity = 1
                        for v, count in enumerate(a["counts"]):
                            multiplicity *= factorial(count) // factorial(count - w.count(v))
                        if probability > 0:
                            entropy -= probability * np.log2(probability / multiplicity)
                    rows.append(dict(t=t, macro_probability=p, unitary_norm_error=norm_error,
                                     assembly_size_probabilities=populations,
                                     full_assembly_probability=populations[-1],
                                     relative_total_energy=energy, conserved_energy_variance=variance,
                                     expected_energy_variance=expected_variance,
                                     mean_geometry_record_work=work,
                                     word_record_entropy_given_root_bits=entropy,
                                     physical_receiver_probability_error_upper=2 * error,
                                     receiver_state_norm_error_upper=error))
                entry = dict(N=n, physical_receiver_modes=64 * n * n,
                             retained_receiver_modes=modes, g=g, nonroot_detunings=nonroot,
                             exact_type_word_dimension=len(a["states"]), rows=rows)
                cases.append(entry)
                lookup[n, bar, g, nonroot] = entry
            for pattern in product([bar - .25, bar + .25], repeat=n):
                p = []
                record_entropies = []
                for k in range(3):
                    total = 0
                    entropy = np.log2(n)
                    for root in range(n):
                        nonroot = tuple(sorted(pattern[v] for v in range(n) if v != root))
                        total += lookup[n, bar, g, nonroot]["rows"][k]["macro_probability"] / n
                        entropy += lookup[n, bar, g, nonroot]["rows"][k]["word_record_entropy_given_root_bits"] / n
                    p.append(total)
                    record_entropies.append(float(entropy))
                ensemble_rows.append(dict(N=n, bar_delta=bar, g=g, pattern=pattern,
                                          root_averaged_probabilities=p,
                                          complete_geometry_record_entropy_bits=record_entropies))
    return dict(cases=cases, all_binary_corner_root_averages=ensemble_rows,
                all_continuous_patterns_proved=False,
                finite_N_not_a_scaling_or_nonshrinking_formation_theorem=True)


def receiver_cut_comparison():
    n, retained, physical = 3, 224, 576
    small = M.vacuum(n, retained, [0, -.25, .25], .9, classes=True)
    large = M.vacuum(n, physical, [0, -.25, .25], .9, classes=True)
    mapping = np.asarray([large["index"][s] for s in small["states"]])
    x = M.unit(len(small["states"]), small["seed"])
    y = M.unit(len(large["states"]), large["seed"])
    rows = []
    previous = 0
    for t in [36, 72]:
        x = M.evolve(small["H"], x, t - previous)
        y = M.evolve(large["H"], y, t - previous)
        previous = t
        embedded = np.zeros_like(y)
        embedded[mapping] = x
        difference = float(np.linalg.norm(y - embedded))
        assert difference < 2e-9
        rows.append(dict(t=t, directly_observed_state_difference=difference,
                         analytic_receiver_cut_error_excluding_roundoff=M.receiver_cut_error(n, retained, t),
                         observed_probability_difference=abs(macro_population(small, x) - macro_population(large, y))))
    return dict(N=n, retained_modes=retained, physical_modes=physical,
                retained_dimension=len(x), physical_dimension=len(y),
                note="observed floating-point discrepancy includes roundoff; the analytic cut error does not",
                rows=rows)


def small_mode_formation():
    rng = np.random.default_rng(20261016)
    rows = []
    for n in [3, 4, 5, 6]:
        patterns = dict(common=[0.0] * n,
                        balanced=[-.25 if i % 2 else .25 for i in range(n)],
                        outlier=[.25] + [-.25] * (n - 1),
                        ramp=np.linspace(-.25, .25, n).tolist(),
                        iid_uniform=rng.uniform(-.25, .25, n).tolist())
        for name, delta in patterns.items():
            means = np.zeros(3)
            entropies = np.zeros(3)
            for root in range(n):
                assignment = [delta[root]] + [delta[v] for v in range(n) if v != root]
                a = M.vacuum(n, 2, assignment, classes=True)
                x = M.unit(len(a["states"]), a["seed"])
                previous = 0
                for k, scale in enumerate([4, 6, 8]):
                    t = scale * n * n
                    x = M.evolve(a["H"], x, t - previous)
                    previous = t
                    assert abs(M.norm2(x) - 1) < 1e-9
                    means[k] += macro_population(a, x) / n
                    probabilities = {}
                    for j, (w, o) in enumerate(a["states"]):
                        probabilities[w] = probabilities.get(w, 0) + abs(x[j]) ** 2
                    for w, p in probabilities.items():
                        multiplicity = 1
                        for v, count in enumerate(a["counts"]):
                            multiplicity *= factorial(count) // factorial(count - w.count(v))
                        if p > 0:
                            entropies[k] -= p * np.log2(p / multiplicity) / n
            entropies += np.log2(n)
            rows.append(dict(N=n, receiver_modes=2, pattern_name=name, delta=delta,
                             g=1.0, seed=20261016 if name == "iid_uniform" else None,
                             root_averaged_macro_probabilities=means.tolist(),
                             complete_word_record_entropy_bits=entropies.tolist(),
                             physical_L_N_comparison_bound=None,
                             interpretation="different finite receiver control; not primary formation evidence"))
    return rows


def conditional_probe(n, modes, delta, word, placement):
    vacuum = M.vacuum(n, modes, delta)
    formed = M.evolve(vacuum["H"], M.unit(len(vacuum["states"]), vacuum["seed"]), 4 * n * n)
    conditional = {o: formed[j] for j, (w, o) in enumerate(vacuum["states"]) if w == word}
    record_probability = sum(abs(z) ** 2 for z in conditional.values())
    assert record_probability > 1e-12
    conditional = {o: z / sqrt(record_probability) for o, z in conditional.items()}
    full = M.full_probe(n, modes, delta)
    p, q = full["P"], full["Q"]
    hp = full["H"][p][:, p]
    b = full["H"][q][:, p]
    hq = full["H"][q][:, q]
    mask = sum(1 << v for v in word)
    positions = (0, 1) if placement == "interior" else (len(word) - 2, len(word) - 1)
    labels = [word[k] for k in positions]
    columns = np.zeros((len(full["states"]), 2), dtype=complex)
    for o, z in conditional.items():
        for k, label in enumerate(labels):
            columns[full["index"][word, mask, o, label], k] = z
    r = max(1, len(word) - 1 - max(positions))
    initial_drive = float(np.linalg.norm(b @ columns[p]))
    if r > 1:
        assert initial_drive < 1e-14
    eig = vectors = None
    if n <= 4:
        eig, vectors = eigh(hq.toarray(), check_finite=False)
    rows = []
    for t in [.005, .05, .25, 1.0]:
        x = M.evolve(full["H"], columns, t)
        xp = M.evolve(hp, columns[p], t)
        qprob = [M.norm2(x[q, k]) for k in range(2)]
        response, consistent_response, inconsistent_response, dec_response = [], [], [], []
        for label in range(n):
            selection = np.asarray([s[3] == label for s in full["states"]])
            response.append(float(np.vdot(x[selection, 1], x[selection, 0]).imag))
            consistent_response.append(float(np.vdot(x[p[selection[p]], 1], x[p[selection[p]], 0]).imag))
            inconsistent_response.append(float(np.vdot(x[q[selection[q]], 1], x[q[selection[q]], 0]).imag))
            dec_response.append(float(np.vdot(xp[selection[p], 1], xp[selection[p], 0]).imag))
        # For theta=+/-pi/2 the identity is exact for this one-particle pulse.
        plus = M.evolve(full["H"], (columns[:, 0] * (-1j) + columns[:, 1]) / sqrt(2), t)
        minus = M.evolve(full["H"], (columns[:, 0] * (1j) + columns[:, 1]) / sqrt(2), t)
        finite_pulse = [float((M.norm2(plus[np.asarray([s[3] == label for s in full["states"]])])
                              - M.norm2(minus[np.asarray([s[3] == label for s in full["states"]])])) / 2)
                        for label in range(n)]
        pulse_error = max(abs(a - z) for a, z in zip(response, finite_pulse))
        assert pulse_error < 1e-10
        tail = dyson_tail(r, M.NU * t)
        feedback_amplitude = min(2.0, M.EPSILON / M.NU * exp(M.NU * t) * gammainc(r, M.NU * t))
        differences = [float(np.linalg.norm(x[:, k] - np.bincount(p, weights=xp[:, k].real, minlength=len(x))
                                          - 1j * np.bincount(p, weights=xp[:, k].imag, minlength=len(x)))) for k in range(2)]
        assert max(differences) <= feedback_amplitude + 1e-10
        assert max(qprob) <= min(1, feedback_amplitude) ** 2 + 1e-10
        assert max(abs(a - z) for a, z in zip(response, dec_response)) <= 2 * feedback_amplitude + 1e-10
        assert max(abs(z) for z in inconsistent_response) <= sqrt(qprob[0] * qprob[1]) + 1e-12
        if t == .005:
            lower = .25 * t - np.expm1(4.2 * t) + 4.2 * t
            assert response[labels[1]] >= lower - 1e-10
        # The measure is a QHQ spectral measure of the actual P-evolved drive.
        # It is not a projection of the prepared state onto full-H eigenstates.
        drive = b @ xp
        spectral = None
        if eig is not None:
            z = vectors.T.conj() @ drive
            spectral = dict(Q_eigenvalue_range=[float(eig[0]), float(eig[-1])],
                            seed_window=[-2.5, 2.5],
                            seed_window_drive_weights=np.sum(abs(z[(eig >= -2.5) & (eig <= 2.5)]) ** 2, axis=0).tolist(),
                            full_complement_drive_weights=np.sum(abs(z) ** 2, axis=0).tolist(),
                            seed_window_eigenvalues=int(np.sum((eig >= -2.5) & (eig <= 2.5))))
        drive_bound = M.EPSILON ** 2 * dyson_tail(r - 1, M.NU * t) ** 2
        assert max(M.norm2(drive[:, k]) for k in range(2)) <= drive_bound + 1e-10
        rows.append(dict(t=t, chi=response, chi_consistent=consistent_response,
                         chi_inconsistent=inconsistent_response, chi_P_only_comparison=dec_response,
                         pulse_identity_error=pulse_error, Q_probabilities=qprob,
                         complete_feedback_state_norm_differences=differences,
                         feedback_state_norm_upper=feedback_amplitude,
                         first_error_Dyson_amplitude_upper=tail,
                         spectral_driving_measure=spectral, all_window_drive_mass_upper=drive_bound))
    return dict(N=n, modes=modes, delta=delta, word=word, source_placement=placement,
                source_labels=labels, number_variance=.25, minimum_first_error_operations=r,
                initial_complement_driving_norm=initial_drive,
                full_dimension=len(full["states"]), conditional_word_probability_given_root=record_probability,
                receiver_conditioned_from_time=4 * n * n, receiver_reset=False,
                note="small-mode physical control; no assertion that its field equals the L_N field",
                rows=rows)


def pole_packet_checks():
    rows = []
    for n in [128, 256]:
        delta = np.linspace(-.75, .75, n)
        for sign in [-1, 1]:
            assignment = sign * delta
            modes, quanta = 64 * n * n, n - 1
            box = modes // 2
            target = -(80 + sum(assignment[1:])) / quanta
            assert abs(target) < 2
            r = np.arange(1, box + 1)
            phi = sqrt(2 / (box + 1)) * np.sin(np.pi * r / (box + 1)) * np.exp(1j * np.arccos(-target / 2) * r)
            hphi = np.zeros(box + 1, dtype=complex)
            hphi[0] = -phi[0]
            hphi[1:-1] = -np.r_[0, phi[:-2]] - phi[1:]
            hphi[-1] = -phi[-2]
            expected = np.r_[0, target * phi]
            residual1 = M.norm2(hphi - expected)
            bias = np.vdot(phi, hphi[1:]).real - target
            residual = sqrt(quanta * residual1 + quanta * (quanta - 1) * bias ** 2)
            bound = 2 * np.pi * quanta / (box + 1)
            assert residual <= bound + 1e-9
            rows.append(dict(N=n, receiver_modes=modes, independent_detuning_pattern="signed_ramp",
                             sign=sign, target_energy=float(target), direct_many_boson_Q_residual=residual,
                             residual_upper=bound, coupled_pole_window_radius=2 * sqrt(n) * bound,
                             inverse_norm_lower_if_invertible=1 / bound,
                             preparation_overlap_claimed=False))
    return rows
