"""Run the preregistered C1 model without modifying archived inputs."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import traceback

import numpy as np
from numpy.polynomial import Polynomial
import scipy

import branch
from branch import ARCHIVE
import spectrum as spec
import tails

HERE = Path(__file__).resolve().parent
PROTO = HERE / "PROTOCOL.json"
POLE = ARCHIVE / "grut_candidate_1_fluctuations_20261019/pole_results.json"
ARCHIVED_TREE = ARCHIVE / "grut_closure_fork_20261019/results.json"


def json_value(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(type(value).__name__)


def save(output, result):
    (output / "RESULTS.json").write_text(json.dumps(result, indent=2, default=json_value, allow_nan=False)+"\n")


def read_json(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def integral(arr, weights):
    return math.fsum((arr*weights).ravel())


def extrapolate(coarse, fine):
    return fine+(fine-coarse)/3


def graph_record(graph):
    return {k: v for k, v in graph.items() if k not in {"k", "vhat", "kernel"}}


def qualify(protocol, result, output):
    settings = protocol["profiles"]
    limits = settings["stationarity_checks"]
    reference = branch.reference_profile()
    profiles = {}
    result["branches"] = {}
    for charge in protocol["charge_order"]:
        name = str(charge)
        coarse, fine = [branch.stationary(reference, *settings[which], settings, charge=charge,
                                          label=f"{name}_{which}") for which in ["coarse", "fine"]]
        coarse.label, fine.label = name+"_coarse", name
        t = branch.tangent(fine, settings)
        controls = [branch.stability(fine, h, settings) for h in [.025, .0125]]
        tree_change = abs(coarse.statistics["energy"]-fine.statistics["energy"])
        good = all(
            p.statistics["success"] and 0 < p.omega < 1 and
            abs(p.statistics["charge"]-charge) < limits["absolute_charge_error_max"] and
            abs(p.statistics["relative_virial"]) < limits["relative_virial_error_max"] and
            p.statistics["min_F"] > -1e-7 and p.statistics["min_S"] > 0 and
            p.statistics["max_active_F_prime"] < 1e-9 and p.statistics["min_active_S_prime"] > -1e-9
            for p in [coarse, fine])
        passes = (good and tree_change < limits["tree_refinement_max"] and t["success"] and
                  t["dQ_dOmega"] < -limits["slope_minimum_negative_magnitude"] and
                  all(r["passes"] for r in controls))
        result["branches"][name] = {"coarse": coarse.statistics, "fine": fine.statistics, "tangent": t,
                                    "tree_refinement": tree_change, "stability": controls, "passes": passes}
        save(output, result)
        print(f"Branch Q={charge}: Omega={fine.omega:.12g}, slope={t['dQ_dOmega']:.7g}, pass={passes}", flush=True)
        if not passes:
            raise ValueError(f"Stationary branch or reduced stability failed at {charge}")
        profiles[name], profiles[coarse.label] = fine, coarse
    parent, neighbor = [profiles[str(q)] for q in protocol["charge_order"][:2]]
    midpoint = (parent.omega+neighbor.omega)/2
    h = protocol["chemical_potential"]["h"]
    dnames = [f"stencil_{i}" for i in range(5)]
    stencil = []
    for name, offset in zip(dnames, np.linspace(-1., 1., 5)):
        prof = branch.stationary(parent, *settings["fine"], settings, omega=midpoint+h*offset, label=name)
        if not prof.statistics["success"] or prof.statistics["min_F"] < -1e-7:
            raise ValueError(f"Auxiliary fixed-frequency branch failed at {name}")
        stencil.append(prof)
        profiles[name] = prof
    result["stencil_profiles"] = [p.statistics for p in stencil]
    return reference, profiles, stencil


def compare_archive(worker, reference, result):
    z, w = spec.quadrature([0., .25, 1., 3., 10.], 8)
    _, hf, hg = worker.settings["fd_h"]
    grids = [worker.high(reference, h, z, "archive_reference", cutoff=40, support=25.)[0] for h in [hf, hg]]
    current = (extrapolate(*grids)*spec.high_weights(w, 40)).sum(axis=1)
    replay = read_json(spec.ENERGY_ARCHIVE / "audit_results.json")["replay_partial_integrals"]
    target = np.array([r["energy"] for r in replay if 2 <= r["ell"] <= 40])
    sample = worker.fd(reference, hf, [(2, .7)], "scalar_vacuum_check", support=0.)
    error = sum(abs(current-target))
    result["independent_reference_check"] = {
        "abs_weighted_difference_by_ell": abs(current-target), "abs_sum": error,
        "sum_signed": sum(current-target), "vacuum_remainder": sample[0],
        "passes": bool(error < .004 and abs(sample[0]) < 1e-12),
    }
    if not result["independent_reference_check"]["passes"]:
        raise ValueError("Independent lattice failed the archived continuum comparator")


def physical_spectrum(worker, profile, protocol):
    settings = worker.settings
    edges = settings["frequency_edges"]
    slow, fast = settings["orders"]
    z, w = spec.quadrature(edges, fast)
    z8, w8 = spec.quadrature(edges, slow)
    wl, wh = spec.low_weights(w), spec.high_weights(w, worker.cutoff)
    wh8 = spec.high_weights(w8, worker.cutoff)
    h0, h1, h2 = settings["fd_h"]
    cross = settings["cross_points"]
    hs, probe = zip(*(worker.high(profile, h, z, "radial_main", extra=cross) for h in [h0, h1, h2]))
    high, high_coarse = extrapolate(*hs[1:]), extrapolate(*hs[:2])
    high8 = extrapolate(*[worker.high(profile, h, z8, "frequency_control")[0] for h in [h1, h2]])
    high_box = worker.high(profile, h1, z, "larger_box_control", box=settings["fd_box_check"])[0]
    high_support = worker.high(profile, h1, z, "profile_support_control", support=settings["fd_support_check"])[0]
    low0 = worker.low(profile.label, z, settings["jost_base"], "radial_base")
    low = worker.low(profile.label, z, settings["jost_refinement"], "radial_refined")
    low8 = worker.low(profile.label, z8, settings["jost_base"], "frequency_control")
    low25 = worker.low(profile.label+"_coarse", z8, settings["jost_support_check"], "profile_support_control")
    jc = worker.jost(profile.label, cross, settings["jost_base"], "cross_method")
    cross_diffs = [r["total_remainder"]-f for r, f in zip(jc, extrapolate(*probe[1:]))]
    ir_z = settings["jost_ir_frequencies"]
    ir = worker.jost(profile.label, [(l, x) for l in [0, 1] for x in ir_z],
                     settings["jost_refinement"], "zero_log_control")
    ir_rows = [[r["logJost_real"]-2*np.log(float(r["zeta"])) for r in ir if r["ell"] == l] for l in [0, 1]]
    ir_ptp = np.max(np.ptp(ir_rows, axis=1))
    graphs = [spec.born_graph(profile, conf) for conf in settings["restored_graph_grids"]]
    energy = spec.low_energy(low, z, w, edges[-1])+integral(high, wh)
    energy8 = spec.low_energy(low8, z8, w8, edges[-1])+integral(high8, wh8)
    graphs_diff = abs(graphs[1]["total"]-graphs[0]["total"])
    convergence = {
        "radial_and_independent_reference": spec.absolute_variation(high, high_coarse, wh),
        "frequency_quadrature": abs(energy-energy8),
        "support_box_and_origin": (
            spec.absolute_variation(low, low0, wl) +
            spec.absolute_variation(low8, low25, spec.low_weights(w8)) +
            spec.absolute_variation(high_box, hs[1], wh) +
            spec.absolute_variation(high_support, hs[1], wh)),
        "graphs": graphs_diff,
    }
    caps = protocol["absolute_allowance"]
    bounds = tails.absolute(profile, worker.support, worker.cutoff, edges[-1], caps["tail_enclosure_caps"])
    acceptance = {
        "controls_each_within_allocation": all(convergence[k] < caps[k] for k in convergence),
        "analytical_tail_caps": bounds["caps_pass"] and bounds["total"] < caps["analytic_total_cap"],
        "finite_grid_cross_method": max(abs(np.array(cross_diffs))) < settings["cross_point_max_abs"],
        "regularized_infrared": ir_ptp < settings["jost_ir_composite_variation_max"],
        "graph_controls": all(g["max_parseval_error"] < settings["graph_parseval_max"] and
                              g["min_bubble_denominator"] > 0 for g in graphs),
    }
    allowance = bounds["total"]+sum(caps[k] for k in [*convergence, "classical"])
    acceptance["total_allowance"] = allowance < caps["total_target_strict"]
    record = {
        "E0": profile.statistics["energy"], "K1_finite": energy+graphs[1]["total"],
        "Born_subtracted_spectrum": energy, "restored_graphs": [graph_record(g) for g in graphs],
        "E0plus1_finite": profile.statistics["energy"]+energy+graphs[1]["total"],
        "infrared_regularized_logJ": ir_rows, "continuum_minus_lattice_cross_points": cross_diffs,
        "retained_spectrum_by_ell_ge_2": (high*wh).sum(axis=1),
        "convergence": convergence, "analytical_tails": bounds, "total_allowance": allowance,
        "checks": acceptance, "passes": all(acceptance.values()),
    }
    raw = {"low": low, "high": high, "low8": low8, "high8": high8, "low0": low0, "low25": low25,
           "high_coarse": high_coarse, "high_at_h1": hs[1], "high_box": high_box, "high_support": high_support,
           "wl": wl, "wh": wh, "wl8": spec.low_weights(w8), "wh8": wh8, "graphs": graphs}
    return record, raw


def emission(protocol, states, spectra, raw, pole):
    p, n = [str(q) for q in protocol["charge_order"][:2]]
    a, b = [raw[q] for q in [p, n]]
    wf = {"low": a["wl"], "high": a["wh"]}
    lohi = {"low": spec.difference(a["low"], b["low"], wf["low"]),
            "high": spec.difference(a["high"], b["high"], wf["high"])}
    spectra_diff = sum(lohi.values())
    graph_diff = spec.graph_difference(a["graphs"][1], b["graphs"][1])
    tree_diff = spectra[p]["E0"]-spectra[n]["E0"]

    def check_same_grid(key, ref, weights):
        return spec.absolute_variation(a[key]-b[key], a[ref]-b[ref], weights)

    quadrature_shift = abs(spectra_diff - spec.difference(a["low8"], b["low8"], a["wl8"]) -
                           spec.difference(a["high8"], b["high8"], a["wh8"]))
    paired_shifts = {
        "empirical_radial": check_same_grid("high", "high_coarse", wf["high"]),
        "empirical_quadrature": quadrature_shift,
        "empirical_support_box_origin": (
            check_same_grid("low", "low0", wf["low"]) + check_same_grid("low8", "low25", a["wl8"]) +
            check_same_grid("high_box", "high_at_h1", wf["high"]) +
            check_same_grid("high_support", "high_at_h1", wf["high"])),
        "empirical_graphs_and_tree": (
            abs(graph_diff-spec.graph_difference(a["graphs"][0], b["graphs"][0])) +
            abs(tree_diff-(states[p+"_coarse"].statistics["energy"]-states[n+"_coarse"].statistics["energy"]))),
    }
    limits = protocol["paired_allowance"]
    det = protocol["determinant"]
    bounds = [spectra[q]["analytical_tails"] for q in [p, n]]
    separate_bound = sum(b["total"] for b in bounds)
    tail_pair = {"total": separate_bound, "method": "sum_of_absolute_tail_bounds"}
    if all(b["caps_pass"] for b in bounds):
        try:
            paired = tails.paired(states[p], states[n], det["fd_support"], det["ell_max"],
                                  det["frequency_edges"][-1], protocol["absolute_allowance"]["tail_enclosure_caps"],
                                  limits["angular_difference_envelope_scale"])
            tail_pair = {**paired, "method": "paired_operator_difference"}
        except ValueError as error:
            tail_pair["fallback_reason"] = str(error)
    allowance = tail_pair["total"]+sum(limits[k] for k in paired_shifts)+limits["pole_numerical"]
    subtraction = math.fsum([tree_diff, spectra_diff, graph_diff, -pole])
    passes = (all(paired_shifts[k] < limits[k] for k in paired_shifts) and
              allowance < limits["total_target_strict"])
    return {"value": subtraction, "allowance": allowance, "signed_interval": [subtraction-allowance, subtraction+allowance],
            "tree_difference": tree_diff, "paired_loop_spectrum_difference": lohi,
            "paired_graph_difference": graph_diff, "pole": pole, "paired_tail": tail_pair,
            "sum_of_individual_tail_bounds": separate_bound, "empirical_refinements": paired_shifts,
            "passes": passes, "aspirational_precision_met": allowance < limits["aspirational_target"]}


def derivative_energy(worker, profile):
    cfg = worker.settings
    edges = cfg["frequency_edges"]
    z, w = spec.quadrature(edges, cfg["orders"][-1])
    hs = [worker.high(profile, h, z, "derivative_stencil")[0] for h in cfg["fd_h"][:2]]
    low = worker.low(profile.label, z, cfg["jost_base"], "derivative_stencil")
    rest = [spec.born_graph(profile, grid) for grid in cfg["restored_graph_grids"]]
    lw = spec.low_energy(low, z, w, edges[-1])
    wh = spec.high_weights(w, worker.cutoff)
    return {
        "K1": lw+integral(extrapolate(*hs), wh)+rest[1]["total"],
        "K1_radial_graph_control": lw+integral(hs[1], wh)+rest[0]["total"],
        "graph_record": [graph_record(g) for g in rest],
    }


def chemical(worker, protocol, states, stencil, direct, output, result, pole):
    p, n = [states[str(i)] for i in protocol["charge_order"][:2]]
    settings = protocol["profiles"]
    cfg = protocol["chemical_potential"]
    omega_middle = (p.omega+n.omega)/2
    half_width = (p.omega-n.omega)/2
    h = cfg["h"]
    x = np.linspace(-1., 1., 5)
    rows = []
    for profile in stencil:
        rows.append(derivative_energy(worker, profile))
        result["chemical_stencil_energies"] = rows
        save(output, result)
        print(f"Derivative stencil {profile.label} complete", flush=True)
    y = np.array([row["K1"] for row in rows])
    y_control = np.array([row["K1_radial_graph_control"] for row in rows])

    def differentiate(y, indices):
        xp = x[indices]
        return Polynomial.fit(xp, y[indices], len(xp)-1).convert().deriv()/h

    all_x, inner, outer = np.arange(5), np.array([1, 2, 3]), np.array([0, 2, 4])
    derivative = differentiate(y, all_x)
    d_inner, d_outer = [differentiate(y, use) for use in [inner, outer]]
    d_radial = differentiate(y_control, all_x)
    nodes = [-1/np.sqrt(3), 1/np.sqrt(3), 0.]
    rows = []
    for i, node in enumerate(nodes):
        omega = omega_middle+half_width*node
        profile = branch.stationary(p, *settings["fine"], settings, omega=omega, label=f"gauss_{i}")
        t = branch.tangent(profile, settings)
        if not (profile.statistics["success"] and t["success"] and
                t["dQ_dOmega"] < -settings["stationarity_checks"]["slope_minimum_negative_magnitude"]):
            raise ValueError("Chemical potential tangent is unstable or singular")
        k_prime = derivative(half_width*node/h)
        rows.append({"omega": omega, "charge": profile.statistics["charge"], **t,
                     "partial_K1_along_branch": k_prime, "mu": omega+k_prime/t["dQ_dOmega"]})
    gauss_x = half_width/h*np.array(nodes[:2])
    omega_width = 2*half_width
    k_part, k_inner, k_outer, k_radial = [omega_width*np.mean(d(gauss_x)) for d in
                                              [derivative, d_inner, d_outer, d_radial]]
    k_midpoint = omega_width*derivative(0.)
    t_part = omega_width/2*sum(r["omega"]*r["dQ_dOmega"] for r in rows[:2])
    t_mid = omega_width*rows[2]["omega"]*rows[2]["dQ_dOmega"]
    q_part = omega_width/2*sum(r["dQ_dOmega"] for r in rows[:2])
    slope = result["branches"][p.label]["tangent"]["dQ_dOmega"]
    mu_p = p.omega+derivative(half_width/h)/slope
    vacuum = math.fsum([t_part, k_part, -pole])
    error = max(abs(k_part-kk) for kk in [k_inner, k_outer, k_radial, k_midpoint])
    empirical = protocol["absolute_allowance"]
    point_allowance = sum(empirical[k] for k in
                          ["radial_and_independent_reference", "frequency_quadrature", "support_box_and_origin", "graphs"])
    divided_weights = omega_width*np.array([np.mean(differentiate(e, all_x)(gauss_x)) for e in np.eye(5)])
    stencil_allowance = point_allowance*sum(abs(divided_weights))
    total = (cfg["loop_stencil_error_allowance"] + cfg["tree_integral_error_allowance"] +
             direct["paired_tail"]["total"] + protocol["paired_allowance"]["pole_numerical"])
    checks = {
        "stencil_refinements_plus_amplification": error+stencil_allowance < cfg["loop_stencil_error_allowance"],
        "tree_gauss_and_midpoint": max(abs(t_part-t_mid), abs(t_part-direct["tree_difference"]),
                                       abs(q_part-1)) < cfg["tree_integral_error_allowance"],
        "direct_agreement": abs(vacuum-direct["value"]) <= total+direct["allowance"],
    }
    return {
        "value": vacuum, "allowance": total, "signed_interval": [vacuum-total, vacuum+total],
        "gauss_and_midpoint_tangents": rows, "tree_integral": t_part, "tree_midpoint": t_mid, "charge_integral": q_part,
        "loop_integral": k_part, "loop_integral_3point_inner": k_inner, "loop_integral_3point_outer": k_outer,
        "loop_midpoint": k_midpoint, "loop_integral_lattice_and_graph_control": k_radial,
        "stencil_refinement_difference_max": error, "propagated_empirical_point_allowance": stencil_allowance,
        "retained_mu_at_parent": mu_p, "retained_parent_curvature_contribution": mu_p-(t_part+k_part),
        "tree_parent_curvature_contribution": p.omega-t_part,
        "classical_dOmega_dQ_at_parent": 1/slope, "difference_from_direct": vacuum-direct["value"],
        "checks": checks, "passes": all(checks.values()),
    }


def classify(row, valid):
    if not valid:
        return "inconclusive"
    low, high = row["signed_interval"]
    if high < 0:
        return "energetically forbidden"
    if low > 0:
        return "energetically accessible"
    return "inconclusive"


def main(args):
    started = time.monotonic()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    protocol = read_json(PROTO)
    result = {"protocol_sha256": sha(PROTO), "versions": {"numpy": np.__version__, "scipy": scipy.__version__,
                                                          "python": sys.version.split()[0]},
              "source_hashes": {p.name: sha(p) for p in [Path(__file__), Path(branch.__file__), Path(spec.__file__),
                                                         Path(tails.__file__), HERE/"radial.cc", POLE, branch.REFERENCE]},
              "energy_unit": "g_X v", "verdict": "D", "completed": False}
    worker = None
    try:
        assert sha(branch.REFERENCE) == protocol["profile_sha256"]
        pole = read_json(POLE)["matching_result"]["strict_pole_mass_over_v"]/protocol["g_X"]
        reference, states, stencil = qualify(protocol, result, output)
        spec.PROFILES.update(states)
        worker = spec.Spectrum(protocol, output, started)
        compare_archive(worker, reference, result)
        spectra, raw = {}, {}
        result["exact_charge_spectra"] = spectra
        p, n, h = [str(q) for q in protocol["charge_order"]]
        for name in [p, n, h]:
            spectra[name], raw[name] = physical_spectrum(worker, states[name], protocol)
            save(output, result)
            print(f"Matched spectrum Q={name} complete, gates={spectra[name]['passes']}", flush=True)
        direct = emission(protocol, states, spectra, raw, pole)
        result["direct_emission"] = direct
        fission = spectra[p]["E0plus1_finite"]-2*spectra[h]["E0plus1_finite"]
        fission_error = spectra[p]["total_allowance"]+2*spectra[h]["total_allowance"]
        result["equal_fission"] = {"value": fission, "allowance": fission_error,
                                  "signed_interval": [fission-fission_error, fission+fission_error],
                                  "loop_difference": spectra[p]["K1_finite"]-2*spectra[h]["K1_finite"]}
        result["archived_classical_comparisons"] = read_json(ARCHIVED_TREE)["classical_channels"]
        save(output, result)
        independent = chemical(worker, protocol, states, stencil, direct, output, result, pole)
        result["chemical_potential_crosscheck"] = independent
        valid_e = direct["passes"] and independent["passes"] and all(spectra[q]["passes"] for q in [p, n])
        valid_f = all(spectra[q]["passes"] for q in [p, h])
        statuses = [classify(direct, valid_e), classify(result["equal_fission"], valid_f)]
        result["channel_status"] = dict(zip(["single_particle_emission", "equal_fission"], statuses))
        forbidden = statuses.count("energetically forbidden")
        result["verdict"] = ["C" if "energetically accessible" in statuses else "D", "B", "A"][forbidden]
        result["completed"] = True
    except Exception as error:
        result["failed_execution"] = {"type": type(error).__name__, "message": str(error)}
        traceback.print_exc()
    finally:
        result["elapsed_seconds"] = time.monotonic()-started
        if worker:
            result["workload"] = {"jost_points_scheduled": worker.jost_points, "fd_batches": worker.fd_batches}
        save(output, result)
    print(json.dumps({k: result[k] for k in ["verdict", "completed", "elapsed_seconds"]}), flush=True)
    return 0 if result["completed"] else 1


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", required=True, type=Path, help="New directory for results and numerical logs")
    raise SystemExit(main(cli.parse_args()))
