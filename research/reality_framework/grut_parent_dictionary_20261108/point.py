"""P1-01: unchanged M1 physics at one further registered parent cutoff."""
import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import eigh

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "grut_detector_matching_20261107"
sys.path.insert(0, str(OLD))
import bosonized as b


def write(path, data):
    with path.open("x") as f:
        json.dump(data, f, indent=2, allow_nan=False)
        f.write("\n")


def inclusion(small, large):
    lookup = {tuple(s): i for i, s in enumerate(large)}
    padding = (0,) * (large.shape[1] - small.shape[1])
    return np.array([lookup[tuple(s) + padding] for s in small])


def field_mean(states, ground, omega, card, meta):
    _, _, fn, _, _ = b.interval_spectrum(card, card["reference_modes"], card["quadrature_epsabs"])
    hs = fn[:len(omega)] / np.sqrt(2*np.pi*omega)
    lookup = {tuple(s): i for i, s in enumerate(states)}
    mean = meta["mean"]
    for i, state in enumerate(states):
        for n, h in enumerate(hs):
            raised = state.copy()
            raised[n] += 1
            j = lookup.get(tuple(raised))
            if j is not None:
                mean += 2*h*np.sqrt(state[n]+1)*ground[i]*ground[j]
    return float(mean)


def run(point):
    start = time.monotonic()
    protocol = json.loads((HERE/"protocol_p1_01.json").read_text())
    card = json.loads((OLD/"protocol_m1_01.json").read_text())
    if not (card["g"] == 1 and card["theta"] == 0 and card["length"] == 5
            and card["boundary_phi_over_sqrt_pi"] == [.25, .25]
            and card["mass_R_over_g"] == .5):
        raise ValueError("This builder only covers the frozen specialized dictionary")
    card["resources"]["dimension_cap"] = protocol["dimension_cap"]
    cutoff, quadrature = point["cutoff"], point["quadrature"]
    states, omega, h = b.build_hamiltonian(card, cutoff, quadrature, card["mass_R_over_g"])
    energies, vectors = eigh(h, driver="evd")
    residual = max(float(np.linalg.norm(h@vectors[:,j]-energies[j]*vectors[:,j])) for j in range(2))
    w, _, _, meta = b.compressed_weyl(card, states, omega)
    defect = np.eye(len(w))-w.conj().T@w
    mineig = float(eigh(defect, eigvals_only=True, subset_by_index=[0,0])[0])
    if mineig < -protocol["gates"]["negative_defect"]:
        raise ValueError(f"Noncontractive Weyl compression: {mineig}")
    del defect
    histories, losses = b.histories(energies, vectors, w, card["protocols"], card["phase"])
    data = {
        "protocol": protocol["id"], "label": point["label"], "cutoff": cutoff,
        "quadrature": quadrature, "dimension": len(states), "occupied_modes": len(omega),
        "energies_low": energies[:6].tolist(), "gap": float(energies[1]-energies[0]),
        "hermiticity_residual": float(np.linalg.norm(h-h.T, ord=np.inf)),
        "eigen_residual": residual, "minimum_defect_eigenvalue": mineig,
        "histories": histories, "overflow_by_step": losses,
        "probability_norm_error": max(abs(sum(p.values())-1) for p in histories.values()),
        "ground_field_mean": field_mean(states, vectors[:,0], omega, card, meta),
        "detector_reference_metadata": meta,
        "reuse": "M1 Hamiltonian, Weyl matrix elements and history algorithm reused without physical changes"
    }
    if cutoff == 14 and quadrature == 128:
        small, small_omega, small_h = b.build_hamiltonian(card, 12, 128, .5)
        small_w, _, _, _ = b.compressed_weyl(card, small, small_omega)
        index = inclusion(small, states)
        projected_w = w[np.ix_(index, index)]
        columns = w[:, index].copy()
        df_block = (np.eye(len(index))-columns.conj().T@columns)/2
        columns[index, :] = 0
        transfer = columns.conj().T@columns/2
        de = (np.eye(len(index))-small_w.conj().T@small_w)/2
        old = np.load(OLD/"raw_m1_01/boson_m05_E12_q128.npz", allow_pickle=False)
        if not np.array_equal(old["basis"], small):
            raise ValueError("Archived occupation dictionary differs")
        embedded = vectors[index, 0]
        data["nested_E12_checks"] = {
            "hamiltonian_max_entry": float(np.max(np.abs(h[np.ix_(index,index)]-small_h))),
            "weyl_max_entry": float(np.max(np.abs(projected_w-small_w))),
            "defect_decomposition_max_entry": float(np.max(np.abs(de-df_block-transfer))),
            "ground_weight_outside_E12": float(1-np.vdot(embedded,embedded).real),
            "embedded_ground_overlap_squared": float(abs(np.vdot(old["ground"],embedded))**2),
            "first_detector_loss_E12_state_into_E12_to_E14": float(np.vdot(old["ground"],transfer@old["ground"]).real),
            "first_detector_loss_E12_state_still_outside_E14": float(np.vdot(old["ground"],df_block@old["ground"]).real),
            "first_detector_loss_E12_state_total": float(np.vdot(old["ground"],de@old["ground"]).real)
        }
    data["seconds"] = time.monotonic()-start
    return data, states, energies, vectors[:,0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    card = json.loads((HERE/"protocol_p1_01.json").read_text())
    point = next(p for p in card["points"] if p["label"] == args.label)
    if args.output.exists() or args.output.with_suffix(".npz").exists():
        raise FileExistsError("Refusing to replace a preserved result")
    data, states, energies, ground = run(point)
    np.savez_compressed(args.output.with_suffix(".npz"), basis=states, energies=energies, ground=ground)
    write(args.output, data)
    print(json.dumps({k:data[k] for k in ("label","dimension","gap","seconds")}), flush=True)
