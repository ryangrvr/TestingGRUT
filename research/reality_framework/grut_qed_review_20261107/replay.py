"""Bounded PR #14 review; reuses its DMRG solver, checks its MPO independently.

Run with TeNPy 1.1.1, NumPy and SciPy installed. Writes only --output.
This is a post-result review, not a newly blinded confirmation campaign.
"""
import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import scipy
from scipy.sparse.linalg import eigsh
import tenpy

PACKET = Path(__file__).resolve().parents[1] / "grut_qed_detector_histories_20261107"
sys.path.insert(0, str(PACKET / "mps"))
sys.path.insert(0, str(PACKET / "ed"))
from schwinger_mps import build, gap
from schwinger import hamiltonian


def mpo_dense(model):
    partial = [np.ones((1, 1))]
    for n in range(model.lat.N_sites):
        w = model.H_MPO.get_W(n).to_ndarray()
        d = partial[0].shape[0]
        nxt = [np.zeros((2*d, 2*d)) for _ in range(w.shape[1])]
        for i, left in enumerate(partial):
            for j in range(w.shape[1]):
                nxt[j] += np.kron(left, w[i, j])
        partial = nxt
    assert len(partial) == 1
    return partial[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = {"source_commit": "8dceba0c5f10b1febc5cc130b7a63da74e2372cc",
           "environment": {"python": platform.python_version(), "numpy": np.__version__,
                           "scipy": scipy.__version__, "tenpy": tenpy.__version__}}
    out["preregistration_sha256"] = {
        str(p.relative_to(PACKET)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [PACKET/"ed/PREREGISTRATION.json", PACKET/"mps/PREREGISTRATION_AMENDMENT_01.json"]}
    out["mpo_algebra"] = []
    for N, x, m in [(4, 2., 0.), (6, 9., .25), (8, 16., .5)]:
        model, _ = build(N, x, m)
        dense = mpo_dense(model)
        h, _, states = hamiltonian(N, x, m)
        indices = np.array([sum(int((s >> n) & 1) << (N-1-n) for n in range(N))
                            for s in states])
        err = float(np.max(np.abs(dense[np.ix_(indices, indices)] - h.toarray())))
        out["mpo_algebra"].append({"N": N, "x": x, "m": m, "max_entry_error": err})
        assert err < 1e-11

    # Fresh ED diagonalization, not a comparison only against an archived cache.
    N, x, m = 12, 9., .25
    h, _, _ = hamiltonian(N, x, m)
    energies = np.sort(eigsh(h, k=2, which="SA", tol=1e-12,
                              v0=np.linspace(1, 2, h.shape[0]), return_eigenvectors=False))
    ed_gap = float(np.diff(energies)[0]/(2*np.sqrt(x)))
    dmrg_gap, info = gap(N, x, m, return_info=True)
    out["fresh_ed"] = {"N": N, "x": x, "m": m, "ed_gap": ed_gap,
                       "dmrg_gap": dmrg_gap, "abs_gap_error": abs(dmrg_gap-ed_gap),
                       "abs_ground_energy_error": abs(info["E0"]-float(energies[0]))}
    assert abs(dmrg_gap-ed_gap) < 1e-8
    print("MPO algebra and fresh ED:", out["mpo_algebra"], out["fresh_ed"], flush=True)

    # One unchanged point of the registered massive production grid.
    N, x, m = 36, 9., .25
    dmrg_gap, info = gap(N, x, m, return_info=True)
    ref = json.loads((PACKET/"mps/production_m0.25.json").read_text())["x9_L12"]
    out["production_point"] = {"N": N, "x": x, "m": m, "gap": dmrg_gap,
                               "archived_gap": ref["omega1"],
                               "abs_gap_difference": abs(dmrg_gap-ref["omega1"]),
                               **{k: float(v) for k, v in info.items()}}
    assert abs(dmrg_gap-ref["omega1"]) < 1e-7

    # A transparent scope check on the quoted soft-photon threshold.
    m_e, energy, hbarc, Ce = .51099895e-3, 2.3e6, 1.973269804e-16, 5/24
    a_natural = m_e/(np.sqrt(2*Ce)*energy**2)
    threshold_shift = m_e**2/(2*energy**2)
    # Untuned Wilson mass term adds m*a/2 to the high-energy velocity.
    a_untuned = 2*threshold_shift / (m_e/2 + np.sqrt((m_e/2)**2 + 4*Ce*energy**2*threshold_shift))
    out["lightcone_scope"] = {
        "body_diagonal_tuned_soft_threshold_m": a_natural*hbarc,
        "omitted_ma_half_at_that_spacing": m_e*a_natural/2,
        "retained_mass_deficit": threshold_shift,
        "ratio_omitted_to_retained": (m_e*a_natural/2)/threshold_shift,
        "untuned_positive_mass_body_diagonal_soft_threshold_m": a_untuned*hbarc,
        "scope": "Leading expansion; illustrative body-diagonal soft threshold, not a full observational exclusion."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
