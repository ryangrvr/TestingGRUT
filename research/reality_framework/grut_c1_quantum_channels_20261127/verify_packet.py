"""Verify the C1 archive and optionally regenerate the saved profile splines."""
import argparse
import gzip
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
ENERGY_ARCHIVE = ARCHIVE / "grut_candidate_1_energy_20261019"
REFERENCE = ARCHIVE / "grut_candidate_1_20261018/code/results_profile.npz"
POLE = ARCHIVE / "grut_candidate_1_fluctuations_20261019/pole_results.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def verify():
    protocol_path = HERE / "PROTOCOL.json"
    protocol, result = read(protocol_path), read(HERE / "RESULTS.json")
    assert result["protocol_sha256"] == digest(protocol_path)
    source_paths = {name: HERE / name for name in ("run.py", "branch.py", "spectrum.py", "tails.py", "radial.cc")}
    source_paths.update({path.name: path for path in (POLE, REFERENCE)})
    assert result["source_hashes"] == {name: digest(path) for name, path in source_paths.items()}
    assert result["completed"] and result["verdict"] == "A"
    assert all(row["passes"] for row in result["branches"].values())
    assert all(row["passes"] for row in result["exact_charge_spectra"].values())
    assert result["independent_reference_check"]["passes"]
    parent, neighbor, half = [result["exact_charge_spectra"][str(q)] for q in protocol["charge_order"]]
    emission, fission, chemical = [
        result[key] for key in ("direct_emission", "equal_fission", "chemical_potential_crosscheck")
    ]
    pole = read(POLE)["matching_result"]["strict_pole_mass_over_v"] / protocol["g_X"]
    assert emission["pole"] == pole
    direct = math.fsum([emission["tree_difference"], *emission["paired_loop_spectrum_difference"].values(),
                        emission["paired_graph_difference"], -pole])
    integrated = math.fsum([chemical["tree_integral"], chemical["loop_integral"], -pole])
    independent_sum = math.fsum([parent["E0plus1_finite"], -neighbor["E0plus1_finite"], -pole])
    total_error = math.fsum([parent["total_allowance"], neighbor["total_allowance"],
                             protocol["paired_allowance"]["pole_numerical"]])
    assert direct == emission["value"] and integrated == chemical["value"]
    assert abs(independent_sum - direct) < 2e-11
    assert independent_sum + total_error < 0
    assert emission["passes"] and chemical["passes"]
    assert fission["value"] == parent["E0plus1_finite"] - 2 * half["E0plus1_finite"]
    assert fission["allowance"] == parent["total_allowance"] + 2 * half["total_allowance"]
    for row in (emission, chemical, fission):
        assert row["signed_interval"] == [row["value"] - row["allowance"], row["value"] + row["allowance"]]
        assert row["signed_interval"][1] < 0
    log = gzip.decompress((HERE / "EXECUTION.jsonl.gz").read_bytes())
    assert hashlib.sha256(log).hexdigest() == "469d1e1923dd32a00bc0af368f4d991bae506f6e09e97e5088a44a8977040e86"
    events = [json.loads(line) for line in log.splitlines()]
    counts = Counter(row["event"] for row in events)
    assert counts["jost_point"] == result["workload"]["jost_points_scheduled"]
    assert counts["fd_start"] == counts["fd_complete"] == result["workload"]["fd_batches"]
    assert all(row["valid"] for row in events if row["event"] == "jost_point")
    assert counts["jost_point"] <= protocol["stopping"]["maximum_jost_points"]
    assert counts["fd_complete"] <= protocol["stopping"]["maximum_fd_batches"]
    assert result["elapsed_seconds"] < protocol["stopping"]["wall_seconds"]
    print(f"Absolute-sum emission interval [{independent_sum-total_error:.12f}, "
          f"{independent_sum+total_error:.12f}] g_X v")
    print(f"Verified {counts['jost_point']} continuum points and {counts['fd_complete']} radial batches")
    return protocol, result


def report_profiles(protocol, result, destination):
    import numpy as np

    import branch

    seed = branch.reference_profile()
    settings = protocol["profiles"]
    data = {}
    for q in protocol["charge_order"]:
        for scale in ("coarse", "fine"):
            label = f"{scale}_{q}"
            profile = branch.stationary(seed, *settings[scale], settings, charge=q, label=label)
            prior = result["branches"][str(q)][scale]
            assert profile.statistics["success"]
            comparisons = {"omega": 2e-12, "charge": 1e-7, "energy": 2e-8}
            for key, allowance in comparisons.items():
                assert abs(profile.statistics[key] - prior[key]) <= allowance, (label, key)
            data[f"{label}_knots"] = profile.spline.x
            data[f"{label}_coefficients"] = profile.spline.c
            data[f"{label}_omega"] = np.array(profile.omega)
    with destination.open("xb") as stream:
        np.savez_compressed(stream, **data)
    print(f"Regenerated six stationary spline records in {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regenerate-profiles", type=Path, metavar="NEW_NPZ")
    options = parser.parse_args()
    settings, results = verify()
    if options.regenerate_profiles:
        report_profiles(settings, results, options.regenerate_profiles)
