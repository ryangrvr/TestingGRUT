"""Post-result spot replays and exploratory matched-variance controls.

No frozen acceptance window or production result is changed.
Controls match Var(c F) by changing detector strength, not the preparation.
They do not establish a continuum limit.
"""
import argparse
import json
import sys
from pathlib import Path
import numpy as np

PACKET = Path(__file__).resolve().parents[1] / "grut_qed_detector_histories_20261107"
sys.path.insert(0, str(PACKET/"ed"))
from histories import run


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    archive = json.loads((PACKET/"ed/histories_results.json").read_text())
    out = {"scope": "Spot replays plus exploratory detector-strength control, not a new holdout", "cases": []}
    for N, x in [(10, 4.), (14, 9.)]:
        baseline = []
        for m in [0., .5]:
            row = run(N, x, m)
            ref = next(r for r in archive if r["N"] == N and r["m_over_g"] == m)
            row["difference_from_archived_deviation"] = row["max_dev_exact_vs_twopoint"] - ref["max_dev_exact_vs_twopoint"]
            assert abs(row["difference_from_archived_deviation"]) < 1e-9
            baseline.append(row)
        c = 2*np.sqrt(baseline[0]["var_cF"]/baseline[1]["var_cF"])
        matched = run(N, x, .5, c=c)
        matched["c"] = float(c)
        assert abs(matched["var_cF"]-baseline[0]["var_cF"]) < 1e-9
        out["cases"].append({"baseline": baseline, "massive_matched_variance": matched})
        print(N, [r["max_dev_exact_vs_twopoint"] for r in baseline],
              "matched", matched["max_dev_exact_vs_twopoint"], flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2)+"\n")


if __name__ == "__main__":
    main()
