"""Replay the original four ED/DMRG checks without modifying sealed data.

Run from any directory. --output must be outside the archived QED packet.
The original acceptance rule, relative gap error < 1e-5, is unchanged.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
sys.path.insert(0, str(HERE))
from schwinger_mps import gap

CASES = ((0.0, 4.0, 16), (0.0, 16.0, 20),
         (0.25, 8.0, 18), (0.25, 32.0, 24))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, default=PACKET / "ed",
                        help="Default is the repository's archived ed directory.")
    parser.add_argument("--case", type=int, choices=range(4), action="append",
                        help="Optional case index for a diagnostic subset; default all four.")
    args = parser.parse_args(argv)
    output = args.output.resolve()
    if output.is_relative_to(PACKET) or output.exists():
        parser.error("Choose a new output file outside the archived QED packet.")
    reference_dir = args.reference_dir.resolve()
    paths = {0.0: reference_dir / "spectrum_cache.json",
             0.25: reference_dir / "spectrum_cache_m025.json"}
    cache = {m: json.loads(p.read_text()) for m, p in paths.items()}
    selected = sorted(set(args.case)) if args.case is not None else list(range(4))
    rows = []
    for index in selected:
        m, x, n = CASES[index]
        start = time.monotonic()
        omega, info = gap(n, x, m, return_info=True)
        archived = cache[m][f"m{m}_x{x}_N{n}"]
        ref = archived["omega1"]
        rel = abs(omega - ref) / abs(ref) if ref else float("inf")
        finite = bool(np.all(np.isfinite([omega, ref, info["E0"], rel])))
        row = {"case_index": index, "m": m, "x": x, "N": n,
               "dmrg": float(omega) if np.isfinite(omega) else None, "ed": ref,
               "E0_per_N_dmrg": float(info["E0"] / n) if np.isfinite(info["E0"]) else None,
               "E0_per_N_ed": archived["E0_per_N"],
               "rel": float(rel) if np.isfinite(rel) else None,
               "pass": bool(finite and rel < 1e-5),
               "sec": round(time.monotonic() - start, 2)}
        rows.append(row)
        print(json.dumps(row, allow_nan=False), flush=True)
    ok = all(row["pass"] for row in rows)
    payload = {"gate": "rel<1e-5 vs exact diagonalization", "pass": ok,
               "full_four_case_replay": selected == list(range(4)),
               "reference_dir": str(reference_dir),
               "reference_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                    for p in paths.values()},
               "cases": rows}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        json.dump(payload, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print("GATE", "PASS" if ok else "FAIL", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
