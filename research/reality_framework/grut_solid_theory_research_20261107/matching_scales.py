"""Arithmetic supporting the GRUT research pass; no physics solver is run.

Read the frozen history parameter ladder without modifying it. The mass
shift is the published bulk staggered-fermion improvement prescription,
not a certified finite-interval mass calibration. Both length conventions
are displayed because their boundary dictionaries must be matched.

Usage, from any directory:
    python /path/to/matching_scales.py
"""
from pathlib import Path
import hashlib
import json
import math
import sys

source = (Path(__file__).resolve().parent.parent /
          "grut_qed_detector_histories_20261107/ed/histories_results.json")
data = json.loads(source.read_text())
ladder = sorted({(int(row["N"]), float(row["x"])) for row in data})
result = {
    "kind": "deterministic_parameter_arithmetic_not_a_new_simulation",
    "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "mass_prescription": "m_lat/g = m_R/g - 1/(8*sqrt(x))",
    "caveat": "Bulk improvement; finite open-boundary matching remains to be established.",
    "frozen_ladder": [
        {
            "N": n,
            "x": x,
            "ag": 1 / math.sqrt(x),
            "reported_Na_g": n / math.sqrt(x),
            "endpoint_span_N_minus_1_a_g": (n - 1) / math.sqrt(x),
            "bulk_improved_m_lat_over_g_for_m_R_0": -1 / (8 * math.sqrt(x)),
            "bulk_improved_m_lat_over_g_for_m_R_half": 0.5 - 1 / (8 * math.sqrt(x)),
            "dimensionless_W_mass_parameter_shift": -0.25,
        }
        for n, x in ladder
    ],
}
json.dump(result, sys.stdout, indent=2, allow_nan=False)
print()
