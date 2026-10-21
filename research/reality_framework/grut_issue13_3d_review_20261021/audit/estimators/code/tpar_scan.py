"""T_par (wall-radiation parity) with an area-count A and with energy-weighted A = R A, R in [1.1236, 1.2625],
using Claude's own cosmo_observable.core.history (imported read-only).
Usage: python3 -I tpar_scan.py <claude_packet_dir> <out_json>"""
import json, sys
sys.path.insert(0, sys.argv[1] + "/modules/cosmo_observable/code")
import core
s = 0.456 ** 2
R = {"area": 1.0, "E_lo": (1 - s) ** -0.5, "E_hi": 1 / (1 - s)}
rows = []
for A in (0.7, 0.8, 0.9):
    for conv in ("t", "H"):
        for br in (False, True):
            row = {"A_area": A, "convention": conv, "backreaction": br}
            for lab, f in R.items():
                row[lab + "_keV"] = core.history(core.SIGMA_BENCH, A=A * f, backreaction=br, convention=conv)["T_switch_GeV"] * 1e6
            rows.append(row)
summary = {lab: [min(r[lab + "_keV"] for r in rows), max(r[lab + "_keV"] for r in rows)] for lab in R}
central = [r for r in rows if r["A_area"] == 0.8 and r["convention"] == "t"]
json.dump({"R": R, "rows": rows, "range_keV": summary, "central_A0.8_t": central}, open(sys.argv[2], "w"), indent=1)
print(json.dumps({"range_keV": summary, "central_A0.8_t": central}, indent=1))
