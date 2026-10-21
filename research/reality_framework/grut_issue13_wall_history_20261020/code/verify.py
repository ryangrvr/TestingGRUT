#!/usr/bin/env python3
"""verify.py -- integrity and cross-module checks for packet grut_issue13_wall_history_20261020.

Usage: python3 code/verify.py      (exit status 0 = all checks pass)
1. sha256 of every file listed in MANIFEST.json.
2. Re-runs code/history.py into a temporary file and compares it with results/history.json.
3. Cross-checks between history.py's own recomputation and the module results it quotes.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

PKT = Path(__file__).resolve().parent.parent
fails = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        fails.append(msg)


def close(a, b, rel):
    return a is not None and b is not None and abs(a - b) <= rel * max(abs(a), abs(b))


def mod(rel):
    return json.loads((PKT / "modules" / rel).read_text())


# 1. manifest
man = json.loads((PKT / "MANIFEST.json").read_text())
for rel, h in man["files"].items():
    p = PKT / rel
    check(p.exists() and hashlib.sha256(p.read_bytes()).hexdigest() == h, f"hash {rel}")

# 2. history.py replay
ref = json.loads((PKT / "results" / "history.json").read_text())
with tempfile.TemporaryDirectory() as td:
    tmp = Path(td) / "pkt"
    tmp.mkdir()
    for sub in ("code", "modules"):
        (tmp / sub).symlink_to(PKT / sub)
    # history.py writes to PKT/results relative to its own location; run a copy that writes to tmp
    src = (PKT / "code" / "history.py").read_text().replace(
        'PKT = Path(__file__).resolve().parent.parent', f'PKT = Path(r"{tmp}")')
    (Path(td) / "h.py").write_text(src)
    subprocess.run([sys.executable, "-I", str(Path(td) / "h.py")], check=True, capture_output=True)
    new = json.loads((tmp / "results" / "history.json").read_text())
check(json.dumps(new, sort_keys=True) == json.dumps(ref, sort_keys=True), "history.py replay identical to results/history.json")

# 3. cross-checks
ph = ref["phase"]
t0 = mod("charge_phase/results/t0_phase.json")
check(close(ph["n_c_over_v3"], t0["n_c_over_v3"], 1e-9), "n_c matches charge_phase")
check(close(ph["n_s1_over_v3"], t0["branch1_spinodal_n_over_v3"], 1e-9), "n_s1 matches charge_phase")
check(close(ph["n_t_over_v3"], t0["homogeneous_global_switch_n_over_v3"], 1e-6), "n_t matches charge_phase")
check(close(ph["n_star_over_v3"], t0["binodal_n_star_over_v3"], 1e-4), "n_* matches charge_phase (grid 5e-6 v^3)")
check(close(ph["Qmatter_E_per_Q_GeV"], t0["mu_c_GeV"], 1e-6), "Q-matter E/Q = m* matches charge_phase")
tp = mod("thermal_potential/results/summary.json")
check(abs(ph["Tc_leading_over_v"] - tp["Tc_leading_highT_over_v"]) < 1e-12, "leading T_c matches thermal_potential")
co = mod("cosmo_observable/results/headline.json")
Tdom_mod = co["T_dom_keV (rho_w=rho_rad, A=0.8, backreaction) [C]"]
Tpar = ref["network"]["T_parity_keV_A"]["0.8"]
check(0.85 < Tpar / Tdom_mod < 1.0, f"simple parity {Tpar:.3f} keV vs cosmo_observable backreaction {Tdom_mod:.3f} keV (expected 0.85-1.0)")
dN_mod = co["stress n=4 (conversion to radiation at T_dom): DeltaNeff at T_rec [C]"]
check(close(ref["network"]["DeltaNeff_conservative_min"], dN_mod, 2e-3), "conservative DeltaN_eff matches cosmo_observable stress n=4")
H1, H2 = ref["network"]["H0_pred_n1_kmsMpc_A0.8_simple"], ref["network"]["H0_pred_n2_kmsMpc_A0.8_simple"]
Hm = co["H0_pred_kmsMpc n=1, n=2 [C]"]
check(0.5 < H1 / Hm[0] < 1.5 and 0.5 < H2 / Hm[1] < 1.5, "simple H0 predictions within 50% of cosmo_observable (matter and backreaction omitted here)")
ce = mod("charged_exit_sim/results/headline.json")
pf = ce["sign_statistics"]["charged_all"]["mean_form"]
check(abs(pf - 0.5) < 0.01, f"charged-exit p(+) at formation {pf:.4f} within 0.01 of 1/2")
ns = mod("network_sim/results/headline.json")
A3 = ns["A_3D_PRS_256_box_safe"]["mean"]
check(0.7 < A3 < 0.9, f"3D PRS A = {A3:.3f} inside the literature band 0.8 +- 0.1")
hs = ref["histories"]["H_star"]
check(hs["omega_sym_estimate"] < 0.12 and hs["maser_product_GeV"] < 63.93, "H-star passes relic and maser product tests")
cs = ref["histories"]["C_star"]
check(cs["n_Q(T_Q)/n_c"] < 1 and cs["below regime-A boundary 145.6 GeV (minority_relic)"], "C-star charge sub-critical and inside regime A")
mr = mod("minority_relic/results/headline.json")
check(close(mr["collapse_energy"]["kappa_E_RD"], 7.786625485748964, 1e-9), "minority_relic kappa_E as quoted")

print(f"\n{len(fails)} failure(s)")
sys.exit(1 if fails else 0)
