#!/usr/bin/env python3
"""Addenda to bcd_envelope.py: flip thresholds for the decisive rows, open-curvature envelope check,
and floor sensitivity to the omega_b / m_b margins. Imports the blind module only."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import brentq

import bcd_envelope as B

out = {}
omb_c, mb_c = B.OMB - 3 * B.SOMB, B.MB * 1.01
Hc = B.hfloor(omb_c, mb_c, B.MSTAR * B.R)[0]
Hn = B.hfloor(B.OMB, B.MB, B.MSTAR * B.R)[0]
for name, D, ep, em, vv, sv, ra, dec in B.TABLE1:
    if name not in ("NGC 5765b", "UGC 3789", "CGCG 074-064"):
        continue
    vobs = vv + 3 * sv
    # peculiar velocity needed to bring D_A,max up to D_low (conservative floor)
    f = lambda vp: B.DA_env(B.zcos(vobs, vp), Hc) - (D - 3 * em)
    vp_flip = brentq(f, -20000, 5000)
    # number of lower-sigma widths at which the row flips, at v_pec = -750 and -1500
    k750 = (D - B.DA_env(B.zcos(vobs, -750.0), Hc)) / em
    k1500 = (D - B.DA_env(B.zcos(vobs, -1500.0), Hc)) / em
    # floor H that would flip the row at -750 / -1500
    out[name] = dict(vpec_needed_to_flip_kms=vp_flip, D_sigma_multiple_to_flip_at_m750=k750,
                     D_sigma_multiple_to_flip_at_m1500=k1500,
                     table4_most_negative_flow_vpec=min(B.TABLE4_PEC[name]))
# open-curvature check: H^2 = Hf^2 (1+z)^3 + K (1+z)^2, K>0, D_M with sinh; must stay <= flat EdS envelope
worst = -np.inf
for K in np.concatenate([np.linspace(0, 1, 11), np.logspace(0, 4, 9)]) * Hc**2:
    for z in (0.001, 0.005, 0.01, 0.03, 0.04, 0.3, 1.0, 3.0):
        H0 = np.sqrt(Hc**2 + K)
        val = B.DA_general(z, H0, Hc**2 / H0**2, K / H0**2, 0.0) / B.DA_env(z, Hc) - 1
        worst = max(worst, val)
out["open_curvature_max_DA_over_env_minus1"] = worst
# alternative conservative choices: nominal floor (no margins), amu mass, 5-sigma omega_b
for lab, omb, mb in [("nominal_floor", B.OMB, B.MB), ("amu_mb", B.OMB, B.MB_AMU), ("omb_minus5sig_mb1.02", B.OMB - 5 * B.SOMB, B.MB * 1.02)]:
    H = B.hfloor(omb, mb, B.MSTAR * B.R)[0]
    rows = B.evaluate(-750.0, H, omb, mb)
    out["floor_variant_" + lab] = dict(Hfloor=H, mismatches=[r["galaxy"] for r in rows if r["mismatch"]],
                                       tightest=min(rows, key=lambda r: r["mX_absr_max_GeV"])["galaxy"],
                                       mX_absr_max=min(r["mX_absr_max_GeV"] for r in rows),
                                       mstar_absr_max=min(r["mstar_absr_max_GeV"] for r in rows))
# the bound evaluated without the 1/mu factor (i.e. if the quantum X-particle mass mX, not m*, were the floor)
out["Hfloor_if_floor_were_mX_cons"] = B.hfloor(omb_c, mb_c, B.MX * B.R)[0]
# all-six-rows-flip-to-consistent condition expressed as a threshold on m*|r| (conservative, -750)
rows = B.evaluate(-750.0, Hc, omb_c, mb_c)
out["mstar_absr_threshold_all_rows_m750"] = min(r["mstar_absr_max_GeV"] for r in rows)
rows = B.evaluate(-1500.0, Hc, omb_c, mb_c)
out["mstar_absr_threshold_all_rows_m1500"] = min(r["mstar_absr_max_GeV"] for r in rows)
json.dump(out, sys.stdout, indent=1, default=float)
print()
