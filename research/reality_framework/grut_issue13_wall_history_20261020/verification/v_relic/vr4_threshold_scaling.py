"""v_relic part 4: velocity scaling of hh -> X Xbar near threshold (is it p-wave, sigma v ~ v^2, or v^4?)."""
import json, numpy as np
import vr_common as C
out = {}
# sigma*v_rel (CM, v_rel = 4p/sqrt(s)) as a function of initial momentum p
rows = []
for p in [0.1, 0.3, 1.0, 3.0, 10.0]:
    s = 4 * (C.mh**2 + p**2)
    sg = C.sigma(np.array([s]), C.M_hh_XX, C.mh, C.mX, 1.0)[0]
    vrel = 4 * p / np.sqrt(s)
    rows.append({"p_GeV": p, "v_rel": vrel, "sigma_v": sg * vrel})
for a, b in zip(rows[:-1], rows[1:]):
    a["log_slope_sigv_vs_v"] = np.log(b["sigma_v"] / a["sigma_v"]) / np.log(b["v_rel"] / a["v_rel"])
out["sigv_vs_vrel"] = rows
# partial-wave content: amplitude even in cos (t<->u), so no l=1; check M(p, cos) ~ p^2 (A + B cos^2)
p = 0.5; s = 4 * (C.mh**2 + p**2)
cs = np.array([0.0, 0.5, 1.0])
out["M_at_p0.5_cos0,0.5,1"] = C.M_hh_XX(s, cs).tolist()
out["M_odd_part_max"] = float(np.max(np.abs(C.M_hh_XX(s, cs) - C.M_hh_XX(s, -cs))))
json.dump(out, open("results/vr4_threshold_scaling.json", "w"), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
