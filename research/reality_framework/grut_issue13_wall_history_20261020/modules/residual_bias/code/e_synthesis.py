"""Synthesis: residual systematic bias m0 = 2p - 1 per restoration route, compared with PLACEHOLDER thresholds
(eps_persist ~ 1e-20, eps_safe ~ 1e-17, from modules/bias_threshold when available; here placeholders only).

Classification of |m0| against the placeholders (labels only, not a result of this module):
  |m0| < eps_persist           -> 'below eps_persist'  (network persists in the bias_threshold sense)
  eps_persist <= |m0| < eps_safe -> 'between'          (late-decaying network)
  |m0| >= eps_safe             -> 'above eps_safe'     (early decay)
For orientation only, the naive extrapolation of network_sim's 2D PRS fit K10 = 0.262 (|m0|/2)^-1.255 [Q via the
task text] is listed; it is not a result of this module.
"""
import json
import math
import os

import numpy as np

from common import RESULTS, dump

EPS_PERSIST = 1e-20
EPS_SAFE = 1e-17


def cls(log10_abs):
    if log10_abs is None:
        return None
    if log10_abs < math.log10(EPS_PERSIST):
        return "below eps_persist"
    if log10_abs < math.log10(EPS_SAFE):
        return "between eps_persist and eps_safe"
    return "above eps_safe"


def k10(log10_abs):
    return math.log10(0.262) - 1.255 * (log10_abs - math.log10(2.0))


def load(name):
    p = os.path.join(RESULTS, name)
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        return json.load(fh)


def d_metrics(r):
    """p(+) definitions from the stored series: first both-signs (>= 5%) output after restoration (as in
    charged_exit_sim), 'settled' = first output with f+ + f- >= 0.5 after restoration, and end."""
    ser = np.array(r["series"])
    fp, fm, fc = ser[:, 4], ser[:, 5], ser[:, 6]
    restored = np.maximum.accumulate(fc >= 0.5)
    out = {}
    i_form = np.where(restored & (fp >= 0.05) & (fm >= 0.05))[0]
    out["p_form"] = float(fp[i_form[0]] / (fp[i_form[0]] + fm[i_form[0]])) if i_form.size else None
    out["tau_form"] = float(ser[i_form[0], 0]) if i_form.size else None
    # decided population: mean of f+/(f+ + f-) over the window [tau_form + 5, tau_form + 25] (after the first
    # two-sign output the coherent swing has to die out; single early outputs can catch a coherent excursion)
    if i_form.size:
        t_f = ser[i_form[0], 0]
        w = (ser[:, 0] >= t_f + 5) & (ser[:, 0] <= t_f + 25) & (fp + fm > 0)
        pw = fp[w] / (fp[w] + fm[w])
        out["p_window"] = float(pw.mean()) if w.any() else None
        out["p_window_sd"] = float(pw.std()) if w.any() else None
        out["coherent_fraction_window"] = float((np.abs(ser[w, 2]) / ser[w, 3]).mean()) if w.any() else None
        out["n_over_nc_form"] = float(ser[i_form[0], 7] / 0.316228)
    else:
        out["p_window"] = out["p_window_sd"] = out["coherent_fraction_window"] = out["n_over_nc_form"] = None
    out["p_end"] = float(fp[-1] / (fp[-1] + fm[-1])) if fp[-1] + fm[-1] > 0 else None
    out["f_vac_end"] = float(fp[-1] + fm[-1])
    # late window: last 100 time units
    tl = ser[:, 0] >= ser[-1, 0] - 100.0
    out["late_mean_S"] = float(ser[tl, 2].mean())
    out["late_p_plus"] = float((fp[tl] / np.maximum(fp[tl] + fm[tl], 1e-12)).mean())
    out["late_f_core"] = float(fc[tl].mean())
    # coherent-envelope decay rate: windows of 10, from the first window below 0.6 to the first below 0.03
    t_w, e_w = [], []
    for t0 in np.arange(ser[0, 0], ser[-1, 0], 10.0):
        w = (ser[:, 0] >= t0) & (ser[:, 0] < t0 + 10.0)
        t_w.append(t0 + 5.0)
        e_w.append(np.abs(ser[w, 2]).max())
    t_w, e_w = np.array(t_w), np.array(e_w)
    i0 = np.where(e_w < 0.6)[0]
    rate = None
    if i0.size:
        i1 = np.where((e_w < 0.03) & (np.arange(len(e_w)) > i0[0]))[0]
        if i1.size and i1[0] - i0[0] >= 2:
            sl = slice(i0[0], i1[0] + 1)
            rate = float(-np.polyfit(t_w[sl], np.log(e_w[sl]), 1)[0])
    out["coherent_envelope_decay_rate_per_ell"] = rate
    out["envelope_floor_late"] = float(e_w[-5:].mean())
    # coherent envelope after settling: max |mean S| over windows of 10
    tau = ser[:, 0]
    env = []
    for t0 in np.arange(tau[0], tau[-1], 10.0):
        w = (tau >= t0) & (tau < t0 + 10.0)
        env.append([float(t0), float(np.abs(ser[w, 2]).max()), float(ser[w, 3].mean())])
    out["envelope_[tau0,max|meanS|,mean|S|]"] = env
    return out


def main():
    A = load("a_damping.json")
    B = load("b_inflation.json")
    out = {"module": "residual_bias", "part": "synthesis", "placeholders": dict(eps_persist=EPS_PERSIST,
                                                                                eps_safe=EPS_SAFE),
           "convention": "m0 = 2 p(+) - 1 at the Z2 exit; + is the previously selected vacuum"}
    routes = {}
    # ---------------- A
    if A:
        s = A["H_star [P inputs]"]
        D = s["D_full"]
        lm = s["decision_full"]["log10_m0_damped"]
        band = [s["decision_full_x2"]["log10_m0_damped"], s["decision_full_x0.5"]["log10_m0_damped"]]
        routes["A_Hstar_thermal [P: xi=5e-4, bath from T_SM=1e8 GeV, phi_0 = v]"] = dict(
            log10_abs_m0_central=lm, log10_abs_m0_band=band, D_amplitude_efolds=D,
            class_central=cls(lm), class_band=[cls(band[0]), cls(band[1])], status="computed [C] with [E] band")
        routes["A_hubble_only (no plasma damping, same history)"] = dict(
            abs_m0=s["hubble_only"]["m0_hubble_only"], r=s["hubble_only"]["r_hubble_only"],
            log10_abs_m0_central=0.0, class_central=cls(0.0),
            meaning="coherent amplitude >> soft rms: the whole patch takes the sign fixed by the oscillation phase "
                    "(single vacuum, no network); m0 = +-1", status="computed [C, E]")
        pm = s["production_map_WKB [Q formula]"]
        routes["A_production_map_WKB_formula_at_Hstar [Q formula, not applicable]"] = dict(
            abs_m0=pm["m0_if_used"], log10_abs_m0_central=math.log10(pm["m0_if_used"]),
            class_central=cls(math.log10(pm["m0_if_used"])), note=pm["applicability"])
        routes["A_xi_min_for_m0_below (production at T_h >= 10 T_c)"] = A["xi_min_for_log10_m0_below [C]"]
        routes["A_xi_min_vs_production_factor"] = A["xi_min_vs_production_factor [C]"]
        routes["A_dilute_gas"] = A["dilute_gas [D, E]"]
        scan = [dict(xi=r["xi"], f=r["T_h_prod_over_Tc"], log10_abs_m0=r["log10_m0"], cls=cls(r["log10_m0"]),
                     kinetic=r["kinetic_bath"]) for r in A["scan_xi_Tprod [C]"]]
        routes["A_scan"] = scan
    # ---------------- B
    if B:
        rows = []
        for r in B["rows [C]"]:
            for N in ("N50", "N60"):
                q = r[N]
                for start, key in (("fresh", "m0_fresh_after_trapping"), ("relaxed", "m0_relaxed_after_trapping")):
                    m = q[key]
                    la = math.log10(abs(m)) if m != 0 else None
                    rows.append(dict(H_GeV=r["H_GeV"], N=N, start=start, m0=m, log10_abs_m0=la, cls=cls(la),
                                     patch_spread_sd_one_efold=q["fresh_one_efold_subpatch_sd"] if start == "fresh"
                                     else None, patch_rms_equilibrium_draw=q["m_patch_rms_equilibrium_draw"],
                                     naive_log10_K10=k10(la) if la is not None else None))
        routes["B_trapping_only [C; P: instant RD, no restoration]"] = rows
        routes["B_module_emulation"] = B["module_emulation_h0.02 [C]"]
        routes["B_grid_convergence"] = B["grid_convergence [C]"]
        cm = B["coherent_mean_fresh_start [C]"]
        prop = []
        if A:
            lA = A["H_star [P inputs]"]["decision_full"]["log10_m0_damped"]
            for c in cm:
                prop.append(dict(H_GeV=c["H_GeV"], mean_phi_end_over_v=c["mean_phi_over_v"],
                                 log10_abs_m0_after_Hstar_thermal_bath=lA + math.log10(abs(c["mean_phi_over_v"])),
                                 note="upper estimate: route-A H-star damping with phi_0 = <phi_end>; extra dephasing "
                                      "of the large-amplitude quartic oscillation only lowers it"))
        routes["B_then_A [C, E]"] = prop
        # B followed by a NON-thermal (free-streaming) restoring gas that exits at T_SM,end (e.g. the heavy-inflaton
        # producer of production_map, exits 1e8-6e9 GeV) [E]: weak damping D = dilute_D(A), A = lambda v^2/H(T_end).
        # Band: lower = trapping-map law 0.36 (v/H)^3 e^-D (oscillating final-sign map), upper = sqrt(2/pi)
        # <phi_end>/phi_rms e^-D (coherent ensemble mean, no dephasing credited).
        from a_damping import dilute_D
        from common import hubble, LPHI, V
        cm_by_H = {c["H_GeV"]: c for c in cm}
        bdil = []
        for H in (1e5, 1e6, 1e8):
            c = cm_by_H[H]
            for T_end in (2.6e7, 1e8, 1e9, 6e9):
                A_ = LPHI * V ** 2 / hubble(T_end)
                Dd = dilute_D(A_)
                lo = math.log10(0.36 * (V / H) ** 3) - Dd / math.log(10)
                hi = math.log10(math.sqrt(2 / math.pi) * abs(c["mean_over_rms"])) - Dd / math.log(10)
                bdil.append(dict(H_inf_GeV=H, T_SM_end_GeV=T_end, A_GeV=A_, D_dilute=Dd,
                                 log10_abs_m0_band=[lo, hi], cls_band=[cls(lo), cls(hi)]))
        routes["B_then_nonthermal_gas [E]"] = bdil
    # ---------------- C
    routes["C_symmetric_thermal_ensemble [D]"] = dict(
        m0=0.0, statement="P rho P^dag = rho (thermal or any Z2-invariant state, Z2 commuting with H and Q_X) => "
                          "Tr(rho O_odd) = 0 => p(+) = p(-) exactly; m0 = 0 identically",
        realization_noise="white-noise patch fluctuation ~ N_dom^-1/2, not a systematic bias")
    # ---------------- D
    rundir = os.path.join(RESULTS, "d_runs")
    drows = []
    if os.path.isdir(rundir):
        for fn in sorted(os.listdir(rundir)):
            with open(os.path.join(rundir, fn)) as fh:
                r = json.load(fh)
            met = d_metrics(r)
            drows.append(dict(name=fn[:-5], kappa=r["kappa"], n_over_ns1=r["n_over_ns1"], mode=r["mode"],
                              tau_i=r["tau_i"], seed=r["seed"], outcome=r["outcome"], **met))
    routes["D_lattice_runs [C]"] = drows
    # aggregate D by (mode, tau_i, kappa)
    agg = {}
    for d in drows:
        key = "%s_t%s_k%.1f" % (d["mode"], d["tau_i"], d["kappa"])
        agg.setdefault(key, []).append(d)
    dsum = {}
    for key, L in agg.items():
        ps = [x["p_window"] for x in L if x["p_window"] is not None]
        pf = [x["p_form"] for x in L if x["p_form"] is not None]
        m = [2 * p - 1 for p in ps]
        dsum[key] = dict(n=len(L), p_window=ps, p_form=pf,
                         n_over_nc_form=[x["n_over_nc_form"] for x in L],
                         m0_window_mean=float(np.mean(m)) if m else None,
                         m0_window_sd=float(np.std(m, ddof=1)) if len(m) > 1 else None,
                         m0_window_sem=float(np.std(m, ddof=1) / math.sqrt(len(m))) if len(m) > 1 else None,
                         p_end=[x["p_end"] for x in L], outcomes=[x["outcome"] for x in L],
                         late_mean_S=[x["late_mean_S"] for x in L], late_p_plus=[x["late_p_plus"] for x in L],
                         late_f_core=[x["late_f_core"] for x in L],
                         envelope_decay_rate_per_ell=[x["coherent_envelope_decay_rate_per_ell"] for x in L],
                         envelope_floor_late=[x["envelope_floor_late"] for x in L])
    routes["D_summary [C]"] = dsum
    out["routes"] = routes
    dump("headline.json", out)
    return out


if __name__ == "__main__":
    o = main()
    for k, v in o["routes"].get("D_summary [C]", {}).items():
        print(k, v["m0_window_mean"], v["m0_window_sd"], [round(p, 3) for p in v["p_window"]],
              [round(p, 2) if p else None for p in v["n_over_nc_form"]], "late p", [round(p, 3) for p in v["late_p_plus"]],
              "fcore", [round(p, 3) for p in v["late_f_core"]], "rate", v["envelope_decay_rate_per_ell"])
