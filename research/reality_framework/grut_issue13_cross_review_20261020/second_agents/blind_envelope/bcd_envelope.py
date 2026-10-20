#!/usr/bin/env python3
"""(b)-(d) Blind recomputation of the flat-FRW expansion floor, the D_A upper envelope, and the
Pesce et al. 2020 (ApJL 891 L1) maser comparison for GRUT Candidate 1.

Written without reading Codex's audit.py/report.py/DERIVATION.md/COMPARISON.json/DISTANCE_ENVELOPE.csv/
maser_table1.json or any Claude script.

Inputs and provenance (transcribed by this script's author):
 * Pesce+2020 Table 1: arXiv:2001.09213v2 LaTeX source main.tex lines 124-129 (identical to pdftotext -layout of
   the v2 PDF, layout lines 211-222). Velocities: CMB frame, optical convention (main.tex line 131).
 * omega_b = 0.02233 +- 0.00036: Mossa et al. 2020 accepted manuscript (Padua repository, mossa_nature_R2.pdf,
   2 Sep 2020) text lines 118-119 and Table 1 row "D+3nu (with new LUNA data)"; Nature 587, 210 (2020).
 * CMB dipole (Planck 2018 I): 369.82 km/s toward (l,b) = (264.021, 48.253) deg.  [from memory]
 * Local Group apex (Karachentsev & Makarov 1996): 316 km/s toward (l,b) = (93, -4) deg. [from memory]
 * Galaxy ICRS coordinates: SIMBAD (queried 2026-10-20), simbad/coords.csv, simbad/cgcg.txt.
"""
import json
import sys

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

C = 299792.458  # km/s

# ---------------- (a) imported result ----------------
v, gX, lp, lX = 1000.0, 0.1, 0.01, 0.001
MX = gX * v
MSTAR = v * min(gX, (lp * lX) ** 0.25)
MU = MSTAR / MX
R = 1.0  # YQ0/YB0

# ---------------- (b) floor ----------------
OMB, SOMB = 0.02233, 0.00036
MB = 0.938
MB_AMU = 0.93149410242


def hfloor(omb, mb, mstar_r):
    omh = mstar_r / mb * omb
    return 100.0 * np.sqrt(omb + omh), omh


# baryon-to-photon cross-check of the m_b convention implied by Omega_b h^2 <-> eta
def mass_per_baryon_from_eta_conversion(eta10_per_ombh2):
    kB = 8.617333262e-5  # eV/K
    hbarc = 197.3269804e-7  # eV cm
    T = 2.7255
    zeta3 = 1.2020569031595942
    ngamma = 2 * zeta3 / np.pi**2 * (kB * T / hbarc) ** 3  # cm^-3
    rhoc_over_h2 = 1.87834e-29 / 1.78266192e-24  # g/cm^3 -> GeV/cm^3
    return rhoc_over_h2 / (eta10_per_ombh2 * 1e-10 * ngamma), ngamma, rhoc_over_h2


# ---------------- (c) envelope ----------------
def DA_env(z, H):
    return 2 * C / (H * (1 + z)) * (1 - (1 + z) ** -0.5)


def DA_numeric(z, Hfun):
    chi = quad(lambda zz: C / Hfun(zz), 0, z, epsabs=0, epsrel=1e-12)[0]
    return chi / (1 + z)


def DA_general(z, H0, Om, Ok, OL, Or=0.0, w_extra=None, Ox=0.0):
    def E(zz):
        e2 = Om * (1 + zz) ** 3 + Ok * (1 + zz) ** 2 + OL + Or * (1 + zz) ** 4
        if w_extra is not None:
            e2 += Ox * (1 + zz) ** (3 * (1 + w_extra))
        return np.sqrt(e2)
    chi = quad(lambda zz: 1 / E(zz), 0, z, epsrel=1e-12)[0] * C / H0
    if abs(Ok) < 1e-14:
        dm = chi
    elif Ok > 0:
        dh = C / H0
        dm = dh / np.sqrt(Ok) * np.sinh(np.sqrt(Ok) * chi / dh)
    else:
        dh = C / H0
        dm = dh / np.sqrt(-Ok) * np.sin(np.sqrt(-Ok) * chi / dh)
    return dm / (1 + z)


# ---------------- (d) data ----------------
# name, Dmed, +err, -err, v_cmb_opt, sigma_v, RA, Dec (ICRS deg)
TABLE1 = [
    ("UGC 3789", 51.5, 4.5, 4.0, 3319.9, 0.8, 109.87871666666668, 59.35513333333334),
    ("NGC 6264", 132.1, 21.0, 17.0, 10192.6, 0.8, 254.31719706930997, 27.84960398059),
    ("NGC 6323", 109.4, 34.0, 23.0, 7801.5, 1.5, 258.32518921447, 43.78243674483999),
    ("NGC 5765b", 112.2, 5.4, 5.1, 8525.7, 0.7, 222.71460599999997, 5.1144870000000004),
    ("CGCG 074-064", 87.6, 7.9, 7.2, 7172.2, 1.9, 210.7685751197300, 8.9475067803600),
    ("NGC 4258", 7.58, 0.11, 0.11, 679.3, 0.4, 184.7400833333333, 47.30371944444444),
]
# Pesce Table 4 (main.tex lines 357-362): flow-model peculiar velocities (Group, 2M++, CF3, M2000)
TABLE4_PEC = {
    "UGC 3789": (-80, -54, 28, -142), "NGC 6264": (-180, 224, 17, -468), "NGC 6323": (-461, 414, -302, -396),
    "NGC 5765b": (-66, 124, 187, -460), "CGCG 074-064": (647, 297, 132, -372), "NGC 4258": (-46, 262, 254, 98),
}


def gal_to_icrs_unit(l_deg, b_deg):
    # IAU galactic->ICRS rotation (Hipparcos/ICRS matrix, transpose of A_G)
    AG = np.array([[-0.0548755604162154, -0.8734370902348850, -0.4838350155487132],
                   [+0.4941094278755837, -0.4448296299600112, +0.7469822444972189],
                   [-0.8676661490190047, -0.1980763734312015, +0.4559837761750669]])
    l, b = np.radians(l_deg), np.radians(b_deg)
    g = np.array([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)])
    return AG.T @ g


def radec_unit(ra, dec):
    a, d = np.radians(ra), np.radians(dec)
    return np.array([np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)])


CMB_APEX = gal_to_icrs_unit(264.021, 48.253)
CMB_V = 369.82
LG_APEX = gal_to_icrs_unit(93.0, -4.0)
LG_V = 316.0


def zcos(v_obs, v_pec, obs_conv="optical", comp="relativistic"):
    """Cosmological redshift from an observed velocity (km/s) and a radial peculiar velocity (km/s)."""
    b = v_obs / C
    if obs_conv == "optical":
        zo = b
    elif obs_conv == "radio":
        zo = b / (1 - b)
    elif obs_conv == "relativistic":
        zo = np.sqrt((1 + b) / (1 - b)) - 1
    else:
        raise ValueError
    bp = v_pec / C
    if comp == "relativistic":
        return (1 + zo) / np.sqrt((1 + bp) / (1 - bp)) - 1
    if comp == "pesce_eq9":  # 1 + z_obs = (1 + z_cos)(1 + v_pec/c)
        return (1 + zo) / (1 + bp) - 1
    if comp == "additive":
        return zo - bp
    raise ValueError


def mstar_r_bound(z, Dlow, omb, mb):
    Hreq = 2 * C / ((1 + z) * Dlow) * (1 - (1 + z) ** -0.5)
    return mb * ((Hreq / 100.0) ** 2 / omb - 1.0), Hreq


def evaluate(vpec, Hcons, omb_c, mb_c, obs_conv="optical", comp="relativistic", frame_shift=None, nsig_D=3, nsig_v=3):
    rows = []
    for name, D, ep, em, vv, sv, ra, dec in TABLE1:
        vobs = vv + nsig_v * sv
        if frame_shift is not None:
            vobs = vobs + frame_shift(ra, dec)
        z = zcos(vobs, vpec, obs_conv, comp)
        Dlow = D - nsig_D * em
        Dmax = DA_env(z, Hcons)
        msr, Hreq = mstar_r_bound(z, Dlow, omb_c, mb_c)
        # extra positive cz (km/s, applied to v_obs before composition) needed to flip the row
        f = lambda dv: DA_env(zcos(vobs + dv, vpec, obs_conv, comp), Hcons) - Dlow
        dv_flip = brentq(f, -5000, 50000) if f(-5000) < 0 < f(50000) else float("nan")
        rows.append(dict(galaxy=name, D_med=D, D_minus=em, D_low=Dlow, v_cmb_opt=vv, sigma_v=sv, v_used=vobs,
                         v_pec=vpec, z_cos=z, DA_max=Dmax, mismatch=bool(Dlow > Dmax),
                         margin_Mpc=Dlow - Dmax, margin_sigma_units=(Dlow - Dmax) / em if em else None,
                         H_required_kms_Mpc=Hreq, mstar_absr_max_GeV=msr, mX_absr_max_GeV=msr / MU,
                         extra_cz_to_flip_kms=dv_flip))
    return rows


def main():
    out = {}
    out["mstar_GeV"], out["mu"] = MSTAR, MU
    # (b)
    Hnom, omh_nom = hfloor(OMB, MB, MSTAR * R)
    Hamu, omh_amu = hfloor(OMB, MB_AMU, MSTAR * R)
    omb_c, mb_c = OMB - 3 * SOMB, MB * 1.01
    Hcons, omh_cons = hfloor(omb_c, mb_c, MSTAR * R)
    out["floor"] = dict(omega_b=OMB, sigma=SOMB, m_b=MB, ratio_mstar_r_over_mb=MSTAR * R / MB,
                        omega_h_min=omh_nom, Hfloor_nominal=Hnom,
                        omega_h_min_amu=omh_amu, Hfloor_amu=Hamu, Hfloor_amu_rel_change=Hamu / Hnom - 1,
                        omega_h_min_rel_change_amu=omh_amu / omh_nom - 1,
                        omega_b_cons=omb_c, m_b_cons=mb_c, omega_h_min_cons=omh_cons, Hfloor_cons=Hcons,
                        Hfloor_with_mX_instead_of_mstar=hfloor(OMB, MB, MX * R)[0])
    for e in (273.279, 273.45, 273.9, 274.0):
        m, ng, rc = mass_per_baryon_from_eta_conversion(e)
        out.setdefault("eta_conversion_mass_per_baryon_GeV", {})[str(e)] = m
    out["n_gamma_cm3"], out["rho_crit_over_h2_GeV_cm3"] = ng, rc

    # (c) envelope checks
    chk = []
    for z in (0.002, 0.01, 0.03, 0.035, 0.1, 1.0):
        num = DA_numeric(z, lambda zz: Hcons * (1 + zz) ** 1.5)
        chk.append(dict(z=z, closed=DA_env(z, Hcons), numeric_EdS=num, rel=num / DA_env(z, Hcons) - 1))
    out["envelope_closed_vs_numeric"] = chk
    rng = np.random.default_rng(7)
    worst = -np.inf
    for _ in range(3000):
        # random flat models whose matter part is >= floor and extra components have rho >= 0
        Hm = Hcons
        extra_L = rng.uniform(0, 3) * Hm**2
        extra_r = rng.uniform(0, 1e-3) * Hm**2
        w = rng.uniform(-1, 1)
        extra_w = rng.uniform(0, 2) * Hm**2
        extra_m = rng.uniform(0, 2) * Hm**2
        Hf = lambda zz: np.sqrt((Hm**2 + extra_m) * (1 + zz) ** 3 + extra_L + extra_r * (1 + zz) ** 4 + extra_w * (1 + zz) ** (3 * (1 + w)))
        z = rng.uniform(0.001, 0.05)
        worst = max(worst, DA_numeric(z, Hf) / DA_env(z, Hm) - 1)
    out["random_flat_positive_component_models_max_DA_over_env_minus1"] = worst

    # (d) main comparison
    out["nominal_vpec_m750"] = evaluate(-750.0, Hcons, omb_c, mb_c)
    out["stress_vpec_m1500"] = evaluate(-1500.0, Hcons, omb_c, mb_c)
    out["central_values_no_margins"] = evaluate(0.0, Hnom, OMB, MB, nsig_D=0, nsig_v=0)

    # convention variants (all with -750 and conservative floor)
    def shift_to_helio(ra, dec):  # v_hel = v_cmb - V cos(theta), additive (NED-style)
        return -CMB_V * float(radec_unit(ra, dec) @ CMB_APEX)

    def shift_to_LG(ra, dec):  # v_LG = v_hel + V_LG cos(theta_LG)
        n = radec_unit(ra, dec)
        return -CMB_V * float(n @ CMB_APEX) + LG_V * float(n @ LG_APEX)

    variants = {
        "pesce_eq9_composition": dict(comp="pesce_eq9"),
        "additive_composition": dict(comp="additive"),
        "table_v_as_relativistic": dict(obs_conv="relativistic"),
        "table_v_as_radio": dict(obs_conv="radio"),
        "heliocentric_frame": dict(frame_shift=shift_to_helio),
        "local_group_frame_KM96": dict(frame_shift=shift_to_LG),
        "heliocentric_frame_radio_additive": dict(frame_shift=shift_to_helio, obs_conv="radio", comp="additive"),
    }
    out["frame_shifts_kms"] = {name: dict(cmb_minus_helio=-shift_to_helio(ra, dec), lg_minus_cmb=shift_to_LG(ra, dec))
                               for name, D, ep, em, vv, sv, ra, dec in TABLE1}
    for k, kw in variants.items():
        for vp in (-750.0, -1500.0):
            out[f"variant_{k}_vpec{int(vp)}"] = evaluate(vp, Hcons, omb_c, mb_c, **kw)
    # largest-favourable combination: radio reading, LG or CMB frame (whichever larger), additive comp
    def best_frame(ra, dec):
        return max(0.0, shift_to_helio(ra, dec), shift_to_LG(ra, dec))
    for vp in (-750.0, -1500.0):
        out[f"variant_most_favourable_vpec{int(vp)}"] = evaluate(vp, Hcons, omb_c, mb_c, obs_conv="radio", comp="additive", frame_shift=best_frame)

    # premise escapes: negative rho_Lambda or closed curvature (flat-FRW / rho>=0 premises)
    def max_mismatch(HL2=0.0, K=0.0, vpec=-750.0):
        # H^2(z) = Hcons^2 (1+z)^3 + HL2 + K (1+z)^2 ; returns max over rows of D_low - D_A (FRW D_A incl. curvature)
        H0sq = Hcons**2 + HL2 + K
        if H0sq <= 0:
            return np.inf
        H0 = np.sqrt(H0sq)
        Om, OL, Ok = Hcons**2 / H0sq, HL2 / H0sq, K / H0sq
        worst = -np.inf
        for name, D, ep, em, vv, sv, ra, dec in TABLE1:
            z = zcos(vv + 3 * sv, vpec)
            if min(Om * (1 + zz) ** 3 + Ok * (1 + zz) ** 2 + OL for zz in np.linspace(0, z, 50)) <= 0:
                return np.inf
            worst = max(worst, (D - 3 * em) - DA_general(z, H0, Om, Ok, OL))
        return worst
    HL2_crit = brentq(lambda x: max_mismatch(HL2=x), -0.99 * Hcons**2, 0.0)
    K_crit = brentq(lambda x: max_mismatch(K=x), -0.99 * Hcons**2, 0.0)
    out["premise_escape_negative_Lambda"] = dict(HL2_crit=HL2_crit, Omega_L_crit=HL2_crit / (Hcons**2 + HL2_crit),
                                                 H0_at_crit=np.sqrt(Hcons**2 + HL2_crit),
                                                 rho_L_over_rho_floor_matter_today=HL2_crit / Hcons**2)
    out["premise_escape_closed_curvature"] = dict(K_crit=K_crit, Omega_k_crit=K_crit / (Hcons**2 + K_crit),
                                                  H0_at_crit=np.sqrt(Hcons**2 + K_crit))
    # opposite corners must also be satisfiable: for the negative-Lambda escape, do the far corners still admit it?
    def far_corner_ok(HL2):
        H0sq = Hcons**2 + HL2
        H0 = np.sqrt(H0sq)
        Om, OL = Hcons**2 / H0sq, HL2 / H0sq
        res = []
        for name, D, ep, em, vv, sv, ra, dec in TABLE1:
            z = zcos(vv - 3 * sv, +750.0)
            res.append((name, D + 3 * ep, DA_general(z, H0, Om, 0.0, OL), bool(D + 3 * ep >= DA_general(z, H0, Om, 0.0, OL))))
        return res
    out["negative_Lambda_far_corner_check_at_crit"] = far_corner_ok(HL2_crit)
    out["negative_Lambda_window"] = None
    # scan HL2 for which every galaxy's 3-sigma box intersects the model curve
    ok = []
    for x in np.linspace(-0.99 * Hcons**2, 0, 4000):
        near = max_mismatch(HL2=x) <= 0
        far = all(r[3] for r in far_corner_ok(x))
        if near and far:
            ok.append(x)
    if ok:
        out["negative_Lambda_window"] = dict(HL2_min=min(ok), HL2_max=max(ok),
                                             H0_range=[float(np.sqrt(Hcons**2 + min(ok))), float(np.sqrt(Hcons**2 + max(ok)))],
                                             OmegaL_range=[min(ok) / (Hcons**2 + min(ok)), max(ok) / (Hcons**2 + max(ok))],
                                             q0_range=[0.5 * (Hcons**2) / (Hcons**2 + x) - x / (Hcons**2 + x) for x in (min(ok), max(ok))])
    json.dump(out, sys.stdout, indent=1, default=float)
    print()


if __name__ == "__main__":
    main()
