"""v_relic part 1: spectrum, decay closure, gravitational h lifetime, hh -> X Xbar threshold amplitude,
thermal <sigma v> of conversions.  Independent code."""
import json, os
import numpy as np
from scipy.special import kve, kn
from scipy.integrate import quad
import vr_common as C

out = {}
hbar = 6.582119569e-25  # GeV s
age = 13.80e9 * 365.25 * 86400.0

# ---------- (1) spectrum and closure ----------
out["masses"] = {"mX": C.mX, "mh": C.mh, "mD": C.mD}
out["closure"] = {
    "h->XXbar_threshold": 2 * C.mX, "h->DD_threshold": 2 * C.mD, "mh": C.mh,
    "h->XXbar_open": bool(C.mh > 2 * C.mX), "h->DD_open": bool(C.mh > 2 * C.mD),
    "D->anything": "lightest Z2_D-odd state (mD=%.3f); D->X Xbar D needs %.1f" % (C.mD, 2 * C.mX + C.mD),
    "X->anything": "lightest U(1)_X-charged state (mX=%.1f < mD)" % C.mX,
}
out["fd_vertices"] = C.check_vertices_fd()

# gravitational lifetimes (orders of magnitude)
G_naive = C.mh**5 / C.MP**4
# loop h R_abcd^2 coefficient (heat-kernel a2 Riem^2/180 per real dof, d ln m^2/d phi)
n_X = 2 * 2.0                                         # 2 real dof * v dln(g^2 phi^2)/dphi = 2
n_D = 2 * C.hD**2 * C.v**2 / C.mD**2                  # 1 dof
n_h = 3.0                                             # v d ln(lam(3phi^2-v^2))/dphi = 6 lam v^2/(2 lam v^2)
c = (n_X + n_D + n_h) / (32 * np.pi**2 * 180)
G_loop = c**2 * C.mh**7 / (4 * np.pi * C.v**2 * C.MP**4)
# trace coupling through running xi: kappa = xi v vH / MP^2 relative to SM Higgs couplings
GSM141 = 8.6e-3  # GeV, SM Higgs total width near 141 GeV [R, recalled]
vH = 246.22
xi_lt = {}
for xi in (1e-3, 3.7e-5, 1e-5):
    kap = xi * C.v * vH / C.MP**2
    xi_lt[str(xi)] = {"tau_s_x1": hbar / (kap**2 * GSM141), "tau_s_x10": hbar / (10 * kap**2 * GSM141)}
# graviton-loop Higgs portal delta-lambda 3.4e-34 (quoted): lambda_p phi^2 |H|^2 -> mass mixing 2 lam_p v vH h H0? (|H|^2=(vH+H0)^2/2)
dl = 3.4e-34
mH = 125.1
mix_conv = {}
for name, dm2 in (("lam_p*phi^2*|H|^2", 2 * dl * C.v * vH), ("(lam_p/2)phi^2|H|^2", dl * C.v * vH),
                  ("(lam_p/4)phi^2 H^2 real", 0.5 * dl * C.v * vH)):
    th = dm2 / (C.mh**2 - mH**2)
    mix_conv[name] = hbar / (th**2 * GSM141)
out["h_lifetime"] = {
    "naive_mh5_over_MP4_GeV": G_naive, "tau_naive_s": hbar / G_naive,
    "loop_c": c, "loop_c_parts": [n_X, n_D, n_h], "tau_loop_s": hbar / G_loop,
    "trace_xi": xi_lt, "portal_3.4e-34_tau_s_by_convention": mix_conv,
    "age_s": age, "min_ratio_tau_over_age_inaction": min(hbar / G_naive, hbar / G_loop) / age,
}

# ---------- (3) threshold amplitude ----------
s_thr = 4 * C.mh**2
Mthr = C.M_hh_XX(s_thr, np.array([0.0, 0.5, -0.7, 1.0]))
formula = -4 * C.gX**2 * (1 - C.gX**2 / C.lam)
scan = []
for a in (0.25, 0.5, 0.9, 1.0, 1.1, 2.0):   # a = g^2/lam, keep lam fixed, mh>2mX requires a<0.5? (only kinematic for exo)
    g2 = a * C.lam
    Mn = C.M_hh_XX(4 * 2 * C.lam * C.v**2, np.array([0.3]), g2=g2, lm=C.lam)[0]
    scan.append({"a": a, "M_numeric": float(Mn), "M_formula": float(-4 * g2 * (1 - g2 / C.lam))})
# diagram decomposition at threshold
mh2 = C.mh**2
contact = -C.V_hhXX
schan = -C.V_hhh * C.V_hXX / (s_thr - mh2)
tu = -2 * C.V_hXX**2 / (-mh2)
out["threshold_amplitude"] = {
    "M_at_threshold_numeric_various_cos": Mthr.tolist(), "formula_-4g2(1-g2/lam)": formula,
    "diagrams": {"contact": contact, "s_channel_h": schan, "t_plus_u": tu, "sum": contact + schan + tu},
    "scan_a=g2/lam": scan,
}
# expansion near threshold: M ~ c1 * (s-4mh^2) + c2 * angular
eps = np.array([1e-4, 1e-3, 1e-2]) * mh2
out["threshold_amplitude"]["M_near_threshold_cos0"] = [float(C.M_hh_XX(s_thr + e, np.array([0.0]))[0]) for e in eps]

# ---------- thermal averages (Gondolo-Gelmini) ----------
def sigv_GG(Mfun, mi, mf, Sf, T):
    """<sigma v_Mol> = 2 R / n_eq^2 with n_eq MB, g=1 per initial particle (GG 1991 eq. 3.8)."""
    smin = max(4 * mi**2, 4 * mf**2) * (1 + 1e-12)
    def integrand(u):  # sqrt(s) = sqrt(smin) + T u
        rs = np.sqrt(smin) + T * u
        s = rs * rs
        sg = C.sigma(np.array([s]), Mfun, mi, mf, Sf)[0]
        # K1(rs/T) = kve * exp(-rs/T); factor out exp(-2 mi/T) consistent with K2(mi/T)^2 scaling
        return sg * (s - 4 * mi**2) * rs * kve(1, rs / T) * np.exp(-(rs - 2 * mi) / T) * 2 * rs * T
    I, err = quad(integrand, 0, 80, limit=200)
    num = T / (32 * np.pi**4) * I  # times exp(-2mi/T)
    neq = mi**2 * T * kve(2, mi / T) / (2 * np.pi**2)  # times exp(-mi/T)
    return num / neq**2


# Rate: R (reactions/volume/time, identical initial counted once) = (1/2) <sigma v> n^2 for identical initial.
# GG: <sigma v> n_eq^2 = g^2 T/(32 pi^4) int ds sigma (s-4m^2) sqrt(s) K1  -> the function above returns <sigma v>.
GeV2_to_cm3s = 1.1673e-17
tab = {}
for lab, Mf, mi, mf, Sf in (("hh->XXbar", C.M_hh_XX, C.mh, C.mX, 1.0), ("hh->DD", C.M_hh_DD, C.mh, C.mD, 0.5),
                            ("DD->XXbar", C.M_DD_XX, C.mD, C.mX, 1.0)):
    rows = []
    for x in (1, 3, 5, 10, 20, 30):
        T = C.mX / x
        sv = sigv_GG(Mf, mi, mf, Sf, T)
        rows.append({"x=mX/T": x, "T": T, "sigv_GeV-2": sv, "sigv_cm3/s": sv * GeV2_to_cm3s})
    tab[lab] = rows
# also hh->XXbar with x defined as mh/T
rows = []
for x in (1, 5, 20):
    T = C.mh / x
    rows.append({"x=mh/T": x, "sigv_GeV-2": sigv_GG(C.M_hh_XX, C.mh, C.mX, 1.0, T)})
tab["hh->XXbar_x_mh"] = rows
# s-wave value if the threshold amplitude were the generic -2g^2 contact only (for scale)
k0 = np.sqrt(C.mh**2 - C.mX**2)
tab["hypothetical_swave_sigv_if_|M|=4g2_GeV-2"] = (4 * C.gX**2)**2 * k0 / (32 * np.pi * C.mh**3)
out["thermal_sigv"] = tab

os.makedirs("results", exist_ok=True)
with open("results/vr1_spectrum_amp.json", "w") as f:
    json.dump(out, f, indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
