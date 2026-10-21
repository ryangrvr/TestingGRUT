"""Item 2 (part): number-changing and conversion processes allowed at tree level in the +v vacuum,
their amplitudes, thermal averages and equilibrium reaction rates R_eq(T) (mu = 0, MB statistics).

Writes results/rates.json. Deterministic (fixed seeds)."""
import math
import time

import numpy as np
from scipy.special import kve

import amplitudes as A
from common import (GX, LAM_PHI, M_D, M_H, M_X, V, dump, numerical_tensors, vertex_tensors)

T0 = time.time()
checks = {}

# ---------------------------------------------------------------- checks of the amplitude engine
H2, N3, N4 = numerical_tensors(step=1.0)
A3, A4 = vertex_tensors()
checks["vertex_tensor_maxdiff_cubic"] = float(np.abs(N3 - A3).max())
checks["vertex_tensor_maxdiff_quartic"] = float(np.abs(N4 - A4).max())
checks["mass2_from_hessian"] = [float(x) for x in np.diag(H2)]

g2, lam = GX ** 2, LAM_PHI
rs = 400.0
s = rs ** 2
pi_ = math.sqrt(s / 4 - M_H ** 2)
p1 = np.array([rs / 2, 0, 0, pi_])
p2 = np.array([rs / 2, 0, 0, -pi_])
rel = []
for c in np.linspace(-0.95, 0.95, 9):
    k1, k2, pf = A.two_body_final(rs, M_X, M_X, np.array([math.sqrt(1 - c * c), 0, c]))
    M = A.amp_process(("h", "h"), ("X", "Xb"), [p1, p2], [k1, k2])
    t = A.mdot(p1 - k1, p1 - k1)
    u = A.mdot(p1 - k2, p1 - k2)
    Man = -(2 * g2 + 12 * lam * g2 * V ** 2 / (s - M_H ** 2) + 4 * g2 ** 2 * V ** 2 * (1 / (t - M_X ** 2) + 1 / (u - M_X ** 2)))
    rel.append(abs(M - Man) / abs(Man))
checks["hh_to_XXb_vs_analytic_maxrel"] = float(max(rel))
checks["charge_violating_amplitude_XX_to_hh"] = float(abs(A.amplitude(["X", "X", "h", "h"], [p1, p2, -k1, -k2])))
checks["Z2D_violating_amplitude_hD_to_XXb"] = float(abs(A.amplitude(["h", "D", "Xb", "X"], [p1, p2, -k1, -k2])))

rng = np.random.default_rng(3)
pin = [A._mb_sampler(M_X, 30.0, rng, 1)[0] for _ in range(3)]
cm, rs3 = A.boost_to_cm(pin)
k1, k2, pf = A.two_body_final(rs3, M_X, M_H, np.array([0.6, 0, 0.8]))
M1 = A.amp_process(("X", "X", "Xb"), ("X", "h"), cm, [k1, k2])


def boost(p, b):
    b = np.asarray(b)
    b2 = b @ b
    gam = 1 / math.sqrt(1 - b2)
    bp = b @ p[1:]
    return np.array([gam * (p[0] + bp), *(p[1:] + ((gam - 1) * bp / b2 + gam * p[0]) * b)])


bb = [0.3, -0.2, 0.5]
M2 = A.amp_process(("X", "X", "Xb"), ("X", "h"), [boost(p, bb) for p in cm], [boost(k1, bb), boost(k2, bb)])
M3 = A.amp_process(("X", "X", "Xb"), ("X", "h"), [cm[1], cm[0], cm[2]], [k1, k2])
checks["lorentz_invariance_rel"] = float(abs(M1 - M2) / abs(M1))
checks["identical_particle_swap_rel"] = float(abs(M1 - M3) / abs(M1))
# nullification check: pure-h threshold amplitudes of the broken lam phi^4 theory vanish (known result)
checks["hhh_to_hh_threshold_M2"] = A.sigma_v2_threshold(("h", "h", "h"), ("h", "h"))["M2"]
checks["XXXb_to_Xh_threshold_M2"] = A.sigma_v2_threshold(("X", "X", "Xb"), ("X", "h"))["M2"]

# ---------------------------------------------------------------- process list
conv, three = A.enumerate_processes()
CC = {"X": "Xb", "Xb": "X", "h": "h", "D": "D"}


def canon(proc):
    ini, fin = proc
    c = (tuple(sorted(CC[t] for t in ini)), tuple(sorted(CC[t] for t in fin)))
    return min((tuple(sorted(ini)), tuple(sorted(fin))), c)


# ---------------------------------------------------------------- 3->2: threshold and MC thermal averages
XGRID = np.array([0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 50.0, 100.0, 300.0])
mc_cache = {}
three_out = []
for proc in three:
    key = canon(proc)
    thr = A.sigma_v2_threshold(*proc)
    if key not in mc_cache:
        heavy = "D" in proc[0] or "D" in proc[1]
        ns = 1000 if heavy else 2500
        vals, errs = [], []
        for i, x in enumerate(XGRID):
            m, e = A.mc_sigma_v2(*proc, M_X / x, nsample=ns, seed=100 + i)
            vals.append(m)
            errs.append(e)
        mc_cache[key] = (vals, errs)
    vals, errs = mc_cache[key]
    three_out.append(dict(initial=list(proc[0]), final=list(proc[1]), sqrt_s_thr=thr["sqrt_s"], p_f_thr=thr["p_f"],
                          M2_thr=thr["M2"], sigma_v2_thr=thr["sigma_v2"],
                          x_grid=XGRID.tolist(), sigma_v2_thermal=vals, sigma_v2_thermal_err=errs,
                          S_i=A.sym_factor(proc[0])))
print("3->2 done", time.time() - T0)

# ---------------------------------------------------------------- 2->2 conversions: sigma(W) and GG thermal average


def sigma_table(proc, npts=140):
    ini, fin = proc
    Wmin = max(sum(A.PMASS[t] for t in ini), sum(A.PMASS[t] for t in fin))
    dW = np.logspace(-4, 4.6, npts)
    sig = np.array([A.sigma_22(ini, fin, Wmin + d, nq=40) for d in dW])
    return Wmin, dW, sig


def lam_kallen(s, a, b):
    return (s - (math.sqrt(a) + math.sqrt(b)) ** 2) * (s - (math.sqrt(a) - math.sqrt(b)) ** 2)


def log_R_eq_22(proc, Wmin, dW, sig, T):
    """log of the mu=0 MB reaction rate per volume: (1/S_i) n1 n2 <sigma v> (Gondolo-Gelmini)."""
    ini, _ = proc
    m1, m2 = A.PMASS[ini[0]], A.PMASS[ini[1]]
    u = np.linspace(0.0, 60.0, 3001)[1:]
    W = Wmin + T * u
    d = W - Wmin
    ls = np.interp(np.log(d), np.log(dW), np.log(np.maximum(sig, 1e-300)))
    sg = np.exp(ls)
    lam_ = np.array([lam_kallen(w * w, m1 * m1, m2 * m2) for w in W])
    # integrand 2 sigma lambda K1(W/T) dW, K1 scaled: K1 = kve * exp(-W/T); factor exp(-Wmin/T) pulled out
    integrand = 2 * sg * lam_ * kve(1, W / T) * np.exp(-u)
    I = np.trapezoid(integrand, W)
    pref = T / (32 * math.pi ** 4) / A.sym_factor(ini)
    return math.log(pref * I) - Wmin / T


conv_out = []
XR = np.logspace(math.log10(0.25), math.log10(3000.0), 90)
for proc in conv:
    Wmin, dW, sig = sigma_table(proc)
    logR = [log_R_eq_22(proc, Wmin, dW, sig, M_X / x) for x in XR]
    # <sigma v> per initial pair, for orientation (mu = 0, MB)
    m1 = A.PMASS[proc[0][0]]
    sv = []
    for x in (1.0, 5.0, 20.0):
        T = M_X / x
        n1 = m1 * m1 * T * kve(2, m1 / T) * math.exp(-m1 / T) / (2 * math.pi ** 2)
        sv.append(math.exp(log_R_eq_22(proc, Wmin, dW, sig, T)) * A.sym_factor(proc[0]) / n1 ** 2)
    conv_out.append(dict(initial=list(proc[0]), final=list(proc[1]), W_threshold=Wmin, dW=dW.tolist(), sigma=sig.tolist(),
                         x_grid=XR.tolist(), log_R_eq=logR, sigma_v_at_x_1_5_20_GeVm2=sv,
                         sigma_v_at_x20_cm3s=sv[2] * 1.1673e-17, S_i=A.sym_factor(proc[0])))
print("2->2 done", time.time() - T0)

# MC thermal average must approach the threshold value as x -> infinity (normalisation check)
mcv, mce = A.mc_sigma_v2(("X", "X", "Xb"), ("X", "h"), M_X / 3000.0, nsample=2000, seed=77)
checks["mc_over_threshold_XXXb_to_Xh_x3000"] = mcv / A.sigma_v2_threshold(("X", "X", "Xb"), ("X", "h"))["sigma_v2"]
checks["mc_over_threshold_XXXb_to_Xh_x3000_stderr"] = mce / A.sigma_v2_threshold(("X", "X", "Xb"), ("X", "h"))["sigma_v2"]
# Gondolo-Gelmini average vs direct Monte Carlo of <sigma v_Mol> for hh -> X Xbar at x = 5
c0 = conv_out[0]
assert c0["initial"] == ["h", "h"] and c0["final"] == ["X", "Xb"]
rngc = np.random.default_rng(7)
Tc5 = M_X / 5.0
P1 = A._mb_sampler(M_H, Tc5, rngc, 200000)
P2 = A._mb_sampler(M_H, Tc5, rngc, 200000)
dot = P1[:, 0] * P2[:, 0] - np.sum(P1[:, 1:] * P2[:, 1:], axis=1)
Wc = np.sqrt(2 * M_H ** 2 + 2 * dot)
vmol = np.sqrt(dot ** 2 - M_H ** 4) / (P1[:, 0] * P2[:, 0])
dWc, sgc = np.array(c0["dW"]), np.array(c0["sigma"])
sg = np.exp(np.interp(np.log(np.maximum(Wc - c0["W_threshold"], 1e-12)), np.log(dWc), np.log(sgc)))
nh5 = M_H ** 2 * Tc5 * kve(2, M_H / Tc5) * math.exp(-M_H / Tc5) / (2 * math.pi ** 2)
gg5 = math.exp(float(np.interp(5.0, c0["x_grid"], c0["log_R_eq"]))) * 2 / nh5 ** 2
checks["GG_over_MC_hh_to_XXb_x5"] = gg5 / float(np.mean(sg * vmol))
checks["runtime_s"] = time.time() - T0
dump("rates.json", dict(
    status="computed (tree level, +v vacuum, MSbar couplings at mu_0, no running); thermal averages with MB statistics at mu=0",
    conventions="sigma_v2 = (1/(8 E1 E2 E3)) int dPhi2 |M|^2 / S_f [GeV^-5]; R_eq = n1 n2 n3 <sigma v2> / S_i; "
                "x = m_X/T; conjugate processes share rates by C",
    checks=checks, three_to_two=three_out, conversions_2to2=conv_out))
print(checks)
