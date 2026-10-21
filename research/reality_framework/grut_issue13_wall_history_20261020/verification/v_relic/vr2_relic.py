"""v_relic part 2: symmetric hot-relic bookkeeping, frozen-number Omega(xi), conversion-only Boltzmann ODE,
xi_max, minimum T_SM for a relic-compatible restoration.  Independent code (numpy/scipy only).

Hidden species: X+Xbar (g=2, common mu), h (g=1), D (g=1); BE thermodynamics by Bessel series;
collision terms MB with fugacities, rates from Gondolo-Gelmini integrals of the hand amplitudes (vr_common).
No 3->2 / 2->4 number-changing (frozen total number) -- this is the 'frozen-number' bookkeeping asked for.
"""
import json, os, sys, time
import numpy as np
from scipy.special import kve, zeta
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
import vr_common as C

# ---------------- cosmological constants (own derivation) ----------------
kB = 8.617333262e-5      # eV/K
hbarc = 1.973269804e-5   # eV cm
T0 = 2.7255 * kB / hbarc  # cm^-1
gs0 = 2 + 7 / 8 * 2 * 3 * 4 / 11
s0 = 2 * np.pi**2 / 45 * gs0 * T0**3      # cm^-3
G = 6.67430e-11; Mpc = 3.0856775814913673e22; c_light = 2.99792458e8
H100 = 1e5 / Mpc                            # s^-1
rhoc_h2_kg_m3 = 3 * H100**2 / (8 * np.pi * G)
GeV_per_kg = 1 / 1.78266192e-27
rhoc_h2 = rhoc_h2_kg_m3 * GeV_per_kg * 1e-6  # GeV / cm^3
OM_FAC = s0 / rhoc_h2                        # Omega h^2 = m[GeV] * Y * OM_FAC
gstar = 106.75
SPECIES = [("X", C.mX, 2.0), ("h", C.mh, 1.0), ("D", C.mD, 1.0)]


def s_SM(T):
    return 2 * np.pi**2 / 45 * gstar * T**3


def rho_SM(T):
    return np.pi**2 / 30 * gstar * T**4


# ---------------- BE thermodynamics ----------------
def thermo(T, psi, m, g):
    """n, rho, P for a BE gas, psi=(mu-m)/T<0. Bessel series sum_k e^{k psi} ..."""
    x = m / T
    kmax = int(min(600, max(4, np.ceil(40.0 / max(-psi, 1e-6)))))
    k = np.arange(1, kmax + 1)
    e = np.exp(k * psi)
    K2 = kve(2, k * x); K1 = kve(1, k * x)
    pre = g / (2 * np.pi**2)
    n = pre * T**3 * np.sum(e * x**2 * K2 / k)
    P = pre * T**4 * np.sum(e * x**2 * K2 / k**2)
    rho = pre * T**4 * np.sum(e * (3 * x**2 * K2 / k**2 + x**3 * K1 / k))
    return n, rho, P


def n_massless(T, g):
    return g * zeta(3) * T**3 / np.pi**2


# ---------------- closed form hot relic ----------------
def closed_form():
    gsum_m = sum(g * m for _, m, g in SPECIES)
    coef = 45 * zeta(3) / (2 * np.pi**4 * gstar)   # Y per boson dof / xi^3
    res = {"Y_per_dof_over_xi3": coef,
           "Omega_coef_all_frozen_massless": coef * gsum_m * OM_FAC,
           "Omega_coef_h_to_X": coef * (2 * C.mX + 1 * C.mX + C.mD) * OM_FAC,
           "Omega_coef_noD": coef * (2 * C.mX + C.mh) * OM_FAC,
           "Omega_coef_noD_h_to_X": coef * (3 * C.mX) * OM_FAC}
    for k_ in list(res.keys()):
        if k_.startswith("Omega_coef"):
            res["xi_max_" + k_[11:]] = (0.12 / res[k_])**(1 / 3)
    # with masses at T_h = 4 mX (thermal reference, mu=0)
    Tr = 4 * C.mX
    nref = {nm: thermo(Tr, -m / Tr, m, g)[0] for nm, m, g in SPECIES}
    nml = {nm: n_massless(Tr, g) for nm, m, g in SPECIES}
    res["n_ref_over_massless"] = {k_: nref[k_] / nml[k_] for k_ in nref}
    # Omega for xi: frozen composition at reference, and with full h->X conversion
    def Om(xi, conv=False, noD=False):
        TSM = Tr / xi
        Y = {k_: nref[k_] / s_SM(TSM) for k_ in nref}
        if noD:
            Y["D"] = 0.0
        if conv:
            return (C.mX * (Y["X"] + Y["h"]) + C.mD * Y["D"]) * OM_FAC
        return (C.mX * Y["X"] + C.mh * Y["h"] + C.mD * Y["D"]) * OM_FAC
    res["Omega_frozen_ref_xi1e-3"] = Om(1e-3)
    res["Omega_fullconv_ref_xi1e-3"] = Om(1e-3, conv=True)
    res["Omega_frozen_ref_xi7.581e-4"] = Om(7.581e-4)
    res["Omega_fullconv_ref_xi7.581e-4"] = Om(7.581e-4, conv=True)
    res["xi_max_frozen_ref"] = brentq(lambda z: Om(z) - 0.12, 1e-5, 1e-2)
    res["xi_max_fullconv_ref"] = brentq(lambda z: Om(z, conv=True) - 0.12, 1e-5, 1e-2)
    res["xi_max_frozen_ref_noD"] = brentq(lambda z: Om(z, noD=True) - 0.12, 1e-5, 1e-2)
    res["xi_max_fullconv_ref_noD"] = brentq(lambda z: Om(z, conv=True, noD=True) - 0.12, 1e-5, 1e-2)
    res["Ytot_ref_per_xi3"] = sum(nref.values()) / s_SM(Tr)
    return res, nref


# ---------------- conversion rates  R~(T) = R_eq(T) exp(2 m_i/T) ----------------
CHANNELS = [  # (name, Mfun, initial species index, final species index, m_i, m_f, S_initial, S_final)
    ("hh->XX", C.M_hh_XX, 1, 0, C.mh, C.mX, 2.0, 1.0),
    ("hh->DD", C.M_hh_DD, 1, 2, C.mh, C.mD, 2.0, 0.5),
    ("DD->XX", C.M_DD_XX, 2, 0, C.mD, C.mX, 2.0, 1.0),
]


def Rtilde(ch, T):
    name, Mf, i, f, mi, mf, Si, Sf = ch
    smin = max(4 * mi**2, 4 * mf**2) * (1 + 1e-12)
    rs0 = np.sqrt(smin)
    def integrand(u):
        rs = rs0 + T * u
        s = rs * rs
        sg = C.sigma(np.array([s]), Mf, mi, mf, Sf)[0]
        return sg * (s - 4 * mi**2) * rs * kve(1, rs / T) * np.exp(-(rs - 2 * mi) / T) * 2 * rs * T
    umax = 60.0
    I = quad(integrand, 0, umax, limit=200, epsrel=1e-8)[0]
    return (1 / Si) * T / (32 * np.pi**4) * I


def build_rate_tables():
    lT = np.linspace(np.log(0.04), np.log(450.0), 160)
    tabs = []
    for ch in CHANNELS:
        vals = np.array([Rtilde(ch, np.exp(l)) for l in lT])
        tabs.append(np.log(np.maximum(vals, 1e-300)))
    return lT, tabs


# ---------------- Boltzmann ODE ----------------
class Model:
    def __init__(self, xi, lT, tabs, collide=True, noD=False):
        self.xi = xi; self.lT = lT; self.tabs = tabs; self.collide = collide; self.noD = noD
        self.Tr = 4 * C.mX
        self.TSMi = self.Tr / xi

    def Rt(self, j, T):
        return np.exp(np.interp(np.log(T), self.lT, self.tabs[j]))

    def state_funcs(self, y):
        lT = y[0]; T = np.exp(lT)
        out = []
        rho_tot = 0.0; P_tot = 0.0
        for (nm, m, g), psi in zip(SPECIES, y[1:]):
            n, r, P = thermo(T, psi, m, g)
            out.append(n); rho_tot += r; P_tot += P
        return np.array(out), rho_tot, P_tot

    def F(self, y):
        n, r, P = self.state_funcs(y)
        return np.concatenate([np.log(n), [np.log(r)]]), n, r, P

    def rhs(self, N, y):
        Fv, n, r, P = self.F(y)
        T = np.exp(y[0])
        TSM = self.TSMi * np.exp(-N)
        H = np.sqrt((rho_SM(TSM) + r) / (3 * C.MP**2))
        Cs = np.zeros(3)
        if self.collide:
            for j, ch in enumerate(CHANNELS):
                name, Mf, i, f, mi, mf, Si, Sf = ch
                rt = self.Rt(j, T)
                net = rt * (np.exp(2 * y[1 + i]) - np.exp(2 * y[1 + f] + 2 * (mf - mi) / T))
                Cs[i] -= 2 * net; Cs[f] += 2 * net
        dF = np.empty(4)
        for s_ in range(3):
            dF[s_] = -3.0 + Cs[s_] / (n[s_] * H)
        dF[3] = -3.0 * (1 + P / r)
        # Jacobian dF/dy by central differences
        J = np.zeros((4, 4))
        for a in range(4):
            hstep = 1e-6 * max(1.0, abs(y[a]))
            yp = y.copy(); ym = y.copy(); yp[a] += hstep; ym[a] -= hstep
            J[:, a] = (self.F(yp)[0] - self.F(ym)[0]) / (2 * hstep)
        dy = np.linalg.solve(J, dF)
        return dy

    def run(self, xend=2000.0):
        y0 = np.array([np.log(self.Tr)] + [-m / self.Tr for _, m, _ in SPECIES])
        if self.noD:
            y0[3] = -40.0  # D essentially unpopulated (n_D/n_D,eq ~ e^-40); hh->DD may still produce D
        Tend = C.mX / xend
        def ev(N, y):
            return y[0] - np.log(Tend)
        ev.terminal = True
        sol = solve_ivp(self.rhs, (0.0, 40.0), y0, method="LSODA", rtol=1e-9, atol=1e-11, events=ev,
                        dense_output=False, max_step=0.05)
        yf = sol.y[:, -1]; Nf = sol.t[-1]
        n, r, P = self.state_funcs(yf)
        TSM = self.TSMi * np.exp(-Nf)
        Y = n / s_SM(TSM)
        Om_s = np.array([m for _, m, _ in SPECIES]) * Y * OM_FAC
        return {"xi": self.xi, "N_end": Nf, "T_h_end": float(np.exp(yf[0])), "Y": Y.tolist(),
                "Ytot": float(Y.sum()), "Omega_species": Om_s.tolist(), "Omega": float(Om_s.sum()),
                "mean_mass": float(Om_s.sum() / Y.sum() / OM_FAC), "nfev": int(sol.nfev), "status": sol.status}


if __name__ == "__main__":
    t0 = time.time()
    out = {"constants": {"s0_cm-3": s0, "rhoc_h2_GeV_cm-3": rhoc_h2, "OM_FAC": OM_FAC, "gstar": gstar}}
    cf, nref = closed_form()
    out["closed_form"] = cf
    print(json.dumps(out, indent=1, default=float)); sys.stdout.flush()
    lT, tabs = build_rate_tables()
    out["rate_build_s"] = time.time() - t0
    runs = {}
    # conservation test: collisionless
    m0 = Model(1e-3, lT, tabs, collide=False)
    r0 = m0.run()
    Yref = np.array([nref[k] for k in ("X", "h", "D")]) / s_SM(4 * C.mX / 1e-3)
    r0["Y_rel_drift"] = (np.array(r0["Y"]) / Yref - 1).tolist()
    runs["collisionless_1e-3"] = r0
    print("collisionless", r0); sys.stdout.flush()
    for xi in (1e-4, 3.16e-4, 7.581e-4, 1e-3):
        rr = Model(xi, lT, tabs).run()
        runs["conv_%g" % xi] = rr
        print(xi, rr); sys.stdout.flush()
    # xi_max with conversions only (frozen total number)
    def f(lx):
        return Model(np.exp(lx), lT, tabs).run()["Omega"] - 0.12
    lx = brentq(f, np.log(6e-4), np.log(9e-4), xtol=1e-6)
    runs["xi_max_conv_only"] = float(np.exp(lx))
    rr = Model(np.exp(lx), lT, tabs).run(); runs["at_xi_max"] = rr
    print("xi_max", np.exp(lx), rr); sys.stdout.flush()
    # no-D variant
    def fnd(lx):
        return Model(np.exp(lx), lT, tabs, noD=True).run()["Omega"] - 0.12
    try:
        lx2 = brentq(fnd, np.log(6e-4), np.log(1.1e-3), xtol=1e-6)
        runs["xi_max_conv_only_noD"] = float(np.exp(lx2))
    except Exception as e:
        runs["xi_max_conv_only_noD"] = "failed: %s" % e
    out["boltzmann_conv_only"] = runs
    out["runtime_s"] = time.time() - t0
    os.makedirs("results", exist_ok=True)
    with open("results/vr2_relic.json", "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))
