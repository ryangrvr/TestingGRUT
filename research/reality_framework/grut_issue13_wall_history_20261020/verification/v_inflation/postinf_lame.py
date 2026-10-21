"""Independent checks (verifier v_inflation): post-inflation zero crossings in RD and the Lame chart.

(a) Zero crossings.  phi'' + 3 H phi' + lam phi (phi^2 - v^2) = 0 with H = 1/(2t) (radiation domination,
    R = 0 so no xi term).  Units Phi = phi/v, tau = sqrt(lam) v t:  Phi'' + (3/(2 tau)) Phi' + Phi^3 - Phi = 0.
    Start: Phi = Phi_end, Phi' = 0 at tau0 = sqrt(lam) v/(2 H_inf) (H_inf = huge => tau0 -> 0; checked).
    Count sign changes of Phi until E = Phi'^2/2 + Phi^4/4 - Phi^2/2 < 0 (trapped; E decreases monotonically).
    Integrated in conformal variables chi = a Phi, a = eta (exact in RD):  chi'' + chi^3 - eta^2 chi = 0,
    which is non-stiff (amplitude of chi constant during the quartic stage).
    Analytic bounds: chi_A^2 = sqrt2 Phi_end (conformal energy), N = chi_A^2/(2 K(1/2) c) with trapping at
    amplitude c v: c = 1 -> 0.3814, c = sqrt2 -> 0.2697 (per Phi_end).

(b) Lame / Floquet chart for x'' + (kappa^2 + q cn^2(T, m = 1/2)) x = 0: monodromy over the period 2K(1/2) of cn^2,
    mu = arccosh(|tr M|/2)/(2K) (growth exponent per unit T, GKLS convention X ~ exp(mu x)).
"""
import json, math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipk, ellipj
from scipy.optimize import minimize_scalar, brentq

K12 = float(ellipk(0.5))


def crossings(Phi_end, eta0=1e-4, rtol=1e-10):
    # conformal: chi = eta * Phi, chi' = Phi + eta dPhi/deta ; with dPhi/dt = 0 at start, dPhi/deta = a dPhi/dt = 0
    y0 = [eta0 * Phi_end, Phi_end]

    def rhs(eta, s):
        chi, dchi = s
        return [dchi, -chi ** 3 + eta ** 2 * chi]

    def zero(eta, s):
        return s[0]

    def trapped(eta, s):
        chi, dchi = s
        Phi = chi / eta
        dPhi_deta = dchi / eta - chi / eta ** 2
        Phidot = dPhi_deta / eta  # d/dt = (1/a) d/deta, a = eta
        E = 0.5 * Phidot ** 2 + 0.25 * Phi ** 4 - 0.5 * Phi ** 2
        return E
    trapped.terminal = True
    trapped.direction = -1
    # scale estimate: crossings ~ 0.31 Phi_end, conformal period 4K/chi_A, chi_A ~ (2 Phi^2)^(1/4)
    chiA = (2 * Phi_end ** 2) ** 0.25
    eta_trap = chiA * 2.0 + 10.0
    sol = solve_ivp(rhs, (eta0, eta_trap * 3), y0, events=[zero, trapped], rtol=rtol, atol=1e-12,
                    method='DOP853', max_step=0.05 / max(chiA, 1.0) * 10)
    nz = len(sol.t_events[0])
    trapped_at = sol.t_events[1][0] if len(sol.t_events[1]) else None
    sign = 1 if sol.y[0, -1] > 0 else -1
    return nz, trapped_at, sign


def lame_mu(kappa2, q):
    P = 2 * K12

    def rhs(T, s):
        sn, cn, dn, ph = ellipj(T, 0.5)
        w2 = kappa2 + q * cn ** 2
        return [s[1], -w2 * s[0], s[3], -w2 * s[2]]
    sol = solve_ivp(rhs, (0, P), [1.0, 0.0, 0.0, 1.0], rtol=1e-12, atol=1e-14, method='DOP853')
    tr = sol.y[0, -1] + sol.y[3, -1]
    if abs(tr) <= 2:
        return 0.0, tr
    return math.acosh(abs(tr) / 2) / P, tr


def band_edges(q, lo, hi, n=400):
    ks = np.linspace(lo, hi, n)
    trs = np.array([abs(lame_mu(k, q)[1]) - 2 for k in ks])
    edges = []
    for i in range(n - 1):
        if trs[i] * trs[i + 1] < 0:
            edges.append(brentq(lambda k: abs(lame_mu(k, q)[1]) - 2, ks[i], ks[i + 1], xtol=1e-10))
    return edges, ks, trs


if __name__ == '__main__':
    out = {}
    # (a) crossings
    rows = []
    for Pe in [1.5, 2.0, 2.5, 2.6, 2.8, 3.0, 3.5, 4, 5, 7, 10, 15, 20, 30, 50, 70, 100, 150, 200, 300, 500, 1000]:
        nz, ta, sg = crossings(Pe)
        rows.append(dict(Phi_end=Pe, N_cross=nz, ratio=nz / Pe, final_sign=sg))
        print(Pe, nz, nz / Pe, sg, flush=True)
    out['crossings'] = rows
    # eta0 sensitivity
    out['eta0_sensitivity_Phi100'] = [crossings(100.0, eta0=e)[0] for e in (1e-5, 1e-4, 1e-3, 1e-2)]
    # fine scan for slope/mean ratio at Phi >= 30
    fine = np.linspace(30, 300, 136)
    ncs = np.array([crossings(p)[0] for p in fine])
    out['fine_scan_30_300'] = dict(mean_ratio=float(np.mean(ncs / fine)), min_ratio=float(np.min(ncs / fine)),
                                   max_ratio=float(np.max(ncs / fine)),
                                   slope_fit=float(np.polyfit(fine, ncs, 1)[0]),
                                   intercept_fit=float(np.polyfit(fine, ncs, 1)[1]))
    out['analytic_bounds'] = dict(c1=math.sqrt(2) / (2 * K12), c_sqrt2=math.sqrt(2) / (2 * K12 * math.sqrt(2)))
    # local slope at small Phi (for the sign-sensitivity 'flips' numbers)
    small = np.linspace(5, 60, 111)
    ncs_s = np.array([crossings(p)[0] for p in small])
    out['small_scan_5_60'] = dict(Phi=small.tolist(), N=ncs_s.tolist())
    # (b) Lame chart
    lame = {}
    for q, (lo, hi) in ((1.0, (0.0005, 1.0)), (3.0, (1.0, 2.0))):
        edges, ks, trs = band_edges(q, lo, hi)
        res = minimize_scalar(lambda k: -lame_mu(k, q)[0], bounds=(edges[0], edges[-1]) if len(edges) >= 2 else (lo, edges[0]),
                              method='bounded', options=dict(xatol=1e-7))
        lame['q=%g' % q] = dict(edges_found_in_scan=edges, kappa2_at_max=res.x, mu_max=-res.fun,
                                mu_at_kappa2_0p0005=lame_mu(0.0005, q)[0] if q == 1 else None)
        print(q, lame['q=%g' % q], flush=True)
    # GKLS table spot values (read in the paper, section VI/VII table)
    spot = {}
    for q, k2 in ((1, 0.1), (1, 0.2), (1, 0.228), (1, 0.3), (1, 0.4), (3, 1.55), (3, 1.6), (3, 1.615), (3, 1.65), (3, 1.7)):
        spot['q=%d,k2=%g' % (q, k2)] = lame_mu(k2, q)[0]
    lame['spot_values'] = spot
    # GKLS analytic formula (34): mu = (2/T) sqrt(2 k2 (1 - 4 k2^2)) ... with T = 4K, I(k) = int_0^{pi/2} sin^{1/2}/(1+2 k2 sin)
    # (equation read in the paper; transcribed as a cross-check, kappa^4 interpreted as (kappa^2)^2)
    from scipy.integrate import quad
    T = 4 * K12
    def mu34(k2):
        I = quad(lambda th: math.sqrt(math.sin(th)) / (1 + 2 * k2 * math.sin(th)), 0, math.pi / 2)[0]
        return (2.0 / T) * math.sqrt(2 * k2 * (1 - 4 * k2 ** 2)) * I if k2 < 0.5 else 0.0
    lame['gkls_eq34_q1'] = {('%g' % k): mu34(k) for k in (0.1, 0.2, 0.228, 0.3, 0.4)}
    out['lame'] = lame
    json.dump(out, open('results_postinf_lame.json', 'w'), indent=1)
