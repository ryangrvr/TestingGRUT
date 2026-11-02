"""mf_gemini: n_s, alpha_s, r for V0-dominated hybrids with the proposal's mild waterfall [C/D].

V = V0 [1 + f(sigma)] during the sigma phase (chi = 0).  Slow roll with V0 dominance (eps ~ 1e-8):
  n_s - 1 = 2 eta - 6 eps + (2/3) eta^2 + ... (second order, Stewart-Lyth, xi^2 terms kept to first order)
  alpha_s = 16 eps eta - 24 eps^2 - 2 xi^2,  r = 16 eps.
Families (model choices):
  quad  : f = +eta0 sigma^2/2 (convex)                  -> blue
  hill  : f = -|eta0| (sigma - sm)^2/2 + const, sm > sigma_c, sigma rolling down to sigma_c -> red, eta const
  log   : f = A ln(sigma/sigma_c) (CW-like)              -> n_s = 1 - 2A/(sigma_c^2 + 2 A N_sigma)
The e-folds from pivot exit to sigma_c: N_sigma = N_* - N_wf, with N_wf from waterfall_dn.json and N_*
from the standard matching (Planck 2018 X eq. 47 form, with the -(1/12) ln g_th term; fast reheating).
A_s = 2.10e-9 fixes eps_* given H_* = H0 (V0 dominance).

Run: python3 -B ns_hybrid.py -> results/ns_hybrid.json
"""
import json
import os
import numpy as np

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
MP = 2.435e18
AS = 2.100e-9
C_SL = -0.7296


def Nstar(H0, Hend_over_H0=0.7126, K=2.0, k=0.05):
    """N_* = 66.89 - ln(k/a0H0) + 1/4 ln(V*^2/(MP^4 rho_end)) + 1/12 ln(rho_th/rho_end) - 1/12 ln g_th.
    rho_th/rho_end from K: a matter-like interval with a^2 H = K (a = 1, H = Hend at the end) ->
    rho_th/rho_end = (H_th/H_end)^2 with a_th = K^{?}: we use the w = 0 interval up to a_th = K
    (H = a^{-3/2}, a^2 H = a^{1/2} = K), i.e. rho_th/rho_end = K^{-6}."""
    hub = H0 / MP
    Vs = 3 * hub**2                       # in MP^4
    rho_end = 3 * (hub * Hend_over_H0)**2
    a0H0 = 67.4 / 2.998e5                 # Mpc^-1 [R]
    lnrat = -6 * np.log(K) if K > 1 else 0.0
    return 66.894 - np.log(k / a0H0) + 0.25 * np.log(Vs**2 / rho_end) + lnrat / 12.0 - np.log(106.75) / 12.0


def ns2(eta, eps, xi2=0.0):
    return 1 + 2 * eta - 6 * eps + (2.0 / 3.0) * eta**2 + 2 * (8 * C_SL - 1) * eps * eta \
        + (2.0 / 3.0 - 2 * C_SL) * 2 * xi2 * 0.0


def main():
    with open(os.path.join(RES, 'waterfall_dn.json')) as f:
        wd = json.load(f)
    Nwf = sorted([c['N_mean_after_sigc'] for c in wd['cases'] if c.get('ok')])
    out = dict(N_wf_range_mean=[Nwf[0], Nwf[-1]])
    for H0 in [1.25e11, 1.3e11, 1.55e11]:
        hub = H0 / MP
        eps = hub**2 / (8 * np.pi**2 * AS)
        d = dict(eps_star=eps, r=16 * eps,
                 N_star_K1=Nstar(H0, K=1.0), N_star_K2=Nstar(H0, K=2.0), N_star_K3=Nstar(H0, K=3.0))
        # convex quadratic: blue
        d['quad_convex'] = {f'eta={e}': float(ns2(e, eps)) for e in [0.002, 0.005, 0.01, 0.02]}
        # hilltop: |eta| needed for n_s window (second order); alpha_s
        lo = [e for e in np.linspace(0.005, 0.03, 25001) if ns2(-e, eps) <= 0.973]
        hi = [e for e in np.linspace(0.005, 0.03, 25001) if ns2(-e, eps) >= 0.967]
        d['hill'] = dict(abs_eta_range_for_ns_0p967_0p973=[float(min(lo)), float(max(hi))],
                         alpha_s=float(16 * eps * (-0.015) - 24 * eps**2), r=float(16 * eps),
                         note='eta constant: n_s independent of N_sigma; slope at sigma_c grows as e^{|eta| N_sigma}')
        # log family: n_s(N_sigma) for sigma_c -> 0 and sigma_c^2 = 2 A N_sigma (comparable)
        logt = {}
        for Ns in [20, 25, 30, 33, 36, 40, 45, 50]:
            # sigma_c << sigma_*: eta = -1/(2N), eps = A/(4N) -> A = 4 N eps
            eta = -1.0 / (2 * Ns)
            A = 4 * Ns * eps
            s2 = 2 * A * Ns                      # V' V'''/V^2 = (A/s)(2A/s^3) = 2A^2/s^4
            xi2 = 2 * A**2 / s2**2
            ns = ns2(eta, A**2 / (2 * s2))
            al = 16 * (A**2 / (2 * s2)) * eta - 24 * (A**2 / (2 * s2))**2 - 2 * xi2
            logt[f'N_sigma={Ns}'] = dict(n_s=float(ns), alpha_s=float(al), kappa=float(np.sqrt(8 * np.pi**2 * A)),
                                          sigma_star=float(np.sqrt(s2)))
        d['log_sigc_small'] = logt
        # N_sigma window for n_s in [0.967, 0.973], log family
        Ngrid = np.linspace(15, 70, 5501)
        nsl = np.array([ns2(-1 / (2 * n), eps) for n in Ngrid])
        win = Ngrid[(nsl >= 0.967) & (nsl <= 0.973)]
        d['log_N_sigma_window'] = [float(win.min()), float(win.max())]
        # with this history
        for K in [1.0, 2.0, 3.0]:
            Nst = Nstar(H0, K=K)
            d[f'log_with_waterfall_K{K:g}'] = dict(
                N_star=float(Nst), N_sigma_range=[float(Nst - Nwf[-1]), float(Nst - Nwf[0])],
                n_s_range=[float(1 - 1 / (Nst - Nwf[-1])), float(1 - 1 / (Nst - Nwf[0]))],
                n_s_at_Nwf_17p4=float(1 - 1 / (Nst - 17.4)), alpha_s_at_Nwf_17p4=float(-1 / (Nst - 17.4)**2))
        out[f'H0={H0:.3g}'] = d
    # retrospective pulls
    out['pulls_ns_0p970'] = dict(Planck2018=float((0.970 - 0.9649) / 0.0042), ACT_DR6_PACTLB=float((0.970 - 0.9743) / 0.0034))
    out['pulls_ns_window'] = {f'{x}': dict(Planck=float((x - 0.9649) / 0.0042), ACT=float((x - 0.9743) / 0.0034))
                              for x in [0.967, 0.973]}
    with open(os.path.join(RES, 'ns_hybrid.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
