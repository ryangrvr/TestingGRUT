"""v_minority: sensitivity checks (verifier's own code; uses only vm_cosmo, vm_fp, vm_collapse of this directory).

1. Where the relaxed N*=50 relic comes from at the threshold H: fraction from domains with wall/radiation ratio
   sigma/(MP^2 H_re) >= 0.1 and >= 0.01 at re-entry (NG without gravity least reliable there).
2. Threshold if only domains with sigma/(MP^2 H_re) < 0.1 (or < 0.01) are counted as hidden-quanta relic.
3. Thresholds when the relic is evaluated with dp/dn = Lambda_1/2 (quasi-stationary) for comparison.
4. r(tau = r0/2) for the RD NG sphere ("frozen while super-horizon").
5. Regime-A boundary for the MD-reheating (T_RH = 10 MeV) and corrected kination N_0 extremes.
6. N_0 max for kination ending at 1 MeV (BBN-edge variant), with the g_* entropy term.
"""
import json, math, os, sys
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vm_cosmo as C
import vm_fp as F

OUT = os.path.dirname(os.path.abspath(__file__))
target = 0.12 * C.rho_c_over_h2 / (2 * math.pi ** 2 / 45 * 3.909 * C.T0 ** 3)
R = {}
Td = C.T_dom_iso()

def T_at_ratio(x):
    """T where sigma/(MP^2 H) = x."""
    return math.exp(brentq(lambda lT: math.log(C.sigma / (C.MP ** 2 * C.H_of_T(math.exp(lT))) / x),
                           math.log(1e-7), math.log(1e2)))

T01, T001 = T_at_ratio(0.1), T_at_ratio(0.01)
R['T_ratio_0.1_keV'] = T01 * 1e6
R['T_ratio_0.01_keV'] = T001 * 1e6

H = 362.35
G, Rr, rel, fr, ng = F.rates_for(H)
tot = F.relic_rho_over_s(H, 50, rel)[0]
c01 = F.relic_rho_over_s(H, 50, rel, Tcut=T01)[0]
c001 = F.relic_rho_over_s(H, 50, rel, Tcut=T001)[0]
R['at_362.35'] = dict(total_over_bound=tot / target, frac_from_ratio_ge_0p1=1 - c01 / tot,
                      frac_from_ratio_ge_0p01=1 - c001 / tot)

def thr(Tcut, Nstar=50, start='rel', lo=340, hi=700):
    def g(h):
        G, Rr, rel, fr, ng = F.rates_for(h)
        fn = rel if start == 'rel' else fr
        return math.log(F.relic_rho_over_s(h, Nstar, fn, Tcut=Tcut)[0] / target)
    return brentq(g, lo, hi, xtol=0.1)

R['threshold_relaxed_N50_f1_cut_ratio0.1'] = thr(T01)
R['threshold_relaxed_N50_f1_cut_ratio0.01'] = thr(T001)
R['threshold_fresh_N50_f1_cut_ratio0.1'] = thr(T01, start='fr', lo=480, hi=800)

# 4. r(tau=r0/2)
def rhs(t, z):
    r, u = z
    return [u / math.sqrt(1 + u * u), -2 * math.sqrt(1 + u * u) / r - 3 * (1 / t) * u]
s = solve_ivp(rhs, [1e-4, 0.5], [1.0, 0.0], method='DOP853', rtol=1e-12, atol=1e-14)
R['RD_r_at_tau_half_r0'] = s.y[0, -1]

# 5./6. regime A at extreme N_0
def NA(h, N0):
    L1, B = F.lambda1_mfpt(h)
    return (L1 / 2) * math.exp(3 * N0) / 3
gcorr = lambda Tx: math.log(C.gs(1e10) / C.gs(Tx)) / 12
def N0_kin(h, Tx):
    return C.N0_instant(h) + math.log(3 * C.MP ** 2 * h ** 2 / C.rho_rad(Tx)) / 12 + gcorr(Tx)
def N0_md(h, Trh):
    return C.N0_instant(h) - math.log(3 * C.MP ** 2 * h ** 2 / C.rho_rad(Trh)) / 12 + gcorr(Trh)
R['regimeA_MD_reheat_10MeV'] = brentq(lambda h: math.log(NA(h, N0_md(h, 1e-2)) / 1e-3), 120, 300, xtol=0.01)
R['N0_reach_corrected'] = {f'{h}': dict(kin_4MeV=N0_kin(h, 4e-3), kin_1MeV=N0_kin(h, 1e-3), md_10MeV=N0_md(h, 1e-2),
                                        g_entropy_term_4MeV=gcorr(4e-3))
                           for h in [150.0, 250.0, 1000.0, 1600.0, 2000.0]}
with open(os.path.join(OUT, 'results_sensitivity.json'), 'w') as f:
    json.dump(R, f, indent=1)
print(json.dumps(R, indent=1))
