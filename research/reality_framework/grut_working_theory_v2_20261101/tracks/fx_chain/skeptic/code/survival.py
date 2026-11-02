"""fx_chain_skeptic (Claude, own code, written fresh): can the condensate / objects survive the cooling of the hidden
plasma at the entropy per charge of the relaxed history?

Ideal Bose gases with vacuum masses (X and Xbar: m = 1, s: m = sqrt 6; m_X = 1), chemical equilibrium
mu_Xbar = -mu_X, mu_s = 0 (the network's fast-s + ss <-> X Xbar closure).  If objects (or a BEC) hold the gas at
mu_X = omega, they carry charge but no entropy, so the gas must carry all the entropy:
    s/n_Q = s_gas(T)/n_Q <= s_gas(T)/q_gas(T) =: sigma(T)
with equality when the objects are gone.  Objects exist at T only if s/n_Q < sigma(T); they survive the whole cooling
only if s/n_Q < min_T sigma(T).  s/n_Q is conserved in equilibrium (no portals).  Compared with the network's own
trajectories (condensed fraction) and with s/n_Q of the pre-onset Hartree plasma (hartree_indep.json).  [C]"""
import json, math, sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar, brentq

def bose(m, mu, T):
    def f(p):
        E = math.sqrt(p * p + m * m); x = (E - mu) / T
        return 0.0 if x > 700 else 1.0 / math.expm1(x)
    pmax = 60 * T + 10 * m
    pts = [0, min(0.3 * math.sqrt(m * T + T * T), pmax), pmax]
    def I(g):
        return sum(quad(lambda p: p * p * g(p) * f(p), a, b, limit=200, epsabs=0, epsrel=1e-11)[0] for a, b in zip(pts[:-1], pts[1:])) / (2 * math.pi ** 2)
    n = I(lambda p: 1.0); rho = I(lambda p: math.sqrt(p * p + m * m)); P = I(lambda p: p * p / (3 * math.sqrt(p * p + m * m)))
    return n, rho, P, (rho + P - mu * n) / T

RS = math.sqrt(6.0)

def sigma(T, om):
    nX, rX, PX, sX = bose(1.0, om, T); nB, rB, PB, sB = bose(1.0, -om, T); nS, rS, PS, sS = bose(RS, 0.0, T)
    return (sX + sB + sS) / (nX - nB)

def main():
    cases = dict(BEC_mu_m=1 - 1e-6, Q_bind_71263=0.9401708054479047, Q_2e5=0.8817204919354501, Q_3e7=0.79378662057509,
                 bulk_Q_inf=0.775865)
    HI = json.load(open(sys.argv[1]))
    sn = {('%.4f' % r['R']): r['S_per_net'] for r in HI['rows'] if r['phi']}
    out = dict(s_per_net_pre_onset=sn, cases={})
    Ts = np.geomspace(0.01, 20, 120)
    for name, om in cases.items():
        sg = [sigma(T, om) for T in Ts]
        i = int(np.argmin(sg))
        r = minimize_scalar(lambda lt: sigma(math.exp(lt), om), bounds=(math.log(Ts[max(i - 1, 0)]), math.log(Ts[min(i + 1, len(Ts) - 1)])),
                            method='bounded', options=dict(xatol=1e-8))
        Tmin = math.exp(r.x); smin = sigma(Tmin, om)
        edge = i == 0                                   # monotone on the grid: the minimum is the T -> 0 limit
        row = dict(omega=om, T_at_min=Tmin, sigma_min=smin, min_at_low_T_edge=edge, windows={})
        for Rk, s in sn.items():
            if s < smin:
                row['windows'][Rk] = 'survives (s/n_Q = %.3f < %.3f)' % (s, smin)
            else:  # temperatures between which the objects are absent
                lo = brentq(lambda lt: sigma(math.exp(lt), om) - s, math.log(1e-3), math.log(Tmin)) if not edge else math.log(1e-3)
                hi = brentq(lambda lt: sigma(math.exp(lt), om) - s, math.log(Tmin), math.log(50.0))
                row['windows'][Rk] = dict(absent_from_T=math.exp(hi), absent_to_T=math.exp(lo), s_per_net=s)
        out['cases'][name] = row
        print(name, json.dumps(row), flush=True)
    json.dump(out, open(sys.argv[2], 'w'), indent=1)

if __name__ == '__main__':
    main()
