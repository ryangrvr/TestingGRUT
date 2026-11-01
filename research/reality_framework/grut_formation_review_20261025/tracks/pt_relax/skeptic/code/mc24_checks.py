"""pt_relax skeptic (Claude, own code): extra validation of mc24.py (thermal 2<->4 collision-integral Monte Carlo).

(a) Constant |M|^2 = 1, massless Boltzmann 2 -> 4: the MC rate must equal the closed form
      R = (1/2)(1/4!) Int dPi1 dPi2 e^{-E1} e^{-E2} Phi_4(s),   Phi_4(s) = (2 pi)^{-8} (pi/2)^3 s^2 / (3! 2!),
      Int dPi1 dPi2 e^{-E1-E2} s^2 = (16/3) [Int dPi e^{-E} E^2]^2 = (16/3) [3/(2 pi^2)]^2 = 12/pi^4   (T = 1).
(b) Massive 3-body RAMBO weight: the mean of the massive/massless weight over RAMBO points must equal
      Phi_3^{massive}(s)/Phi_3^{massless}(s), with Phi_3^{massive} from the standard Dalitz integral
      Phi_3 = (1/(2 pi)^5) (pi^2/4) (1/s) Int dm12^2 dm23^2 / ... -> numerically: (1/(2pi)^5)(pi^2/(4 s)) Int dm12^2 dm23^2 over the Dalitz region
      hmm: Phi_3 = (2 pi)^{-5} (pi^2 / 4) (1/s) Int dm_12^2 dm_23^2 ... see below (lambda-function form).
(c) Convergence of the Bose c_eff at m/T = 0.16 and 0.36 with 4x the events and a different seed.
Writes results/mc24_checks.json.
"""
import math, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mc24 as M

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'mc24_checks.json')
M.rng = np.random.default_rng(777)


def const_m2_massless(nev=2_000_000, batch=250_000):
    draw, n, chi = M.make_sampler(0.0, 'boltzmann')
    acc, acc2, ntot = 0.0, 0.0, 0
    while ntot < nev:
        p1, p2 = draw(batch), draw(batch); P = p1 + p2; s = M.mdot(P, P)
        ph = M.phi_massless(4, s)                       # |M|^2 = 1, no Bose factors, massless final state
        val = (1.0 / 24) * ph / (4 * p1[:, 0] * p2[:, 0])
        acc += val.sum(); acc2 += (val * val).sum(); ntot += batch
    mean = acc / ntot; err = math.sqrt(max(acc2 / ntot - mean * mean, 0) / ntot)
    R = 0.5 * n * n * mean
    # analytic: n = 1/pi^2 ; (1/2)(1/24) (2pi)^-8 (pi/2)^3 /12 * 12/pi^4
    Ran = 0.5 * (1 / 24) * (2 * math.pi) ** -8 * (math.pi / 2) ** 3 / 12 * 12 / math.pi ** 4
    return dict(mc=R, mc_err=0.5 * n * n * err, analytic=Ran, ratio=R / Ran)


def three_body_massive(m=1.0, rs_list=(3.5, 5.0, 10.0), k=400_000):
    """mean massive weight from RAMBO vs Phi3(massive)/Phi3(massless) from the Dalitz integral."""
    out = []
    for rs in rs_list:
        s = rs * rs
        q = M.rambo(np.full(k, rs), 3, k); _, w = M.massive(q, np.full(k, rs), m)
        # Dalitz: Phi_3 = (2pi)^-5 (pi^2/4) ... use the standard form  dPhi3 = (1/(2pi)^5) (1/16 s) dm12^2 dm23^2 (pi^2) ... do it by direct 1D integral:
        # Phi_3(s) = (1/(2pi)^5) (pi^2/(4 s)) Int_{(2m)^2}^{(rs-m)^2} dm12^2  [m23max^2 - m23min^2]
        def lam(a, b, c): return a * a + b * b + c * c - 2 * a * b - 2 * b * c - 2 * c * a
        m2 = m * m
        x = np.linspace((2 * m) ** 2, (rs - m) ** 2, 200001)[1:-1]
        E2s = (x - m2 + m2) / (2 * np.sqrt(x)); E3s = (s - x - m2) / (2 * np.sqrt(x))
        p2s = np.sqrt(np.maximum(E2s ** 2 - m2, 0)); p3s = np.sqrt(np.maximum(E3s ** 2 - m2, 0))
        width = 4 * p2s * p3s                                      # m23max^2 - m23min^2
        I = np.trapezoid(width, x)
        phi3_massive = (2 * math.pi) ** -5 * math.pi ** 2 / (4 * s) * I
        phi3_massless = M.phi_massless(3, s)
        out.append(dict(rs=rs, mc_mean_weight=float(w.mean()), mc_err=float(w.std() / math.sqrt(k)),
                        dalitz_ratio=phi3_massive / phi3_massless))
    return out


def main():
    t0 = time.time()
    res = {}
    res['const_M2_massless_2to4'] = const_m2_massless()
    print('CHECK a', res['const_M2_massless_2to4'], flush=True)
    res['three_body_massive_weight'] = three_body_massive()
    print('CHECK b', res['three_body_massive_weight'], flush=True)
    conv = []
    for m in [0.36, 0.16]:
        r = M.rate(m, 4, 'bose', nev=4_000_000, batch=200_000); conv.append(r)
        print('CONV m/T=%.2f c_eff=%.4e +- %.1e  (t=%.0fs)' % (m, r['c_eff'], r['Gamma_chem_err'], time.time() - t0), flush=True)
    res['convergence_bose_4M_events_seed777'] = conv
    json.dump(res, open(OUT, 'w'), indent=1)


if __name__ == '__main__':
    main()
