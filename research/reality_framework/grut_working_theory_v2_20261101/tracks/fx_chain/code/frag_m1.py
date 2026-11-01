"""fx_chain item 2b (Claude): slow passage of the relaxed (circular) condensate through the M1 spinodal [E].

The linear cold-fluid dispersion of the homogeneous rotating condensate + phi in the M1 potential is copied from
tracks/m1_vacuum/code/coldfluid_m1.py (functions Uphi_p, Uphi_pp, S_of_F, cubic_coeffs, growth; unchanged).  In the
relaxed history the condensate density falls as n = n_sp exp(-3 H t) through n_sp = 233.2 m_X^3 (spinodal), so the
growth rate of the fastest mode rises from zero.  Amplification A(t) = int Gamma_peak dt reaches the nonlinear level
A_* (20-60) at depth delta_* = (n_sp - n)/n_sp; the fastest wavelength there sets the first lumps,
Q_frag = n (2 pi / k_peak)^3 (one cube per wavelength; (pi/k)^3 is quoted as the lower end).  T = 0 dispersion;
thermal damping, thermal masses, noise level and nucleation in the metastable band 233-449 m_X^3 are not included.
"""
import json, math, os, sys
import numpy as np
from scipy.optimize import brentq, minimize_scalar

a, b, c = 0.1, 3.0, 0.2
G2 = 0.01
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'results')


def Uphi_p(S): return b * (S * S - 1) * S - 3 * c * (S - 1) ** 2
def Uphi_pp(S): return b * (3 * S * S - 1) - 6 * c * (S - 1)
def S_of_F(F): return brentq(lambda s: Uphi_p(s) + F * F * s, 1e-14, 1.0, xtol=1e-15)


def cubic_coeffs(k2, F, S, Om, D):
    p1 = np.polymul([1, k2 + 2 * a * F * F], [1, k2])
    p = np.polymul(p1, [1, k2 + D])
    p = np.polyadd(p, 4 * Om * Om * np.polymul([1, 0], [1, k2 + D]))
    p = np.polysub(p, 4 * F * F * S * S * np.array([1, k2]))
    return p


def growth(k, F, S, Om, D):
    z = np.roots(cubic_coeffs(k * k, F, S, Om, D))
    nu = np.sqrt(z.astype(complex))
    return float(max(np.max(nu.real), np.max((-nu).real), 0.0))


def state(F):
    S = S_of_F(F); Om = math.sqrt(S * S + a * F * F); D = Uphi_pp(S) + F * F
    return S, Om, D, 100 * Om * F * F


def peak(F):
    S, Om, D, n = state(F)
    ks = np.geomspace(1e-7, 2.0, 400)
    g = np.array([growth(k, F, S, Om, D) for k in ks])
    i = int(np.argmax(g))
    if g[i] <= 0:
        return dict(F=F, n=n, Gamma=0.0, k=0.0)
    r = minimize_scalar(lambda lk: -growth(math.exp(lk), F, S, Om, D),
                        bracket=(math.log(ks[max(i - 1, 0)]), math.log(ks[i]), math.log(ks[min(i + 1, len(ks) - 1)])), tol=1e-10)
    k = math.exp(r.x)
    return dict(F=F, n=n, S=S, Omega=Om, Gamma=growth(k, F, S, Om, D), k=k)


def main():
    # spinodal F (min Omega along the family)
    Fsp = minimize_scalar(lambda F: state(F)[1], bracket=(1.7, 1.85, 2.0), tol=1e-13).x
    nsp = state(Fsp)[3]
    rows = []
    for dl in np.geomspace(1e-7, 0.4, 36):
        F = brentq(lambda F: state(F)[3] - nsp * (1 - dl), 0.3, Fsp, xtol=1e-15)
        p = peak(F); p['delta'] = dl; rows.append(p)
    d = np.array([r['delta'] for r in rows]); G = np.array([r['Gamma'] for r in rows]); k = np.array([r['k'] for r in rows])
    small = d < 1e-4
    c1 = float(np.exp(np.mean(np.log(G[small] / d[small]))))           # Gamma_peak ~ c1 delta
    c2 = float(np.exp(np.mean(np.log(k[small] / np.sqrt(d[small])))))  # k_peak ~ c2 sqrt(delta)
    # Hubble rate at the spinodal crossing (relaxed history, e = 1e-6 band): T_vis from a_ref scaling, g* = 106.75
    HO = json.load(open(os.path.join(OUT, 'heat_onset.json')))
    hist = {h['label']: h for h in HO['history_M1']}
    chain = {r['label']: r for r in HO['chain'] if r['phi_in_equilibrium']}
    mx = 0.0924089637375362
    out = dict(spinodal=dict(F=Fsp, n=nsp), table=rows, near_spinodal_fit=dict(c1_Gamma_over_delta=c1, c2_k_over_sqrt_delta=c2), cases=[])
    for lab in ('heavy_release:R0', 'heavy_release:tot_e1e6_band_4.97', 'heavy_release:tot_e1e7_band_41.7'):
        a_sp = chain[lab]['onset_untilted']['a_over_aref'] * hist[lab]['crossings']['spinodal_233']['a_over_aon']
        Tvis = 45.66997652455368 / a_sp * (86.25 / 106.75) ** (1 / 3)      # entropy conservation g_s 86.25 -> 106.75 [E]
        H = math.sqrt(math.pi ** 2 * 106.75 / 90) * Tvis ** 2 / 2.435e18
        h = H / mx
        for A in (20.0, 40.0, 60.0):
            # A = (1/3h) int_0^delta Gamma(delta') d delta'  (numerical, with the table; linear near zero)
            from scipy.interpolate import interp1d
            lg = interp1d(np.log(d), np.log(G), fill_value='extrapolate')
            def Aof(dl):
                xs = np.geomspace(1e-14, dl, 400)
                return np.trapezoid(np.exp(lg(np.log(xs))), xs) / (3 * h)
            dstar = brentq(lambda ld: Aof(math.exp(ld)) - A, math.log(1e-14), math.log(0.39))
            dstar = math.exp(dstar)
            F = brentq(lambda F: state(F)[3] - nsp * (1 - dstar), 0.3, Fsp, xtol=1e-15)
            pk = peak(F)
            Qw = pk['n'] * (2 * math.pi / pk['k']) ** 3
            out['cases'].append(dict(label=lab, a_sp_over_aref=a_sp, T_vis_GeV=Tvis, H_over_mX=h, A=A, delta_star=dstar,
                                     k_peak=pk['k'], Gamma_peak=pk['Gamma'], Gamma_over_H=pk['Gamma'] / h,
                                     Q_wavelength_cube=Qw, Q_half_wavelength_cube=Qw / 8, Q_bind=71262.76, Q_min=56709.8))
            print('%-34s Tvis=%.1f H/m=%.2e A=%g delta*=%.2e k=%.3e Q(2pi/k)^3=%.2e Q(pi/k)^3=%.2e' % (
                lab, Tvis, h, A, dstar, pk['k'], Qw, Qw / 8), flush=True)
    # sudden-onset comparison: the fastest lump deep in the band (where the eccentric/old picture put the onset)
    for n in (164.3, 100.0, 50.0, 24.0):
        F = brentq(lambda F: state(F)[3] - n, 0.05, Fsp, xtol=1e-15)
        pk = peak(F)
        out.setdefault('deep_band', []).append(dict(n=n, k_peak=pk['k'], Gamma=pk['Gamma'], Q_wavelength_cube=n * (2 * math.pi / pk['k']) ** 3))
        print('deep band n=%.1f k=%.4f Gamma=%.4f Q=%.3e' % (n, pk['k'], pk['Gamma'], n * (2 * math.pi / pk['k']) ** 3))
    print('spinodal F=%.6f n=%.4f; fit Gamma ~ %.4f delta, k ~ %.4f sqrt(delta)' % (Fsp, nsp, c1, c2))
    json.dump(out, open(os.path.join(OUT, 'frag_m1.json'), 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
