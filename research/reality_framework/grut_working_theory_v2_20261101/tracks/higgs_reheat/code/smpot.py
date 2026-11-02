"""SM Higgs effective potential with the curvature coupling xi_H R h^2/2, for track higgs_reheat (Claude, issue #13,
decisive pass 2 Nov 2026).  Own code.

Reused (copied unchanged into this directory, stated):
  * sm_rge.py      -- track `vacua`'s 3-loop SM RGEs, Higgs anomalous dimension and the Landau-gauge 2-loop lambda_eff
                      (transcribed from Buttazzo et al. 1307.3536 App. A/B), as used by cf_vacuum.
  * c1_fit_grid.json -- cf_vacuum's 11 SMDR (Martin-Robertson 1907.02500) full-precision fits at PDG 2026
                      (M_h = 125.13, alpha_s = 0.1180, M_t = 170.6 ... 174.6).
  * build_matching()  -- copied from cf_vacuum/code/c2_sm_relation.py (quadratic in M_t, linear in M_h, alpha_s).
cf_vacuum's own instability scales at M_t = 172.60 are the validation target: Lambda_I(2-loop) = 3.49e11,
h_max(2-loop) = 2.69e11 GeV (results/sm_relation.json there).

New here:
  * one-loop running of xi:  16 pi^2 beta_xi = (xi - 1/6)(12 lambda + 2 Y2 - 3 g'^2/2 - 9 g^2/2),
    Y2 = 3 yt^2 + 3 yb^2 + ytau^2  [Q: Markkanen-Nurmi-Rajantie-Stopyra 1804.02020 eq. (4.21); Herranen et al.
    1407.3141].  Solved exactly along the 3-loop trajectory: (xi - 1/6)(mu) = (xi0 - 1/6) exp(int gamma_xi dln mu).
  * the de Sitter one-loop potential of 1804.02020 eq. (5.2)-(5.3), Tables 1-2 (Landau gauge, zeta = 0): masses
        W:  mW^2 + H^2 (n = 2, d = 3/2, n' = -34/15), (6, 5/6, -34/5), mW^2 - 2H^2 (-2, 3/2, 4/15)
        Z:  mZ^2 + H^2 (1, 3/2, -17/15), (3, 5/6, -17/5), mZ^2 - 2H^2 (-1, 3/2, 2/15)
        t:  mt^2 + H^2 (-12, 3/2, 38/5)      (b and tau analogously, -12 and -4, small)
        h:  3 lam phi^2 + 12(xi - 1/6)H^2 (1, 3/2, -2/15);  chi_W: lam phi^2 + 12(xi-1/6)H^2 (2, 3/2, -4/15);
        chi_Z: (1, 3/2, -2/15)
    with phi = e^{Gamma(mu)} h, couplings at mu, 3-loop running (as Mantziris-Markkanen-Rajantie 2011.03763 do:
    "three-loop RG improvement with one-loop curvature corrections" [Q abstract]).  Field-independent pieces
    (ghosts, photon, gluons, light fermions, alpha H^4, kappa H^2, V_Lambda) are removed by subtracting V(0) at the
    SAME mu (our prescription, stated; it makes Delta V one-loop mu-consistent).
  * scale choices mu^2 = h^2 + c R (c = 1 main; c = 0, 1/12, 1/4 bracket the scheme error).
"""
import json
import os

import numpy as np
from scipy.interpolate import CubicSpline

import sm_rge as R

HERE = os.path.dirname(os.path.abspath(__file__))
K16 = 16 * np.pi ** 2
MPL_RED = 2.435e18
MH0, AS0 = 125.13, 0.1180
MT0, SMT = 172.60, 0.27
KEYS = ['g3', 'g', 'gp', 'yt', 'yb', 'ytau', 'k']
_grid = json.load(open(os.path.join(HERE, 'c1_fit_grid.json')))['fits']


def build_matching(fits, mh0, as0, mt_ref):
    """copied from cf_vacuum/code/c2_sm_relation.py (unchanged logic)."""
    base = [f for f in fits if abs(f['inputs'][1] - mh0) < 1e-9 and abs(f['inputs'][2] - as0) < 1e-9]
    mts = np.array([f['inputs'][0] for f in base])
    coef = {kk: np.polyfit(mts, [f[kk] for f in base], 2) for kk in KEYS}
    get = lambda p: [f for f in fits if np.allclose(f['inputs'], p)][0]
    dmh = {kk: (get((mt_ref, mh0 + 0.5, as0))[kk] - get((mt_ref, mh0 - 0.5, as0))[kk]) / 1.0 for kk in KEYS}
    das = {kk: (get((mt_ref, mh0, as0 + 0.0009))[kk] - get((mt_ref, mh0, as0 - 0.0009))[kk]) / 0.0018 for kk in KEYS}

    def m(Mt, Mh, als):
        return {kk: float(np.polyval(coef[kk], Mt) + dmh[kk] * (Mh - mh0) + das[kk] * (als - as0)) for kk in KEYS}
    return m


match = build_matching(_grid, MH0, AS0, 171.0)


class SM:
    """3-loop running from M_t to 1.3e19 GeV; xi running at one loop; potentials."""

    def __init__(self, Mt=MT0, Mh=MH0, als=AS0, n=1600):
        c = match(Mt, Mh, als)
        c0 = dict(g1=np.sqrt(5 / 3) * c['gp'], g2=c['g'], g3=c['g3'], yt=c['yt'], yb=c['yb'], ytau=c['ytau'],
                  lam=c['k'])
        self.Mt = Mt
        self.run = R.run(c0, Mt, mu_max=1.3e19, loops=3, n=n)
        self.lnmu = np.log(self.run['mu'])
        Y = self.run['Y']
        self.Y = Y
        # gamma_xi = (12 lam + 2Y2 - 3/2 g'^2 - 9/2 g^2)/16pi^2, g'^2 = 3/5 G1
        gxi = (12 * Y[6] + 2 * (3 * Y[3] + 3 * Y[4] + Y[5]) - 1.5 * 0.6 * Y[0] - 4.5 * Y[1]) / K16
        Ixi = np.concatenate([[0.0], np.cumsum(0.5 * (gxi[1:] + gxi[:-1]) * np.diff(self.lnmu))])
        self.sp = {nm: CubicSpline(self.lnmu, Y[i]) for i, nm in enumerate(['G1', 'G2', 'G3', 'T', 'B', 'L', 'lam', 'Gam'])}
        self.sp['Ixi'] = CubicSpline(self.lnmu, Ixi)
        self.gxi = gxi

    # ------------------------------------------------------------------ couplings
    def cpl(self, mu):
        lm = np.clip(np.log(mu), self.lnmu[0], self.lnmu[-1])
        return {k: self.sp[k](lm) for k in self.sp}

    def xi(self, mu, xi0, mu0=None):
        """xi(mu) given xi(mu0) = xi0 (mu0 = M_t by default)."""
        mu0 = self.Mt if mu0 is None else mu0
        c, c0 = self.cpl(mu), self.cpl(mu0)
        return 1 / 6 + (xi0 - 1 / 6) * np.exp(c['Ixi'] - c0['Ixi'])

    # ------------------------------------------------------------------ flat space (cf_vacuum's lambda_eff)
    def lam_eff_flat(self, h, two_loop=True):
        """Landau-gauge lambda_eff(h) at mu = h, INCLUDING e^{4 Gamma}: V = lam_eff h^4/4 (flat space)."""
        h = np.atleast_1d(h)
        out = np.empty(len(h))
        for i, hh in enumerate(h):
            c = self.cpl(hh)
            y = np.array([c['G1'], c['G2'], c['G3'], c['T'], c['B'], c['L'], c['lam'], c['Gam']])
            out[i] = R.lambda_eff(y, two_loop=two_loop) * np.exp(4 * c['Gam'])
        return out

    def flat_scales(self):
        hs = np.exp(np.linspace(np.log(1e9), np.log(1e16), 1401))
        le = self.lam_eff_flat(hs)
        V = le * hs ** 4 / 4
        j = np.where(le < 0)[0][0]
        LamI = float(np.exp(np.interp(0.0, [le[j], le[j - 1]], [np.log(hs[j]), np.log(hs[j - 1])])))
        k = int(np.argmax(np.where(hs < 1e15, V, -np.inf)))
        return dict(Lambda_I=LamI, h_max=float(hs[k]), V_max=float(V[k]), lam_eff_min=float(le.min()),
                    h_at_lam_eff_min=float(hs[np.argmin(le)]))

    # ------------------------------------------------------------------ de Sitter one loop (1804.02020)
    def V_dS(self, h, H, xi0, cR=1.0, mu0=None, loop=True, mu_fixed=None):
        """Delta V(h) = V(h) - V(0) at the same mu, de Sitter with Hubble H (R = 12 H^2), Landau gauge.
        xi0 = xi(mu0), mu0 = M_t default.  mu^2 = h^2 + cR*R unless mu_fixed."""
        h = np.atleast_1d(np.asarray(h, float))
        Rc = 12 * H ** 2
        mu = np.sqrt(h ** 2 + cR * Rc) if mu_fixed is None else np.full_like(h, mu_fixed)
        c = self.cpl(mu)
        xi = self.xi(mu, xi0, mu0)
        phi = np.exp(c['Gam']) * h
        lam = c['lam']
        g2, gp2, yt2, yb2, yl2 = c['G2'], 0.6 * c['G1'], c['T'], c['B'], c['L']
        V = 0.5 * xi * Rc * phi ** 2 + 0.25 * lam * phi ** 4
        if not loop:
            return V
        H2, H4 = H ** 2, H ** 4
        mW, mZ, mt, mb, ml = g2 * phi ** 2 / 4, (g2 + gp2) * phi ** 2 / 4, yt2 * phi ** 2 / 2, yb2 * phi ** 2 / 2, \
            yl2 * phi ** 2 / 2
        cs = 12 * (xi - 1 / 6) * H2
        terms = [(2, 1.5, -34 / 15, mW + H2), (6, 5 / 6, -34 / 5, mW + H2), (-2, 1.5, 4 / 15, mW - 2 * H2),
                 (1, 1.5, -17 / 15, mZ + H2), (3, 5 / 6, -17 / 5, mZ + H2), (-1, 1.5, 2 / 15, mZ - 2 * H2),
                 (-12, 1.5, 38 / 5, mt + H2), (-12, 1.5, 38 / 5, mb + H2), (-4, 1.5, 38 / 15, ml + H2),
                 (1, 1.5, -2 / 15, 3 * lam * phi ** 2 + cs), (3, 1.5, -6 / 15, lam * phi ** 2 + cs)]
        terms0 = [(2, 1.5, -34 / 15, H2 + 0 * h), (6, 5 / 6, -34 / 5, H2 + 0 * h), (-2, 1.5, 4 / 15, -2 * H2 + 0 * h),
                  (1, 1.5, -17 / 15, H2 + 0 * h), (3, 5 / 6, -17 / 5, H2 + 0 * h), (-1, 1.5, 2 / 15, -2 * H2 + 0 * h),
                  (-12, 1.5, 38 / 5, H2 + 0 * h), (-12, 1.5, 38 / 5, H2 + 0 * h), (-4, 1.5, 38 / 15, H2 + 0 * h),
                  (1, 1.5, -2 / 15, cs + 0 * h), (3, 1.5, -6 / 15, cs + 0 * h)]

        def cw(tt):
            s = 0.0
            for n, d, n1, M2 in tt:
                L = np.log(np.abs(M2) / mu ** 2 + 1e-300)
                s = s + n * M2 ** 2 * (L - d) + n1 * H4 * L
            return s / (64 * np.pi ** 2)
        return V + cw(terms) - cw(terms0)


def barrier(sm, H, xi0, cR=1.0, mu0=None, hmin=None, hmax_s=1e17, n=1500, loop=True):
    """Barrier of Delta V_dS: (h_bar, Delta V_bar).  Returns (None, None) if no barrier (unstable at small h)."""
    hmin = 1e-3 * H if hmin is None else hmin
    hs = np.exp(np.linspace(np.log(hmin), np.log(hmax_s), n))
    V = sm.V_dS(hs, H, xi0, cR=cR, mu0=mu0, loop=loop)
    k = int(np.argmax(V))
    if k == 0 or k == n - 1 or V[k] <= 0:
        return None, None, hs, V
    # refine
    a, b = hs[max(k - 2, 0)], hs[min(k + 2, n - 1)]
    hh = np.exp(np.linspace(np.log(a), np.log(b), 401))
    VV = sm.V_dS(hh, H, xi0, cR=cR, mu0=mu0, loop=loop)
    kk = int(np.argmax(VV))
    return float(hh[kk]), float(VV[kk]), hs, V


if __name__ == '__main__':
    for Mt in (172.33, 172.60, 172.87):
        sm = SM(Mt)
        print(Mt, sm.flat_scales())
