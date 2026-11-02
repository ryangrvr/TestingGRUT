"""mf_hybrid common module (Claude, multi-field inflaton pass; GitHub issue #13 context).

Model (every term is a model choice, none follows from closure):
  V = (lam/4)(chi^2 - v^2)^2 + (1/2) g^2 sigma^2 chi^2 - c chi                      (waterfall chi, real; tadpole bias c)
      + V_SR(sigma) = c1 (sigma - sigma_c) - (1/2) m2 (sigma - sigma_c)^2            (tree clock potential, red tilt)
      + (1/2) g_s^2 chi^2 s^2   [mechanism a]   or   (1/2)(m_s0^2 + xi_s R) s^2  [mechanism b]   (spectator holder s)
      + lam_chiH (chi^2 - v^2)|H|^2 (lam_chiH < 0) - (y/2)(chi - v + M1/y) Nbar N     (SM side; N -> L H)
      + one-loop Coleman-Weinberg (MS-bar) of every field-dependent mass.
Units: H_ref = 1 (the inflationary Hubble rate at the waterfall), mP = M_P/H_ref (reduced Planck mass).
"""
import math
import numpy as np

MP_GEV = 2.435e18           # reduced Planck mass
AS = 2.100e-9               # Planck 2018 scalar amplitude (target)
NS_PLANCK, NS_PLANCK_ERR = 0.9649, 0.0042      # arXiv:1807.06211
NS_ACT, NS_ACT_ERR = 0.9743, 0.0034            # ACT DR6 P-ACT-LB, arXiv:2503.14452
R_BOUND = 0.036                                 # BICEP/Keck 2021, arXiv:2110.00483
KSTAR_GEV = 0.05 * 1.9733e-14 / 3.0857e24       # 0.05/Mpc in GeV
T0_GEV = 2.7255 * 8.617333e-14
GS0 = 3.931                                     # same convention as the inf_chain skeptic
GSTAR = 106.75
OMEGA_YL = 22.25                                # M1r Y_L frequency (H_i units)


def derived(H_GeV, m_chi, m_sig, y_N=0.0, lam_chiH=0.0):
    """Waterfall parameters from (H, m_chi, m_sigma) in H units, with lam v^4/4 = V0 = 3 mP^2 (V0-dominated).

    m_chi^2 = 2 lam v^2 (chi mass at the vacuum), m_sig = g v (sigma mass at the vacuum), sigma_c = sqrt(lam) v / g.
    """
    mP = MP_GEV / H_GeV
    V0 = 3.0 * mP ** 2
    lam = m_chi ** 4 / (48.0 * mP ** 2)
    v = m_chi / math.sqrt(2 * lam)
    g = m_sig / v
    sc = math.sqrt(lam) * v / g
    X = g * g * sc                                   # slope of M_chi^2(sigma) at sigma_c is 2X
    return dict(H_GeV=H_GeV, mP=mP, V0=V0, lam=lam, v=v, g=g, sigma_c=sc, X=X, m_chi=m_chi, m_sig=m_sig,
                y_N=y_N, lam_chiH=lam_chiH, v_MP=v / mP, sigma_c_MP=sc / mP)


def cw(M2, mu2, ndof=1.0, fermion=False):
    """MS-bar one-loop term n M^4 (ln|M^2|/mu^2 - 3/2)/(64 pi^2); fermions with a minus sign (ndof = 2 per Majorana)."""
    M2 = np.asarray(M2, dtype=float)
    a = np.abs(M2)
    with np.errstate(divide='ignore', invalid='ignore'):
        t = np.where(a > 0, M2 ** 2 * (np.log(a / mu2) - 1.5), 0.0)
    return (-1.0 if fermion else 1.0) * ndof * t / (64 * math.pi ** 2)


class Inflation:
    """Effective single-field potential along the inflationary valley chi = chi_min(sigma) ~ 0, with the chi loop."""

    def __init__(self, p, c1, m2, c_bias=0.0, mu2=None, loop=True):
        self.p, self.c1, self.m2, self.cb, self.loop = p, c1, m2, c_bias, loop
        self.mu2 = mu2 if mu2 is not None else p['g'] ** 2 * (2 * p['sigma_c']) ** 2

    def M2(self, s):
        p = self.p
        return p['g'] ** 2 * s * s - p['lam'] * p['v'] ** 2

    def V(self, s):
        p = self.p
        x = s - p['sigma_c']
        out = p['V0'] + self.c1 * x - 0.5 * self.m2 * x * x
        if self.loop:
            out = out + cw(self.M2(s), self.mu2)
        if self.cb:
            out = out - 0.5 * self.cb ** 2 / np.maximum(self.M2(s), 1e-300)
        return out

    def _cwd(self, s):
        """analytic chi-loop derivatives (first, second, third)."""
        p = self.p
        M2 = self.M2(s)
        a1 = 2 * p['g'] ** 2 * s
        a2 = 2 * p['g'] ** 2
        if M2 == 0:
            return 0.0, 0.0, 0.0
        L = math.log(abs(M2) / self.mu2)
        k = 32 * math.pi ** 2
        return (M2 * a1 * (L - 1.0) / k, (a1 * a1 * L + M2 * a2 * (L - 1.0)) / k, (3 * a1 * a2 * L + a1 ** 3 / M2) / k)

    def derivs(self, s):
        p = self.p
        s = float(s)
        x = s - p['sigma_c']
        d1, d2, d3 = self.c1 - self.m2 * x, -self.m2, 0.0
        if self.loop:
            l1, l2, l3 = self._cwd(s)
            d1, d2, d3 = d1 + l1, d2 + l2, d3 + l3
        if self.cb:
            M2 = self.M2(s)
            a1, a2 = 2 * p['g'] ** 2 * s, 2 * p['g'] ** 2
            c2 = self.cb ** 2
            d1 += 0.5 * c2 * a1 / M2 ** 2
            d2 += 0.5 * c2 * (a2 / M2 ** 2 - 2 * a1 * a1 / M2 ** 3)
            d3 += 0.5 * c2 * (-6 * a1 * a2 / M2 ** 3 + 6 * a1 ** 3 / M2 ** 4)
        return d1, d2, d3

    def dV(self, s, h=None):
        return self.derivs(s)[0]

    def d2V(self, s, h=None):
        return self.derivs(s)[1]

    def parts(self, s):
        """Tree and loop pieces of V' and V'' (diagnostics of the tuning)."""
        p = self.p
        x = s - p['sigma_c']
        M2 = self.M2(s)
        dM2 = 2 * p['g'] ** 2 * s
        L = math.log(abs(M2) / self.mu2)
        cw1 = M2 * dM2 * (L - 1.0) / (32 * math.pi ** 2)
        cw2 = (dM2 ** 2 * (L - 1.0) + M2 * dM2 ** 2 / M2 + M2 * 2 * p['g'] ** 2 * (L - 1.0)) / (32 * math.pi ** 2)
        return dict(tree_dV=self.c1 - self.m2 * x, loop_dV=cw1, tree_d2V=-self.m2, loop_d2V=cw2)


def Nstar_from_K(Hstar_GeV, Hi_GeV, K, g=GSTAR):
    """e-folds between pivot exit and the slice (a=1, H=H_i), for radiation rho a^4 = 3 K^2 M_P^2 H_i^2 afterwards."""
    Ta = (90.0 * K * K * MP_GEV ** 2 * Hi_GeV ** 2 / (math.pi ** 2 * g)) ** 0.25    # T a in GeV (a = 1 at slice)
    a0 = (g / GS0) ** (1.0 / 3.0) * Ta / T0_GEV
    return math.log(Hstar_GeV / (KSTAR_GEV * a0))


def T_of_rho(rho_GeV4, g=GSTAR):
    return (30.0 * rho_GeV4 / (math.pi ** 2 * g)) ** 0.25
