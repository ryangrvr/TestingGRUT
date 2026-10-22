"""Standard Model equation of state for the scenario track.

Source (credited; read level: ancillary data tables + TeX sections 3 intro, 3.3, 3.4):
  K. Saikawa and S. Shirai, "Primordial gravitational waves, precisely: The role of thermodynamics in the Standard
  Model", JCAP 05 (2018) 035, arXiv:1803.01038. Tables: standardmodel2018.dat (T, g_*rho, err, g_*s, err; arXiv
  ancillary file) and eos_w.dat (T >= 10 MeV: w = p/rho with error; authors' web page). They renormalise the pressure
  so that it vanishes at T = 0 (their section 3.3) and build it by integrating the trace anomaly
  Delta = (rho - 3p)/T^4 (their eq. for Delta(T)); lattice QCD input Borsanyi et al. 1606.07494, EW crossover
  D'Onofrio-Rummukainen and Laine-Meyer 1503.04935 (as cited by them).
The p(T = 0) = 0 normalisation is exactly the one P_q dictates: the T = 0 matter vacuum energy belongs to eps_inf.

Trace S = rho - 3p:
  T >= 10 MeV : S/rho = 1 - 3 w (SS w table; quoted error carried as 'S_err').
  T <  10 MeV : neutrinos decoupled (two temperatures), so 1 - 3w from g's is not a trace; use the ideal-gas
                e+- and mu+- contributions at the photon temperature (computed here, QED corrections neglected,
                [C/E]); photons and (massless) neutrinos are traceless.
Downloads live in the scratchpad (never in the repository); `fetch(dirname)` gets them.
"""
import math
import os
import urllib.request

import numpy as np

URLS = {
    'standardmodel2018.dat': 'https://member.ipmu.jp/satoshi.shirai/standardmodel2018.dat',
    'eos_w.dat': 'https://member.ipmu.jp/satoshi.shirai/eos_w.dat',
}
DEFAULT_DIR = ('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/qtheory/'
               'scenario_dl/ss')


def fetch(d=DEFAULT_DIR):
    os.makedirs(d, exist_ok=True)
    for fn, url in URLS.items():
        p = os.path.join(d, fn)
        if not os.path.exists(p) or os.path.getsize(p) < 10000:
            alt = os.path.join(d, 'src', 'anc', fn)
            if os.path.exists(alt):
                continue
            urllib.request.urlretrieve(url, p)
    return d


def _path(d, fn):
    for p in (os.path.join(d, fn), os.path.join(d, 'src', 'anc', fn)):
        if os.path.exists(p) and os.path.getsize(p) > 10000:
            return p
    raise FileNotFoundError(fn)


def fermion_trace_over_T4(mu, g):
    """(rho - 3p)/T^4 for an ideal Fermi gas, g dof, mass/T = mu, zero chemical potential."""
    x = np.linspace(0.0, 60.0, 6001)
    E = np.sqrt(x * x + mu * mu)
    with np.errstate(over='ignore'):
        integrand = x * x / E / (np.exp(E) + 1.0)
    return g / (2 * math.pi ** 2) * mu * mu * np.trapezoid(integrand, x)


class SMEoS:
    ME, MMU = 0.51099895e-3, 0.1056583755   # GeV [R: PDG]

    def __init__(self, d=DEFAULT_DIR):
        a = np.loadtxt(_path(d, 'standardmodel2018.dat'))
        b = np.loadtxt(_path(d, 'eos_w.dat'))
        self.T = a[:, 0]
        self.grho, self.grho_err, self.gs, self.gs_err = a[:, 1], a[:, 2], a[:, 3], a[:, 4]
        self.lnT = np.log(self.T)
        self.lnTw, self.w_tab, self.w_err = np.log(b[:, 0]), b[:, 1], b[:, 2]
        # trace fraction on the g-table grid
        Sfrac = np.empty_like(self.T)
        Serr = np.empty_like(self.T)
        hi = self.T >= 0.010
        Sfrac[hi] = 1 - 3 * np.interp(self.lnT[hi], self.lnTw, self.w_tab)
        Serr[hi] = 3 * np.interp(self.lnT[hi], self.lnTw, self.w_err)
        for i in np.where(~hi)[0]:
            Tk = self.T[i]
            tr = fermion_trace_over_T4(self.ME / Tk, 4) + fermion_trace_over_T4(self.MMU / Tk, 4)
            rho_T4 = math.pi ** 2 / 30 * self.grho[i]
            Sfrac[i] = tr / rho_T4
            Serr[i] = 0.1 * Sfrac[i]
        self.Sfrac, self.Serr = Sfrac, Serr
        # smooth N-derivative support: N decreases as lnT decreases; dlnT/dN = -1/(1 + dln g_s/(3 dlnT))
        lg = np.log(self.gs)
        dlg = np.gradient(lg, self.lnT)
        self.dlnT_dN = -1.0 / (1.0 + dlg / 3.0)

    # --- thermodynamics in GeV units
    def rho(self, T):
        return math.pi ** 2 / 30 * np.interp(np.log(T), self.lnT, self.grho) * T ** 4

    def s(self, T):
        return 2 * math.pi ** 2 / 45 * np.interp(np.log(T), self.lnT, self.gs) * T ** 3

    def g_rho(self, T):
        return np.interp(np.log(T), self.lnT, self.grho)

    def g_s(self, T):
        return np.interp(np.log(T), self.lnT, self.gs)

    def S_over_rho(self, T):
        return np.interp(np.log(T), self.lnT, self.Sfrac)

    def S_over_rho_err(self, T):
        return np.interp(np.log(T), self.lnT, self.Serr)

    def dlnT_dN_at(self, T):
        return np.interp(np.log(T), self.lnT, self.dlnT_dN)

    def S(self, T):
        return self.S_over_rho(T) * self.rho(T)

    def S_N(self, T, h=2e-3):
        """dS/dN along the adiabatic plasma history (entropy conservation)."""
        lT = math.log(T)
        Sp, Sm = self.S(math.exp(lT + h)), self.S(math.exp(lT - h))
        return (Sp - Sm) / (2 * h) * self.dlnT_dN_at(T)

    # --- w(ln rho) for the full-equation integrations (rho in GeV^4 times `unit`)
    def w_of_lnrho(self, rho_unit, Tmin=1e-3, Tmax=1e7):
        sel = (self.T >= Tmin) & (self.T <= Tmax)
        T = self.T[sel]
        rho = math.pi ** 2 / 30 * self.grho[sel] * T ** 4 / rho_unit
        w = (1 - self.Sfrac[sel]) / 3.0
        x = np.log(rho)
        order = np.argsort(x)
        return x[order], w[order]
