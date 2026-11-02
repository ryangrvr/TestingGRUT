"""qm_transfer: LCDM background for the linear Boltzmann code (own code, Claude, 2 Nov 2026).  [C]

Planck 2018 TT,TE,EE+lowE+lensing base-LCDM parameters [R] (Planck 2018 VI, arXiv:1807.06209):
h = 0.6736, omega_b = 0.02237, omega_c = 0.1200 (here: ALL dark matter = objects + vapour), N_eff = 3.046 (massless),
T_cmb = 2.7255 K, Y_p = 0.2454, n_s = 0.9649.
Provides: conformal Hubble rate calH(a) [1/Mpc], conformal time tau(a) [Mpc], recombination x_e(a) (Saha He +
Peebles H with the RECFAST fudge 1.14, Compton-coupled matter temperature), Thomson rate 1/tau_c(a) = a n_e sigma_T,
baryon sound speed c_s^2(a), R(a) = 4 rho_g / 3 rho_b.  An optional dark-matter energy factor e(a)/e(1) (used only by
the thermal-WDM validation, where the relic is relativistic early) multiplies Omega_dm a^-3.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid

C_KMS = 299792.458
MPC_M = 3.0856775814913673e22
GEV_TO_INVMPC = MPC_M / 1.973269804e-16          # 1 GeV = 1.5637e38 / Mpc
KB_EV = 8.617333262e-5
SIGMA_T = 6.6524587321e-29                      # m^2
M_E_EV = 510998.95
M_H_KG = 1.6735575e-27
RHO_C_KG_H2 = 1.87834e-26                       # kg/m^3 per h^2


class Cosmo:
    def __init__(self, h=0.6736, omb=0.02237, omc=0.1200, Neff=3.046, T0=2.7255, Yp=0.2454, ns=0.9649,
                 dm_efac=None, recomb=True):
        self.h, self.omb, self.omc, self.Neff, self.T0, self.Yp, self.ns = h, omb, omc, Neff, T0, Yp, ns
        self.H0 = h * 100.0 / C_KMS                       # 1/Mpc
        self.Og = 2.4728e-5 * (T0 / 2.7255)**4 / h**2
        self.Onu = self.Og * Neff * (7 / 8) * (4 / 11)**(4 / 3)
        self.Or = self.Og + self.Onu
        self.Ob = omb / h**2
        self.Odm = omc / h**2
        self.Om = self.Ob + self.Odm
        self.OL = 1 - self.Om - self.Or
        self.dm_efac = dm_efac                            # callable a -> e(a)/e(1) or None
        self.a_eq = self.Or / self.Om
        # conformal time table
        self.lna_tab = np.linspace(math.log(1e-14), 0.0, 40001)
        a = np.exp(self.lna_tab)
        cH = np.array([self.calH(x) for x in a])
        tau = cumulative_trapezoid(1.0 / cH, self.lna_tab, initial=0.0)
        tau += a[0] / cH[0]                              # RD: tau = a/calH at the first point
        self.tau_tab = tau
        self.tau0 = float(tau[-1])
        if recomb:
            self._recomb()

    # --- background ---
    def efac(self, a):
        return 1.0 if self.dm_efac is None else self.dm_efac(a)

    def calH(self, a):
        """conformal Hubble rate a H [1/Mpc]"""
        return self.H0 * math.sqrt(self.Or / a**2 + (self.Ob + self.Odm * self.efac(a)) / a + self.OL * a * a)

    def calH_v(self, a):
        a = np.asarray(a, float)
        ef = np.ones_like(a) if self.dm_efac is None else np.vectorize(self.dm_efac)(a)
        return self.H0 * np.sqrt(self.Or / a**2 + (self.Ob + self.Odm * ef) / a + self.OL * a * a)

    def tau(self, a):
        return float(np.interp(math.log(a), self.lna_tab, self.tau_tab))

    def a_of_tau(self, t):
        return math.exp(float(np.interp(t, self.tau_tab, self.lna_tab)))

    def Omega_m_a(self, a):
        """matter fraction of the critical density at a (radiation, Lambda included)"""
        m = (self.Ob + self.Odm * self.efac(a)) / a
        return m / (self.Or / a**2 + m + self.OL * a * a)

    def A_species(self, a):
        """4 pi G a^2 rho_i [1/Mpc^2] for (dm, b, g, nu)"""
        f = 1.5 * self.H0**2
        return f * self.Odm * self.efac(a) / a, f * self.Ob / a, f * self.Og / a**2, f * self.Onu / a**2

    # --- recombination ---
    def _recomb(self):
        Yp = self.Yp
        nH0 = (1 - Yp) * self.omb * RHO_C_KG_H2 / M_H_KG             # m^-3
        fHe = Yp / (3.9715 * (1 - Yp))
        self.nH0, self.fHe = nH0, fHe
        T0 = self.T0
        hbar_c_eVm = 1.973269804e-7
        lam_a = 121.5670e-9
        Lam2s = 8.2245809
        aR = 7.5657e-16                                               # J m^-3 K^-4
        me_kg, c = 9.1093837e-31, 2.99792458e8

        def saha_fac(TeV):  # (m_e T / 2 pi)^{3/2} in m^-3
            return (M_E_EV * TeV / (2 * math.pi))**1.5 / hbar_c_eVm**3

        def He_saha(z, xH):
            # electrons from H fraction xH; solve HeII/HeI and HeIII/HeII Saha with n_e self-consistent (iterate)
            TeV = KB_EV * T0 * (1 + z)
            nH = nH0 * (1 + z)**3
            S2 = 4 * saha_fac(TeV) * math.exp(-24.587 / TeV) / nH        # x_e x_HeII / x_HeI
            S3 = 1 * saha_fac(TeV) * math.exp(-54.418 / TeV) / nH        # x_e x_HeIII / x_HeII
            xe = xH + 2 * fHe
            for _ in range(60):
                r2 = S2 / xe; r3 = S3 / xe
                fI = 1.0 / (1 + r2 + r2 * r3)
                xe_new = xH + fHe * (r2 * fI + 2 * r2 * r3 * fI)
                if abs(xe_new - xe) < 1e-12:
                    break
                xe = 0.5 * (xe + xe_new)
            return xe

        def H_si(z):
            a = 1 / (1 + z)
            return self.calH(a) / a * C_KMS * 1e3 / MPC_M                  # 1/s

        def alphaB(Tm):
            t = Tm / 1e4
            return 1e-19 * 1.14 * 4.309 * t**(-0.6166) / (1 + 0.6703 * t**0.5300)   # m^3/s

        # Saha for H down to where x_H ~ 0.99, then Peebles
        zs = np.linspace(9000, 1700, 3000)
        xHs = []
        for z in zs:
            TeV = KB_EV * T0 * (1 + z)
            nH = nH0 * (1 + z)**3
            S = saha_fac(TeV) * math.exp(-13.605693 / TeV) / nH
            xH = (-S + math.sqrt(S * S + 4 * S)) / 2 if S < 1e8 else 1.0      # x^2/(1-x) = S (He electrons ignored here)
            xHs.append(min(xH, 1.0))
        z_start = 1700.0

        def rhs(lna, y):
            xp, Tm = y
            a = math.exp(lna); z = 1 / a - 1
            Tg = T0 * (1 + z)
            TeV = KB_EV * Tm
            nH = nH0 * (1 + z)**3
            xe = He_saha(z, xp) if z > 1500 else xp + (He_saha(z, xp) - xp)
            aB = alphaB(Tm)
            bB = aB * saha_fac(KB_EV * Tg) * math.exp(-3.4014 / (KB_EV * Tg))
            Hz = H_si(z)
            K = lam_a**3 / (8 * math.pi * Hz)
            n1s = nH * max(1 - xp, 1e-30)
            Cfac = (1 + K * Lam2s * n1s) / (1 + K * (Lam2s + bB) * n1s)
            dxp_dt = -Cfac * (aB * nH * xe * xp - bB * (1 - xp) * math.exp(-10.2043 / (KB_EV * Tg)))
            GC = 8 * SIGMA_T * aR * Tg**4 / (3 * me_kg * c) * xe / (1 + fHe + xe)
            dTm = -2 * Tm + GC / Hz * (Tg - Tm)
            return [dxp_dt / Hz, dTm]

        a0 = 1 / (1 + z_start)
        TeV0 = KB_EV * T0 * (1 + z_start)
        S0 = saha_fac(TeV0) * math.exp(-13.605693 / TeV0) / (nH0 * (1 + z_start)**3)
        xp0 = (-S0 + math.sqrt(S0 * S0 + 4 * S0)) / 2
        sol = solve_ivp(rhs, (math.log(a0), 0.0), [xp0, T0 * (1 + z_start)], method='LSODA', rtol=1e-7, atol=[1e-12, 1e-8],
                        dense_output=True, max_step=0.01)
        lg = np.linspace(math.log(1e-8), 0.0, 6000)
        xe = np.zeros_like(lg); Tm = np.zeros_like(lg)
        for i, l in enumerate(lg):
            a = math.exp(l); z = 1 / a - 1
            if l < math.log(a0):
                xH = float(np.interp(-z, -zs, xHs)) if z < 9000 else 1.0
                xe[i] = He_saha(z, xH) if z < 2e4 else 1 + 2 * fHe
                Tm[i] = T0 * (1 + z)
            else:
                xp, T = sol.sol(l)
                xe[i] = He_saha(z, xp) if z > 1500 else xp + max(He_saha(z, xp) - xp, 0.0)
                Tm[i] = T
        # reionisation (tanh in y=(1+z)^1.5, z_re = 7.7, dz = 0.5), He II at z = 3.5 [model choice; irrelevant for P(k)]
        z = 1 / np.exp(lg) - 1
        y = (1 + z)**1.5; yre = (1 + 7.7)**1.5; dy = 1.5 * (1 + 7.7)**0.5 * 0.5
        xre = (1 + fHe) * 0.5 * (1 + np.tanh((yre - y) / dy)) + fHe * 0.5 * (1 + np.tanh((3.5 - z) / 0.5))
        xe = np.maximum(xe, xre)
        self.rec_lna, self.rec_xe, self.rec_Tm = lg, xe, Tm
        dlT = np.gradient(np.log(Tm), lg)
        mu = 1.0 / ((1 - Yp) * (1 + xe) + Yp / 4.0)       # mean molecular weight (approx.)
        self.rec_cs2 = KB_EV * Tm / (mu * 938.272e6) * (1 - dlT / 3.0)
        # Thomson: a n_e sigma_T [1/Mpc]
        a = np.exp(lg)
        self.rec_dk = xe * nH0 * SIGMA_T / a**2 * MPC_M
        # z_* (visibility peak proxy) and z_drag-like: where tau_c^-1 = calH and R/tau_c = calH
        self.z_dec = None

    def xe(self, a):
        return float(np.interp(math.log(a), self.rec_lna, self.rec_xe, left=self.rec_xe[0]))

    def dkappa(self, a):
        """a n_e sigma_T [1/Mpc] (fully ionised scaling before the table)"""
        la = math.log(a)
        if la < self.rec_lna[0]:
            return self.rec_dk[0] * math.exp(2 * (self.rec_lna[0] - la))
        return float(np.interp(la, self.rec_lna, self.rec_dk))

    def cs2(self, a):
        la = math.log(a)
        if la < self.rec_lna[0]:
            return self.rec_cs2[0] * math.exp(self.rec_lna[0] - la)
        return float(np.interp(la, self.rec_lna, self.rec_cs2))

    def Rgb(self, a):
        """R = (4/3) rho_g / rho_b"""
        return 4.0 / 3.0 * self.Og / (self.Ob * a)

    def Tm(self, a):
        return float(np.interp(math.log(a), self.rec_lna, self.rec_Tm))


if __name__ == '__main__':
    c = Cosmo()
    print('Om=%.5f Or=%.3e a_eq=%.4e tau0=%.1f Mpc' % (c.Om, c.Or, c.a_eq, c.tau0))
    for z in (3000, 1500, 1200, 1100, 1090, 1000, 800, 500, 200, 50, 10, 6, 0):
        a = 1 / (1 + z)
        print('z=%5d xe=%.4e Tm/Tg=%.3f 1/tau_c/calH=%.3e cs2=%.3e' % (z, c.xe(a), c.Tm(a) / (c.T0 * (1 + z)), c.dkappa(a) / c.calH(a), c.cs2(a)))
    # z_* where optical depth = 1
    from scipy.integrate import quad
    lnas = np.linspace(math.log(1 / 3000), 0, 4000)
    dk = np.array([c.dkappa(math.exp(l)) / c.calH(math.exp(l)) for l in lnas])
    kap = cumulative_trapezoid(dk[::-1], -lnas[::-1], initial=0)[::-1]
    i = np.argmin(abs(kap - 1)); print('z_* (kappa=1) ~ %.1f, tau(z*)=%.1f Mpc' % (1 / math.exp(lnas[i]) - 1, c.tau(math.exp(lnas[i]))))
