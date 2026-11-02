"""qm_transfer: hidden-sector background for the Q-matter + vapour transfer calculation (own code, Claude, 2 Nov 2026). [C]

Physical model (M1r Q-matter branch, tracks/fx_formation, tracks/fx_observables):
  objects of charge Q (bulk spheres R = (3Q/4 pi n_B)^(1/3), n_B = 449.31 m_X^3) hold the charge fraction C_o = 1 - N_v;
  an ideal NR vapour X (mass m_X = 1) at chemical potential mu = (1 - Delta) m_X sits in two-phase equilibrium with
  them, n_sat(T) = g (T/2 pi)^(3/2) exp(-Delta/T) (g = 1: ideal vapour of fx_formation; g != 1 parametrises a different
  saturated density, e.g. from tracks/qm_micro).
  Quasi-static equilibrium conserves the adiabat N_v (5/2 + Delta/T) = S/Q (derived in NOTES.md for any g).
Kinetic model (used from a_kbg, where the capture rate Gamma_a(vbar)/H = G_BG, onward): the vapour phase-space
  distribution on a fixed comoving-momentum grid q_j (per-charge normalisation N_j = dN_v/dq), evolving by
      dN_j/dln a = (Gamma_a(q_j)/H) [ N_eq(q_j; T_o) - N_j ],
  Gamma_a = s pi R^2 n_obj v(q) = zeta K n_Q C_o^(2/3) v(q),  K = pi (3/4 pi n_B)^(2/3),  zeta = s Q^(-1/3),
  with the object temperature T_o fixed at every instant by the objects' energy balance (zero heat capacity: every
  absorbed charge's latent heat Delta + E_k must be re-emitted as evaporated X):
      sum_j w_j Gamma_a,j (Delta + E_k,j) [N_j - N_eq,j(T_o)] = 0.
  'thermo'  : this energy-balance thermostat throughout (re-emission continues after capture freeze-out);
  'sudden'  : the same until capture freeze-out (Gamma_a(vbar) = H), then absorption AND emission switched off
              (the fx_formation / fx_observables picture: vapour number and comoving momenta frozen);
  'fx'      : the exact history_B arrays of tracks/fx_observables/code/vapour_transfer.py (validation only);
  'frozen'  : no collisions at all (thermal-WDM validation; Fermi-Dirac or Maxwellian distribution).
Elastic scattering: sigma_el = c_el pi R^2 (default c_el = 1 - s, so that absorption + elastic = pi R^2 is the
  momentum-transfer cross-section for any sticking; c_el = 1 is the alternative).
Outputs tables on a uniform ln a grid for the Boltzmann code (qt_boltz.Mode, DM kind 'kin'/'tca').
"""
import math
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

from qt_bg import GEV_TO_INVMPC

MX = 0.0924089637375362          # GeV
NB = 449.31336788
KGEO = math.pi * (3 / (4 * math.pi * NB))**(2 / 3)
RHOC_GEV4_H2 = 8.0992e-47        # critical density / h^2 [GeV^4]


def D_moment(N, q, w, D=None):
    """D = q dN/dq - 2N on the ln q grid: cubic-spline derivative (if D is None) plus a correction (a + b x + c x^2) N
    (x = q/q_N) that makes the discrete moments exact: sum w q^n D = -(n+3) sum w q^n N for n = 0, 1, 2
    (number, momentum, energy identities of int q^n (q N' - 2N) dq).  N: (nt, nq) or (nq,)."""
    from scipy.interpolate import CubicSpline
    one = (N.ndim == 1)
    N2 = np.atleast_2d(N)
    lq = np.log(q)
    if D is None:
        D2 = CubicSpline(lq, N2, axis=1)(lq, 1) - 2 * N2
    else:
        D2 = np.atleast_2d(D).copy()
    out = D2.copy()
    for i in range(N2.shape[0]):
        Ni = N2[i]
        tot = float(w @ Ni)
        if tot <= 0:
            continue
        qN = float(w @ (Ni * q)) / tot
        x = q / qN
        M = np.array([[float(w @ (x**(n + m) * Ni)) for m in range(3)] for n in range(3)])
        r = np.array([-(n + 3) * float(w @ (x**n * Ni)) - float(w @ (x**n * D2[i])) for n in range(3)])
        cc = np.linalg.solve(M, r)
        out[i] = D2[i] + (cc[0] + cc[1] * x + cc[2] * x * x) * Ni
    return out[0] if one else out


def Omega_fit(Q):
    return 0.775865 + (2 / 3) * 8.1429 * Q**(-1 / 3) + (1 / 3) * 125.79 * Q**(-2 / 3)


def EQ_fit(Q):
    return 0.775865 + 8.1429 * Q**(-1 / 3) + 125.79 * Q**(-2 / 3) - 2847.0 / Q


class Vapour:
    """hidden-sector background + per-q tables.  Parameters:
       SQ (S/Q adiabat), Q, s (sticking), g (saturation factor), mode, c_el, NQ0 (comoving n_Q in m_X^3; default:
       from omega_dm and the mass per charge), sigma overrides for small Q (validation), Nq, dlnq, Lv."""
    kind = 'vapour'

    def __init__(self, cos, SQ=1.235, Q=1e12, s=1.0, g=1.0, mode='thermo', c_el=None, NQ0=None, Delta=None,
                 sig_a=None, sig_el=None, Nq=28, Lv=10, G_BG=60.0, G_KIN=30.0, eps_TC=0.03, eps_coll=0.01,
                 qlo_fac=0.15, qhi_fac=25.0, fx_hist=None, a_lo=1e-14, verbose=False, s_of_E=None, s_scale=1.0):
        self.c = cos
        self.SQ, self.Q, self.s, self.g, self.mode = SQ, Q, s, g, mode
        self.c_el = (1 - s) if c_el is None else c_el
        self.Delta = (1 - Omega_fit(Q)) if Delta is None else Delta
        self.eps_o = 1 - self.Delta
        self.Lv = Lv
        self.G_BG, self.G_KIN, self.eps_TC, self.eps_coll = G_BG, G_KIN, eps_TC, eps_coll
        # energy-dependent sticking s(E_k) = s_scale * s_of_E(E_k) (e.g. tracks/qm_micro threshold-law table); then
        # zeta is quoted at the freeze-out thermal energy and nsig_a/nsig_el below are GEOMETRIC (s = 1)
        self.s_of_E, self.s_scale = s_of_E, s_scale
        if s_of_E is not None:
            self.s = 1.0; self.c_el = 0.0 if c_el is None else c_el
        self.zeta = self.s * Q**(-1 / 3)
        # comoving charge density: rho_dm = n_Q m_X e_DM  (e_DM ~ eps_o + N_v Delta; iterate once)
        rho_dm = cos.omc * RHOC_GEV4_H2
        if NQ0 is None:
            self.NQ0 = rho_dm / (MX**4 * (self.eps_o + 0.04 * self.Delta))
            self._auto_nq = True
        else:
            self.NQ0 = NQ0; self._auto_nq = False
        self.sig_a, self.sig_el = sig_a, sig_el          # explicit cross sections [m_X^-2] (with n_obj = C_o n_Q / Q)
        self.fx_hist = fx_hist
        self.verbose = verbose
        self.Nq_req, self.qlo_fac, self.qhi_fac = Nq, qlo_fac, qhi_fac
        self._equilibrium(a_lo)
        if self._auto_nq:     # one refinement with the actual late vapour fraction
            self.NQ0 = rho_dm / (MX**4 * (self.eps_o + self.Nv_f_eq * self.Delta))
            self._equilibrium(a_lo)
        if mode == 'fx':
            self._fx_tables()
        else:
            self._kinetic()

    # ---------------- validation mode: the exact fx_observables history (history_B of vapour_transfer.py)
    def _fx_tables(self):
        """fx_hist: dict(a, Nv, Tv, af, am, sig_c, sig_m, Q) from tracks/fx_observables (two_fluid_check.history_full).
        Before a_f the vapour is in chemical + thermal equilibrium at T_v (N = N_eq); absorption/emission at
        n_obj sig_c v(q) with n_obj = (1 - N_v) n_Q / Q; elastic n_obj (sig_m - sig_c) v(q).  After a_f: no exchange,
        elastic n_obj sig_m v(q), T_v = T_f (a_f/a)^2 (frozen comoving Maxwellian)."""
        H = self.fx_hist
        self.NQ0 = H['NQ0']
        af = H['af']
        Tf = float(np.interp(math.log(af), np.log(H['a']), H['Tv']))
        a_lo = H['a'][0]
        T0 = H['Tv'][0]
        q_lo = self.qlo_fac * a_lo * math.sqrt(T0) * 3
        q_hi = 8.0 * af * math.sqrt(Tf)
        lq = np.linspace(math.log(q_lo), math.log(q_hi), self.Nq_req)
        self.q = np.exp(lq); self.dlnq = lq[1] - lq[0]; self.wq = self.q * self.dlnq; self.Nq = len(self.q)
        lg = np.linspace(math.log(a_lo), 0.0, 3001)
        a = np.exp(lg)
        Nv = np.interp(lg, np.log(H['a']), H['Nv'], right=H['Nv'][-1])
        Tv = np.where(a < af, np.exp(np.interp(lg, np.log(H['a']), np.log(H['Tv']))), Tf * (af / a)**2)
        sig2 = Tv[:, None] * a[:, None]**2
        shape = self.q[None, :]**2 * np.exp(-self.q[None, :]**2 / (2 * sig2))
        N = shape / (shape @ self.wq)[:, None] * Nv[:, None]
        self.k_lna, self.k_dl = lg, lg[1] - lg[0]
        self.k_N, self.k_T = N, Tv
        self.af_kin = af
        q = self.q[None, :]
        E = np.sqrt(1 + (q / a[:, None])**2); Ek = E - 1; v = (q / a[:, None]) / E
        nQ = self.NQ0 / a**3
        nobj = (1 - Nv) * nQ / H['Q']
        Hh = np.array([self.H_mX(x) for x in a]); cH = self.c.calH_v(a)
        pre = a < af
        Ga_H = np.where(pre[:, None], (nobj * H['sig_c'] / Hh)[:, None] * v, 0.0)
        Ge_H = np.where(pre[:, None], (nobj * (H['sig_m'] - H['sig_c']) / Hh)[:, None] * v, (nobj * H['sig_m'] / Hh)[:, None] * v)
        ga, ge = Ga_H * cH[:, None], Ge_H * cH[:, None]
        Neq = N.copy()
        D = -(q * q / (a[:, None]**2 * E * Tv[:, None])) * N     # analytic Maxwellian (relativistic E_k)
        D = D_moment(N, self.q, self.wq, D)
        self.eps_o = 1.0          # the fx background carries no latent-heat bookkeeping: keep rho_dm a^3 fixed
        T = Tv[:, None]
        dTN = Neq * (self.Delta + Ek) / T**2
        dTden = (self.wq * ga * (self.Delta + Ek) * dTN).sum(axis=1)
        Co = 1 - Nv
        eDM = self.eps_o * Co + (N * E) @ self.wq
        self.T = dict(N=N, Neq=Neq, E=E, Ek=Ek, v=v, ga=ga, ge=ge, D=D, Deq=D.copy(), dTNeq=dTN, dTden=np.maximum(dTden, 1e-300),
                      jc=np.zeros(len(lg)), Co=Co, eo=self.eps_o * Co, eDM=eDM, dNdlna=np.gradient(N, lg, axis=0), GaH=Ga_H, GeH=Ge_H,
                      dlnGa=1 / E**2)
        self.Nv_tab = Nv
        Gmean = ((Ga_H + Ge_H) * N) @ self.wq / Nv
        self.Gmean_H = Gmean
        jj = np.where(Gmean > self.eps_coll)[0]
        self.a_coll_end = float(a[min(jj[-1] + 1, len(lg) - 1)]) if len(jj) else float(a[0])
        self.a_kbg = float(a[0])
        self.f_v_charge = float(Nv[-1]); self.f_v_mass = float((N[-1] * E[-1]) @ self.wq / eDM[-1])
        self.sigma_p = math.sqrt(float((N[-1] * self.q**2) @ self.wq / Nv[-1]) / 3)
        self.a_f = af

    # ---------------- rates
    def H_mX(self, a):
        return self.c.calH(a) / a / (MX * GEV_TO_INVMPC)

    def nsig_a(self, a, Co):
        """n_obj sigma_a [m_X]"""
        nQ = self.NQ0 / a**3
        if self.sig_a is not None:
            return Co * nQ / self.Q * self.sig_a
        return self.zeta * KGEO * nQ * Co**(2 / 3)

    def nsig_el(self, a, Co):
        nQ = self.NQ0 / a**3
        if self.sig_el is not None:
            return Co * nQ / self.Q * self.sig_el
        return self.c_el * (self.zeta / self.s) * KGEO * nQ * Co**(2 / 3) if self.s > 0 else 0.0

    def sE(self, Ek):
        """per-particle sticking factor (multiplies the geometric absorption); 1 for constant-s models (s in zeta)"""
        if self.s_of_E is None:
            return np.ones_like(Ek)
        return np.clip(self.s_scale * self.s_of_E(Ek), 0.0, 1.0)

    # ---------------- quasi-static equilibrium
    def eq_T(self, a):
        nQ = self.NQ0 / a**3
        SQ, D, g = self.SQ, self.Delta, self.g
        f = lambda lT: (math.log(g) + 1.5 * (lT - math.log(2 * math.pi)) - D / math.exp(lT) - math.log(nQ)
                        - math.log(SQ / (2.5 + D / math.exp(lT))))
        lo, hi = math.log(1e-5), math.log(0.3)
        if f(hi) < 0:
            return None
        lT = brentq(f, lo, hi, xtol=1e-14)
        T = math.exp(lT)
        return T, SQ / (2.5 + D / T)

    def _equilibrium(self, a_lo):
        """equilibrium history on a ln a grid; locate the capture freeze-out estimate a_f (Gamma_a(vbar)/H = 1)."""
        lg = np.linspace(math.log(a_lo), 0.0, 6001)
        rows = []
        for l in lg:
            a = math.exp(l)
            r = self.eq_T(a)
            if r is None or r[0] > 0.05 or r[1] > 0.9:      # two-phase NR model only (vapour-dominated states excluded)
                continue
            T, Nv = r
            vb = math.sqrt(8 * T / math.pi)
            sf = float(self.sE(np.array([1.5 * T]))[0])
            G = self.nsig_a(a, 1 - Nv) * sf * vb / self.H_mX(a)
            Gm = (self.nsig_a(a, 1 - Nv) * (sf if self.s_of_E is None else 1.0) + self.nsig_el(a, 1 - Nv)) * vb / self.H_mX(a)
            rows.append((l, T, Nv, G, Gm))
        R = np.array(rows)
        self.eq_lna, self.eq_T_arr, self.eq_Nv, self.eq_G, self.eq_Gm = R[:, 0], R[:, 1], R[:, 2], R[:, 3], R[:, 4]
        j = np.argmax(R[:, 3] < 1.0) if (R[:, 3] < 1.0).any() else len(R) - 1
        self.a_f_eq = math.exp(R[j, 0]); self.T_f_eq = R[j, 1]; self.Nv_f_eq = R[j, 2]
        jb = np.argmax(R[:, 3] < self.G_BG) if (R[:, 3] < self.G_BG).any() else 0
        self.a_kbg = math.exp(R[max(jb, 0), 0])
        self.a_start_eq = math.exp(R[0, 0])

    def eq_at(self, a):
        l = math.log(a)
        T = float(np.exp(np.interp(l, self.eq_lna, np.log(self.eq_T_arr))))
        Nv = float(np.interp(l, self.eq_lna, self.eq_Nv))
        return T, Nv

    # ---------------- kinetic background on the q grid
    def _grid(self):
        T0 = self.eq_at(self.a_kbg)[0]
        q_lo = self.qlo_fac * self.a_kbg * math.sqrt(T0)
        q_hi = self.qhi_fac * self.a_f_eq * math.sqrt(self.T_f_eq)
        lq = np.linspace(math.log(q_lo), math.log(q_hi), self.Nq_req)
        self.q = np.exp(lq)
        self.dlnq = lq[1] - lq[0]
        self.wq = self.q * self.dlnq                     # trapezoid in ln q (ends are negligible)
        self.Nq = len(self.q)

    def Neq_q(self, a, T):
        q = self.q
        E = np.sqrt(1 + (q / a)**2)
        Ek = E - 1
        pref = self.g * 4 * math.pi * q * q / ((2 * math.pi)**3 * self.NQ0)
        return pref * np.exp(-(self.Delta + Ek) / T), E, Ek

    def solve_To(self, a, N, ga, T_guess):
        """energy balance for T_o: sum w ga (Delta+Ek) (N - Neq(T)) = 0"""
        q = self.q
        E = np.sqrt(1 + (q / a)**2); Ek = E - 1
        D = self.Delta
        wgt = self.wq * ga * (D + Ek)
        B = float(wgt @ N)
        if B <= 0:
            return T_guess
        lpref = np.log(self.g * 4 * math.pi * q * q / ((2 * math.pi)**3 * self.NQ0) * wgt + 1e-300)
        x = D + Ek
        lnB = math.log(B)

        def F(lT):
            T = math.exp(lT)
            z = lpref - x / T
            zm = z.max()
            return zm + math.log(np.exp(z - zm).sum()) - lnB

        lT = math.log(T_guess)
        for _ in range(60):
            T = math.exp(lT)
            z = lpref - x / T
            zm = z.max(); e = np.exp(z - zm)
            Fv = zm + math.log(e.sum()) - lnB
            dF = float((e * x).sum() / e.sum()) / T            # dF/dlnT
            step = -Fv / dF
            step = max(min(step, 0.5), -0.5)
            lT += step
            if abs(step) < 1e-13:
                break
        return math.exp(lT)

    def _kinetic(self):
        self._grid()
        c = self.c
        a0 = self.a_kbg
        T0, Nv0 = self.eq_at(a0)
        N0, _, _ = self.Neq_q(a0, T0)
        N0 *= Nv0 / float(self.wq @ N0)
        self.af_kin = None
        mode = self.mode
        lna_end = 0.0
        lg = np.linspace(math.log(a0), lna_end, 3001)
        state = {'T': T0, 'frozen': False}

        def rates(a, N):
            Co = 1 - float(self.wq @ N)
            E = np.sqrt(1 + (self.q / a)**2)
            v = (self.q / a) / E
            Hh = self.H_mX(a)
            Ga = self.nsig_a(a, Co) * v * self.sE(E - 1) / Hh
            return Ga, Co

        def rhs(l, N):
            a = math.exp(l)
            if state['frozen']:
                return np.zeros_like(N)
            Ga, Co = rates(a, N)
            To = self.solve_To(a, N, Ga, state['T'])
            state['T'] = To
            Neq, E, Ek = self.Neq_q(a, To)
            return Ga * (Neq - N)

        def jac(l, N):
            a = math.exp(l)
            Ga, Co = rates(a, N)
            To = state['T']
            Neq, E, Ek = self.Neq_q(a, To)
            x = self.Delta + Ek
            dT = Neq * x / To**2
            wgt = self.wq * Ga * x
            den = float(wgt @ dT)
            J = -np.diag(Ga) + np.outer(Ga * dT, wgt / den)
            return J

        # integrate piecewise, record table
        tab_N = np.zeros((len(lg), self.Nq)); tab_T = np.zeros(len(lg))
        tab_N[0] = N0; tab_T[0] = T0
        N = N0.copy()
        for i in range(1, len(lg)):
            a = math.exp(lg[i])
            if not state['frozen']:
                sol = solve_ivp(rhs, (lg[i - 1], lg[i]), N, method='BDF', jac=jac, rtol=1e-9, atol=1e-16 + 1e-10 * float(N.max()))
                N = np.maximum(sol.y[:, -1], 0.0)
                Ga, Co = rates(a, N)
                state['T'] = self.solve_To(a, N, Ga, state['T'])
                Gbar = float(self.wq @ (N * Ga)) / float(self.wq @ N)
                if self.af_kin is None and Gbar < 1.0:
                    self.af_kin = a
                    if mode == 'sudden':
                        state['frozen'] = True
                if mode == 'thermo' and Gbar < 1e-6:
                    state['frozen'] = True
            tab_N[i] = N; tab_T[i] = state['T']
        self.k_lna = lg
        self.k_dl = lg[1] - lg[0]
        self.k_N = tab_N
        self.k_T = tab_T
        self._derive_tables()

    def _derive_tables(self):
        lg = self.k_lna
        a = np.exp(lg)[:, None]
        q = self.q[None, :]
        E = np.sqrt(1 + (q / a)**2)
        Ek = E - 1
        v = (q / a) / E
        N = self.k_N
        T = self.k_T[:, None]
        Nv = N @ self.wq
        Co = 1 - Nv
        pref = self.g * 4 * math.pi * q * q / ((2 * math.pi)**3 * self.NQ0)
        Neq = pref * np.exp(-(self.Delta + Ek) / T)
        Hh = np.array([self.H_mX(x) for x in a[:, 0]])
        cH = self.c.calH_v(a[:, 0])
        nsa = np.array([self.nsig_a(x, y) for x, y in zip(a[:, 0], Co)])
        nse = np.array([self.nsig_el(x, y) for x, y in zip(a[:, 0], Co)])
        sEq = self.sE(Ek)
        Ga_H = nsa[:, None] * v * sEq / Hh[:, None]
        Ge_H = nse[:, None] * v / Hh[:, None] if self.s_of_E is None else nsa[:, None] * v * (1 - sEq) / Hh[:, None]
        # q d ln Gamma_a / d q (enters the boost term x (N_eq - N))
        dlnGa = 1 / E**2 + (np.gradient(np.log(np.maximum(sEq, 1e-300)), np.log(self.q), axis=1) if self.s_of_E is not None else 0.0)
        froz = np.zeros(len(lg), bool)
        if self.mode == 'sudden' and self.af_kin is not None:
            froz = a[:, 0] > self.af_kin
        Ga_H = np.where(froz[:, None], 0.0, Ga_H)
        if self.mode == 'sudden' and self.af_kin is not None:
            # fx picture: after a_f the momentum coupling is pi R^2 (elastic) only
            Ge_H = np.where(froz[:, None], (nse + nsa)[:, None] * v / Hh[:, None], Ge_H)
        ga = Ga_H * cH[:, None]          # conformal [1/Mpc]
        ge = Ge_H * cH[:, None]
        # derivatives in ln q
        D = D_moment(N, self.q, self.wq)
        Deq = -(q * q / (a * a * E * T)) * Neq            # analytic: N_eq ~ q^2 exp(-E_k/T)
        dTN = Neq * (self.Delta + Ek) / T**2
        dTden = ((self.wq * ga * (self.Delta + Ek) * dTN)).sum(axis=1)
        jc = (self.wq * ga * (N - Neq)).sum(axis=1)       # conformal charge gain of the objects per charge
        dNdl = np.gradient(N, lg, axis=0)
        eDM = self.eps_o * Co + (N * E) @ self.wq
        self.T = dict(N=N, Neq=Neq, E=E, Ek=Ek, v=v, ga=ga, ge=ge, D=D, Deq=Deq, dTNeq=dTN, dTden=np.maximum(dTden, 1e-300),
                      jc=jc, Co=Co, eo=self.eps_o * Co, eDM=eDM, dNdlna=dNdl, GaH=Ga_H, GeH=Ge_H, dlnGa=dlnGa * np.ones_like(N))
        self.Nv_tab = Nv
        # end of collisions: max over q of (Ga+Ge)/H < eps_coll after freeze-out
        Gmean = ((Ga_H + Ge_H) * N) @ self.wq / np.maximum(Nv, 1e-300)      # N-weighted collision rate / H
        self.Gmean_H = Gmean
        jj = np.where(Gmean > self.eps_coll)[0]
        self.a_coll_end = float(a[min(jj[-1] + 1, len(lg) - 1), 0]) if len(jj) else float(a[0, 0])
        # summary
        late = len(lg) - 1
        self.f_v_charge = float(Nv[late])
        self.f_v_mass = float((N[late] * E[late]) @ self.wq / eDM[late])
        q2 = float((N[late] * self.q**2) @ self.wq / Nv[late])
        self.sigma_p = math.sqrt(q2 / 3.0)                 # 1D comoving velocity dispersion today [c]
        # mean comoving momentum etc.
        self.a_f = self.af_kin if self.af_kin is not None else self.a_f_eq

    # ---------------- interface for qt_boltz
    def kin_at(self, a):
        l = math.log(a)
        x = (l - self.k_lna[0]) / self.k_dl
        if x <= 0:
            i, w = 0, 0.0
        elif x >= len(self.k_lna) - 1:
            i, w = len(self.k_lna) - 2, 1.0
        else:
            i = int(x); w = x - i
        T = self.T
        out = {}
        for key in ('N', 'Neq', 'D', 'Deq', 'dTNeq', 'ga', 'ge', 'dNdlna', 'dlnGa'):
            out[key] = (1 - w) * T[key][i] + w * T[key][i + 1]
        for key in ('Co', 'eo', 'eDM', 'jc', 'dTden'):
            out[key] = (1 - w) * T[key][i] + w * T[key][i + 1]
        q = self.q
        E = np.sqrt(1 + (q / a)**2)
        out['E'] = E; out['Ek'] = E - 1; out['v'] = (q / a) / E
        out['Delta'] = self.Delta; out['eps_o'] = self.eps_o
        # conformal rates scale with a between table points ~ smooth; keep linear interpolation
        return out

    def ceff2(self, a):
        T, Nv = self.eq_at(a)
        x = self.Delta / T
        gam = (2.5 + x) / (1.5 + x - x / (2.5 + x))
        return gam * Nv * T / (self.eps_o * (1 - Nv) + Nv * (1 + 1.5 * T))

    def tca_local(self, a):
        T, Nv = self.eq_at(a)
        x = self.Delta / T
        Dfac = 1.5 + x - x / (2.5 + x)
        r = (1.5 + x) / Dfac
        return r, Dfac, T

    def a_kin_for_k(self, k, cos):
        """TCA -> kinetic switch for mode k: the earliest of (Gamma_m/H = G_KIN, eps_TC(k) = eps_TC), not before a_kbg."""
        best = None
        for l, T, Nv, G, Gm in zip(self.eq_lna, self.eq_T_arr, self.eq_Nv, self.eq_G, self.eq_Gm):
            a = math.exp(l)
            if a < self.a_kbg:
                continue
            vb = math.sqrt(8 * T / math.pi)
            eps = k * vb / (Gm * cos.calH(a))
            if Gm < self.G_KIN or eps > self.eps_TC:
                best = a; break
        if best is None:
            best = self.a_kbg
        return max(best, self.a_kbg)

    def summary(self):
        zeff = self.zeta
        if self.s_of_E is not None and self.a_f:
            j = int(np.argmin(abs(self.k_lna - math.log(self.a_f))))
            N = self.k_N[j]; Ek = np.sqrt(1 + (self.q / self.a_f)**2) - 1
            v = np.sqrt(1 - 1 / (1 + Ek)**2)
            zeff = float((self.wq * N * v * self.sE(Ek)).sum() / (self.wq * N * v).sum()) * self.Q**(-1 / 3)
        return dict(mode=self.mode, SQ=self.SQ, Q=self.Q, s=self.s, zeta=self.zeta, zeta_eff_f=zeff, g=self.g, c_el=self.c_el,
                    Delta=self.Delta, NQ0=self.NQ0, a_kbg=self.a_kbg, a_f_eq=self.a_f_eq, T_f_eq=self.T_f_eq,
                    Nv_f_eq=self.Nv_f_eq, a_f_kin=self.af_kin, a_coll_end=self.a_coll_end, f_v_charge=self.f_v_charge,
                    f_v_mass=self.f_v_mass, sigma_p=self.sigma_p, Nq=self.Nq, dlnq=self.dlnq, q_lo=float(self.q[0]),
                    q_hi=float(self.q[-1]))


class Frozen(Vapour):
    """collisionless relic with a fixed comoving-momentum distribution (all DM or a fraction): thermal-WDM validation.
       dist = 'FD' (Fermi-Dirac, comoving temperature Tq in m units) or 'MB' (Maxwellian, sigma = 1D dispersion)."""
    kind = 'frozen'

    def __init__(self, cos, fraction=1.0, dist='FD', Tq=None, sigma=None, Nq=24, Lv=12, qrange=(0.05, 14.0), a_lo=1e-10):
        self.c = cos
        self.Lv = Lv
        self.Delta = 0.0
        self.eps_o = 1.0
        self.a_coll_end = 0.0
        scale = Tq if dist == 'FD' else sigma
        lq = np.linspace(math.log(qrange[0] * scale), math.log(qrange[1] * scale), Nq)
        self.q = np.exp(lq); self.dlnq = lq[1] - lq[0]; self.wq = self.q * self.dlnq; self.Nq = Nq
        if dist == 'FD':
            shape = self.q**2 / (np.exp(self.q / Tq) + 1)
        else:
            shape = self.q**2 * np.exp(-self.q**2 / (2 * sigma**2))
        Nn = shape / float(self.wq @ shape) * fraction            # number fraction with the same mass per particle
        self.k_lna = np.linspace(math.log(a_lo), 0.0, 4001); self.k_dl = self.k_lna[1] - self.k_lna[0]
        nT = len(self.k_lna)
        a = np.exp(self.k_lna)[:, None]
        E = np.sqrt(1 + (self.q[None, :] / a)**2)
        N = np.repeat(Nn[None, :], nT, axis=0)
        if dist == 'FD':
            D0 = -Nn * (self.q / Tq) / (1 + np.exp(-self.q / Tq))
        else:
            D0 = -Nn * self.q**2 / sigma**2
        D0 = D_moment(Nn, self.q, self.wq, D0)
        D = np.repeat(D0[None, :], nT, axis=0)
        Z = np.zeros_like(N)
        Co = np.full(nT, 1 - fraction)
        eDM = Co + (N * E) @ self.wq
        self.T = dict(N=N, Neq=Z, D=D, Deq=Z, dTNeq=Z, ga=Z, ge=Z, dNdlna=Z, Co=Co, eo=Co, eDM=eDM,
                      jc=np.zeros(nT), dTden=np.ones(nT), dlnGa=Z)
        self.eDM_late = float(eDM[-1])
        self.f_v_mass = fraction
        self.sigma_p = math.sqrt(float(Nn @ (self.wq * self.q**2)) / float(Nn @ self.wq) / 3)

    def efac(self, a):
        l = math.log(a)
        return float(np.interp(l, self.k_lna, self.T['eDM'])) / self.eDM_late

    def a_kin_for_k(self, k, cos):
        return 0.0
