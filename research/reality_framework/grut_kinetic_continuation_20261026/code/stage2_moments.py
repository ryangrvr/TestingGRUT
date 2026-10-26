# Continued by Codex from Claude partial code at e71730c24fe6a0b28ddd72cd30982c5b39c9aedf.
# Modifications and limitations are recorded in this packet.
#!/usr/bin/env python3
"""Port of Claude's unfinished stage 2. Production imports only its network,
Case, Solver.rates and capture helpers. The legacy inversion integrator is
retained for attribution/context; it is NOT the production solver and no
stage1_relax.py was supplied or run. Use moment_controls.py for the labelled
controls. No full Boltzmann or formation solution is asserted.

State (per physical volume, converted to comoving by (a/a_ref)^3):
    n_Q(a)  net charge, comoving by construction (= Y_Q s_vis)
    P       pair inventory n_-^gas + n_s + n_cb   (ODE variable ln P)
    rho_tot gas energy + e_cap n_c + e_cb n_cb    (ODE variable E = rho_tot a^3)
plus ledgers: Y_+^tot (independent integration of the X total for the charge check), the P dV work, the hard-X
energy emitted by capture, event counts per channel.  Closure: one kinetic temperature (elastic rates >> H,
reported), fast s (mu_s = mu_+ + mu_-; closure remains conditional and rates are reported), Bose state functions (kin_thermo thermo.py), caps
(objects at omega, or Bose saturation).  Collision terms: MB rate densities of kin_rates with the Bose/MB factors
and exact MB chemical-potential factors; each net term vanishes at its own equilibrium relation (stationarity).
    ss -> X Xb (+ reverse)          gamma_ss(T) B_ss(T) [e^{2mu_s/T} - e^{(mu_+ + mu_-)/T}]            dP = -1, dN_+ = +1
    n -> 2 subtracted remainders    gamma_sub(T) [e^{sum mu_in/T} - e^{sum mu_out/T}] (resolved only)  stoichiometric
    capture Xb + L(Q) -> X + L(Q-2) n_L n_-^gas <sv>_f(T) (1 + f_+(E_out)) [1 - e^{(mu_+ - mu_- - 2 omega)/T}]  dP = -1, dN_+^tot = -1
    condensate 4->2 (rule a)        A_eff(T) [(n_+^tot)^2 - (n_+^gas)^2] n_-^2/4 [1 - e^{-(mu_+ + mu_-)/T}] etc.
Expansion: d(rho a^3) = -P_gas d(a^3); objects pressureless; no portal; the sources are complete before the start
(kin_source: 220-7000x before a_ref) and enter through the initial energy.  Object energy: e_cap per condensed charge
= omega (linear comparator) or E(Q)/Q (finite-Q booking at fixed Q; dE_L = omega d(n_L Q) + (E - Q omega) dn_L with
Q fixed is exactly (E/Q) d n_c).
"""
import json
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import ks_common as K
from thermo_state import Inverter, equilibrium_inventory

MX = K.MX
RS = K.RS
RATES = K.Rates()
OBJ = K.Objects()
_CONV = {}


def conv_tables(extrap):
    if extrap not in _CONV:
        ts, sv, Ein = K.conversion_table(extrap)
        _CONV[extrap] = (np.log(ts), np.log(sv), Ein)
    return _CONV[extrap]


def capture_kernel(t, extrap):
    lts, lsv, Ein = conv_tables(extrap)
    lt = np.clip(np.log(t), lts[0], lts[-1])
    return float(np.exp(np.interp(lt, lts, lsv))), float(np.interp(lt, lts, Ein))


# n <-> 2 channels used (P-changing), with (dP, dN_+^tot) per forward event and the exponents of e^{mu/T}:
# forward exponent (a_+, a_-, a_s) in units of (mu_+, mu_-, mu_s); reverse likewise
N2 = {
    'A_XXXbXb_XXb': dict(dP=-1, dNp=-1, fwd=(2, 2, 0), rev=(1, 1, 0)),
    'B_XXXXb_XX':   dict(dP=-1, dNp=-1, fwd=(3, 1, 0), rev=(2, 0, 0)),
    'E_XXbs_XXb':   dict(dP=-1, dNp=0,  fwd=(1, 1, 1), rev=(1, 1, 0)),
    'E2_XXs_XX':    dict(dP=-1, dNp=0,  fwd=(2, 0, 1), rev=(2, 0, 0)),
    'F_sss_ss':     dict(dP=-1, dNp=0,  fwd=(0, 0, 3), rev=(0, 0, 2)),
    'G2_Xss_Xs':    dict(dP=-1, dNp=0,  fwd=(1, 0, 2), rev=(1, 0, 1)),
    'L_sss_XXb':    dict(dP=-2, dNp=+1, fwd=(0, 0, 3), rev=(1, 1, 0)),
    'N_ssss_ss':    dict(dP=-2, dNp=0,  fwd=(0, 0, 4), rev=(0, 0, 2)),
    'H_XXbss_XXb':  dict(dP=-2, dNp=0,  fwd=(1, 1, 2), rev=(1, 1, 0)),
    'J_XXXbs_Xs':   dict(dP=-1, dNp=-1, fwd=(2, 1, 1), rev=(1, 0, 1)),
    'K_XXbss_ss':   dict(dP=-1, dNp=-1, fwd=(1, 1, 2), rev=(0, 0, 2)),
}


class Case:
    def __init__(self, label, R_total, rule='b', with_V0=True, e_cap_variant='omega', capture='hold', start='alpha',
                 P0_mode='equilibrium', bose=True, n2=True, catalysis=True, Q_L=None, background='table',
                 x_end=9.5, rtol=1e-8, atol=1e-12, t_stop=0.01, Pfloor=1e-28, capture_Q_scaling=False, note=''):
        self.__dict__.update(locals()); del self.__dict__['self']
        if rule == 'b':
            self.Q_L = K.Q_REF if Q_L is None else Q_L
            self.omega_L = K.OMEGA_REF if self.Q_L == K.Q_REF else OBJ.omega(self.Q_L)
            self.eq_L = K.E_OVER_QM_REF if self.Q_L == K.Q_REF else OBJ.e_over_q(self.Q_L)
            self.cap_p = self.omega_L
            self.e_cap = self.omega_L if e_cap_variant == 'omega' else self.eq_L
            self.cap_m = self.omega_L
            self.e_cm = self.e_cap
        elif rule == 'c':
            self.Q_L = 1.521e5 if Q_L is None else Q_L
            self.omega_L = OBJ.omega(self.Q_L)
            self.eq_L = OBJ.e_over_q(self.Q_L)
            self.cap_p = 1.0 - 1e-7
            self.e_cap = self.eq_L                 # saturated gas condenses at m, fragments at once to E(Q_L)/Q_L
            self.cap_m = 1.0 - 1e-7
            self.e_cm = 1.0
        else:                                      # 'a' gas only: homogeneous condensate at m_X
            self.Q_L = None; self.omega_L = None; self.eq_L = 1.0
            self.cap_p = 1.0 - 1e-7; self.e_cap = 1.0; self.cap_m = 1.0 - 1e-7; self.e_cm = 1.0
            self.capture = None

    def to_dict(self):
        d = {k: v for k, v in self.__dict__.items() if not isinstance(v, (np.ndarray,))}
        return d


class Solver:
    def __init__(self, case):
        self.c = case
        self.B = K.Background(case.background)
        self.inv = Inverter(cap_plus=case.cap_p, e_cap_plus=case.e_cap, cap_minus=case.cap_m, e_cap_minus=case.e_cm)
        self.nfev = 0
        self.diag = []
        self.last_x = None

    # ---------------------------------------------------------------- initial state
    def initial(self):
        c = self.c
        rho0 = K.rho_onset(c.R_total, c.with_V0)
        if c.start == 'alpha':
            x0 = 0.0
            nQ, rho = K.N_Q0, rho0
        else:   # 'beta': start at the thermal transition, T = T_c (relativistic gas before it: rho - V0 ∝ a^-4)
            def tc_resid(x):
                nQ = K.N_Q0 * np.exp(-3 * x); rho = (rho0 - (K.V0 if c.with_V0 else 0.0)) * np.exp(-4 * x) + (K.V0 if c.with_V0 else 0.0)
                inv = Inverter(cap_plus=c.cap_p, e_cap_plus=c.e_cap, cap_minus=c.cap_m, e_cap_minus=c.e_cm)
                P, st = equilibrium_inventory(inv, nQ / MX ** 3, rho / MX ** 4, tguess=10.0)
                return st['t'] - RATES.T_c
            x0 = brentq(tc_resid, -2.0, 0.5, xtol=1e-6)
            nQ = K.N_Q0 * np.exp(-3 * x0)
            rho = (rho0 - (K.V0 if c.with_V0 else 0.0)) * np.exp(-4 * x0) + (K.V0 if c.with_V0 else 0.0)
        Peq, st = equilibrium_inventory(self.inv, nQ / MX ** 3, rho / MX ** 4, tguess=7.0)
        Peq *= MX ** 3
        if c.P0_mode == 'equilibrium':
            P0 = Peq
        elif c.P0_mode == 'R1':
            P0 = 1e-6 * nQ            # coherent release: no free antiparticles (a floor for ln P)
        elif c.P0_mode == 'R2':
            P0 = 2.585998632082164 * nQ   # adiabatic orbit count [kin_source R2]
        else:
            P0 = float(c.P0_mode) * nQ
        a3 = np.exp(3 * x0)
        y0 = np.array([np.log(P0 * a3), rho * a3, (nQ + P0) * a3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.x0, self.P_eq0, self.rho0, self.nQ0 = x0, Peq, rho, nQ
        self.inv.last = dict(t=st['t'], up=st['up'], um=st['um'], n_c=st['n_c'], n_cb=st['n_cb'])
        return x0, y0

    # ---------------------------------------------------------------- state and rates at (x, y)
    def state(self, x, y):
        a3 = np.exp(3 * x)
        nQ = K.N_Q0 / a3
        P = np.exp(y[0]) / a3
        rho = y[1] / a3
        st = self.inv.invert(nQ / MX ** 3, P / MX ** 3, rho / MX ** 4)
        st['nQ'] = nQ; st['P'] = P; st['rho_tot'] = rho; st['a3'] = a3
        return st

    def rates(self, st):
        """All rate densities (GeV^4) and their stoichiometry; returns dict of channel -> (gamma_net, dP, dN_+^tot)."""
        c = self.c
        t, up, um, us = st['t'], st['up'], st['um'], st['us']
        R = {}
        def net(logg,ef,er):
            if logg is None:return 0.
            if ef>=er:return np.exp(logg+ef)*(-np.expm1(er-ef))
            return np.exp(logg+er)*np.expm1(ef-er)
        # ss <-> X Xb
        lg=RATES.log_gamma_ss(t)+np.log(RATES.bose_ss(t) if c.bose else 1.)
        R['ss_XXb'] = (net(lg,2*us/t,us/t),-1,+1)
        # the above: gamma e^{2 mu_s/T} - gamma e^{(mu_+ + mu_-)/T} with mu_s = mu_+ + mu_- -> e^{2u_s}(1 - e^{-u_s}); written stably
        # n <-> 2 remainders
        if c.n2:
            expanded=[]
            for ch,d in N2.items():
                expanded.append((ch,ch,d))
                if ch in ('B_XXXXb_XX','E2_XXs_XX','G2_Xss_Xs','J_XXXbs_Xs'):
                    cc=dict(d);cc['fwd']=(d['fwd'][1],d['fwd'][0],d['fwd'][2]);cc['rev']=(d['rev'][1],d['rev'][0],d['rev'][2])
                    expanded.append((ch+'_cc',ch,cc))
            for name,ch,d in expanded:
                logg0 = RATES.log_gamma_n2(ch, t)
                if logg0 is None:
                    continue
                ef = (d['fwd'][0] * up + d['fwd'][1] * um + d['fwd'][2] * us) / t
                er = (d['rev'][0] * up + d['rev'][1] * um + d['rev'][2] * us) / t
                m = max(ef, er)
                if m < -700:
                    continue
                R[name] = (net(logg0,ef,er), d['dP'], d['dNp'])
        # capture on objects (rules b, c)
        if c.capture and c.rule in ('b', 'c') and st['n_c'] > 0:
            sv, Ein = capture_kernel(t, c.capture)
            if c.capture_Q_scaling and c.Q_L != K.Q_REF:
                sv *= (c.Q_L / K.Q_REF) ** (2.0 / 3.0)
            n_L = st['n_c'] * MX ** 3 / c.Q_L
            n_m = st['n_minus'] * MX ** 3
            E_out = Ein + 2 * c.omega_L
            xo = (E_out - up) / t
            bose_out = 1.0 + (1.0 / np.expm1(xo) if xo < 600 else 0.0)
            ex = (up - um - 2 * c.omega_L) / t
            net = 1.0 - np.exp(ex) if ex < 700 else -np.inf
            g = n_L * n_m * sv * (bose_out if c.bose else 1.0) * net
            R['capture'] = (g, -1, -1)
            self._cap_energy = g * E_out * MX           # GeV^5: energy carried by the emitted X per volume per time
        else:
            self._cap_energy = 0.0
        # condensate-catalysed 4 -> 2 (rule a; homogeneous condensate at mu = m)
        if c.rule == 'a' and c.catalysis and st['n_c'] > 0:
            n_g = st['n_plus'] * MX ** 3; n_c = st['n_c'] * MX ** 3; n_m = st['n_minus'] * MX ** 3
            fac = 1.0 - np.exp(-(up + um) / t) if (up + um) / t > -700 else 1.0
            rA = RATES.A_eff_ratio('A_XXXbXb_XXb', t); rB = RATES.A_eff_ratio('B_XXXXb_XX', t)
            gA = RATES.A_old * rA / 4 * ((n_g + n_c) ** 2 - n_g ** 2) * n_m ** 2 * fac
            gB = RATES.B_old * rB / 6 * ((n_g + n_c) ** 3 - n_g ** 3) * n_m * fac
            R['A_condensate'] = (gA, -1, -1)
            R['B_condensate'] = (gB, -1, -1)
        return R

    def rhs(self, x, y):
        self.nfev += 1
        st = self.state(x, y)
        R = self.rates(st)
        H = self.B.H_of_x(x)
        a3 = st['a3']
        dP = sum(g * dp for g, dp, dn in R.values())
        # The redundant ledger is I_plus=N_plus(total)+N_s=n_Q+P.
        # Fast s closure eliminates its separate decay term; its change
        # per reaction is dP, not the gas-X multiplicity alone.
        dNp = sum(g * dp for g, dp, dn in R.values())
        P_gas = st['P_gas'] * MX ** 4
        dy = np.zeros_like(y)
        dy[0] = dP / (H * st['P'])
        dy[1] = -3 * P_gas * a3
        dy[2] = a3 * dNp / H
        dy[3] = 3 * P_gas * a3
        dy[4] = a3 * self._cap_energy / H
        dy[5] = a3 * (R['capture'][0] if 'capture' in R else 0.0) / H
        dy[6] = a3 * R['ss_XXb'][0] / H
        dy[7] = a3 * sum(g for ch, (g, dp, dn) in R.items() if ch in N2 or ch.endswith('_cc')) / H
        dy[8] = a3 * sum(g for ch, (g, dp, dn) in R.items() if ch.endswith('condensate')) / H
        return dy

    # ---------------------------------------------------------------- diagnostics along the solution
    def diagnostics(self, x, y):
        st = self.state(x, y)
        R = self.rates(st)
        H = self.B.H_of_x(x)
        t = st['t']
        nQ = st['nQ']
        P = st['P']
        rate_P = sum(abs(g * dp) for g, dp, dn in R.values())
        # closures: s decay (net Weldon width) and elastic kinetic closure
        gam_s_net = K.GAMMA_S_VAC * RATES.mean_m_over_E(t) * (RATES.bose_dec_net(t) if self.c.bose else 1.0)
        n_tot_gas = (st['n_plus'] + st['n_minus'] + st['n_s']) * MX ** 3
        el = {}
        for ch, partner in (('XXb_XXb', st['n_plus']), ('XX_XX', st['n_plus']), ('Xs_Xs', st['n_s']), ('ss_ss', st['n_s'])):
            f, lo, hi = RATES.elastic[ch]
            sv = float(np.exp(f(np.clip(np.log(t), np.log(lo), np.log(hi)))))
            el[ch] = sv * partner * MX ** 3 / H
        out = dict(x=float(x), a_over_aref=float(np.exp(x)), T_vis_GeV=self.B.T_of_x(x), H_GeV=H,
                   t=t, T_GeV=t * MX, mu_plus=st['up'], mu_minus=st['um'], mu_s=st['us'],
                   n_Q=nQ, P=P, P_over_nQ=P / nQ,
                   Xbar_free_per_net=st['n_minus'] * MX ** 3 / nQ, s_per_net=st['n_s'] * MX ** 3 / nQ,
                   anti_condensed_per_net=st['n_cb'] * MX ** 3 / nQ, condensed_fraction=st['n_c'] * MX ** 3 / nQ,
                   X_gas_per_net=st['n_plus'] * MX ** 3 / nQ,
                   rho_tot_GeV4=st['rho_tot'], rho_gas_GeV4=st['rho_gas'] * MX ** 4, P_gas_GeV4=st['P_gas'] * MX ** 4,
                   scalar_energy_fraction_of_gas=st['rho_s'] / st['rho_gas'] if st['rho_gas'] > 0 else 0.0,
                   energy_per_net_mX=st['rho_tot'] / (nQ * MX),
                   f_bench=K.f_bench(P / nQ), f_old=K.f_old(P / nQ),
                   rates_over_H={ch: g / H / max(P, 1e-300) * 0 + g / (H * max(nQ, 1e-300)) for ch, (g, dp, dn) in R.items()},
                   P_destruction_rate_over_H=rate_P / (H * P) if P > 0 else 0.0,
                   Gamma_s_net_over_H=gam_s_net / H, elastic_over_H=el,
                   active_cap=bool(st['active_p']), inversion_iterations=st['iterations'])
        out['channel_share'] = {ch: abs(g * dp) / rate_P if rate_P > 0 else 0.0 for ch, (g, dp, dn) in R.items()}
        return out

    # ---------------------------------------------------------------- run
    def run(self, verbose=False):
        c = self.c
        x0, y0 = self.initial()
        t_start = time.time()

        def ev_T(x, y):
            st = self.state(x, y)
            return st['t'] - c.t_stop
        ev_T.terminal = True; ev_T.direction = -1

        def ev_P(x, y):
            return y[0] - np.log(K.N_Q0 * c.Pfloor)
        ev_P.terminal = True; ev_P.direction = -1

        sparsity = np.zeros((9, 9), int); sparsity[:, :2] = 1
        sol = solve_ivp(self.rhs, (x0, c.x_end), y0, method='Radau', rtol=c.rtol, atol=c.atol, dense_output=True,
                        events=[ev_T, ev_P], jac_sparsity=sparsity, max_step=0.25, first_step=1e-6)
        if not sol.success:
            raise RuntimeError('solve_ivp failed: %s' % sol.message)
        xs = sol.t
        # diagnostics on a log-spaced sample plus the end point
        xd = np.unique(np.concatenate([np.linspace(xs[0], xs[-1], 161), xs[-1:]]))
        traj = []
        for x in xd:
            try:
                traj.append(self.diagnostics(x, sol.sol(x)))
            except Exception as e:   # pragma: no cover
                traj.append(dict(x=float(x), error=str(e)))
        fin = traj[-1]
        yF = sol.y[:, -1]
        stF = self.state(xs[-1], yF)
        a3 = stF['a3']
        # checks
        E0 = y0[1]
        energy_resid = (yF[1] + yF[3] - E0) / E0
        # state energy vs integrated energy
        rho_state = (stF['rho_gas'] + c.e_cap * stF['n_c'] + c.e_cm * stF['n_cb']) * MX ** 4
        energy_state_resid = (rho_state - yF[1] / a3) / (yF[1] / a3)
        charge_resid = (yF[2] / a3 - (stF['nQ'] + (stF['n_minus'] + stF['n_s'] + stF['n_cb']) * MX ** 3)) / stF['nQ']
        # freeze-out: last x where the P-destruction rate / H crossed 1
        xf, tf = None, None
        for i in range(1, len(traj)):
            if 'P_destruction_rate_over_H' in traj[i] and traj[i - 1].get('P_destruction_rate_over_H', 0) > 1 >= traj[i]['P_destruction_rate_over_H']:
                xf, tf = traj[i]['x'], traj[i]['t']
        pairs = fin['P_over_nQ']
        res = dict(case=c.to_dict(), x0=x0, x_end=float(xs[-1]), stop_reason=('T<t_stop' if sol.t_events[0].size else ('P floor' if sol.t_events[1].size else 'x_end')),
                   initial=dict(nQ_GeV3=self.nQ0, rho_GeV4=self.rho0, P0_GeV3=float(np.exp(y0[0]) / np.exp(3 * x0)), P_eq_GeV3=self.P_eq0,
                                energy_per_net_mX=self.rho0 / (self.nQ0 * MX), T_over_m=traj[0].get('t'), condensed_fraction=traj[0].get('condensed_fraction')),
                   final=dict(pairs_per_net=pairs, Xbar_free_per_net=fin['Xbar_free_per_net'], s_per_net=fin['s_per_net'],
                              anti_condensed_per_net=fin['anti_condensed_per_net'], condensed_fraction=fin['condensed_fraction'],
                              f_sym_fixed_DM=K.f_old(pairs),f_reference_Q200000=K.f_bench(pairs),f_bench=K.f_bench(pairs), f_old=K.f_old(pairs), t=fin['t'], T_vis_GeV=fin['T_vis_GeV'], a_over_aref=fin['a_over_aref'],
                              energy_per_net_mX=fin['energy_per_net_mX']),
                   freeze_out=dict(x=xf, t=tf),
                   ledgers=dict(hard_X_energy_from_capture_over_E0=float(yF[4] / E0), capture_events_per_net=float(yF[5] / K.N_Q0),
                                ss_events_per_net=float(yF[6] / K.N_Q0), n2_events_per_net=float(yF[7] / K.N_Q0),
                                condensate_42_events_per_net=float(yF[8] / K.N_Q0), pdv_work_over_E0=float(yF[3] / E0)),
                   checks=dict(energy_integration_resid=float(energy_resid), energy_state_vs_integrated_resid=float(energy_state_resid),
                               charge_resid=float(charge_resid), min_P_over_nQ=float(min(d.get('P_over_nQ', np.inf) for d in traj)),
                               min_Gamma_s_over_H=float(min(d.get('Gamma_s_net_over_H', np.inf) for d in traj)),
                               min_elastic_over_H=float(min(min(d['elastic_over_H'].values()) for d in traj if 'elastic_over_H' in d)),
                               nonnegative=bool(all(d.get('P_over_nQ', 1) >= 0 and d.get('condensed_fraction', 0) >= -1e-12 for d in traj)),
                               nfev=self.nfev, nsteps=len(xs), seconds=time.time() - t_start),
                   trajectory=traj)
        if verbose:
            print('%-28s R=%7.3f rule=%s: pairs/net=%.4e f_bench=%.4e f_old=%.4e  end t=%.4f x=%.2f (%s)  freeze x=%s  cond=%.3f  dE=%.1e dQ=%.1e  %d fev %.0fs' % (
                c.label, c.R_total, c.rule, pairs, K.f_bench(pairs), K.f_old(pairs), fin['t'], xs[-1], res['stop_reason'], xf, fin['condensed_fraction'],
                energy_resid, charge_resid, self.nfev, res['checks']['seconds']))
        return res


def run_case(case, verbose=True):
    return Solver(case).run(verbose=verbose)


if __name__ == '__main__':
    raise SystemExit('Legacy partial solver is not the production method. Use production.py for explicitly conditional controls.')
