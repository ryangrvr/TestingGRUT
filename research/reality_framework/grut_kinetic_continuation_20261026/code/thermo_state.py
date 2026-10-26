# Continued by Codex from Claude partial code at e71730c24fe6a0b28ddd72cd30982c5b39c9aedf.
# Modifications and limitations are recorded in this packet.
#!/usr/bin/env python3
"""kin_solve (Claude, own code): the thermodynamic inversion the moments solver calls at every step.

Given the physical densities of the hidden sector
    n_Q     = n_+^gas + n_c - n_-^gas - n_cb   (net charge; charge conservation is built in)
    P       = n_-^gas + n_s + n_cb             (pair inventory: free Xbar, s quanta, anti-condensed charge)
    rho_tot = rho_gas + e_cap n_c + e_cb n_cb
find (T, mu_+, mu_-) with the SUPPLIED fast-s closure mu_s = mu_+ + mu_-
and the caps mu_+ <= cap_+ (objects at omega, or a Bose saturation regulator for the gas-only and
fragmentation rules), mu_- <= cap_- (anti-objects).  Active-set Newton with the analytic Bose susceptibilities of
kin_thermo/code/thermo.py (reuse stated; its series was checked there against an independent quadrature to 2e-15).
Units: m_X = 1 (densities in m_X^3, energies in m_X^4) inside; callers convert.
"""
import numpy as np

import ks_common as K
import thermo as TH          # kin_thermo state functions (Bose series + quadrature)
TH.KMAX=20000  # switch to the existing direct quadrature before million-term arrays

RS = K.RS


class SpeciesSat:
    """One X degree at the saturation cap: tabulated quadrature (ks_common.saturated_table)."""

    def __init__(self):
        S = K.saturated_table()
        self.lt = np.log(S['t'])
        self.tab = {k: S[k] for k in ('n', 'P', 'rho', 'dn_dt', 'drho_dt')}

    def __call__(self, t):
        lt = np.log(t)
        out = {}
        for k, v in self.tab.items():
            # n, P, rho positive: log interpolation; derivatives positive too
            out[k] = float(np.exp(np.interp(lt, self.lt, np.log(np.maximum(v, 1e-300)))))
        out['dn_du'] = np.inf; out['drho_du'] = np.inf
        out['s'] = (out['rho'] + out['P'] - (1.0 - 1e-7) * out['n']) / t
        return out


SAT = None


def species(t, r, u):
    """Bose state functions with derivatives; the saturated X table when u is at the saturation cap."""
    global SAT
    if r == 1.0 and u >= 1.0 - 2e-7:
        if SAT is None:
            SAT = SpeciesSat()
        return SAT(t)
    return TH.species(t, r, u, derivs=True)


def initial_equilibrium(nq,rho,cap,e_cap):
    """Match the initial state using the same densities as the ODE.

    TH.net_charge uses a truncated series even when TH.species switches to
    quadrature. With the reduced KMAX that is inaccurate near the gas cap.
    Keep the pinned input untouched and perform this local matching instead.
    """
    from scipy.optimize import brentq
    def at(t):
        def net(u):return species(t,1.,u)['n']-species(t,1.,-u)['n']
        if net(cap)<=nq:
            u=cap;nc=nq-net(cap)
        else:
            u=brentq(lambda z:net(z)-nq,0.,cap,xtol=1e-15,rtol=1e-14);nc=0.
        xp=species(t,1.,u);xm=species(t,1.,-u);s=species(t,RS,0.)
        return dict(t=t,up=u,um=-u,n_cond=nc,rho_total=xp['rho']+xm['rho']+s['rho']+e_cap*nc)
    lt=brentq(lambda z:at(np.exp(z))['rho_total']-rho,np.log(1e-4),np.log(3e3),xtol=1e-13,rtol=1e-13)
    return at(np.exp(lt))


class Inverter:
    def __init__(self, cap_plus=K.OMEGA_REF, e_cap_plus=None, cap_minus=K.OMEGA_REF, e_cap_minus=None, with_s=True):
        self.cap_p = cap_plus
        self.e_p = cap_plus if e_cap_plus is None else e_cap_plus
        self.cap_m = cap_minus
        self.e_m = cap_minus if e_cap_minus is None else e_cap_minus
        self.with_s = with_s
        self.last = None

    # ---- residuals and Jacobian for a given active set
    def _eval(self, lt, up, um, nc, ncb, act_p, act_m, targets):
        t = np.exp(lt)
        X = species(t, 1.0, up); Xb = species(t, 1.0, um)
        us = up + um
        if self.with_s:
            S = species(t, RS, us) if us < RS - 1e-9 else species(t, RS, RS - 1e-9)
        else:
            S = dict(n=0.0, P=0.0, rho=0.0, s=0.0, dn_du=0.0, drho_du=0.0, dn_dt=0.0, drho_dt=0.0)
        n_Q, P_tot, rho_tot = targets
        F = np.array([X['n'] + nc - Xb['n'] - ncb - n_Q,
                      Xb['n'] + S['n'] + ncb - P_tot,
                      X['rho'] + Xb['rho'] + S['rho'] + self.e_p * nc + self.e_m * ncb - rho_tot])
        # columns: d/dln t, d/d(second), d/d(third) where second = n_c if act_p else u_+; third = n_cb if act_m else u_-
        J = np.zeros((3, 3))
        J[0, 0] = t * (X['dn_dt'] - Xb['dn_dt'])
        J[1, 0] = t * (Xb['dn_dt'] + S['dn_dt'])
        J[2, 0] = t * (X['drho_dt'] + Xb['drho_dt'] + S['drho_dt'])
        if act_p:
            J[0, 1] = 1.0; J[1, 1] = 0.0; J[2, 1] = self.e_p
        else:
            J[0, 1] = X['dn_du']; J[1, 1] = S['dn_du']; J[2, 1] = X['drho_du'] + S['drho_du']
        if act_m:
            J[0, 2] = -1.0; J[1, 2] = 1.0; J[2, 2] = self.e_m
        else:
            J[0, 2] = -Xb['dn_du']; J[1, 2] = Xb['dn_du'] + S['dn_du']; J[2, 2] = Xb['drho_du'] + S['drho_du']
        state = dict(t=t, up=up, um=um, us=us, X=X, Xb=Xb, S=S, n_c=nc, n_cb=ncb)
        return F, J, state

    def invert(self, n_Q, P_tot, rho_tot, guess=None, tol=1e-11, itmax=120):
        """Returns the state dict (m_X units) with keys t, up, um, us, n_plus, n_minus, n_s, n_c, n_cb, rho_gas,
        P_gas, s_gas, active_p, active_m, iterations."""
        if guess is None:
            guess = self.last if self.last is not None else dict(t=3.0, up=0.3, um=-0.3, n_c=0.0, n_cb=0.0)
        lt = np.log(guess['t']); up = min(guess['up'], self.cap_p); um = min(guess['um'], self.cap_m)
        nc = max(guess.get('n_c', 0.0), 0.0); ncb = max(guess.get('n_cb', 0.0), 0.0)
        targets = (n_Q, P_tot, rho_tot)
        scale = np.array([max(abs(n_Q), 1e-300), max(P_tot, 1e-300), max(rho_tot, 1e-300)])
        # initial active set from the capacity at the guess temperature
        t0 = np.exp(lt)
        Xb0 = species(t0, 1.0, um)
        cap_n = species(t0, 1.0, self.cap_p)['n']
        act_p = (nc > 0) or (n_Q + Xb0['n'] + ncb > cap_n)
        if act_p:
            up = self.cap_p; nc = max(n_Q + Xb0['n'] + ncb - cap_n, 0.0)
        act_m = ncb > 0
        if act_m:
            um = self.cap_m
        switches = 0
        converged = False
        for it in range(itmax):
            F, J, st = self._eval(lt, up, um, nc, ncb, act_p, act_m, targets)
            if not np.all(np.isfinite(F)):
                raise RuntimeError('non-finite residual in inversion at t=%.3g up=%.3g um=%.3g' % (np.exp(lt), up, um))
            if np.max(np.abs(F) / scale) < tol:
                # consistency of the active set at convergence
                if act_p and nc < 0:
                    act_p = False; nc = 0.0; up = self.cap_p - 1e-4; switches += 1; continue
                if act_m and ncb < 0:
                    act_m = False; ncb = 0.0; um = self.cap_m - 1e-4; switches += 1; continue
                converged = True
                break
            try:
                dx = np.linalg.solve(J, -F)
            except np.linalg.LinAlgError:
                dx = np.linalg.lstsq(J, -F, rcond=None)[0]
            if not np.all(np.isfinite(dx)):
                dx = np.sign(-F) * 1e-3
            t_now = np.exp(lt)
            cap_step = min(1.0, 0.7 / max(abs(dx[0]), 1e-300))
            if not act_p:
                cap_step = min(cap_step, max(0.5, t_now) / max(abs(dx[1]), 1e-300))
            if not act_m:
                cap_step = min(cap_step, max(0.5, t_now) / max(abs(dx[2]), 1e-300))
            dx = dx * cap_step
            # condensates may not step below zero: stop at the floor; deactivate only if the floor is pushed again
            floor_p = floor_m = False
            if act_p and dx[1] < 0 and nc + dx[1] < 0:
                if nc <= 0:
                    floor_p = True
                else:
                    dx = dx * (nc / (-dx[1]))
            if act_m and dx[2] < 0 and ncb + dx[2] < 0:
                if ncb <= 0:
                    floor_m = True
                else:
                    dx = dx * (ncb / (-dx[2]))
            if floor_p:
                act_p = False; nc = 0.0; up = self.cap_p - 1e-4; switches += 1; continue
            if floor_m:
                act_m = False; ncb = 0.0; um = self.cap_m - 1e-4; switches += 1; continue
            lam = 1.0
            for _ in range(40):
                nlt = lt + lam * dx[0]
                nup = up if act_p else up + lam * dx[1]
                nnc = nc + lam * dx[1] if act_p else nc
                num = um if act_m else um + lam * dx[2]
                nncb = ncb + lam * dx[2] if act_m else ncb
                if (nup < 1.0 - 1e-9) and (num < 1.0 - 1e-9) and (nup + num < RS - 1e-9):
                    break
                lam *= 0.5
            lt, up, um, nc, ncb = nlt, nup, num, nnc, nncb
            # on-the-fly switching
            if (not act_p) and up > self.cap_p:
                act_p = True; up = self.cap_p
                tn = np.exp(lt)
                nc = max(n_Q + species(tn, 1.0, um)['n'] + ncb - species(tn, 1.0, self.cap_p)['n'], 0.0)
                switches += 1
            if (not act_m) and um > self.cap_m:
                act_m = True; um = self.cap_m
                tn = np.exp(lt)
                ncb = max(P_tot - species(tn, 1.0, self.cap_m)['n'] - species(tn, RS, min(up + um, RS - 1e-9))['n'], 0.0)
                switches += 1
        F, J, st = self._eval(lt, up, um, nc, ncb, act_p, act_m, targets)
        resid = float(np.max(np.abs(F) / scale))
        if resid > 1e-6 or nc < -1e-12 * max(n_Q, 1e-300):
            raise RuntimeError('inversion did not converge: resid=%.2e targets=%r state t=%.4g up=%.4g um=%.4g nc=%.3g it=%d' % (
                resid, targets, st['t'], up, um, nc, it))
        X, Xb, S = st['X'], st['Xb'], st['S']
        out = dict(t=st['t'], up=up, um=um, us=st['us'], n_plus=X['n'], n_minus=Xb['n'], n_s=S['n'], n_c=max(nc, 0.0), n_cb=max(ncb, 0.0),
                   rho_gas=X['rho'] + Xb['rho'] + S['rho'], P_gas=X['P'] + Xb['P'] + S['P'],
                   rho_X=X['rho'], rho_Xb=Xb['rho'], rho_s=S['rho'],
                   s_gas=(X['rho'] + Xb['rho'] + S['rho'] + X['P'] + Xb['P'] + S['P'] - up * X['n'] - um * Xb['n'] - st['us'] * S['n']) / st['t'],
                   active_p=act_p, active_m=act_m, iterations=it + 1, switches=switches, residual=resid)
        self.last = dict(t=out['t'], up=up, um=um, n_c=max(nc, 0.0), n_cb=max(ncb, 0.0))
        return out


def equilibrium_inventory(inv, n_Q, rho_tot, tguess=5.0):
    """Full chemical equilibrium (mu_- = -mu_+): the pair inventory P such that the inversion gives up + um = 0.
    Used for the forced-thermal initial state and for the stationarity check.  (up + um is monotone decreasing
    in P at fixed n_Q and rho: fewer pairs means a more negative mu_-.)"""
    from scipy.optimize import brentq

    def resid(lP):
        st = inv.invert(n_Q, np.exp(lP), rho_tot, guess=dict(t=tguess, up=0.3, um=-0.3, n_c=0.0, n_cb=0.0))
        return st['up'] + st['um']
    lo, hi = np.log(n_Q * 1e-4), np.log(n_Q * 30.0)
    while resid(lo) > 0:                 # small P: mu_- very negative, up + um < 0
        lo -= 2
    while resid(hi) < 0:                 # large P: Xbar-rich, up + um > 0
        hi += 1
    lP = brentq(resid, lo, hi, xtol=1e-13, rtol=1e-13)
    st = inv.invert(n_Q, np.exp(lP), rho_tot, guess=dict(t=tguess, up=0.3, um=-0.3, n_c=0.0, n_cb=0.0))
    return np.exp(lP), st


if __name__ == '__main__':
    import time
    MX = K.MX
    nq = K.N_Q0 / MX ** 3
    rho = K.rho_onset(40.56262288525937) / MX ** 4
    inv = Inverter(cap_plus=K.OMEGA_REF, cap_minus=K.OMEGA_REF)
    # compare with kin_thermo state_two_potential at Codex's with-s state (P = 0.06850712 GeV^3)
    P = 0.06850712020373081 / MX ** 3
    t0 = time.time()
    st = inv.invert(nq, P, rho, guess=dict(t=7.0, up=0.5, um=-0.5, n_c=0.0, n_cb=0.0))
    dt = time.time() - t0
    print('own inversion: t=%.9f up=%.6f um=%.6f n_c=%.9e GeV^3  (%d it, %d switches, %.1f ms)' % (
        st['t'], st['up'], st['um'], st['n_c'] * MX ** 3, st['iterations'], st['switches'], dt * 1e3))
    ref = TH.state_two_potential(nq, P, rho, K.OMEGA_REF, tguess=7.3)
    print('kin_thermo     : t=%.9f up=%.6f um=%.6f n_c=%.9e GeV^3' % (ref['t'], ref['up'], ref['um'], ref['n_cond'] * MX ** 3))
    print('Codex control  : t=7.343904195 n_c=7.072206e-03 (equilibrium P)')
    # warm-start timing
    t0 = time.time()
    for i in range(50):
        st2 = inv.invert(nq * (1 + 1e-4 * i), P * (1 - 1e-4 * i), rho * (1 - 2e-4 * i))
    print('warm-start inversion: %.2f ms each, %d it last' % ((time.time() - t0) / 50 * 1e3, st2['iterations']))
    # full equilibrium inventory
    Peq, steq = equilibrium_inventory(inv, nq, rho, tguess=7.3)
    print('equilibrium P = %.9e GeV^3 (Codex 6.850712e-02), t=%.9f, n_c=%.9e, up+um=%.1e' % (Peq * MX ** 3, steq['t'], steq['n_c'] * MX ** 3, steq['up'] + steq['um']))
    # saturation-cap inverter (gas_only), undersaturated and saturated cases
    inv2 = Inverter(cap_plus=1.0 - 1e-7, e_cap_plus=1.0, cap_minus=1.0 - 1e-7, e_cap_minus=1.0)
    st = inv2.invert(nq, P, rho, guess=dict(t=7.0, up=0.5, um=-0.5, n_c=0.0, n_cb=0.0))
    print('saturation cap at onset: t=%.6f up=%.7f um=%.6f n_c=%.5e active=%s' % (st['t'], st['up'], st['um'], st['n_c'] * MX ** 3, st['active_p']))
    # low temperature, NR: n_tot small
    st = inv.invert(1e-4 / MX ** 3, 1e-9 / MX ** 3, (1.0e-4 * 0.95 + 1e-9 * 1.5) * MX / MX ** 4, guess=dict(t=0.2, up=0.8, um=-0.8, n_c=0.0, n_cb=0.0))
    print('   n_-=%.3e n_s=%.3e rho_gas=%.4e' % (st['n_minus'] * MX ** 3, st['n_s'] * MX ** 3, st['rho_gas'] * MX ** 4))
    # very late: tiny inventory
    st = inv.invert(1e-8 / MX ** 3, 1e-40 / MX ** 3, (1.0e-8 * 0.93) * MX / MX ** 4, guess=dict(t=0.02, up=0.86, um=-0.86, n_c=0.0, n_cb=0.0))
    print('very late state: t=%.5f up=%.6f um=%.4f n_c/n_Q=%.6f active=%s it=%d' % (st['t'], st['up'], st['um'], st['n_c'] / (1e-8 / MX ** 3), st['active_p'], st['iterations']))
    print('late NR state: t=%.5f up=%.6f um=%.6f n_c/n_tot=%.4f active=%s it=%d' % (st['t'], st['up'], st['um'], st['n_c'] / (1e-4 / MX ** 3), st['active_p'], st['iterations']))
