#!/usr/bin/env python3
"""inf_reheat SKEPTIC (Claude, own code, independent of inf_reheat/code): T-model benchmarks recomputed by different methods.

Differences from the track (benchmarks.py: tmodel, Nstar_formula; fluid.py):
  * n_s, r, A_s from numerical Mukhanov-Sasaki mode integration (scalar zeta and tensor h), not first-order slow roll at phi_*.
  * Background in cosmic time (variables phi, phidot, ln a), V0 = 1 units, M_P = 1.
  * N_* from an explicit e-fold bookkeeping that uses the radiation constant K = a^2 H directly (a = 1 at epsilon_H = 1),
    with T_0, g_s0, k_* = 0.05/Mpc written out, instead of the Planck 2018 X eq. 47 constant.
  * The reheating constant K on the ROLLING T-model slice: the full T-model oscillator continued from epsilon_H = 1 with the
    engine's decay term Gamma phidot (I'' + (3H + Gamma) I' + V' = 0, rho_r' = -4 H rho_r + Gamma I'^2), then the cycle-averaged
    fluid; compared with the track's static-slice K (released at rest, quadratic, K = c_K(b) Gamma^{-1/3}).
  * The first turning-point amplitude Phi_1 after epsilon_H = 1 (enters the lambda_IH window and the S1 blocking).
Usage: sk_tmodel.py -> ../results/sk_tmodel.json.  Labels: [C].
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
TRK = HERE.parent.parent / 'results'
MP_GEV = 2.435e18
AS = 2.10e-9
MPC_GEV = 1.97327e-14 / 3.08568e24          # 1/Mpc in GeV
T0_GEV = 2.7255 * 8.617333e-14
GS0 = 3.931
GSTAR = 106.75


class TModel:
    def __init__(self, alpha):
        self.alpha = alpha
        self.w = math.sqrt(6 * alpha)

    def V(self, p):
        return math.tanh(p / self.w) ** 2

    def dV(self, p):
        t = math.tanh(p / self.w)
        return 2 * t * (1 - t * t) / self.w

    def background(self, N_before=75.0):
        """cosmic-time integration from slow-roll attractor far up the plateau to epsilon_H = 1; returns dense arrays in N."""
        w = self.w
        # start: phi such that the slow-roll e-fold count to the end is ~N_before (rough), then integrate
        xs = 0.5 * math.acosh(math.cosh(2 * 0.5 * math.asinh(math.sqrt(8) / w)) + 8 * N_before / w ** 2)
        p0 = w * xs
        H0 = math.sqrt(self.V(p0) / 3)
        pd0 = -self.dV(p0) / (3 * H0)

        def rhs(t, y):
            p, pd, la = y
            H = math.sqrt((0.5 * pd * pd + self.V(p)) / 3)
            return [pd, -3 * H * pd - self.dV(p), H]

        def ev(t, y):
            p, pd, la = y
            H2 = (0.5 * pd * pd + self.V(p)) / 3
            return 0.5 * pd * pd / H2 - 1.0
        ev.terminal = True; ev.direction = 1
        s = solve_ivp(rhs, (0, 1e7), [p0, pd0, 0.0], method='LSODA', rtol=1e-12, atol=1e-14, events=ev, dense_output=True)
        te = s.t_events[0][0]
        pe, pde, lae = s.y_events[0][0]
        self.sol, self.te, self.lae = s, te, lae
        self.pe, self.pde = pe, pde
        self.He = math.sqrt((0.5 * pde * pde + self.V(pe)) / 3)
        return s

    def at_N(self, Nb):
        """state Nb e-folds before the end (a = 1 at the end)."""
        f = lambda t: self.sol.sol(t)[2] - (self.lae - Nb)
        t = brentq(f, 0.0, self.te)
        p, pd, la = self.sol.sol(t)
        H = math.sqrt((0.5 * pd * pd + self.V(p)) / 3)
        return t, p, pd, H

    def ms_spectra(self, Nstar, dlnk=(-0.4, -0.2, 0.0, 0.2, 0.4)):
        """numerical Mukhanov-Sasaki: zeta and h modes for k = k_* e^{d}, k_* = a H at Nstar before the end."""
        t_s, p_s, pd_s, H_s = self.at_N(Nstar)
        kst = math.exp(-Nstar) * H_s
        out = []
        for d in dlnk:
            k = kst * math.exp(d)
            # start where k/(aH) = 300, stop where k/(aH) = 1e-3
            def kaH(t):
                p, pd, la = self.sol.sol(t)
                H = math.sqrt((0.5 * pd * pd + self.V(p)) / 3)
                return k / (math.exp(la - self.lae) * H)
            t0 = brentq(lambda t: kaH(t) - 300.0, 0.0, t_s + 50.0 / H_s) if kaH(0.0) > 300 else 0.0
            t1 = brentq(lambda t: kaH(t) - 1e-3, t_s, self.te)

            # mode equations in cosmic time for zeta and h, complex split; background from dense output
            def bg(t):
                p, pd, la = self.sol.sol(t)
                H = math.sqrt((0.5 * pd * pd + self.V(p)) / 3)
                a = math.exp(la - self.lae)
                eps = 0.5 * pd * pd / H ** 2
                # d ln eps / dt = 2 pdd/pd - 2 Hdot/H, Hdot = -pd^2/2
                pdd = -3 * H * pd - self.dV(p)
                deps = 2 * pdd / pd + pd * pd / H
                return a, H, eps, deps

            def rhs(t, y):
                a, H, eps, deps = bg(t)
                zr, zi, vr, vi, hr, hi, ur, ui = y
                # zeta'' + (3H + dln eps/dt) zeta' + (k/a)^2 zeta = 0 ; h'' + 3H h' + (k/a)^2 h = 0
                fz = -(3 * H + deps)
                q = (k / a) ** 2
                return [vr, vi, fz * vr - q * zr, fz * vi - q * zi, ur, ui, -3 * H * ur - q * hr, -3 * H * ui - q * hi]
            a, H, eps, deps = bg(t0)
            z = a * math.sqrt(2 * eps)
            # BD: u = e^{-ik tau}/sqrt(2k); zeta = u/z; dzeta/dt = (1/a)(-ik u/z - u z'/z^2), z'/z = aH(1 + deps/(2H))
            u = 1 / math.sqrt(2 * k)
            zeta = u / z
            dzeta = (1 / a) * (-1j * k * zeta) - zeta * H * (1 + deps / (2 * H))
            h = 2 * u / a
            dh = (1 / a) * (-1j * k * h) - h * H
            y0 = [zeta.real, zeta.imag, dzeta.real, dzeta.imag, h.real, h.imag, dh.real, dh.imag]
            s = solve_ivp(rhs, (t0, t1), y0, method='DOP853', rtol=1e-10, atol=1e-30)
            zr, zi, vr, vi, hr, hi, ur, ui = s.y[:, -1]
            Pz = k ** 3 * (zr * zr + zi * zi) / (2 * math.pi ** 2)
            Pt = 2 * k ** 3 * (hr * hr + hi * hi) / (2 * math.pi ** 2)
            out.append(dict(dlnk=d, Pz=Pz, Pt=Pt))
        lk = np.array([o['dlnk'] for o in out]); lP = np.log([o['Pz'] for o in out])
        c = np.polyfit(lk, lP, 2)
        Pz0 = math.exp(np.polyval(c, 0.0))
        Pt0 = [o['Pt'] for o in out if o['dlnk'] == 0.0][0]
        lPt = np.log([o['Pt'] for o in out]); ct = np.polyfit(lk, lPt, 2)
        # slow-roll comparison (first order, as the track)
        e1 = 0.5 * (self.dV(p_s) / self.V(p_s)) ** 2
        w = self.w
        d2V = (2 / w ** 2) * (1 / math.cosh(p_s / w) ** 4 - 2 * math.tanh(p_s / w) ** 2 / math.cosh(p_s / w) ** 2)
        eta = d2V / self.V(p_s)
        V0 = AS / Pz0
        return dict(phi_star=p_s, H_star_V0=H_s, ns_MS=1 + c[1], running_MS=2 * c[0], r_MS=Pt0 / Pz0, nt_MS=ct[1],
                    ns_SR1=1 - 6 * e1 + 2 * eta, r_SR1=16 * e1, V0=V0, H_star_GeV=H_s * math.sqrt(V0) * MP_GEV,
                    H_end_GeV=self.He * math.sqrt(V0) * MP_GEV, m_GeV=math.sqrt(2 * V0 / w ** 2) * MP_GEV,
                    b=math.sqrt(2 / w ** 2) / self.He, V_star=V0 * self.V(p_s), modes=out)

    def first_turn(self):
        """amplitude at the first turning point after epsilon_H = 1 (no decay)."""
        def rhs(t, y):
            p, pd, la = y
            H = math.sqrt((0.5 * pd * pd + self.V(p)) / 3)
            return [pd, -3 * H * pd - self.dV(p), H]
        ev = lambda t, y: y[1]
        ev.terminal = True; ev.direction = 1
        s = solve_ivp(rhs, (0, 200.0 / self.He), [self.pe, self.pde, 0.0], method='DOP853', rtol=1e-11, atol=1e-14, events=ev)
        return abs(s.y_events[0][0][0]), s.t_events[0][0] * self.He, s.y_events[0][0][2]

    def K_rolling(self, Gam_over_He, t_match_He=60.0):
        """K = a^2 H / H_end in radiation, a = 1 at epsilon_H = 1: full oscillator with decay to t_match, then averaged fluid."""
        He = self.He; G = Gam_over_He * He

        def rhs(t, y):
            p, pd, la, rr = y
            H = math.sqrt((0.5 * pd * pd + self.V(p) + rr) / 3)
            return [pd, -(3 * H + G) * pd - self.dV(p), H, -4 * H * rr + G * pd * pd]
        tm = t_match_He / He
        s = solve_ivp(rhs, (0, tm), [self.pe, self.pde, 0.0, 0.0], method='DOP853', rtol=1e-11, atol=1e-16,
                      max_step=0.01 * self.w)
        p, pd, la, rr = s.y[:, -1]
        a = math.exp(la)
        rhoI = 0.5 * pd * pd + self.V(p)
        # near the minimum V ~ (1/w^2) p^2: quadratic, w_I -> 0 (anharmonic corrections ~ (Phi/w)^2, small here)
        x0, y0 = rhoI * a ** 3, rr * a ** 4

        def frhs(N, z):
            x, y = z
            aa = math.exp(N)
            H = math.sqrt(max((x / aa ** 3 + y / aa ** 4) / 3, 1e-300))
            return [-(G / H) * x, (G / H) * aa * x]

        def done(N, z):
            aa = math.exp(N)
            return z[0] / aa ** 3 - 1e-13 * z[1] / aa ** 4
        done.terminal = True
        f = solve_ivp(frhs, (la, la + 90.0), [x0, y0], method='LSODA', rtol=1e-11, atol=1e-40, events=done)
        y = f.y[1, -1]
        return math.sqrt(y / 3) / He, x0 / He ** 2


def K_static(b, Gam, t_match=60.0):
    """the track's static slice recomputed here with my own integrator: quadratic oscillator released at rest, H = 1."""
    I0 = math.sqrt(6.0) / b

    def rhs(t, y):
        p, pd, la, rr = y
        H = math.sqrt((0.5 * pd * pd + 0.5 * b * b * p * p + rr) / 3)
        return [pd, -(3 * H + Gam) * pd - b * b * p, H, -4 * H * rr + Gam * pd * pd]
    s = solve_ivp(rhs, (0, t_match), [I0, 0.0, 0.0, 0.0], method='DOP853', rtol=1e-11, atol=1e-16, max_step=0.01)
    p, pd, la, rr = s.y[:, -1]
    a = math.exp(la)
    x0, y0 = (0.5 * pd * pd + 0.5 * b * b * p * p) * a ** 3, rr * a ** 4

    def frhs(N, z):
        x, y = z
        aa = math.exp(N)
        H = math.sqrt(max((x / aa ** 3 + y / aa ** 4) / 3, 1e-300))
        return [-(Gam / H) * x, (Gam / H) * aa * x]

    def done(N, z):
        aa = math.exp(N)
        return z[0] / aa ** 3 - 1e-13 * z[1] / aa ** 4
    done.terminal = True
    f = solve_ivp(frhs, (la, la + 90.0), [x0, y0], method='LSODA', rtol=1e-11, atol=1e-40, events=done)
    return math.sqrt(f.y[1, -1] / 3), x0


def universal_cK():
    """pure-matter start (a -> 0) analytic limit: K Gamma^{1/3} = C x0^{2/3}; C from the dimensionless fluid."""
    # units Gamma = 1, x0 = 1: start deep in matter era at a = a0 with y = 0 and the matter-era radiation already sourced
    a0 = 1e-6

    def frhs(N, z):
        x, y = z
        aa = math.exp(N)
        H = math.sqrt(max((x / aa ** 3 + y / aa ** 4) / 3, 1e-300))
        return [-(1 / H) * x, (1 / H) * aa * x]
    # leading matter-era solution: H = sqrt(x/3) a^{-3/2}; y = int a x dt = int a x da/(aH) ~ (2/5) sqrt(3) a^{5/2}
    y0 = 0.4 * math.sqrt(3.0) * a0 ** 2.5
    f = solve_ivp(frhs, (math.log(a0), math.log(a0) + 60.0), [1.0, y0], method='LSODA', rtol=1e-12, atol=1e-60)
    return math.sqrt(f.y[1, -1] / 3)


def N_star_from_K(H_star_GeV, H_end_GeV, K, kpiv=0.05):
    """e-folds between the exit of k_* and epsilon_H = 1, from K = a^2 H/H_end in radiation (a = 1 at the end)."""
    Ta = (90 * MP_GEV ** 2 * K ** 2 * H_end_GeV ** 2 / (math.pi ** 2 * GSTAR)) ** 0.25   # T a in radiation (a_end = 1)
    aend_over_a0 = T0_GEV * (GS0 / GSTAR) ** (1 / 3) / Ta
    k = kpiv * MPC_GEV
    return math.log(aend_over_a0 * H_star_GeV / k)


def main():
    trk = json.load(open(TRK / 'benchmarks.json'))
    tb = {(bm['family'], bm['channel']): bm for bm in trk['benchmarks']}
    res = dict(method=__doc__, rows=[])
    # universal constant and the static-slice check against the track's c_K
    C = universal_cK()
    res['universal_C'] = C
    flu = json.load(open(TRK / 'fluid.json'))
    st = []
    for b in (2.0, 2.5, 3.0):
        for G in (1e-4, 1e-6):
            K, x0 = K_static(b, G)
            row = [g for g in flu['grid'] if g['b'] == b and abs(g['Gamma'] / G - 1) < 1e-9][0]
            st.append(dict(b=b, Gamma=G, K_mine=K, K_track=row['K'], rel=K / row['K'] - 1, x0=x0, C_x0_23=C * x0 ** (2 / 3) * G ** (-1 / 3)))
            print('static b=%.1f G=%.0e K mine %.5f track %.5f rel %.2e  x0 %.4f  C x0^(2/3) G^(-1/3) = %.5f' % (
                b, G, K, row['K'], K / row['K'] - 1, x0, C * x0 ** (2 / 3) * G ** (-1 / 3)), flush=True)
    res['static_check'] = st
    for alpha in (0.3, 1.0, 3.0):
        fam = tb[('T(alpha=%g)' % alpha, 'S1')]
        s = fam['solution']; f = fam['fam']
        tm = TModel(alpha); tm.background()
        sp = tm.ms_spectra(s['N_star'])
        Phi1, tturn, laturn = tm.first_turn()
        Gh = s['Gamma_over_Hi']
        Kr, x0r = tm.K_rolling(Gh)
        Ks = s['K']
        Nmine = N_star_from_K(sp['H_star_GeV'], sp['H_end_GeV'], Ks)
        Nmine_roll = N_star_from_K(sp['H_star_GeV'], sp['H_end_GeV'], Kr)
        # shift of the self-consistent Gamma if the rolling K replaces the static one: H_req ~ K at fixed N_X, K ~ Gamma^{-1/3}
        G_roll = Gh * (Kr / Ks) ** 3
        row = dict(alpha=alpha, N_star_track=s['N_star'], N_star_mine_from_trackK=Nmine, N_star_mine_from_rollingK=Nmine_roll,
                   ns_track=f['ns'], ns_MS=sp['ns_MS'], ns_SR1=sp['ns_SR1'], running_MS=sp['running_MS'], r_track=f['r'], r_MS=sp['r_MS'],
                   nt_MS=sp['nt_MS'], H_end_track=f['H_end'], H_end_mine=sp['H_end_GeV'], H_star_track=f['H_star'], H_star_mine=sp['H_star_GeV'],
                   m_track=f['m'], m_mine=sp['m_GeV'], b_track=f['b'], b_mine=sp['b'], phi_end_track=f['phi_end'], phi_end_mine=tm.pe,
                   Phi1_track=f['Phi1'], Phi1_mine=Phi1, t_turn_He=tturn, Gamma_over_Hi=Gh, K_static_track=Ks, K_rolling=Kr,
                   K_ratio=Kr / Ks, x0_rolling_over_He2=x0r, Gamma_selfconsistent_if_rolling=G_roll, Gamma_shift=G_roll / Gh)
        res['rows'].append(row)
        print(json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in row.items()}), flush=True)
    json.dump(res, open(RES / 'sk_tmodel.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
