"""Item 3c of track higgs_reheat: thermal effects once the (3, 0.5) reheating radiation appears.  Own code.
Writes ../results/thermal.json.

[C] pieces:
  * thermal-mass coefficient kappa_T = m_T^2/T^2 = 3g^2/16 + g'^2/16 + yt^2/4 + lam/2 (h << T), couplings at mu = pi T;
  * the one-loop finite-T potential V(h,T) = V0(h) + T^4/(2 pi^2) [6 J_B(mW^2/T^2) + 3 J_B(mZ^2/T^2) + 12 J_F(mt^2/T^2)]
    - (same at h = 0), J_B(y) = int x^2 ln(1 - e^{-sqrt(x^2+y)}), J_F(y) = -int x^2 ln(1 + e^{-sqrt(x^2+y)}),
    V0 = lam_eff(max(h, T)) h^4/4 (flat 2-loop lambda_eff of smpot; using lam_eff at T below h = T is conservative,
    since lam_eff(T) < lam_eff(h) there).  Higgs/Goldstone loops and daisy resummation omitted (stated; both small:
    lam/2 ~ -0.006 against kappa_T ~ 0.11, and the daisy correction lowers the gauge-boson thermal contribution by
    O(g T) terms) [E].
  * the thermal barrier h_bar(T) and the O(3) bounce action S3/T (overshoot/undershoot shooting), the
    high-temperature estimate S3/T ~ 18.9 sqrt(kappa_T)/|lambda| for comparison [Q: the O(3) bounce of
    m^2 h^2/2 - |lambda| h^4/4 has S3 = 18.9 m/|lambda| (Arnold-Espinosa hep-ph/9212235 normalisation) [R]];
  * expected number of thermally nucleated bubbles in our past light cone along the engine's actual T(tau):
    N_b = e^{3N} int dtau a(tau)^3 Gamma(T)/H_i^4, Gamma = T^4 (S3/2 pi T)^{3/2} e^{-S3/T}.
"""
import json
import os

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import interp1d

import smpot as S
from background import T_GeV

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')


def _J(y, sgn):
    f = (lambda x: x * x * np.log1p(-np.exp(-np.sqrt(x * x + y)))) if sgn > 0 else \
        (lambda x: -x * x * np.log1p(np.exp(-np.sqrt(x * x + y))))
    return quad(f, 0, 60 + np.sqrt(y), limit=200)[0]


_ly = np.linspace(-8, np.log(3000.0), 500)
_JB = interp1d(_ly, [_J(np.exp(v), +1) for v in _ly], bounds_error=False, fill_value=(_J(0.0, 1), 0.0))
_JF = interp1d(_ly, [_J(np.exp(v), -1) for v in _ly], bounds_error=False, fill_value=(_J(0.0, -1), 0.0))
JB0, JF0 = _J(0.0, 1), _J(0.0, -1)


def JB(y):
    y = np.asarray(y, float)
    return np.where(y < np.exp(-8), JB0 + np.pi ** 2 * y / 12, _JB(np.log(np.maximum(y, 1e-300))))


def JF(y):
    y = np.asarray(y, float)
    return np.where(y < np.exp(-8), JF0 + np.pi ** 2 * y / 24, _JF(np.log(np.maximum(y, 1e-300))))


class ThermalV:
    def __init__(self, sm, T):
        self.sm, self.T = sm, T
        c = sm.cpl(np.pi * T)
        self.g2, self.gp2, self.yt2 = c['G2'], 0.6 * c['G1'], c['T']
        self.kT = 3 * c['G2'] / 16 + 0.6 * c['G1'] / 16 + c['T'] / 4 + c['lam'] / 2
        hs = np.exp(np.linspace(np.log(1e8), np.log(1e19), 600))
        self.lam = interp1d(np.log(hs), sm.lam_eff_flat(hs), bounds_error=False, fill_value='extrapolate')

    def V(self, u):
        """V/T^4 at h = u T."""
        u = np.asarray(u, float)
        h = u * self.T
        lam = self.lam(np.log(np.maximum(h, self.T)))
        V0 = lam * u ** 4 / 4
        yW, yZ, yt = self.g2 * u ** 2 / 4, (self.g2 + self.gp2) * u ** 2 / 4, self.yt2 * u ** 2 / 2
        VT = (6 * (JB(yW) - JB0) + 3 * (JB(yZ) - JB0) + 12 * (JF(yt) - JF0)) / (2 * np.pi ** 2)
        return V0 + VT

    def dV(self, u, d=1e-5):
        return (self.V(u * (1 + d)) - self.V(u * (1 - d))) / (2 * u * d)

    def splined(self, umax=400.0, n=6000):
        """fast copy: V and dV/du from a cubic spline of V/T^4 on [0, umax] (used by the shooting)."""
        from scipy.interpolate import CubicSpline
        us = np.concatenate([[0.0], np.exp(np.linspace(np.log(1e-4), np.log(umax), n))])
        cs = CubicSpline(us, self.V(np.maximum(us, 1e-12)))
        d1 = cs.derivative()

        class _F:
            V = staticmethod(lambda u: cs(np.asarray(u, float)))
            dV = staticmethod(lambda u, d=None: float(d1(u)))
        return _F()


def bounce(tv, umax=400.0, rmax=200.0):
    us = np.exp(np.linspace(np.log(1e-3), np.log(umax), 4000))
    V = tv.V(us)
    k = int(np.argmax(V))
    neg = np.where((us > us[k]) & (V < 0))[0]
    if not len(neg):
        return dict(barrier_u=float(us[k]), Vbar=float(V[k]), S3_over_T=np.inf)
    u_esc = us[neg[0]]

    def shoot(u0):
        def rhs(r, y):
            return [y[1], tv.dV(max(y[0], 1e-12)) - (2 / r) * y[1] if r > 0 else tv.dV(y[0]) / 3]

        def over(r, y):
            return y[0] - 1e-6 * u0
        over.terminal = True

        def under(r, y):
            return y[1]
        under.terminal = True
        under.direction = 1
        r0 = 1e-6
        a0 = tv.dV(u0) / 3
        sol = solve_ivp(rhs, (r0, rmax), [u0 + 0.5 * a0 * r0 ** 2, a0 * r0], events=[over, under], rtol=1e-10,
                        atol=1e-12, dense_output=True, max_step=1.0)
        if len(sol.t_events[0]):
            return +1, sol
        return -1, sol
    lo, hi = u_esc, u_esc
    s, _ = shoot(hi)
    while s < 0:
        hi *= 1.3
        s, _ = shoot(hi)
        if hi > 1e5:
            return dict(barrier_u=float(us[k]), Vbar=float(V[k]), S3_over_T=None, note='no overshoot found')
    for _ in range(36):
        mid = 0.5 * (lo + hi)
        s, sol = shoot(mid)
        if s > 0:
            hi = mid
        else:
            lo = mid
        if (hi - lo) / hi < 1e-11:
            break
    s, sol = shoot(lo)
    rr = np.linspace(sol.t[0], sol.t[-1], 6000)
    uu, up = sol.sol(rr)
    S3 = 4 * np.pi * np.trapezoid(rr ** 2 * (0.5 * up ** 2 + tv.V(np.maximum(uu, 1e-9))), rr)
    # Derrick check: kinetic = 3 * (-potential integral)... for O(3): T = -3 V_int, S3 = T + V = (2/3)T
    K = 4 * np.pi * np.trapezoid(rr ** 2 * 0.5 * up ** 2, rr)
    return dict(barrier_u=float(us[k]), barrier_h_over_T=float(us[k]), Vbar=float(V[k]), u_escape=float(u_esc),
                u0=float(lo), S3_over_T=float(S3), derrick_ratio_S3_over_2K3=float(S3 / (2 * K / 3)))


def main(Mts=(172.33, 172.60, 172.87), Tlist=(1e13, 6e13, 1.6e14, 3e14)):
    path = os.path.join(RES, 'thermal.json')
    prev = json.load(open(path)) if os.path.exists(path) else {}
    bg = json.load(open(os.path.join(RES, 'background.json')))
    tab = bg['table_downsampled']
    tau, a, rr = np.array(tab['tau']), np.array(tab['a']), np.array(tab['rho_r'])
    out = prev
    for Mt in Mts:
        sm = S.SM(Mt)
        rows = {}
        for T in Tlist:
            tv = ThermalV(sm, T)
            b = bounce(tv.splined())
            lamT = float(tv.lam(np.log(T * (b.get('u0') or 3.0))))
            b.update(kappa_T=float(tv.kT), lam_eff_at_bounce_center=lamT,
                     highT_estimate=(float(18.9 * np.sqrt(tv.kT) / abs(lamT)) if lamT < 0 else None))
            rows[f'{T:.3g}'] = b
            print(Mt, '%.3g' % T, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in b.items()}, flush=True)
        # bubbles along the engine history (H_i = 1.29e11 and 1.51e11)
        Ts = np.array([float(k) for k in rows])
        Ss = np.array([rows[k]['S3_over_T'] if rows[k]['S3_over_T'] is not None else np.inf for k in rows])
        fin = np.isfinite(Ss)
        Sint = interp1d(np.log(Ts[fin]), Ss[fin], bounds_error=False, fill_value=(Ss[fin][0], Ss[fin][-1]))
        nb = {}
        for H_GeV in (1.29e11, 1.51e11):
            T = T_GeV(rr, H_GeV)
            m = T > 1e12
            S3 = Sint(np.log(T[m]))
            lnG = 4 * np.log(T[m] / H_GeV) + 1.5 * np.log(S3 / (2 * np.pi)) - S3     # ln(Gamma/H_i^4)
            integ = np.log(np.trapezoid(np.exp(lnG - lnG.max()) * a[m] ** 3, tau[m])) + lnG.max()
            nb[f'H={H_GeV:.3g}'] = {f'N={N}': float((integ + 3 * N) / np.log(10)) for N in (50, 55, 60)}
            nb[f'H={H_GeV:.3g}']['min_S3_over_T_along_history'] = float(S3.min())
            nb[f'H={H_GeV:.3g}']['T_max_GeV'] = float(T.max())
        out[f'Mt={Mt}'] = dict(bounce=rows, log10_bubbles_in_past_light_cone=nb)
        print(Mt, nb, flush=True)
        json.dump(out, open(path, 'w'), indent=1, default=float)


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        main(Mts=(float(sys.argv[1]),), Tlist=tuple(float(v) for v in sys.argv[2:]) or (1e13, 6e13, 1.6e14, 3e14))
    else:
        main()
