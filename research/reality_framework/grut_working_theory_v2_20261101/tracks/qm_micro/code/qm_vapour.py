"""qm_micro item 2 (Claude, 2 Nov 2026): saturated vapour next to a finite object (Gibbs-Thomson), the sticking-limited
exchange fluxes and the ripening coefficient; plus the equilibrium two-phase histories at fixed S/Q used by item 3.

Inputs (read-only):  Omega(Q) from fx_formation/results/family.json (solver table to 1e8, the first-law derivative of
the three-term fit beyond; via thermo_copy.Family), the vapour by the Bessel series of fx_formation's skeptic
(bessel_vapour_copy.py), the sticking tables of results/stick_thermal.json (this track), the two-phase equilibrium
state of fx_formation (eq_history_copy.state, coexistence with the one-loop liquid thermal pressure).
Definitions (m_X = 1):
  n_sat(T, Q) = n_vap(mu = Omega(Q), T)  for X, Xbar (phi is neutral, mu = 0);
  sticking-weighted one-sided flux  Phi_s(mu, T) = (1/8 pi^2) int p^2 f(E; mu, T) s_ang(E - 1) dE,
     s_ang = int_0^1 s(E, mu) 2 mu dmu (angle-averaged planar sticking, X) or d_ang (Xbar destruction);
  by detailed balance the evaporation flux of an object at chemical potential Omega(Q) equals Phi_s(Omega(Q), T);
  net charge flux per area j = [Phi_X(mu_v) - Phi_X(Omega)] - [Phi_Xbar(-mu_v) - Phi_Xbar(-Omega)]
     ~ L (mu_v - Omega),  L = d/dmu [Phi_X(mu) + Phi_Xbar(-mu)] (kinetic coefficient per area);
  Gibbs-Thomson: Omega(R) = mu_c(T) + 2 sigma/(n_B R) (sigma = 98.8 m_X^3 from the E(Q) fit; planar value 99.04 here);
  interface-limited (Wagner) ripening: dR/dt = K_W (1/R_c - 1/R), K_W = 2 sigma L / n_B^2.
"""
import json, math, os, sys, time
import numpy as np
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bessel_vapour_copy as bv
import thermo_copy as th
import eq_history_copy as eh

OMB, NB = th.OM_B, th.N_B
SIGMA = 98.8
MX = 0.0924089637375362
RES = os.path.join(HERE, '..', 'results')


class Stick:
    def __init__(self, path=os.path.join(RES, 'stick_thermal.json')):
        d = json.load(open(path))
        self.Tk = np.array(d['T_kin']); self.sX = np.array(d['s_X_angle_avg']); self.dB = np.array(d['d_Xbar_angle_avg'])

    def s(self, Tkin, which='X'):
        y = self.sX if which == 'X' else self.dB
        Tkin = np.asarray(Tkin, float)
        out = np.exp(np.interp(np.log(np.clip(Tkin, self.Tk[0], self.Tk[-1])), np.log(self.Tk), np.log(y)))
        lo = Tkin < self.Tk[0]
        return np.where(lo, y[0] * np.sqrt(np.maximum(Tkin, 0) / self.Tk[0]), out)


ST = None


def flux(mu, T, which='X', stick=True, deriv=False):
    """One-sided flux of X (chemical potential mu) or Xbar (chemical potential -mu is passed as mu) per area,
    sticking-weighted if stick; deriv returns d/dmu as well."""
    global ST
    if ST is None:
        ST = Stick()
    u = np.linspace(0, 1, 6001)[1:] * 9.0
    Tkin = T * u * u
    E = 1 + Tkin
    x = (E - mu) / T
    f = 1 / np.expm1(np.minimum(x, 700))
    w = (E * E - 1) * 2 * T * u / (8 * math.pi ** 2)
    s = ST.s(Tkin, which) if stick else 1.0
    F = np.trapezoid(w * f * s, u)
    if not deriv:
        return float(F)
    dF = np.trapezoid(w * f * (1 + f) * s, u) / T
    return float(F), float(dF)


def omega_of_Q(fam, Q):
    """Omega(Q): solver table below 1e8 (interpolated from family.json rows), first-law derivative of the fit above."""
    d = json.load(open(os.path.join(th.RES, 'family.json')))
    rows = sorted(d['rows'], key=lambda r: r['Q'])
    Qs = np.array([r['Q'] for r in rows]); Os = np.array([r['Omega'] for r in rows])
    f = d['fit3']
    out = []
    for q in np.atleast_1d(Q):
        if q <= Qs[-1] and q >= Qs[0]:
            out.append(float(np.interp(math.log(q), np.log(Qs), Os)))
        elif q < Qs[0]:
            out.append(float(d['Omega_at_Qmin']))
        else:
            out.append(OMB + 2 / 3 * f['A'] * q ** (-1 / 3) + 1 / 3 * f['B'] * q ** (-2 / 3))
    return np.array(out)


def gibbs_thomson_table():
    fam = None
    Qs = [5.6709e4, 7.1263e4, 2e5, 1e6, 1e8, 1e10, 1e12, 1e14, 1e16, 1e18, 1e20, float('inf')]
    Ts = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.336, 0.5, 1.0, 2.0, 3.0, 5.0]
    Om = [float(omega_of_Q(fam, q)[0]) if math.isfinite(q) else OMB for q in Qs]
    rows = []
    for T in Ts:
        mc = bv.mu_coex(T) if T <= 1.6 else None
        for q, o in zip(Qs, Om):
            v = bv.vapour(o, T)
            r = dict(T=T, Q=q if math.isfinite(q) else 'bulk', Omega=o, n_X=v['nX'], n_Xbar=v['nXb'], n_net=v['n'], s_vap=v['s'],
                     n_over_s=v['n'] / v['s'], fugacity=math.exp((o - 1) / T))
            vb = bv.vapour(OMB, T)
            r['GT_ratio_nX'] = v['nX'] / vb['nX']
            r['GT_ratio_NR'] = math.exp((o - OMB) / T)
            if mc is not None:
                r['mu_c_vapour_shifted'] = mc
            rows.append(r)
    return dict(Q=[q if math.isfinite(q) else 'bulk' for q in Qs], Omega=Om, rows=rows)


def kinetic_coefficients():
    rows = []
    for T in [0.003, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0]:
        mc = bv.mu_coex(T) if T < 0.5 else th.mu_coex(T)      # bulk coexistence; the liquid thermal pressure matters above ~0.5
        FX, dFX = flux(mc, T, 'X', True, True)
        FB, dFB = flux(-mc, T, 'Xbar', True, True)
        FX0 = flux(mc, T, 'X', False)
        v = bv.vapour(mc, T)
        L = dFX + dFB
        KW = 2 * SIGMA * L / NB ** 2
        rows.append(dict(T=T, mu_c=mc, n_X=v['nX'], n_Xbar=v['nXb'], Phi_X_bare=FX0, Phi_X_stick=FX, Phi_Xbar_stick=FB,
                         s_flux_avg=FX / FX0, vbar_eff=4 * FX0 / v['nX'], L=L, K_W=KW,
                         ripening_time_R1e5_over_R2=1.0 / KW))
        print('T=%.4g mu_c=%.6f n_X=%.3e Phi_s/Phi=%.4f L=%.3e K_W=%.3e' % (T, mc, v['nX'], FX / FX0, L, KW), flush=True)
    return rows


def history(SQ, YQ, Tgrid):
    """Two-phase equilibrium history at fixed Y_h = 1/(S/Q) (fx_formation's state(), copied) with the visible clock
    for the given Y_Q, plus the sticking-weighted fluxes at mu_c(T)."""
    Yh = 1.0 / SQ
    rows = []
    for T in Tgrid:
        st = eh.state(float(T), Yh)
        if st['phase'] != 'two-phase':
            rows.append(dict(T=float(T), phase=st['phase'])); continue
        mc = st['mu']
        v = bv.vapour(mc, T) if T <= 1.6 else None
        if v is None:
            vq = th.vapour(mc, T); v = dict(nX=vq['n_X'], nXb=vq['n_Xbar'], n=vq['n_net'], s=vq['s'], P=vq['P'], rho=vq['e'])
        FX, dFX = flux(mc, T, 'X', True, True)
        FB, dFB = flux(-mc, T, 'Xbar', True, True)
        FX0 = flux(mc, T, 'X', False)
        try:
            Tv = eh.T_vis_of(st['nQ'], YQ)
        except ValueError:
            break                                  # the hidden T_h falls only logarithmically: past today
        if Tv < 2.35e-13:
            break
        H = eh.hubble(Tv) / MX
        rows.append(dict(T=float(T), phase='two-phase', eta=st['eta'], nQ=st['nQ'], mu_c=mc, n_X=v['nX'], n_Xbar=v['nXb'],
                         n_vap_net=v['n'], s_vap=v['s'], eP_vap=v['rho'] + v['P'], Phi_X=FX, Phi_Xbar=FB, L=dFX + dFB,
                         Phi_X_bare=FX0, T_vis_GeV=Tv, H=H, phi_vol=st['eta'] * st['nQ'] / NB, Xbar_per_charge=st['Xbar_per_charge']))
    return dict(S_over_Q=SQ, Y_h=Yh, Y_Q=YQ, rows=rows)


def main(out, part):
    t0 = time.time()
    if part == 'gt':
        res = dict(gibbs_thomson=gibbs_thomson_table(), kinetic=kinetic_coefficients(), sigma=SIGMA, n_B=NB)
        json.dump(res, open(out, 'w'), indent=1, default=float)
    else:
        SQ = float(part)
        M = float(os.environ.get('QM_MASS', '0.785'))
        YQ = 5.07941423e-9 * 0.9374 / M
        Tg = np.unique(np.r_[np.geomspace(0.0015, 0.05, 70), np.geomspace(0.05, 1.0, 70), np.geomspace(1.0, 12.0, 50)])[::-1]
        res = history(SQ, YQ, Tg)
        res['mass_per_charge_assumed'] = M
        json.dump(res, open(out, 'w'), indent=1, default=float)
        tw = [r for r in res['rows'] if r['phase'] == 'two-phase']
        print('S/Q=%.3f: %d two-phase points; eta min %.4f at T=%.4f; T range %.4g..%.4g; %.0fs' % (
            SQ, len(tw), min(r['eta'] for r in tw), min(tw, key=lambda r: r['eta'])['T'], tw[-1]['T'], tw[0]['T'], time.time() - t0))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
