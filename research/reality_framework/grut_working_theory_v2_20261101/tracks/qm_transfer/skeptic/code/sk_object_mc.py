"""qm_transfer SKEPTIC: can a finite Q-matter object act as the thermostat at T_o that the track's 'thermo' kinetics
assumes?  Own code (Claude, skeptic pass, 2 Nov 2026).  [C on E inputs]

The track sets the object temperature T_o by instantaneous energy balance and re-emits every absorbed charge as an
evaporated X drawn from the saturated-vapour flux at T_o (thermal, <E_k> = 2 T_o).  That requires the object to HOLD
a thermal excitation energy U(T_o) well above the binding energy Delta = 1 - Omega(Q) of one X; otherwise it cannot
emit spontaneously and every emission is triggered by, and limited to the energy of, the previous absorption.

Model (microcanonical, Weisskopf evaporation theory [R]): internal modes = bulk phonons (c_s = 0.5467, qm_micro [Q])
+ capillary ripplons on the sphere (omega^2 = sigma k^3 / w, sigma = 98.8 m_X^3, w = e + P = 348.6 m_X^4 [Q]);
S(U) from the canonical mode sums, rho(U) = exp S(U); perfect internal thermalisation (the most favourable case for
the track: qm_micro's Beliaev trapping would make the object even more elastic).
  absorption: rate pi R^2 n_v <v>, kinetic energy from the incident flux (Maxwellian T_v), U -> U + Delta + E;
  emission:   dW/dE' = (R^2/pi) E' rho(U - Delta - E')/rho(U) (geometric, m_X = 1, one internal state), U -> U - Delta - E'.
The canonical limit of this kernel is exactly the track's thermostat: total rate (R^2/pi) T^2 e^{-Delta/T} = pi R^2
n_sat(T) <v>_T (detailed balance) and <E'> = 2T.  Gillespie simulation of one object in a fixed vapour (n_v, T_v)
taken from the track's own post-freeze-out background (qt_hidden 'sudden'/'thermo' runs, probe_bg.json).
Continuous U (no energy grid; an earlier grid version stalled at the threshold U = Delta and
produced spurious heating for small objects).
Output: per absorbed X, the energy given to the vapour, (Sum E' - Sum E_in)/N_abs, compared with the canonical
thermostat value; N_emit/N_abs; mean U/Delta.
usage: sk_object_mc.py probe_bg.json OUT.json [n_events]"""
import json, math, sys
import numpy as np
from scipy.special import gamma as Gf, zeta as Zf
from scipy.optimize import brentq

NB, CS, SIG, WEP = 449.31336788, 0.5467, 98.8, 348.6
C_S = (1 / (3 * math.pi)) * (SIG / WEP)**(-2 / 3) * Gf(7 / 3) * Zf(7 / 3)      # ripplon energy per area / T^(7/3)
C_B = math.pi**2 / 30 / CS**3                                                 # phonon energy per volume / T^4


def geom(Q):
    R = (3 * Q / (4 * math.pi * NB))**(1 / 3)
    return R, 4 * math.pi * R * R, Q / NB


def U_of_T(T, Q):
    R, A, V = geom(Q)
    return C_B * V * T**4 + C_S * A * T**(7 / 3)


def S_of_T(T, Q):
    R, A, V = geom(Q)
    return (4 / 3) * C_B * V * T**3 + (7 / 4) * C_S * A * T**(4 / 3)


class Obj:
    """continuous-U microcanonical object: S(U) by interpolation in ln U; W(U) tabulated in eps = U - Delta on a log
    grid (no threshold stall: W > 0 for every U > Delta); E' sampled continuously."""
    def __init__(self, Q, Delta, nx=1500):
        self.Q, self.D = Q, Delta
        Ts = np.exp(np.linspace(math.log(1e-8), math.log(0.3), 6000))
        self.lU = np.log(U_of_T(Ts, Q)); self.Stab = S_of_T(Ts, Q)
        self.nx = nx
        self.leps = np.linspace(math.log(1e-12), math.log(max(5 * Delta, 3 * U_of_T(0.02, Q))), 700)
        self.lW = np.array([math.log(max(self._Wint(math.exp(l))[0], 1e-300)) for l in self.leps])

    def S(self, U):
        U = np.asarray(U, float)
        return np.where(U > 0, np.interp(np.log(np.maximum(U, 1e-300)), self.lU, self.Stab, left=0.0), 0.0)

    def _Wint(self, eps):
        # E' = x in [0, eps], residual eps - x; integrand x exp(S(eps - x) - S(D + eps)); grid dense near both ends
        u = (1 - np.cos(np.linspace(0, math.pi, self.nx))) / 2
        x = eps * u
        f = x * np.exp(self.S(eps - x) - float(self.S(self.D + eps)))
        return np.trapezoid(f, x) / math.pi, x, f

    def W(self, U):
        e = U - self.D
        if e <= 0:
            return 0.0
        le = math.log(e)
        if le < self.leps[0]:
            return math.exp(self.lW[0]) * (e / math.exp(self.leps[0]))**2
        if le > self.leps[-1]:
            return self._Wint(e)[0]
        return math.exp(float(np.interp(le, self.leps, self.lW)))

    def emit(self, U, rng):
        e = U - self.D
        _, x, f = self._Wint(e)
        c = np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(x))))
        c /= c[-1]
        return float(np.interp(rng.random(), c, x))


def canonical_To(n_v, T_v, D):
    # (1/pi) T^2 e^{-D/T} (D + 2T) = pi n_v <v> (D + 2 T_v)   (flux-weighted mean incident energy 2 T_v)
    rhs = math.pi * n_v * math.sqrt(8 * T_v / math.pi) * (D + 2 * T_v)
    f = lambda lT: math.log(math.exp(2 * lT) / math.pi * (D + 2 * math.exp(lT))) - D / math.exp(lT) - math.log(rhs)
    return math.exp(brentq(f, math.log(1e-4), math.log(0.1)))


def run(Q, D, n_v, T_v, nev, seed=1, U_start=None):
    rng = np.random.default_rng(seed)
    To = canonical_To(n_v, T_v, D)
    U0 = U_of_T(To, Q)
    ob = Obj(Q, D)
    a_r = math.pi * n_v * math.sqrt(8 * T_v / math.pi)          # absorption rate per R^2
    U = U0 if U_start is None else U_start
    Na = Ne = 0; Ein = Eout = 0.0; Usum = 0.0; burn = nev // 4
    for ev in range(nev):
        W = ob.W(U)
        if rng.random() < a_r / (a_r + W):
            E = rng.gamma(2.0, T_v)                             # flux-weighted Maxwellian kinetic energy
            U += D + E
            if ev >= burn: Na += 1; Ein += E
        else:
            Ep = ob.emit(U, rng)
            U = max(U - D - Ep, 0.0)
            if ev >= burn: Ne += 1; Eout += Ep
        if ev >= burn: Usum += U
    can_heat = (D + 2 * T_v) * 2 * To / (D + 2 * To) - 2 * T_v           # per absorption, canonical thermostat
    mc_heat = (Eout - Ein) / max(Na, 1)
    SD = float(ob.S(D))
    return dict(Q=Q, Delta=D, n_v=n_v, T_v=T_v, T_o_can=To, U_of_To=U0, U_of_To_over_Delta=U0 / D,
                S_of_Delta=SD, Delta_over_To=D / To, T_hold=D / SD,
                N_abs=Na, N_emit=Ne, Ne_over_Na=Ne / max(Na, 1), Ne_over_Na_can=(D + 2 * T_v) / (D + 2 * To),
                mean_Ein=Ein / max(Na, 1), mean_Eout=Eout / max(Ne, 1), mean_Eout_can=2 * To,
                heat_per_abs_mc=mc_heat, heat_per_abs_can=can_heat, heating_efficiency=mc_heat / can_heat,
                mean_U_over_Delta=Usum / max(nev - burn, 1) / D)


if __name__ == '__main__':
    pb = json.load(open(sys.argv[1])); out_fn = sys.argv[2]; nev = int(sys.argv[3]) if len(sys.argv) > 3 else 40000
    res = []
    # (1) the track's own models: object Q = zeta^-3 with the vapour of the same model at a = 2, 5, 20 a_f
    for z in (1e-2, 3e-3, 1e-3, 3e-4, 1e-4):
        th = pb['z%g_thermo' % z]; su = pb['z%g_sudden' % z]
        for r_th, r_su in zip(th['rows'], su['rows']):
            if r_th['fac'] not in (2, 5, 20):
                continue
            r = run(th['Q'], th['Delta'], r_su['n_v'], r_su['Tv_rms'], nev)
            r.update(zeta=z, fac=r_th['fac'], T_o_track=r_th['T_o'])
            res.append(r)
            print('zeta=%.0e Q=%.1e a/a_f=%3d  T_o can %.4f (track %.4f)  S(Delta)=%6.1f Delta/T_o=%5.1f U(T_o)/Delta=%7.2f  <U>/Delta=%6.2f  Ne/Na=%.4f (can %.4f)  <E_out>/2T_o=%.3f  heat eff=%.3f' % (
                z, th['Q'], r['fac'], r['T_o_can'], r['T_o_track'], r['S_of_Delta'], r['Delta_over_To'], r['U_of_To_over_Delta'], r['mean_U_over_Delta'],
                r['Ne_over_Na'], r['Ne_over_Na_can'], r['mean_Eout'] / r['mean_Eout_can'], r['heating_efficiency']), flush=True)
    # (2) crossover scan in Q at fixed vapour (zeta = 1e-3 model at a = 5 a_f)
    su = pb['z0.001_sudden']['rows'][2]; D = pb['z0.001_sudden']['Delta']
    for Q in (1e7, 3e7, 1e8, 3e8, 6e8, 1e9, 2e9, 3e9, 1e10, 3e10, 1e11, 1e12):
        r = run(Q, D, su['n_v'], su['Tv_rms'], nev)
        r.update(scan='Q at fixed vapour (zeta=1e-3 model, a=5a_f)')
        res.append(r)
        print('scan Q=%.0e  S(Delta)=%6.1f Delta/T_o=%5.1f U(T_o)/Delta=%7.2f  heat eff=%.3f  Ne/Na=%.4f  <E_out>/2T_o=%.3f' % (
            Q, r['S_of_Delta'], r['Delta_over_To'], r['U_of_To_over_Delta'], r['heating_efficiency'], r['Ne_over_Na'], r['mean_Eout'] / r['mean_Eout_can']), flush=True)
    json.dump(dict(scope=__doc__.strip().splitlines()[0], constants=dict(c_s=CS, sigma=SIG, w=WEP, C_ripplon=C_S, C_phonon=C_B),
                   runs=res), open(out_fn, 'w'), indent=1)
