"""qm_micro SKEPTIC (Claude, 2 Nov 2026): energy-conserving two-moment model of the vapour freeze-out.  Own code.

Why: qm_micro's coarsen.py imposes the equilibrium T_h(a) on the vapour and lets only its chemical potential relax;
after Gamma_c = H its variant (b) lets the vapour cool as a^-2 with CAPTURE ONLY.  Neither conserves energy on the
objects: Q-matter has no light neutral state and a negligible heat capacity at T < 0.05, so the heat Delta + eps of
every absorbed X (Delta = 1 - Omega_b = 0.224) must leave with re-emitted X.  Here the vapour is a Maxwellian
(n_v, T_v) with exact moment exchange, and the objects are zero-heat-capacity bodies whose temperature T_o solves
the heat balance  Phi_a (Delta + <eps>_a) = Phi_e(T_o) (Delta + <eps>_e)  at every instant.

Per comoving charge (y = n_v/n_Q, K = vapour kinetic energy per charge = 1.5 y T_v; alpha = object area per charge):
   dy/dt = -alpha (Phi_a - Phi_e),          dK/dt = -2 H K - alpha (E_a - E_e),
   Phi_a = n_v sqrt(T_v/2pi) S0(T_v),  E_a = n_v sqrt(T_v/2pi) S1(T_v),
   Phi_e = n_sat(T_o) sqrt(T_o/2pi) S0(T_o),  E_e = n_sat(T_o) sqrt(T_o/2pi) S1(T_o)   (detailed balance),
   S_n(T) = int s(eps) eps^n (eps/T^2) e^{-eps/T} d eps  (angle-averaged planar sticking, qm_micro stick_thermal.json),
   n_sat(T) = (T/2pi)^{3/2} e^{-Delta/T} (1 + 15T/8)   (NR Boltzmann with the leading relativistic correction).
Objects: fixed number after ripening (taken from a coarsening run at T = 0.05), area alpha = alpha_0 (eta/eta_0)^{2/3}.
Background H(a), n_Q(a) and the start point from the same two-phase history as qm_micro (hist_<S/Q>.json).
Modes: 'energy' (above); 'qm_a' (T_v = T_eq(a) imposed, T_o = T_v: qm_micro's exchange, run to the end);
'qm_b' (qm_a until Gamma_c = H, then capture only with T_v ~ a^-2: qm_micro's variant b).  The last two validate the
set-up against qm_micro's own numbers.
"""
import json, math, os, sys, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import coarsen_heat as ch           # History class (copy of qm_micro's), radius, constants

QMR = os.path.join(HERE, '..', '..', 'results')
SKR = os.path.join(HERE, '..', 'results')
SCR = os.environ.get('SK_SCR', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/qm_micro_skeptic')
DELTA = 1 - ch.OMB
NB = ch.NB

_st = json.load(open(os.path.join(QMR, 'stick_thermal.json')))
_Tk = np.array(_st['T_kin']); _sa = np.array(_st['s_X_angle_avg'])


def s_of(eps):
    eps = np.asarray(eps, float)
    out = np.exp(np.interp(np.log(np.clip(eps, _Tk[0], _Tk[-1])), np.log(_Tk), np.log(_sa)))
    return np.where(eps < _Tk[0], _sa[0] * np.sqrt(np.maximum(eps, 0) / _Tk[0]), out)


_U = np.linspace(0, 1, 1601)[1:] * 7.0


TRAP_R = None          # set to an object radius to multiply the sticking by qm_micro's Beliaev trapping factor


def Smom(T, stick=True):
    eps = T * _U * _U
    w = (eps / T ** 2) * np.exp(-eps / T) * 2 * T * _U
    s = s_of(eps) if stick else 1.0
    if TRAP_R is not None:
        s = s * float(ch.trap_factor(np.array([TRAP_R]), T)[0])     # qm_micro's [E] factor, evaluated at the quantum's T
    return float(np.trapezoid(w * s, _U)), float(np.trapezoid(w * s * eps, _U))


def nsat(T):
    return (T / (2 * math.pi)) ** 1.5 * math.exp(-DELTA / T) * (1 + 15 * T / 8)


def T_object(Pa_heat, To_guess):
    """solve Phi_e(T_o)(Delta + <eps>_e) = Pa_heat (heat in per area)"""
    def g(lt):
        To = math.exp(lt)
        S0, S1 = Smom(To)
        pref = nsat(To) * math.sqrt(To / (2 * math.pi))
        return math.log(pref * (DELTA * S0 + S1)) - math.log(Pa_heat)
    lo, hi = math.log(To_guess) - 1.5, math.log(To_guess) + 1.5
    while g(lo) > 0:
        lo -= 1.0
    while g(hi) < 0:
        hi += 1.0
    return math.exp(brentq(g, lo, hi, xtol=1e-13))


class Model:
    def __init__(self, SQ, coarsen_json, case, T0=0.05):
        self.H = ch.History(os.path.join(SCR, 'hist_%s.json' % SQ))
        d = json.load(open(coarsen_json))
        if 'runs' in d:
            run = [r for r in d['runs'] if r['case'] == case][0]
        else:                                   # qm_micro's compact results/coarsen_runs.json
            run = d[SQ][case]
        rows = run['rows']
        # object population at T0: number per charge and capture-effective charge
        i0 = int(np.argmin([abs(math.log(r['T'] / T0)) for r in rows]))
        r0 = rows[i0]
        self.T0 = r0['T']
        self.la0 = self.H.la_of_T(self.T0)
        self.eta0 = r0['eta']
        self.alpha0 = 4 * math.pi * (3 / (4 * math.pi * NB)) ** (2 / 3) * r0['eta'] / r0['Q_capture_eff'] ** (1 / 3)
        self.Qeff0 = r0['Q_capture_eff']
        self.case = case
        self.To = self.T0

    def alpha(self, y):
        return self.alpha0 * (max(1 - y, 1e-6) / self.eta0) ** (2 / 3)

    def rhs_energy(self, la, z):
        y = math.exp(z[0]); Tv = math.exp(z[1])
        s = self.H.at(la)
        nQ, H = s['nQ'], s['H']
        nv = y * nQ
        S0a, S1a = Smom(Tv)
        pa = nv * math.sqrt(Tv / (2 * math.pi))
        Pa, Ea = pa * S0a, pa * S1a
        To = T_object(DELTA * Pa + Ea, self.To)
        self.To = To
        S0e, S1e = Smom(To)
        pe = nsat(To) * math.sqrt(To / (2 * math.pi))
        Pe, Ee = pe * S0e, pe * S1e
        al = self.alpha(y)
        dy = -al * (Pa - Pe) / H
        K = 1.5 * y * Tv
        dK = -2 * K - al * (Ea - Ee) / H
        # d ln y, d ln Tv  (K = 1.5 y Tv -> dlnK = dlny + dlnTv)
        dlny = dy / y
        dlnT = dK / K - dlny
        return [dlny, dlnT]

    def run_energy(self, la_end=None):
        s0 = self.H.at(self.la0)
        y0 = 1 - s0['eta']
        la_end = self.H.lna[-1] if la_end is None else la_end
        rows = []
        def rec(la, z):
            y = math.exp(z[0]); Tv = math.exp(z[1]); s = self.H.at(la)
            Gc = self.alpha(y) * s['nQ'] * math.sqrt(Tv / (2 * math.pi)) * Smom(Tv)[0] / s['H']   # gross absorption per vapour X / H
            return dict(lna=la, T_eq=s['T'], T_vis_GeV=s['T_vis_GeV'], f_v=y, f_v_eq=1 - s['eta'], T_v=Tv, T_o=self.To,
                        Gamma_a_over_H=Gc)
        sol = solve_ivp(self.rhs_energy, (self.la0, la_end), [math.log(y0), math.log(s0['T'])], method='LSODA',
                        rtol=1e-8, atol=1e-10, dense_output=True)
        xs = np.linspace(self.la0, sol.t[-1], 300)
        for la in xs:
            z = sol.sol(la)
            self.rhs_energy(la, z)       # sets T_o
            rows.append(rec(la, z))
        return dict(success=bool(sol.success), rows=rows, final=rows[-1])

    def run_qm(self, variant='b'):
        """qm_micro's model on the same object area: T_v = T_eq imposed and isothermal exchange; implicit balance
        y = f_v solved each step; variant b: after Gamma_c = H, capture only with T_v ~ a^-2."""
        la = self.la0
        s0 = self.H.at(la)
        y = 1 - s0['eta']
        la_end = self.H.lna[-1]
        dla = 2e-3
        frozen = None
        rows = []
        while la < la_end:
            s = self.H.at(la + dla)
            T, nQ, H = s['T'], s['nQ'], s['H']
            dt = dla / H
            if frozen is None:
                # implicit: dy/dt = -alpha nQ sqrt(T/2pi) S0 (y - y_sat), y_sat = nsat(T)/nQ
                S0 = Smom(T)[0]
                rate = self.alpha(y) * nQ * math.sqrt(T / (2 * math.pi)) * S0
                ysat = nsat(T) / nQ
                y = (y + rate * dt * ysat) / (1 + rate * dt)
                if rate / H < 1.0 and T < 0.1:
                    frozen = dict(lna=la + dla, T=T, f_v=y, Gamma_over_H=rate / H, T_vis_GeV=s['T_vis_GeV'])
                    Tf, laf = T, la + dla
            else:
                if variant == 'a':
                    S0 = Smom(T)[0]
                    rate = self.alpha(y) * nQ * math.sqrt(T / (2 * math.pi)) * S0
                    ysat = nsat(T) / nQ
                    y = (y + rate * dt * ysat) / (1 + rate * dt)
                else:
                    Tv = Tf * math.exp(-2 * (la + dla - laf))
                    S0 = Smom(Tv)[0]
                    rate = self.alpha(y) * nQ * math.sqrt(Tv / (2 * math.pi)) * S0
                    y *= math.exp(-rate * dt)
            la += dla
            dla = min(dla * 1.02, 0.05)
            rows.append(dict(lna=la, T=T, f_v=y))
        return dict(freeze=frozen, final_f_v=y, rows=rows[::20])


def main(out):
    global TRAP_R
    t0 = time.time()
    res = dict(notes=__doc__, runs=[])
    jobs = []
    for SQ in ('1.5', '2.5', '3.0', '3.2'):
        bh = os.path.join(SCR, ('cheat2_%s.json' if SQ in ('2.5', '3.0') else 'cheat_%s.json') % SQ)
        jobs += [(SQ, os.path.join(SCR, 'cheat_%s.json' % SQ), 'control_computed', False),
                 (SQ, bh, 'bose_heat_sphi1', False),
                 (SQ, os.path.join(QMR, 'coarsen_runs.json'), 'Q0_4e11_trap_beliaev', True)]
        if SQ in ('2.5', '3.0'):
            jobs.append((SQ, os.path.join(SCR, 'cheat_%s.json' % SQ), 'heat_sphi1', False))
    for SQ, cj, case, trap in jobs:
        if True:
            TRAP_R = None
            M = Model(SQ, cj, case)
            if trap:
                TRAP_R = float(ch.radius(M.Qeff0))
            rq_b = M.run_qm('b'); rq_a = M.run_qm('a')
            re = M.run_energy()
            fin = re['final']
            # kinetic freeze-out of the energy run: Gamma_a/H = 1 crossing and where f_v leaves f_v_eq by 10%
            rows = re['rows']
            def first(cond):
                for r in rows:
                    if cond(r):
                        return r
                return None
            fo = first(lambda r: r['Gamma_a_over_H'] < 1)
            dev = first(lambda r: r['f_v'] > 1.1 * r['f_v_eq'])
            rec = dict(S_over_Q=float(SQ), objects_from=case, beliaev_trap=trap, T0=M.T0, alpha0=M.alpha0, Qeff0=M.Qeff0,
                       qm_b=dict(freeze=rq_b['freeze'], final_f_v=rq_b['final_f_v']), qm_a=dict(final_f_v=rq_a['final_f_v']),
                       energy=dict(final_f_v=fin['f_v'], final_T_v=fin['T_v'], final_T_o=fin['T_o'], Gamma_a_eq_H=fo,
                                   f_v_10pct_above_eq=dev, rows=rows[::3]))
            res['runs'].append(rec)
            print('S/Q=%s objects=%s: alpha0=%.3e Qeff0=%.3e | qm-like: freeze T=%.4g f_v=%.4f -> b %.4f / a %.4f | energy: f_v final %.4f (at Gamma_a=H: T_eq=%.4g T_v=%.4g T_o=%.4g f_v=%.4f, f_v_eq=%.4f; 10%% excess at T_eq=%.4g Tvis=%.3g MeV, Gamma/H=%.0f) | %.0fs' % (
                SQ, case, M.alpha0, M.Qeff0, rq_b['freeze']['T'], rq_b['freeze']['f_v'], rq_b['final_f_v'], rq_a['final_f_v'],
                fin['f_v'], fo['T_eq'], fo['T_v'], fo['T_o'], fo['f_v'], fo['f_v_eq'], dev['T_eq'] if dev else float('nan'),
                dev['T_vis_GeV'] * 1e3 if dev else float('nan'), dev['Gamma_a_over_H'] if dev else float('nan'), time.time() - t0), flush=True)
            json.dump(res, open(out, 'w'), indent=None, default=float, separators=(',', ':'))


if __name__ == '__main__':
    main(sys.argv[1])
