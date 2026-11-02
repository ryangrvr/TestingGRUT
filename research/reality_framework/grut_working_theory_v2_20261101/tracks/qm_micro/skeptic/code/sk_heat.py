"""qm_micro SKEPTIC (Claude, 2 Nov 2026): the non-isothermal (latent-heat-limited) exchange factor F_h(T).  Own code.

Physics [D]: Q-matter has no light neutral state (m_phi = sqrt6, the vector and Y, Sigma are heavy) and a negligible
heat capacity at T <~ 1 (phonon gas).  An object therefore cannot keep the temperature of the vapour while it
grows or shrinks: the heat h = E - q Omega released per absorbed quantum (E the energy, q = +1 X, -1 Xbar, 0 phi)
must leave again with emitted quanta.  The object temperature T_o adjusts so that the NET HEAT FLUX VANISHES
(zero-heat-capacity steady state; the thermal relaxation time of an object is shorter than its ripening time by
~(lambda^2+2)/(c_o F_h g) >> 1, g = 2 sigma/(n_B R T)).

Linear response (detailed balance; Bose weights; one-sided flux measure p^2 dE/(8 pi^2)): with per-species weights
  w_i(E) = s_i(E) f_i (1 + f_i) p^2/(8 pi^2),   f_i at chemical potential q_i mu_c and temperature T,
and the moments  W0 = sum q_i^2 int w_i,  W1 = sum q_i int w_i h_i/T,  W2 = sum int w_i h_i^2/T^2,
the net charge flux into an object at affinity Delta a = (mu_v - Omega(R))/T is
  J = W0 (1 - W1^2/(W0 W2)) Delta a    (heat-balanced)   vs   J_iso = W0 Delta a (isothermal, qm_micro).
So F_h = 1 - W1^2/(W0 W2).  For a pure Maxwellian X vapour with s = 1: F_h = 2T^2/((Delta + 2T)^2 + 2T^2),
Delta = 1 - Omega_b; this is the Feder-Russell-Lothe-Pound (1966) non-isothermal factor b^2/(b^2 + q^2) [R] with no
carrier gas.  The ripening coefficient becomes K_W F_h; the growth rate of every object is multiplied by F_h.
phi (neutral) can carry heat without charge; its absorption probability on Q-matter is not computed (qm_micro
computed X and Xbar only), so it is bracketed by s_phi = 0 and s_phi = 1.
Sticking tables: qm_micro results/stick_thermal.json (angle-averaged; verified point-wise by sk_scatter.py).
"""
import json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
QM = os.path.join(HERE, '..', '..')
OMB = 0.7758647074595397
MPHI = math.sqrt(6.0)


class Stick:
    def __init__(self):
        d = json.load(open(os.path.join(QM, 'results', 'stick_thermal.json')))
        self.Tk = np.array(d['T_kin'])
        self.sX = np.array(d['s_X_angle_avg']); self.dB = np.array(d['d_Xbar_angle_avg'])

    def s(self, Tkin, which):
        y = self.sX if which == 'X' else self.dB
        Tkin = np.asarray(Tkin, float)
        out = np.exp(np.interp(np.log(np.clip(Tkin, self.Tk[0], self.Tk[-1])), np.log(self.Tk), np.log(y)))
        return np.where(Tkin < self.Tk[0], y[0] * np.sqrt(np.maximum(Tkin, 0) / self.Tk[0]), out)


ST = Stick()


def moments(T, mu=OMB, s_phi=1.0, include_Xbar=True, stat='bose', stick=True):
    """W0, W1, W2 and per-species pieces at temperature T, coexistence mu (object Omega = mu)."""
    out = {}
    W0 = W1 = W2 = 0.0
    specs = [('X', 1, 1.0, mu)]
    if include_Xbar:
        specs.append(('Xbar', -1, 1.0, -mu))
    if s_phi > 0:
        specs.append(('phi', 0, MPHI, 0.0))
    for name, q, m, muq in specs:
        u = np.linspace(0, 1, 20001)[1:] * 10.0
        Tkin = T * u * u
        E = m + Tkin
        p2 = E * E - m * m
        x = (E - muq) / T
        f = 1 / np.expm1(x) if stat == 'bose' else np.exp(-x)
        ff = f * (1 + f) if stat == 'bose' else f
        if name == 'phi':
            s = s_phi * np.ones_like(E)
        else:
            s = ST.s(Tkin, name) if stick else 1.0
        w = s * ff * p2 / (8 * math.pi ** 2) * 2 * T * u        # dE = 2 T u du
        h = E - q * mu
        w0 = np.trapezoid(w, u) * q * q
        w1 = np.trapezoid(w * h / T, u) * q
        w2 = np.trapezoid(w * h * h / T / T, u)
        out[name] = dict(W0=float(w0), W1=float(w1), W2=float(w2))
        W0 += w0; W1 += w1; W2 += w2
    Fh = 1 - W1 * W1 / (W0 * W2)
    return dict(T=T, W0=W0, W1=W1, W2=W2, F_h=Fh, species=out, dTo_over_T_per_unit_affinity=-W1 / W2)


def F_h_NR(T, Delta=1 - OMB):
    return 2 * T * T / ((Delta + 2 * T) ** 2 + 2 * T * T)


def main(out):
    rows = []
    Ts = [0.003, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0]
    for T in Ts:
        r_phi1 = moments(T, s_phi=1.0)
        r_phi0 = moments(T, s_phi=0.0)
        r_nr = moments(T, s_phi=0.0, include_Xbar=False, stat='boltz', stick=False)
        rows.append(dict(T=T, F_h_sphi1=r_phi1['F_h'], F_h_sphi0=r_phi0['F_h'], F_h_X_only_Boltzmann_s1=r_nr['F_h'],
                         F_h_closed_form_NR=F_h_NR(T), W0=r_phi0['W0'],
                         species_sphi1=r_phi1['species']))
        print('T=%.4g  F_h: s_phi=1 %.4f  s_phi=0 %.4f  | X only Boltzmann s=1 %.4f (closed form %.4f)' % (
            T, r_phi1['F_h'], r_phi0['F_h'], r_nr['F_h'], F_h_NR(T)), flush=True)
    json.dump(dict(notes=__doc__, rows=rows), open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1])
