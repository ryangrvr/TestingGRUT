"""qm_micro item 1c (Claude, 2 Nov 2026): sticking s(T, mu = cos theta) tables for X and Xbar on planar bulk
Q-matter and their thermal-flux averages.  Own code; uses scatter1d.solve (Numerov, exact at linear order).

Flux of a vapour species hitting a plane: dPhi ~ f(E) p^2 dE * mu dmu (mu = cos theta in (0, 1]).
  s_bar(T_h) = int s(E, mu) f p^2 mu dmu dE / int f p^2 mu dmu dE,
with f = Boltzmann e^{-E/T} or Bose 1/(e^{(E - mu_c)/T} - 1) at the coexistence chemical potential mu_c = Omega_b
(the two agree to the digits quoted for T <= 0.5, since e^{(mu_c - m)/T} <= 0.64 enters only through the
energy dependence of s).  The table runs over kinetic energies 1e-6 .. 30 m_X; s is extended below 1e-6 by its
sqrt(T) threshold law and set to its T = 30 value above (it is 1 to 1e-9 there).
For Xbar the table holds the destruction probability d = 1 - R(Xbar -> Xbar) and its channels.
"""
import json, math, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import scatter1d as sc

MU_NODES, MU_W = np.polynomial.legendre.leggauss(16)
MU_NODES = 0.5 * (MU_NODES + 1); MU_W = 0.5 * MU_W


def tables(I):
    Om = I.Om
    Tk = np.unique(np.r_[np.geomspace(1e-6, 1e-2, 14), np.geomspace(1e-2, 1.0, 30), np.geomspace(1.0, 30.0, 12)])
    sX = np.zeros((Tk.size, MU_NODES.size)); dB = np.zeros_like(sX)
    chB = {k: np.zeros_like(sX) for k in ('out_X', 'out_phi', 'bulk')}
    worst = 0.0
    for i, T in enumerate(Tk):
        E = 1 + T; p = math.sqrt(E * E - 1)
        for j, mu in enumerate(MU_NODES):
            kp = p * math.sqrt(max(0.0, 1 - mu * mu))
            II = I      # h = 0.01: Numerov error (k h)^4/240 <= 5e-5 even at T = 30 (where s = 1 anyway)
            r = sc.solve(II, E - Om, kp, 'X')
            sX[i, j] = 1 - r['reflected']
            q = sc.solve(II, E + Om, kp, 'Xbar')
            dB[i, j] = 1 - q['reflected']
            chB['out_X'][i, j] = q['out'].get('out_X', 0.0)
            chB['out_phi'][i, j] = q['out'].get('out_phi', 0.0)
            chB['bulk'][i, j] = q['trans_total']
            worst = max(worst, abs(r['unitarity']), abs(q['unitarity']))
    return Tk, sX, dB, chB, worst


def interp_T(Tk, tab, T):
    """log-log interpolation in kinetic energy; sqrt(T) extension below the table, constant above."""
    out = np.empty((T.size, tab.shape[1]))
    for j in range(tab.shape[1]):
        y = np.log(np.maximum(tab[:, j], 1e-300))
        lt = np.log(np.clip(T, Tk[0], Tk[-1]))
        out[:, j] = np.exp(np.interp(lt, np.log(Tk), y))
        lo = T < Tk[0]
        out[lo, j] = tab[0, j] * np.sqrt(T[lo] / Tk[0])
    return out


def thermal_average(Tk, tab, Th, stat='boltz', mu_c=0.7758647074595397, mass=1.0):
    # kinetic-energy quadrature: substitution t = T_h u^2 (dense near threshold), u in (0, 8)
    u = np.linspace(0, 1, 4001)[1:] * 8.0
    T = Th * u * u
    E = mass + T
    p2 = E * E - mass * mass
    if stat == 'boltz':
        f = np.exp(-(E - mass) / Th)
    else:
        f = 1 / np.expm1((E - mu_c) / Th)
    jac = 2 * Th * u
    w = f * p2 * jac
    sT = interp_T(Tk, tab, T)
    num = np.trapezoid(w * (sT @ (MU_W * MU_NODES)), u)
    den = np.trapezoid(w * np.sum(MU_W * MU_NODES), u)
    # also the normal-incidence and angle-averaged-at-mean-energy values
    return float(num / den)


def main(out):
    t0 = time.time()
    I = sc.Interface(h=0.01)
    Tk, sX, dB, chB, worst = tables(I)
    print('tables done (%d x %d), max |unitarity| %.1e, %.0fs' % (Tk.size, MU_NODES.size, worst, time.time() - t0), flush=True)
    # angle-averaged (flux, 2 mu dmu) sticking per energy: the geometric-optics capture efficiency of a large sphere
    s_ang = sX @ (2 * MU_W * MU_NODES)
    d_ang = dB @ (2 * MU_W * MU_NODES)
    res = dict(T_kin=Tk.tolist(), mu_nodes=MU_NODES.tolist(), mu_weights=MU_W.tolist(), s_X=sX.tolist(), d_Xbar=dB.tolist(),
               Xbar_channels={k: v.tolist() for k, v in chB.items()}, s_X_angle_avg=s_ang.tolist(), d_Xbar_angle_avg=d_ang.tolist(),
               max_unitarity=worst)
    Ths = [0.003, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]
    avg = []
    for Th in Ths:
        row = dict(T_h=Th, s_X_boltz=thermal_average(Tk, sX, Th, 'boltz'), s_X_bose_coex=thermal_average(Tk, sX, Th, 'bose'),
                   d_Xbar_boltz=thermal_average(Tk, dB, Th, 'boltz'),
                   Xbar_to_X=thermal_average(Tk, chB['out_X'], Th, 'boltz'), Xbar_to_phi=thermal_average(Tk, chB['out_phi'], Th, 'boltz'),
                   Xbar_to_bulk=thermal_average(Tk, chB['bulk'], Th, 'boltz'))
        avg.append(row)
        print('T_h=%.4g  s_X: Boltzmann %.5f Bose(mu_c) %.5f | Xbar destroyed %.5f (to X %.2e, phi %.2e, bulk %.5f)' % (
            Th, row['s_X_boltz'], row['s_X_bose_coex'], row['d_Xbar_boltz'], row['Xbar_to_X'], row['Xbar_to_phi'], row['Xbar_to_bulk']), flush=True)
    res['thermal'] = avg
    # low-energy law: s ~ c sqrt(T) at normal incidence (mu -> 1) and its coefficient
    j1 = int(np.argmax(MU_NODES))
    res['threshold_law'] = dict(c_normal=float(sX[0, j1] / math.sqrt(Tk[0])), c_angle_avg=float(s_ang[0] / math.sqrt(Tk[0])),
                                T_ref=float(Tk[0]))
    res['runtime_s'] = time.time() - t0
    json.dump(res, open(out, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main(sys.argv[1])
