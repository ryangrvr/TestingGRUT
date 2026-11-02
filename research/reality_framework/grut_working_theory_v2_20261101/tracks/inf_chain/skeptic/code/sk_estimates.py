#!/usr/bin/env python3
"""inf_chain SKEPTIC (Claude, own code): order-of-magnitude estimates [E] of effects the track omitted for the quiet
T-models (TB1 alpha = 0.3; numbers for TB2 alongside).  Units H_i = H(S1) = 1, a = 1 at S1, M_P = 2.435e18 GeV.

 (a) Fermion preheating at the zero crossings while phi -> NN is kinematically blocked (y Phi > m_I/2):
     Greene-Kofman n_psi = (y m Phi)^{3/2}/(4 pi^3) per crossing (2 spin states, harmonic phidot = m Phi), energy y Phi
     each at the turning point: fraction per crossing f_c = y^{5/2} (Phi/m)^{1/2}/(2 pi^3); compared with the
     perturbative fraction per period, 2 pi Gamma/m.  Amplitude Phi from my own T-model oscillation (first turning
     point after eps_H = 1), Phi ∝ a^{-3/2} afterwards.
 (b) Gravitational (oscillating-R) production during the matter era: R_osc = 9 H^2 cos 2mt (quadratic oscillator, rho = 3H^2);
     X (xi_X = -5/6 in the engine's sign, m^2 = -xi R, conformal at -1/6) couples to (2/3) R_osc; the Goldstone pi
     (minimal) to R_osc/6.  Pair rate A^2/(32 pi) (complex) or A^2/(64 pi) (real) per volume; injected X energy
     E_inj = int m A^2/(16 pi) a^4 dt; number of pi = int A^2 a^3/(64 pi) dt (each pi = one long-lived A_L, deposit as
     the relic vectors).
 (c) The S1-channel R-odd induced A phi_I |H|^2 = y y_nu^2 M_N/8 pi^2 x [1, ln(M_P^2/M_N^2)] (inf_reheat skeptic form) at
     the TB couplings, seesaw y_nu^2 = m_nu M_N/v^2 (m_nu = 0.05 eV, v = 174 GeV, M_N = m_I/4), against inf_higgs's
     inflationary bound |A phi_I| <~ 0.75 H_*^2 on the plateau.
Output ../results/sk_estimates.json.
"""
import json, math, sys, pathlib
import numpy as np
from scipy.integrate import solve_ivp
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sk_cmb_fluid as sk

MP = 2.435e18


def first_amplitude(m, sol, Nend):
    """first turning point |phi| after eps_H = 1 (H_i units, a = 1 at S1), and a there."""
    f, fN = sol.sol(Nend - 1.0); Hs, _ = sk.H_of(m, f, fN); sc = 1 / Hs ** 2
    def rhs(t, y):
        f, p, lna = y
        rho = 0.5 * p * p + sc * m.V(f)
        return [p, -3 * math.sqrt(rho / 3) * p - sc * m.dV(f), math.sqrt(rho / 3)]
    ev = lambda t, y: y[1]
    ev.direction = 0
    s = solve_ivp(rhs, (0, 60), [f, fN, 0.0], method='DOP853', rtol=1e-11, atol=1e-14, events=ev, dense_output=True)
    tp = list(s.t_events[0])                         # phidot = 0: the first is the first turning point after inflation
    y = s.sol(tp[0])
    return abs(y[0]), math.exp(y[2]), tp[0]


def main():
    cmb = json.load(open(HERE.parent / 'results' / 'sk_cmb_fluid.json'))
    out = {}
    for key, alpha in (('TB1', 0.3), ('TB2', 1.0)):
        m = sk.TModel(alpha)
        sol, Nend = sk.background(m, m.phi_start(80.0))
        Phi1, a1, t1 = first_amplitude(m, sol, Nend)
        r = cmb[key]
        Hi = r['H_i_GeV']; mI = r['m_I_over_Hi']; gam = r['Gamma_over_Hi']; x0 = r['x0']; y = r['y_majorana']
        # (a)
        fc = y ** 2.5 * (Phi1 * MP / (mI * Hi)) ** 0.5 / (2 * math.pi ** 3)
        fpert = 2 * math.pi * gam / mI
        block1 = y * Phi1 * MP / (0.5 * mI * Hi)
        a_unblock = a1 * block1 ** (2 / 3)
        # integrated NP fraction until unblocking: crossings 2 per period, f_c ∝ a^{-3/4}, dt = da/(aH), H = (x0/3)^1/2 a^-3/2
        aa = np.geomspace(a1, a_unblock, 4000)
        H = math.sqrt(x0 / 3) * aa ** -1.5
        rate = (2 * mI / (2 * math.pi)) * fc * (aa / a1) ** -0.75      # per unit time
        F_np = float(np.trapezoid(rate / (aa * H), aa))
        F_pert = float(np.trapezoid(gam / (aa * H), aa))
        # (b) gravitational production, matter era from a = e
        ae = math.e
        H1 = math.sqrt(x0 / 3)
        # X energy: dE/da = m A^2/(16 pi) a^4/(a H), A = (2/3) 9 H^2
        EX_inj = mI * (6.0) ** 2 / (16 * math.pi) * H1 ** 3 * 2 / ae ** 0.5          # int a^{-3/2} da from ae
        # pi number: dN/da = (1.5 H^2)^2/(64 pi) a^3/(a H) = 2.25 H^3 a^2/(64 pi)
        Npi = 2.25 / (64 * math.pi) * H1 ** 3 * (2 / 3) * ae ** -1.5
        rel = json.load(open(HERE.parent.parent / 'results' / ('release_%s.json' % key)))['result']
        vec = json.load(open(HERE.parent.parent / 'results' / ('vectors_%s_g%.3e.json' % (key, {'TB1': 7.512e-9, 'TB2': 2.113e-9}[key]))))
        v6 = [x for x in vec['by_e'] if abs(x['e'] - 1e-6) < 1e-12][0]
        R_grav_X = EX_inj / rel['E_X']
        R_grav_pi = Npi * v6['deposit_per_vector'] / rel['E_X']
        # (c) induced A-term
        MN = 0.25 * mI * Hi
        ynu2 = 0.05e-9 * MN / 174.0 ** 2
        A = y * ynu2 * MN / (8 * math.pi ** 2)
        A_ln = A * math.log(MP ** 2 / MN ** 2)
        f55, _ = sol.sol(Nend - 55.0)
        Hst = r['H_star_GeV']
        bound = 0.75 * Hst ** 2 / (f55 * MP)
        out[key] = dict(Phi_first_MP=Phi1, a_first=a1, block_factor_first=block1, a_unblock=a_unblock,
                        N_unblock_after_S1=math.log(a_unblock), frac_per_crossing_np=fc, frac_per_period_pert=fpert,
                        F_np_until_unblock=F_np, F_pert_until_unblock=F_pert,
                        EX_grav_injected=EX_inj, R_grav_X=R_grav_X, N_pi_grav=Npi, R_grav_pi_e1e6=R_grav_pi, N_A_relic=v6['N_A'],
                        M_N_GeV=MN, y_nu2=ynu2, A_ind_GeV=[A, A_ln], A_bound_GeV=bound, phi_plateau_N55_MP=f55)
        print(key, out[key], flush=True)
    json.dump(out, open(HERE.parent / 'results' / 'sk_estimates.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
