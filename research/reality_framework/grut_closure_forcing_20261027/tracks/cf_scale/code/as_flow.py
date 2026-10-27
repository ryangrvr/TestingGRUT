#!/usr/bin/env python3
"""cf_scale P4: asymptotic-safety (AS) boundary conditions applied to the GRUT hidden sector (own code).

Model of the gravity contributions (standard parametrisation, e.g. Eichhorn-Held 1707.01107, de Brito et al.
2312.02086): above M_P (sharp switch) beta_lambda_i += A_lam lambda_i (universal for every quartic monomial; A_lam > 0
favoured: Shaposhnikov-Wetterich 0912.0208, Eichhorn-Pauly 2009.13543), beta_e += -f_g e (universal over gauge groups).
Below M_P gravity is off. One loop throughout (betas from beta_derive.py).

Parts
 A  gravity inputs: f_g calibrated on hypercharge (g_Y at M_Pl from Buttazzo et al. 1307.3536 eq. 61a) and the
    Eichhorn-Pauly leading-order critical exponent of the potential couplings at three quoted fixed points.
 B  P4-S (fully predictive: e at its interacting fixed point e*): real fixed points of the hidden quartics.
 C  P4-M (standard AS: quartics irrelevant, e on the Gaussian (relevant) trajectory with free e(M_P) <= e*):
    critical-surface values at M_P and at M_Y = 4.0e12 GeV, from a long transplanckian flow with two different
    initial conditions (irrelevance check).
 D  the GRUT benchmarks run up to M_P without gravity (how far off the critical surface they sit), and the
    transplanckian scale at which they would leave perturbation theory for A_lam > 0.
 E  the Standard-Model sector under the same postulate: lambda_H(M_P) from 1307.3536 eq. (61e) with PDG 2026 inputs.
Usage: python3 -B as_flow.py <out.json>
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import rg_common as R

NAMES, NC, K16 = R.NAMES, R.NC, R.K16
M_Y, H_SIG = 4.0e12, 5.33e12
OUT = {}


# ------------------------------------------------------------------ A. gravity inputs
def part_A():
    g1 = 0.6154                      # [Q] 1307.3536 eq. (61a), g1 = sqrt(5/3) gY at M_Pl
    gY = g1 * np.sqrt(3 / 5)
    bY = 41 / 6
    fg_min = gY**2 * bY / K16
    estar = lambda fg: np.sqrt(K16 * fg / R.B_E)
    th = lambda g, l: -2 * (51 - 144 * l + 104 * l**2) * g / (3 * np.pi * (3 - 4 * l)**2 * (1 - 2 * l)**2)
    fps = {'pure gravity (EP 2009.13543)': (1.32, 0.13), 'SM content (EP 2009.13543)': (7.63, -9.97),
           'U(1)_D content (de Brito et al. 2312.02086)': (2.78, -3.58)}
    A_vals = {k: -th(*v) for k, v in fps.items()}
    res = dict(gY_MPl=gY, f_g_min_from_hypercharge=fg_min, e_star_hidden_at_fg_min=estar(fg_min),
               e_star_over_gY=np.sqrt(bY / R.B_E), f_g_EV_bound_note='Eichhorn-Versteegen 1709.07252: f_g such that the IR upper '
               'bound is ~35% above the SM hypercharge; larger f_g raises e* as sqrt(f_g)',
               e_star_for_fg_x2=estar(2 * fg_min), A_lam_leading_order=A_vals,
               A_lam_formula='A = -theta, theta = -2(51-144L+104L^2) g / [3 pi (3-4L)^2 (1-2L)^2] (EP 2009.13543, quoted)')
    print('A:', {k: (v if not isinstance(v, dict) else {kk: round(vv, 4) for kk, vv in v.items()}) for k, v in res.items() if 'note' not in k}, flush=True)
    return res


# ------------------------------------------------------------------ B. P4-S fixed points (e = e*)
SUBNAMES = ['lam_X', 'lam_Y', 'lam_S', 'l_XY', 'l_XS', 'l_YS']
SUBIDX = [NAMES.index(n) for n in SUBNAMES]


def fp_eqs(z, A, e):
    y = np.zeros(NC + 1); y[SUBIDX] = z; y[NC] = e
    return R.beta(y, A_lam=A, f_g=0.0)[SUBIDX]


def stability(z, A, e, h=1e-7):
    J = np.zeros((6, 6))
    for j in range(6):
        dz = np.zeros(6); dz[j] = h
        J[:, j] = (fp_eqs(z + dz, A, e) - fp_eqs(z - dz, A, e)) / (2 * h)
    return np.linalg.eigvals(J)


def part_B(estar):
    rng = np.random.default_rng(11)
    rows = []
    for A in (0.0237, 0.1, 0.5, 0.78, 1.0, 2.83, 10.0):
        sols = []
        for _ in range(700):
            z0 = rng.normal(size=6) * 10**rng.uniform(-3, 1.5, size=6)
            z, info, ier, _ = fsolve(fp_eqs, z0, args=(A, estar), full_output=True, xtol=1e-14, maxfev=400)
            if ier == 1 and np.max(np.abs(fp_eqs(z, A, estar))) < 1e-12 and not any(np.allclose(z, s, rtol=1e-6, atol=1e-9) for s in sols):
                sols.append(z)
        out = []
        for z in sols:
            y = np.zeros(NC); y[SUBIDX] = z
            vmin, _ = R.V4_unit_min(y, n=30, fields='hidden_charged')
            ev = stability(z, A, estar)
            out.append(dict(couplings=dict(zip(SUBNAMES, z.tolist())), V4_min_unit=vmin,
                            n_irrelevant=int(np.sum(ev.real > 0)), eig_real=sorted(ev.real.tolist())))
        # analytic single-field Sigma criterion: 20 L^2 + (A K16 - 432 e^2) L + 7776 e^4 = 0
        disc = (A * K16 - 432 * estar**2)**2 - 80 * 7776 * estar**4
        rows.append(dict(A_lam=A, n_real_fixed_points=len(out), n_bounded_below=sum(o['V4_min_unit'] > 0 for o in out),
                         fixed_points=out, lamS_alone_discriminant=disc))
        print('B: A=%.4g real FPs=%d bounded-below=%d  max V4min=%s  disc_S=%.3g' % (
            A, len(out), rows[-1]['n_bounded_below'], max([o['V4_min_unit'] for o in out]) if out else None, disc), flush=True)
    Acrit = (432 * estar**2 + np.sqrt(80 * 7776) * estar**2) / K16
    return dict(e_star=estar, A_crit_real_lamS_FP=Acrit, scan=rows)


# ------------------------------------------------------------------ C. P4-M critical-surface values
def e_of_t(t, eP, fg):
    """Exact one-loop+gravity e(t), t = ln(k/M_P) >= 0, with e(0) = eP."""
    c = R.B_E / (K16 * fg)
    u = c + (1 / eP**2 - c) * np.exp(2 * fg * t)
    return 1 / np.sqrt(u)


def transplanckian(eP, A, fg, init, N=40.0):
    tUV = N / A
    y0 = np.zeros(NC + 1); y0[:NC] = init; y0[NC] = e_of_t(tUV, eP, fg)
    ev = lambda t, y: np.max(np.abs(y[:NC])) - 50.0
    ev.terminal = True
    sol = solve_ivp(lambda t, y: R.beta(y, A_lam=A, f_g=fg), (tUV, 0.0), y0, method='LSODA', rtol=1e-10,
                    atol=1e-40, max_step=max(tUV / 400, 0.5), events=ev)
    if len(sol.t_events[0]):
        return None, float(sol.t_events[0][0])
    return sol.y[:, -1], None


def run_down(yP, mu_to, mu_from=R.MP_RED):
    """Returns (couplings at mu_to, None) or (couplings at the scale where max|lambda| = 50, that scale in GeV)."""
    ev = lambda t, y: np.max(np.abs(y[:NC])) - 50.0
    ev.terminal = True
    sol = solve_ivp(lambda t, y: R.beta(y), (0.0, np.log(mu_to / mu_from)), yP, method='LSODA', rtol=1e-10, atol=1e-40,
                    events=ev)
    if len(sol.t_events[0]):
        return sol.y[:, -1], float(mu_from * np.exp(sol.t_events[0][0]))
    return sol.y[:, -1], None


def part_C(fg):
    estar = np.sqrt(K16 * fg / R.B_E)
    rows = []
    init_b = np.array([0.05, 0.3, 0.3, 0.01, 0.01, 0.01, 0.01, 0.02, 0.01, 0.01, 0.05])   # arbitrary UV start (nonzero kap, g2)
    for A in (0.0237, 0.1, 1.0, 2.83):
        for eP in (1e-7, 1e-4, 1e-2, 0.1, 0.9 * estar):
            ya, ta = transplanckian(eP, A, fg, np.zeros(NC))
            yb, tb = transplanckian(eP, A, fg, init_b)
            if ya is None or yb is None:
                rows.append(dict(A_lam=A, e_MP=float(eP), diverged=True, t_diverge_ln_k_over_MP=ta,
                                 note='quartics run to |lambda| = 50 in the transplanckian flow: no real critical surface'))
                print('C: A=%.4g e_P=%.3g DIVERGES at ln(k/M_P) = %s (no real critical surface)' % (A, eP, ta), flush=True)
                continue
            rel = float(np.max(np.abs(ya[:NC] - yb[:NC])) / max(np.max(np.abs(ya[:NC])), 1e-300))
            yM, mu_div = run_down(ya, M_Y)
            vP, nP = R.V4_unit_min(ya[:NC], n=20, fields='hidden_charged')
            vM, nM = R.V4_unit_min(yM[:NC], n=20, fields='hidden_charged')
            rows.append(dict(A_lam=A, e_MP=float(ya[NC]), at_MP=R.asdict(ya), at_MY=R.asdict(yM), V4min_MP=vP, V4min_MY=vM,
                             below_MP_lambda_reaches_50_at_GeV=mu_div,
                             init_insensitivity_maxdiff=float(np.max(np.abs(ya[:NC] - yb[:NC]))),
                             init_insensitivity_rel=rel))
            print('C: A=%.4g e_P=%.3g | MP: lamX=%.3e lamS=%.3e kap=%.1e g2=%.1e | MY: lamX=%.3e lamY=%.3e lamS=%.3e V4min=%.3e | init diff %.1e | div %s' % (
                A, ya[NC], ya[0], ya[2], ya[3], ya[7], yM[0], yM[1], yM[2], vM, rows[-1]['init_insensitivity_maxdiff'], mu_div), flush=True)
    return dict(e_star=estar, rows=rows)


# ------------------------------------------------------------------ D. benchmarks run up (no gravity)
def part_D():
    out = {}
    spot = R.vec(dict(lam_X=1e-3, lam_Y=1.0, lam_S=1.0, kap=0.0034586847211042523, g2=0.01, lam_phi=0.03, e=1e-7))
    ev = lambda t, y: np.max(np.abs(y[:NC])) - 4 * np.pi
    ev.terminal = True
    s = solve_ivp(lambda t, y: R.beta(y), (0, np.log(R.MP_RED / M_Y)), spot, method='LSODA', rtol=1e-10, atol=1e-30, events=ev)
    pole = float(M_Y * np.exp(s.t_events[0][0])) if len(s.t_events[0]) else None
    out['heavy_spot'] = dict(start_scale=M_Y, pole_4pi_GeV=pole, end=R.asdict(s.y[:, -1]), end_scale=float(M_Y * np.exp(s.t[-1])))
    # Candidate 1: X-phi sector from m_X (U(1) broken at h: e off below H_SIG; e = 1e-7 above, negligible)
    c1 = R.vec(dict(lam_X=1e-3, g2=0.01, lam_phi=0.03))
    mX = 0.0924
    s1 = solve_ivp(lambda t, y: R.beta(y, gauge_active=False), (0, np.log(R.MP_RED / mX)), c1, method='LSODA', rtol=1e-11, atol=1e-30)
    yP = s1.y[:, -1]
    out['candidate1'] = dict(start_scale=mX, at_MP=R.asdict(yP))
    # transplanckian departure for A > 0: lambda grows ~ (k/M_P)^A; scale where the largest reaches 4 pi (linear estimate)
    for tag, y in (('candidate1', yP), ('heavy_spot_kappa', np.array([0, 0, 0, 0.0034587] + [0] * 8))):
        lmax = np.max(np.abs(y[:NC]))
        out[tag + '_exit'] = {('A=%g' % A): dict(log10_k_over_MP=float(np.log10(4 * np.pi / lmax) / A)) for A in (0.0237, 0.1, 1.0, 2.83)}
    print('D:', json.dumps(out, default=float)[:900], flush=True)
    return out


# ------------------------------------------------------------------ E. SM sector (Shaposhnikov-Wetterich condition)
def part_E():
    def lamMP(mt, als, mh):        # [Q] 1307.3536 eq. (61e)
        return -0.0143 - 0.0066 * (mt - 173.34) + 0.0018 * (als - 0.1184) / 0.0007 + 0.0029 * (mh - 125.15)
    inputs = dict(mh=(125.13, 0.11, 'PDG 2026 (pdglive, fetched 2026-10-27)'), als=(0.1180, 0.0009, 'PDG 2024/2025 average as cited'),
                  mt_direct=(172.60, 0.27, 'PDG 2026 direct average (MC mass); +-0.3 GeV MC-pole ambiguity added, as 1307.3536'),
                  mt_xsec=(172.1, 0.6, 'PDG 2026 pole mass from cross sections'))
    rows = {}
    for tag, (mt, smt) in (('mt_direct', inputs['mt_direct'][:2]), ('mt_xsec', inputs['mt_xsec'][:2])):
        smt_eff = np.hypot(smt, 0.3) if tag == 'mt_direct' else smt
        l = lamMP(mt, inputs['als'][0], inputs['mh'][0])
        sig = np.sqrt((0.0066 * smt_eff)**2 + (0.0018 * inputs['als'][1] / 0.0007)**2 + (0.0029 * inputs['mh'][1])**2 + 0.0003**2)
        mh_pred = 125.15 + (0.0143 + 0.0066 * (mt - 173.34) - 0.0018 * (inputs['als'][0] - 0.1184) / 0.0007) / 0.0029
        sig_mh = np.sqrt((0.0066 * smt_eff)**2 + (0.0018 * inputs['als'][1] / 0.0007)**2 + 0.0003**2) / 0.0029
        mt_pred = 173.34 + (-0.0143 + 0.0018 * (inputs['als'][0] - 0.1184) / 0.0007 + 0.0029 * (inputs['mh'][0] - 125.15)) / 0.0066
        rows[tag] = dict(lambda_MPl=l, sigma=sig, pull_vs_zero=l / sig, pull_vs_deBrito_FP_5p56e_3=(l - 5.56e-3) / sig,
                         mh_predicted_by_lambdaMPl0=mh_pred, sigma_mh_pred=sig_mh,
                         mt_predicted_by_lambdaMPl0=mt_pred)
    print('E:', rows, flush=True)
    return dict(inputs=inputs, results=rows,
                note='retrospective; lambda(M_Pl) = 0 is the Shaposhnikov-Wetterich A_lam > 0 condition; the same near-criticality '
                     'relation as the q-vacuum P_q stability relation (HANDOFF, unreviewed)')


def main(out):
    t0 = time.time()
    OUT['method'] = __doc__
    OUT['A'] = part_A()
    fg = OUT['A']['f_g_min_from_hypercharge']
    OUT['B_P4S_fixed_points'] = part_B(OUT['A']['e_star_hidden_at_fg_min'])
    print('t after B', time.time() - t0, flush=True)
    OUT['C_P4M_critical_surface'] = part_C(fg)
    print('t after C', time.time() - t0, flush=True)
    OUT['D_benchmarks_run_up'] = part_D()
    OUT['E_SM_SW'] = part_E()
    OUT['runtime_s'] = time.time() - t0
    Path(out).write_text(json.dumps(OUT, indent=1, default=float))
    print('done', OUT['runtime_s'])


if __name__ == '__main__':
    main(sys.argv[1])
