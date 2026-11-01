#!/usr/bin/env python3
"""fx_holder track (Claude, own code) -- items 1, 2 and the static parts of item 4 of the repair brief:
the renormalizable inflaton-trilinear phase holder  L ⊃ -(g_I e^{i alpha_I} phi_I Sigma* Y^2 + h.c.),
mu_eff = mu + G phihat, G = g_I e^{i alpha_I}, phihat = phi_I/H_i.

  1. coupling: G from mu_eff(phi_I(0)) = 10 i H_i, phi_I(0) from the Friedmann constraint (engine convention)      [D, C]
  2. stationary point and Hessian at the preparation slice (lambda_Y = lambda_Sigma = 1 and 0.5); 7x7 Hessian with the
     inflaton; multi-start global search at fixed phi_I                                                              [C]
  3. one-loop Coleman-Weinberg corrections from the hidden sector: inflaton mass, tadpole, phi_I-dependent hidden masses [C]
  4. tree-level boundedness of the joint (phi_I, hidden) potential; kappa bound; hidden tree bound; large-Sigma margin  [D, C]
  5. decay kinematics and the induced six-X channel                                                                 [D, E]
  6. CP: rephasing invariant, conjugate and opposite-inflaton-branch preparations; dependence of the prepared state on
     the inflaton field value (s = phi_I/phi_I(0)) at fixed H                                                        [D, C]
Output: ../results/tri_point.json
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fx_common as fc

rng = np.random.default_rng(20261101)
MPL_GEV = 2.435e18
REF_1A = dict(amps=[87.19556213287582, 5.754247963285605, 28.38718194436162],
              hess=[13.091408829452917, 14.589249978667548, 538.1466434275553, 1546.659721973108, 3239.8852511586965],
              V=-32480.58531140619)          # [Q own m1_holder/results/holder_point.json, stationary_point.json]


def gev(x):
    return float(x * fc.H_GEV)


def physical(ev):
    return [float(v) for v in ev if abs(v) > 1e-6]


def global_search(u, p, R=12.0, nstart=300, box=200.0):
    best = []
    for _ in range(nstart):
        amps = np.exp(rng.uniform(np.log(0.3), np.log(box), 3)); ph = rng.uniform(-np.pi, np.pi, 3)
        r0 = fc.r_of_z(amps * np.exp(1j * ph))
        m = minimize(lambda r: fc.Vr(r, u, p, R), r0, jac=lambda r: fc.gradr(r, u, p, R), method='L-BFGS-B',
                     options=dict(maxiter=5000, ftol=1e-15, gtol=1e-9))
        z = fc.z_of_r(m.x); th = np.angle(z[0]); z = z * np.exp(-1j * th * np.array([1, 3, 6]))
        best.append((float(m.fun), np.abs(z).tolist(), float(fc.wrap(np.angle(z[2]) - 6 * np.angle(z[0])))))
    best.sort(key=lambda t: t[0])
    distinct = []
    for bb in best:
        if not any(abs(bb[0] - d[0]) <= 1e-6 * max(1, abs(d[0])) for d in distinct):
            distinct.append(bb)
    return [dict(V=d[0], amps=d[1], delta=d[2]) for d in distinct[:6]]


# ---------------------------------------------------------------------------------------------- one-loop CW
def V_cw(r, u, p, R, Q2):
    ev = np.linalg.eigvalsh(fc.hessr(r, u, p, R))
    ev = ev[np.abs(ev) > 1e-9]
    m4 = ev ** 2
    return float(np.sum(m4 * (np.log(np.abs(ev) / Q2) - 1.5)) / (64 * np.pi ** 2))


def cw_inflaton(r, G, ph, p, R, Q2, dph):
    """d V_CW/d phihat and d^2 V_CW/d phihat^2 at fixed hidden fields (central differences in phihat)."""
    f = lambda x: V_cw(r, p['mu'] + G * x, p, R, Q2)
    f0, fp, fm = f(ph), f(ph + dph), f(ph - dph)
    fp2, fm2 = f(ph + 2 * dph), f(ph - 2 * dph)
    d1 = (-fp2 + 8 * fp - 8 * fm + fm2) / (12 * dph)
    d2 = (-fp2 + 16 * fp - 30 * f0 + 16 * fm - fm2) / (12 * dph ** 2)
    return float(d1), float(d2)


def cw_mass_shift(r, u, p, R, Q2, step=1e-3):
    """one-loop correction to the hidden mass matrix: Hessian of V_CW in the six real fields (finite differences)."""
    Hm = np.zeros((6, 6))
    f = lambda x: V_cw(x, u, p, R, Q2)
    for i in range(6):
        for j in range(i, 6):
            ei = np.zeros(6); ei[i] = step; ej = np.zeros(6); ej[j] = step
            if i == j:
                Hm[i, i] = (-f(r + 2 * ei) + 16 * f(r + ei) - 30 * f(r) + 16 * f(r - ei) - f(r - 2 * ei)) / (12 * step ** 2)
            else:
                Hm[i, j] = Hm[j, i] = (f(r + ei + ej) - f(r + ei - ej) - f(r - ei + ej) + f(r - ei - ej)) / (4 * step ** 2)
    return Hm


def main():
    out = {}
    # ------------------------------------------------------------------ 1. coupling
    coup = {}
    for lam in (1.0, 0.5):
        p = fc.params(lam)
        st = fc.stationary(fc.U_PREP, p)
        for b in (3.0, 5.0):
            G, I0, ph0 = fc.trilinear_G(st['z'], p, b)
            coup[f'lam{lam:g}_b{b:g}'] = dict(G=[G.real, G.imag], g_I=float(abs(G)), alpha_I=float(np.angle(G)), I0_over_Mr=I0, phihat0=float(ph0),
                                              phi_I0_GeV=gev(ph0), g_I_phi0_over_H=float(abs(G) * ph0),
                                              check_u0=[float((p['mu'] + G * ph0).real), float((p['mu'] + G * ph0).imag)])
    out['coupling'] = coup
    out['coupling']['alpha_1A'] = float(np.angle(fc.LOCK_1A))
    out['coupling']['note'] = ('alpha_I = arg(10i - mu) = alpha of Option 1A because R(0) = 12 H^2 > 0 and phi_I(0) > 0 enter the same way; '
                               'g_I phi_I(0) = |10i - mu| = 12.298 H_i = |lock| * 12 of Option 1A')

    # ------------------------------------------------------------------ 2. stationary points
    sp_ = {}
    for lam in (1.0, 0.5):
        p = fc.params(lam)
        st = fc.stationary(fc.U_PREP, p)
        ev = st['eig']; phys = physical(ev)
        G, I0, ph0 = fc.trilinear_G(st['z'], p, 3.0)
        # 7x7 with the inflaton (canonical phihat, tree mass^2 b^2); mixing d^2V/(dphihat dr) = Re G dV_ur + Im G dV_ui
        mix = fc.mixing_u(st['r'], fc.U_PREP, p)
        mvec = G.real * mix[0] + G.imag * mix[1]
        H7 = np.zeros((7, 7)); H7[:6, :6] = fc.hessr(st['r'], fc.U_PREP, p); H7[6, 6] = 9.0; H7[:6, 6] = H7[6, :6] = mvec
        ev7 = np.linalg.eigvalsh(H7)
        rec = dict(amplitudes=st['amps'].tolist(), delta=st['delta'], delta_y=st['delta_y'], V=st['V'], gradient_norm=st['grad'],
                   hessian_eigenvalues=ev.tolist(), physical=phys, lightest_physical=min(phys), all_positive=bool(min(phys) > 0),
                   mixing_vector_norm=float(np.linalg.norm(mvec)), eig7=ev7.tolist(),
                   eig7_shift_of_inflaton=float(ev7[np.argmin(abs(ev7 - 9.0))] - 9.0),
                   eig7_max_shift_hidden=float(np.max(np.abs(np.sort(np.delete(ev7, np.argmin(abs(ev7 - 9.0)))) - np.sort(ev)))),
                   inflaton_force_from_hidden_over_tree=float(abs(2 * np.real(G * np.conj(st['z'][2]) * st['z'][1] ** 2)) / (9.0 * ph0)),
                   X_over_MY=float(st['amps'][0] / p['M']), Y_over_MY=float(st['amps'][1] / p['M']), S_over_MY=float(st['amps'][2] / p['M']))
        if lam == 1.0:
            rec['vs_option_1A'] = dict(amp_rel_diff=(st['amps'] / np.array(REF_1A['amps']) - 1).tolist(),
                                       hess_rel_diff=(np.array(phys) / np.array(REF_1A['hess']) - 1).tolist(), V_rel_diff=st['V'] / REF_1A['V'] - 1)
        rec['global_search_fixed_phiI0'] = global_search(fc.U_PREP, p, nstart=300)
        # held-phase (original, unheld) preparation for this lambda
        zo, info = fc.prepare_held(p, 3.0)
        rec['held_phase_preparation'] = dict(z=[[c.real, c.imag] for c in zo], amps=np.abs(zo).tolist(), residual=info['residual'], R0=info['R0'])
        sp_[f'lam{lam:g}'] = rec
    out['stationary_points'] = sp_

    # ------------------------------------------------------------------ 3. one-loop CW (hidden-sector loops, Q = M_Y)
    cw = {}
    for lam in (1.0, 0.5):
        p = fc.params(lam); Q2 = p['M'] ** 2
        st = fc.stationary(fc.U_PREP, p)
        G, I0, ph0 = fc.trilinear_G(st['z'], p, 3.0)
        rvac = fc.r_of_z(np.array([0.0, 0.0, p['h'] / np.sqrt(2)], complex))
        d1v, d2v = cw_inflaton(rvac, G, 0.0, p, 0.0, Q2, 0.01 * ph0)
        d1p, d2p = cw_inflaton(st['r'], G, ph0, p, 12.0, Q2, 0.01 * ph0)
        # phi_I-dependent one-loop hidden masses at the late vacuum: tree + CW at phihat = 0 and phihat = phi0 (R = 0)
        shifts = {}
        for nm, x in (('phihat_0', 0.0), ('phihat_phi0', ph0), ('phihat_minus_phi0', -ph0)):
            u = p['mu'] + G * x
            Ht = fc.hessr(rvac, u, p, 0.0); Hc = cw_mass_shift(rvac, u, p, 0.0, Q2)
            et = np.linalg.eigvalsh(Ht); e1 = np.linalg.eigvalsh(Ht + Hc)
            shifts[nm] = dict(u=[u.real, u.imag], tree=et.tolist(), one_loop=e1.tolist(), shift=(e1 - et).tolist())
        # induced operators: Tr M^4 (field-dependent hidden mass matrix squared, real basis) is quadratic in phihat;
        # its phihat^2 and phihat coefficients at Sigma-only and Y-only backgrounds give the one-loop running
        # dc/dlnQ = coef/(32 pi^2) of the portal counterterms lambda_IS phihat^2|Sigma|^2, lambda_IY phihat^2|Y|^2, A_IS phihat|Sigma|^2
        def trm4(z, x):
            Hm = fc.hessr(fc.r_of_z(np.asarray(z, complex)), p['mu'] + G * x, p, 0.0); return float(np.trace(Hm @ Hm))
        def coefs(z):
            sc = 1e6; f0, fp, fm = trm4(z, 0.0), trm4(z, sc), trm4(z, -sc)
            return (fp - fm) / (2 * sc), (fp - 2 * f0 + fm) / sc ** 2
        c1S, c2S = coefs([0.0, 0.0, 1.0]); c1S0, c2S0 = coefs([0.0, 0.0, 0.0])
        c1Y, c2Y = coefs([0.0, 1.0, 0.0])
        lamIS = (c2S - c2S0) / 2 / (32 * np.pi ** 2); lamIY = (c2Y - c2S0) / 2 / (32 * np.pi ** 2); AIS = (c1S - c1S0) / (32 * np.pi ** 2)
        induced = dict(TrM4_phi2_Sigma2_coef=float((c2S - c2S0) / 2), TrM4_phi2_Y2_coef=float((c2Y - c2S0) / 2), TrM4_phi_Sigma2_coef=float(c1S - c1S0),
                       expected_8_gI2=float(8 * abs(G) ** 2), expected_16_Re_mu_Gconj=float(16 * np.real(p['mu'] * np.conj(G))),
                       beta_lambda_IS=float(lamIS), beta_lambda_IY=float(lamIY), beta_A_IS=float(AIS),
                       delta_mS2_at_phi0_per_unit_lnQ=float(lamIS * ph0 ** 2), delta_mY2_at_phi0_per_unit_lnQ=float(lamIY * ph0 ** 2),
                       delta_mS2_trilinear_at_phi0_per_unit_lnQ=float(AIS * ph0),
                       ln_MPl_over_MY=float(np.log(MPL_GEV / (p['M'] * fc.H_GEV))),
                       note='|Sigma|^2 and |Y|^2 normalised to the complex fields; coefficients per unit ln Q (MSbar running), masses in H_i^2')
        # analytic Y-loop check at the vacuum: V1 = (1/64pi^2) sum_pm m_pm^4 (ln m_pm^2/Q^2 - 3/2), m^2 = M^2 + 2 lam_Y... (Y = 0) +- 2|u| s
        s = p['h'] / np.sqrt(2)
        def V1(x):
            uu = abs(p['mu'] + G * x)
            return sum((p['M'] ** 2 + sg * 2 * uu * s) ** 2 * (np.log((p['M'] ** 2 + sg * 2 * uu * s) / Q2) - 1.5) for sg in (1, -1)) / (64 * np.pi ** 2)
        dph = 0.01 * ph0
        d2an = (-V1(2 * dph) + 16 * V1(dph) - 30 * V1(0) + 16 * V1(-dph) - V1(-2 * dph)) / (12 * dph ** 2)
        g = abs(G)
        cw[f'lam{lam:g}'] = dict(
            Q='M_Y', delta_mI2_late_vacuum=d2v, delta_mI2_late_vacuum_Yloop_analytic=float(d2an), tadpole_late_vacuum=d1v,
            inflaton_shift_from_tadpole=float(-d1v / 9.0), delta_mI2_prep=d2p, tadpole_prep=d1p,
            ratio_delta_mI2_to_mI2=float(abs(d2v) / 9.0),
            estimate_gI2_MY2_over_16pi2=float(g ** 2 * p['M'] ** 2 / (16 * np.pi ** 2)),
            estimate_gI2_S2_over_4pi2=float(g ** 2 * s ** 2 / (4 * np.pi ** 2)),
            two_loop_sunset_cutoff_1e15GeV=float(g ** 2 * (1e15 / fc.H_GEV) ** 2 / (16 * np.pi ** 2) ** 2),
            induced_portal_lambda_IS_natural=float(g ** 2 / (8 * np.pi ** 2)), induced_portal_times_phi0sq=float(g ** 2 * ph0 ** 2 / (8 * np.pi ** 2)),
            induced_trilinear_phiS2_natural=float(g * p['mu'] / (8 * np.pi ** 2)),
            hidden_masses_late_vacuum=shifts, induced_operators=induced,
            portal_bound_lambda_IS_for_1H2=float(1.0 / ph0 ** 2),
            note='hidden-sector loops only (the inflaton-hidden mixing is O(1e-11)); MSbar at Q = M_Y; masses in H_i^2')
    out['one_loop'] = cw

    # ------------------------------------------------------------------ 4. bounds
    bd = {}
    for lam in (1.0, 0.5):
        p = fc.params(lam)
        kmax = (16 * p['lx'] ** 3 * p['ly'] / 27) ** 0.25
        lamY_eff = p['ly'] - 27 * p['kap'] ** 4 / (16 * p['lx'] ** 3)
        def margin(uabs):
            s_c = p['M'] ** 2 / (2 * uabs); r2 = p['h'] ** 2 / (2 * s_c ** 2)
            return float(p['ls'] * (1 - r2) ** 2 * s_c ** 2 - uabs ** 2 / lamY_eff), float(s_c), float(r2)
        st = fc.stationary(fc.U_PREP, p)
        G, I0, ph0 = fc.trilinear_G(st['z'], p, 3.0)
        # Codex chain at the preparation slice, best reserve
        chain = {}
        for res in (0.05, 0.1, 0.15, 0.2, 0.3):
            A = res * p['lx']; c = (1 - res) * p['lx']; B = p['ly'] - 27 * p['kap'] ** 4 / (16 * c ** 3)
            if B <= 0:
                chain[str(res)] = None; continue
            ax = 10.0
            chain[str(res)] = float(p['ls'] * p['h'] ** 4 / 4 - ax ** 2 / (4 * A) - (p['ls'] * p['h'] ** 2 + 100.0 / B) ** 2 / (4 * p['ls']))
        # late vacuum global search at u = mu, R = 0 (is (0, 0, h/sqrt2) the global minimum?)
        late = global_search(p['mu'] + 0j, p, R=0.0, nstart=200)
        # joint boundedness with the inflaton (tree): lambda_I >= g^4/(lambda_S lambda_Y^2); runaway |Sigma| for lambda_I = 0
        g = abs(G)
        bd[f'lam{lam:g}'] = dict(kappa_over_kappa_max_tree=float(p['kap'] / kmax), kappa_max=float(kmax), lambda_Y_eff=float(lamY_eff),
                                 large_Sigma_margin_late_mu=margin(p['mu']), large_Sigma_margin_prep_abs_u_10=margin(10.0),
                                 codex_chain_bound_prep=chain, late_vacuum_global_search=late,
                                 inflaton_quartic_needed=float(g ** 4 / (p['ls'] * p['ly'] ** 2)),
                                 runaway_Sigma_without_inflaton_quartic_over_H=float(3.0 * np.sqrt(p['ly'] / 2) / g),
                                 runaway_Sigma_GeV=gev(3.0 * np.sqrt(p['ly'] / 2) / g),
                                 one_loop_inflaton_quartic_estimate=float(g ** 4 / (16 * np.pi ** 2)))
    bd['note'] = ('kappa/kappa_max = kappa/(16 lx^3 ly/27)^(1/4) at M_Y (tree); cf_scale reports 0.83 at lambda = 0.5 including the running; '
                  'large-Sigma margin = s^2-coefficient lower bound of V for |Sigma| >= s_c = M_Y^2/(2|u|) (sufficient, not necessary)')
    out['bounds'] = bd

    # ------------------------------------------------------------------ 5. decays
    p = fc.params(1.0)
    s = p['h'] / np.sqrt(2); mL = np.sqrt(p['M'] ** 2 - 2 * p['mu'] * s); mH = np.sqrt(p['M'] ** 2 + 2 * p['mu'] * s)
    mS = np.sqrt(4 * p['ls'] * s ** 2)
    G3 = complex(*coup['lam1_b3']['G'])
    A_I = abs(G3) * s                                  # effective phi_I Y Y trilinear in the Higgs phase
    c6 = A_I * p['kap'] ** 2 / mL ** 4                 # phi_I X^6 after integrating out two Y (H_i^-3)
    def gamma6(mI):
        # massless 6-body phase space: Phi_6 = M^8 / (2 (4 pi)^9 5! 4!); |amp|^2 ~ (6! c6)^2 / 6! (identical X) ~ 720 c6^2
        Phi6 = mI ** 8 / (2 * (4 * np.pi) ** 9 * 120 * 24)
        return 720 * c6 ** 2 * Phi6 / (2 * mI)
    out['decays'] = dict(mYL=float(mL), mYH=float(mH), m_Sigma_radial=float(mS), A_I_phiYY_over_H=float(A_I),
                         phiI_to_YY_open_b3=bool(3.0 > 2 * mL), phiI_to_YY_open_b5=bool(5.0 > 2 * mL), threshold_2mYL=float(2 * mL),
                         phiI_to_SigmaYY_threshold=float(mS + 2 * mL), c6_over_H_minus3=float(c6),
                         Gamma_phiI_to_6X_b3_over_H=float(gamma6(3.0)), Gamma_phiI_to_6X_b5_over_H=float(gamma6(5.0)),
                         note='phi_I -> YY and phi_I -> Sigma Y Y closed (m_I << 2 m_YL); the open phi_I -> 6X via two virtual Y is ~1e-38 H_i [E]')

    # ------------------------------------------------------------------ 6. CP and the inflaton field value
    cp = {}
    p = fc.params(1.0)
    st = fc.stationary(fc.U_PREP, p)
    G, I0, ph0 = fc.trilinear_G(st['z'], p, 3.0)
    # rephasing check: rotate X, Y, Sigma by random phases, transform couplings, re-make kappa and mu real, read arg G
    inv = []
    for _ in range(6):
        a1, a2, a3 = rng.uniform(-np.pi, np.pi, 3)
        kap_c = p['kap'] * np.exp(1j * (3 * a1 - a2)); mu_c = p['mu'] * np.exp(1j * (2 * a2 - a3)); G_c = G * np.exp(1j * (2 * a2 - a3))
        # undo: choose b2 = a2' to make kappa real (Y rotation), then b3 to make mu real (Sigma rotation)
        b2 = np.angle(kap_c); kap_r = kap_c * np.exp(-1j * b2); mu_t = mu_c * np.exp(2j * b2); G_t = G_c * np.exp(2j * b2)
        b3 = np.angle(mu_t); mu_r = mu_t * np.exp(-1j * b3); G_r = G_t * np.exp(-1j * b3)
        inv.append(float(fc.wrap(np.angle(G_r) - np.angle(G))))
    cp['rephasing_invariant_residuals'] = inv
    cp['invariant'] = 'arg(G/mu) with kappa, mu made real = alpha_I; under the inflaton reflection phi_I -> -phi_I (a symmetry of 1/2 m^2 phi^2 and of the decay) alpha_I ~ alpha_I + pi'
    # conjugate preparation u = -10 i: stationary point is the complex conjugate (delta -> +pi/2)
    stc = fc.stationary(-10j, p)
    cp['conjugate_preparation'] = dict(amps=stc['amps'].tolist(), delta=stc['delta'], lightest=min(physical(stc['eig'])))
    # opposite inflaton branch: phi_I(0) -> -phi_I(0) with the same G: u = mu - (10i - mu) = 2 mu - 10 i
    uo = 2 * p['mu'] - 10j
    sto = fc.stationary(uo, p)
    s_ = p['h'] / np.sqrt(2)
    cp['opposite_branch'] = dict(u=[uo.real, uo.imag], abs_u=float(abs(uo)), amps=sto['amps'].tolist(), delta=sto['delta'], V=sto['V'],
                                 eig=sto['eig'].tolist(), lightest=min(physical(sto['eig'])), forced_light_Y_mass2_at_vev=float(p['M'] ** 2 - 2 * abs(uo) * s_),
                                 global_search=global_search(uo, p, nstart=200))
    # scan of the prepared state in s = phi_I/phi_I(0) at fixed H (R = 12): u(s) = mu + (10i - mu) s
    scan = []
    guess = st['z']
    for s in np.r_[np.linspace(1.0, 2.0, 21), np.linspace(1.0, 0.0, 21)[1:]]:
        if abs(s - 1.0) < 1e-12:
            guess = st['z']
        u = p['mu'] + (10j - p['mu']) * s
        r = fc.stationary(u, p, guess=guess)
        gs = global_search(u, p, nstart=40)
        guess = r['z']
        scan.append(dict(s=float(s), u=[u.real, u.imag], abs_u=float(abs(u)), amps=r['amps'].tolist(), delta=r['delta'], V=r['V'],
                         lightest=min(physical(r['eig'])), global_min_V=gs[0]['V'], global_min_amps=gs[0]['amps'], global_min_delta=gs[0]['delta'],
                         is_global=bool(abs(r['V'] - gs[0]['V']) <= 1e-6 * abs(gs[0]['V'])),
                         delta_eq_formula=float(fc.wrap(np.angle(u) + np.pi))))
    scan.sort(key=lambda d: d['s'])
    cp['phiI_scan_fixed_H'] = scan
    cp['d_delta_eq_d_ln_phiI_at_prep'] = float(np.imag((10j - p['mu']) / 10j))
    cp['m2phi2_relation'] = 'in m^2 phi^2 slow roll R = 12 H^2 ∝ phi^2, so Option 1A tracks s^2 while the trilinear tracks s'
    out['cp_and_clock'] = cp

    (fc.RES).mkdir(parents=True, exist_ok=True)
    json.dump(out, open(fc.RES / 'tri_point.json', 'w'), indent=1, default=float)
    # short summary on stdout
    print(json.dumps(out['coupling'], indent=1, default=float))
    for k, v in out['stationary_points'].items():
        print(k, v['amplitudes'], v['delta'], v['physical'], v['eig7_shift_of_inflaton'], v['eig7_max_shift_hidden'], v.get('vs_option_1A'), v['global_search_fixed_phiI0'][:2])
    for k, v in out['one_loop'].items():
        print(k, {a: b for a, b in v.items() if not isinstance(b, (dict, list))})
        print('   induced', v['induced_operators'])
        print('   vac', {n: (w['tree'], w['shift']) for n, w in v['hidden_masses_late_vacuum'].items()})
    print(json.dumps(out['bounds'], indent=1, default=float)[:4000])
    print(json.dumps(out['decays'], indent=1, default=float))
    print('CP', cp['rephasing_invariant_residuals'], cp['conjugate_preparation'], {k: v for k, v in cp['opposite_branch'].items() if k != 'global_search'}, cp['opposite_branch']['global_search'][:3])
    for d in scan:
        print('  s=%.2f |u|=%.2f amps=%s delta=%.4f eq=%.4f lightest=%.3f global=%s gV=%.2f V=%.2f' % (d['s'], d['abs_u'], np.round(d['amps'], 3), d['delta'], d['delta_eq_formula'], d['lightest'], d['is_global'], d['global_min_V'], d['V']))


if __name__ == '__main__':
    main()
