#!/usr/bin/env python3
"""inf_reheat SKEPTIC (Claude, own code): channel-level checks by other methods, and the omitted effects.

Sections (keys of ../results/sk_channels.json):
  gk            Greene-Kofman production after one zero crossing by an exact 2x2 matrix-exponential propagator (vectorised
                over k; piecewise-constant mass at the midpoint), not an RK integration as in fermion_modes.py [C]
  width         Majorana width from an explicit two-component spinor sum for -(y/2) phi Nbar N (vertex y, identical 1/2) [C]
  s1_odd        S1 soft breaking: R-odd one-loop inflaton terms (tadpole + cubic from (M + y phi)^4 ln) against the tree slope
                at phi_*, and the log-divergent induced A phi |H|^2 ~ y y_nu^2 M_N ln/(8 pi^2), at B1-B3 and F0 [C/E]
  portal        lambda_IH conventions (mass^2 lambda phi^2 vs (lambda/2) phi^2 |H|^2), the window, and the condensate annihilation
                phi phi -> H H* through the oscillating portal mass as an early Higgs-coupled plasma [C/E]
  s2_eft        S2 axial EFT: mu5_max/(4 pi f), phidot/f^2, rho^{1/4}/(4 pi f), winding, xi_psi during inflation [C]
  optionB       Option B numbers recomputed; the stability bounds at the benchmark H_i [C]
  cmb           n_s of B1-B3 (track and Mukhanov-Sasaki) against Planck 2018, CMB-SPA (SPT-3G D1) and P-ACT-LB (ACT DR6) [Q/C]
  p_feedback    the omitted symmetric component: with H_i = H_end fixed by the inflaton, the DM calibration moves Gamma by
                (1+2p)^{-2}, not H_i by (1+2p)^{-2/3}; consequences for T_RH, y, K, S/Q [D/E]
  sq_rolling    S/Q(1e-6) with the rolling-slice x0 (sk_tmodel) in the track's own rescaling [E]
Usage: sk_channels.py (after sk_tmodel.py) -> ../results/sk_channels.json
"""
import json, math, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
TRK = HERE.parent.parent / 'results'
MP = 2.435e18
VEW = 174.1
MNU = 0.05e-9


# ------------------------------------------------------------------------------------------------ GK by matrix exponentials
def gk(yPhi, m, nk=60, steps_per=0.004):
    kgk = math.sqrt(yPhi * m / math.pi)
    ks = np.linspace(0.0, 3.0 * kgk, nk)
    T = math.pi / m
    n = int(math.ceil(T / (steps_per / yPhi)))
    dt = T / n
    # start in the + eigenvector of [[M, k], [k, -M]] at t = 0 (M = yPhi > 0)
    w0 = np.hypot(ks, yPhi)
    u = np.stack([w0 + yPhi, ks], axis=1).astype(complex)
    u /= np.linalg.norm(u, axis=1)[:, None]
    for j in range(n):
        Mt = yPhi * math.cos(m * (j + 0.5) * dt)
        w = np.sqrt(ks * ks + Mt * Mt)
        c = np.cos(w * dt); s = np.sin(w * dt) / np.where(w > 0, w, 1.0)
        u0 = c * u[:, 0] - 1j * s * (Mt * u[:, 0] + ks * u[:, 1])
        u1 = c * u[:, 1] - 1j * s * (ks * u[:, 0] - Mt * u[:, 1])
        u = np.stack([u0, u1], axis=1)
    Mt = yPhi * math.cos(m * T)                       # = -yPhi
    w = np.sqrt(ks * ks + Mt * Mt)
    # negative-frequency eigenvector of [[Mt, k], [k, -Mt]]: (k, -(w + Mt))... use (Mt - w, k)
    vm = np.stack([Mt - w, ks], axis=1).astype(complex)
    vm /= np.linalg.norm(vm, axis=1)[:, None]
    nk_num = np.abs(np.sum(np.conj(vm) * u, axis=1)) ** 2
    ngk = np.exp(-math.pi * ks * ks / (yPhi * m))
    N = 2 * np.trapezoid(ks * ks * nk_num, ks) / (2 * math.pi ** 2)
    Nan = (yPhi * m) ** 1.5 / (4 * math.pi ** 3)
    return dict(yPhi=yPhi, m=m, steps=n, max_abs_dev=float(np.max(np.abs(nk_num - ngk))), N_ratio=float(N / Nan),
                n_at_k0=float(nk_num[0]))


# ------------------------------------------------------------------------------------------------ Majorana width by spinors
def width_spinor(mphi=3.0, m=0.75):
    """-(y/2) phi Nbar N, N Majorana: vertex -i y; |M|^2 summed over the two spins = y^2 sum |ubar(p1,s) v(p2,s')|^2.
    Explicit Dirac-basis spinors, spin basis along z.  Two-body: Gamma = (1/2)(identical) sum|M|^2 p / (8 pi mphi^2)."""
    E = mphi / 2; p = math.sqrt(E * E - m * m)
    g0 = np.diag([1, 1, -1, -1]).astype(complex)
    def u(pv, s):
        chi = np.array([1, 0], complex) if s == 0 else np.array([0, 1], complex)
        sig = np.array([[pv[2], pv[0] - 1j * pv[1]], [pv[0] + 1j * pv[1], -pv[2]]], complex)
        En = math.sqrt(m * m + pv @ pv)
        return math.sqrt(En + m) * np.concatenate([chi, sig @ chi / (En + m)])
    def v(pv, s):
        chi = np.array([0, 1], complex) if s == 0 else np.array([-1, 0], complex)
        sig = np.array([[pv[2], pv[0] - 1j * pv[1]], [pv[0] + 1j * pv[1], -pv[2]]], complex)
        En = math.sqrt(m * m + pv @ pv)
        return math.sqrt(En + m) * np.concatenate([sig @ chi / (En + m), chi])
    p1 = np.array([0, 0, p]); p2 = -p1
    S = 0.0
    for s1 in (0, 1):
        for s2 in (0, 1):
            ub = np.conj(u(p1, s1)) @ g0
            S += abs(ub @ v(p2, s2)) ** 2
    beta = p / E
    G_over_y2 = 0.5 * S * p / (8 * math.pi * mphi ** 2)    # identical-particle 1/2; |M|^2 = y^2 S
    return dict(sum_spins=S, expected_2mphi2beta2=2 * mphi ** 2 * beta ** 2, Gamma_over_y2m=G_over_y2 / mphi,
                formula_beta3_over_16pi=beta ** 3 / (16 * math.pi))


# ------------------------------------------------------------------------------------------------ S1 R-odd one-loop terms
def F(mN, mu):
    return mN ** 4 * (math.log(mN * mN / (mu * mu)) - 1.5)


def s1_odd(b, y, phis, Vs, eps, MPh, MN, mus):
    """V_CW = -(1/32 pi^2) F(M + y phi) (Majorana, 2 dof); odd part slope at phi_* against the tree slope sqrt(2 eps) V/M_P."""
    tree = math.sqrt(2 * eps) * Vs / MPh
    out = {}
    for name, mu in mus.items():
        h = 1e-4 * phis
        Vodd = lambda ph: -(F(MN + y * ph, mu) - F(abs(MN - y * ph), mu)) / (2 * 32 * math.pi ** 2)
        sl = (Vodd(phis + h) - Vodd(phis - h)) / (2 * h)
        out[name] = dict(mu=mu, odd_slope=sl, ratio_to_tree=sl / tree)
    even = lambda ph: -(F(MN + y * ph, phis * y) + F(abs(MN - y * ph), phis * y)) / (2 * 32 * math.pi ** 2)
    out['V_CW_even_over_Vtree_at_phistar_mu_yphi'] = abs(even(phis)) / Vs
    out['odd_over_even_size'] = 4 * MN / (y * phis)
    out['tree_slope'] = tree
    return out


def a_induced(y, MN, Hi_GeV):
    ynu2 = MNU * MN * Hi_GeV / VEW ** 2
    base = y * ynu2 * MN / (8 * math.pi ** 2)
    L = math.log((MP / Hi_GeV) ** 2 / MN ** 2)
    return dict(ynu2=ynu2, A_no_log=base, A_log_MP=base * L, log_MP=L)


# ------------------------------------------------------------------------------------------------ portal: condensate annihilation
def portal(fam, sol, MPh):
    m = fam['b']; Phi1 = fam['Phi1'] * MPh
    lmin = sol['lambda_IH']['lam_min']; lmax = sol['lambda_IH']['lam_max_q1']
    rows = []
    for lam in (lmin, 2 * lmin, 4 * lmin, 0.5 * lmax, lmax):
        q = lam * Phi1 ** 2 / (4 * m * m)                    # track's convention: Higgs mass^2 = lam phi^2
        beta = math.sqrt(max(0.0, 1 - 2 * q))               # <m_h^2> = lam Phi^2/2 against (omega/2)^2 = m^2
        Geff = lam ** 2 * Phi1 ** 2 * beta / (16 * math.pi * m)     # rho-loss rate / rho_phi, 4 real Higgs dof
        # SM energy injected over the first Hubble time and the hard-gas induced mass^2 (track's coefficient 0.018 rho/E^2, E ~ m)
        rho_phi = 0.5 * m * m * Phi1 ** 2
        dE = Geff * rho_phi * 1.0
        m2 = 0.018 * dE / m ** 2
        T = (30 * dE / (math.pi ** 2 * 106.75)) ** 0.25
        rows.append(dict(lam_track=lam, lam_lagrangian_half=2 * lam, q=q, beta=beta, Gamma_eff_over_Hi=Geff,
                         Gamma_eff_over_GammaI=Geff / sol['Gamma_over_Hi'], SM_energy_first_Hubble_time=dE, m2_induced=m2, T_inst=T))
    return dict(rows=rows, lam_open_max=2 * m * m / Phi1 ** 2,
                note='Track convention: (lambda/2) phi^2 |H|^2 is said to give m_H^2 = lambda phi^2; for a canonical complex doublet '
                     'it gives lambda phi^2/2. The window in the Lagrangian coefficient of (lambda/2) phi^2 |H|^2 is twice the quoted one.')


def s2_eft(fam, sol, MPh, Hstar_over_Hi):
    f = sol['f_over_Hi']; m = fam['b']; Phi1 = fam['Phi1'] * MPh
    mu5 = sol['axial_nonperturbative']['mu5_max']
    phidot_end = math.sqrt(2.0) * MPh                       # epsilon_H = 1: phidot = sqrt(2) H M_P, H = H_i = 1
    rho_end = 3 * MPh ** 2
    xi_cmb = math.sqrt(2 * fam['eps_star']) * MPh * Hstar_over_Hi / (2 * f * Hstar_over_Hi)
    xi_end = phidot_end / (2 * f)
    # description-B anomaly coupling (rotating L by phi/2f, 3 doublets): 1/f_A = 3 alpha_2/(4 pi f); xi_W = phidot/(2 f_A H)
    a2 = 1 / 40.0
    fA = f / (3 * a2 / (4 * math.pi))
    return dict(f=f, four_pi_f=4 * math.pi * f, mu5_max=mu5, mu5_over_4pif=mu5 / (4 * math.pi * f),
                phidot_end_over_f2=phidot_end / f ** 2, phidot_end_over_4pif_sq=phidot_end / (4 * math.pi * f) ** 2,
                rho_end_quarter_over_4pif=rho_end ** 0.25 / (4 * math.pi * f), Phi1_over_f=Phi1 / f,
                winding_cmb_to_end=(fam['phi_star'] - fam['phi_end']) * MPh / (2 * math.pi * f), xi_psi_cmb=xi_cmb, xi_psi_end=xi_end,
                xi_W_descB_cmb=math.sqrt(2 * fam['eps_star']) * MPh / (2 * fA), xi_W_descB_end=phidot_end / (2 * fA))


def optionB(b=3.0, Gam=0.5, MPh=MP / 1.291e11):
    A = math.sqrt(8 * math.pi * b * Gam)
    Phi = math.sqrt(6) / b * MPh
    from scipy.integrate import quad
    I = quad(lambda x: math.sqrt(math.cos(x)), -math.pi / 2, math.pi / 2)[0]
    return dict(A=A, mH2=A * Phi, growth=I * math.sqrt(A * Phi) / b, dlam=A * A / (2 * b * b), int_sqrt_cos=I)


def main():
    trk = json.load(open(TRK / 'benchmarks.json'))
    skt = json.load(open(RES / 'sk_tmodel.json'))
    res = dict(method=__doc__)
    res['gk'] = [gk(60.0, 3.0), gk(300.0, 3.0)]
    print('GK', res['gk'], flush=True)
    res['width'] = width_spinor()
    print('width', res['width'])
    res['optionB'] = optionB()
    print('optionB', res['optionB'])
    bms = {(bm['family'], bm['channel']): bm for bm in trk['benchmarks']}
    res['s1_odd'] = []; res['portal'] = []; res['s2_eft'] = []; res['cmb'] = []; res['p_feedback'] = []; res['sq_rolling'] = []
    data = dict(planck2018=(0.9649, 0.0042), cmb_spa_spt3g=(0.9684, 0.0030), p_act_lb=(0.9743, 0.0034))
    for alpha, rr in zip((0.3, 1.0, 3.0), skt['rows']):
        bm = bms[('T(alpha=%g)' % alpha, 'S1')]; s = bm['solution']; fam = bm['fam']
        bm2 = bms[('T(alpha=%g)' % alpha, 'S2')]; s2 = bm2['solution']
        Hi = fam['H_end']; MPh = MP / Hi
        b = fam['b']; MN = b / 4
        phis = fam['phi_star'] * MPh; Vs = fam['V_star'] * MPh ** 4
        o = s1_odd(b, s['y'], phis, Vs, fam['eps_star'], MPh, MN, dict(mu_MN=MN, mu_yphistar=s['y'] * phis))
        o.update(alpha=alpha, A_ind=a_induced(s['y'], MN, Hi), A_max_inflation_branch_no_portal=0.5 / (fam['phi_end'] * MPh),
                 lam_min_phi_slice=s['lambda_IH']['lam_min'] * fam['phi_end'] * MPh)
        res['s1_odd'].append(o)
        print('S1 odd alpha=%g: ratio(mu=M)=%.3g ratio(mu=y phi*)=%.3g  V_CW/V=%.3g  A_ind=%.2e..%.2e  A_max(no portal)=%.2e  lam phi_slice=%.2e' % (
            alpha, o['mu_MN']['ratio_to_tree'], o['mu_yphistar']['ratio_to_tree'], o['V_CW_even_over_Vtree_at_phistar_mu_yphi'],
            o['A_ind']['A_no_log'], o['A_ind']['A_log_MP'], o['A_max_inflation_branch_no_portal'], o['lam_min_phi_slice']), flush=True)
        p = portal(fam, s, MPh); p['alpha'] = alpha
        res['portal'].append(p)
        for r_ in p['rows']:
            print('  portal lam=%.2e q=%.3f beta=%.2f Geff=%.2e (x Gamma_I %.2e) dE=%.2e m2=%.2e T=%.2f' % (
                r_['lam_track'], r_['q'], r_['beta'], r_['Gamma_eff_over_Hi'], r_['Gamma_eff_over_GammaI'], r_['SM_energy_first_Hubble_time'],
                r_['m2_induced'], r_['T_inst']))
        e = s2_eft(fam, s2, MPh, fam['H_star'] / Hi); e['alpha'] = alpha
        res['s2_eft'].append(e)
        print('  S2 EFT', {k: round(v, 4) if isinstance(v, float) else v for k, v in e.items()}, flush=True)
        c = dict(alpha=alpha, ns_track=fam['ns'], ns_MS=rr['ns_MS'])
        # e-fold correction: track N_* is too large by 0.39 (missing -(1/12) ln g) + sudden-decay rho_th; use mine (rolling K)
        # self-consistently K is fixed by the DM calibration on either slice, so N_* follows from that K (my bookkeeping)
        dN = rr['N_star_mine_from_trackK'] - rr['N_star_track']
        c['dN_bookkeeping'] = dN
        c['ns_corrected'] = rr['ns_MS'] + 2 * dN / rr['N_star_track'] ** 2       # n_s ~ 1 - 2/N
        for k, (mu, sg) in data.items():
            c['pull_' + k] = (c['ns_corrected'] - mu) / sg
        res['cmb'].append(c)
        print('  CMB', c, flush=True)
        # p feedback: Gamma -> Gamma (1+2p)^-2 at fixed H_i = H_end; K -> K (1+2p)^{2/3}; T_RH ~ Gamma^{1/2}; y_S1 ~ Gamma^{1/2}
        for pp in (0.01, 0.1, 1.0, 10.0, 30.0):
            fct = (1 + 2 * pp)
            G = s['Gamma_over_Hi'] / fct ** 2
            K = s['K'] * fct ** (2 / 3)
            Rv = s['SQ_e1e-6']['R_vec'] * math.sqrt(fct ** (2 / 3))
            W = s['SQ_e1e-6']['E_over_Ec'] * (1 + Rv)
            SQ = 0.190255 * (W - 0.97494) ** 0.75013 - 0.00428 + trk['SQ_offset']
            res['p_feedback'].append(dict(alpha=alpha, p=pp, Gamma_over_m=G / b, T_RH_GeV=s['T_RH_GeV'] / fct, K=K, y_S1=s['y'] / fct,
                                          y_over_ynat=s['y_over_ynat'] / fct, z_RH=s['z_RH'] * fct, SQ=SQ,
                                          N_matter_era=s['N_matter_era'] + (2 / 3) * math.log(fct ** 2)))
        # S/Q with the rolling-slice x0 (matter-era vector crossing a_*^3 = x0/(3 m_A^2)); K fixed by the calibration
        x0s = trk['x0']; kb = sorted(float(k) for k in x0s)
        x0_static = float(np.exp(np.interp(math.log(b), np.log(kb), np.log([x0s[str(k) if str(k) in x0s else '%g' % k] for k in kb]))))
        fac = rr['x0_rolling_over_He2'] / x0_static
        Rv = s['SQ_e1e-6']['R_vec'] * fac
        W = s['SQ_e1e-6']['E_over_Ec'] * (1 + Rv + s['SQ_e1e-6']['R_Y'] + s['SQ_e1e-6']['R_casc'])
        SQ = 0.190255 * (W - 0.97494) ** 0.75013 - 0.00428 + trk['SQ_offset']
        Rv2 = Rv * 0.40     # with the decisive pass's rolling horizon-exit factor f_eff ~ 0.35-0.44 / 0.47-equivalent (illustrative)
        W2 = s['SQ_e1e-6']['E_over_Ec'] * (1 + Rv2 + s['SQ_e1e-6']['R_Y'] + s['SQ_e1e-6']['R_casc'])
        SQ2 = 0.190255 * (W2 - 0.97494) ** 0.75013 - 0.00428 + trk['SQ_offset']
        # induced A_ISigma phi |Sigma|^2 (natural size, Y loop with the holder and mu): phi -> A_L A_L through sigma mixing
        gI = abs(10j - 7.159456) / (fam['phi_end'] * MPh)
        ReGmu = gI * 7.159456 * math.cos(2.192144074)
        A0 = 2 * abs(ReGmu) / (16 * math.pi ** 2)
        msig = 40.0; MY = 30.0
        EXb = float(np.interp(b, [2.0434, 2.1735, 2.477, 3.2142], [297896.5, 300068.5, 285787.1, 184116.4]))
        leak = dict(g_I=gI, ReGmu=ReGmu)
        for nm, L in (('vanish_at_MY_log_MY_over_mI', math.log(MY ** 2 / b ** 2)), ('vanish_at_MP_log_MP_over_MY', math.log(MPh ** 2 / MY ** 2))):
            A = A0 * L
            G = A * A * b ** 3 / (32 * math.pi * msig ** 4)
            leak[nm] = dict(A_ISigma=A, Gamma_hid=G, BR=G / s['Gamma_over_Hi'], R_inf=3 * MPh ** 2 * s['K'] ** 2 * (G / s['Gamma_over_Hi']) / EXb)
        res.setdefault('induced_sigma_leak', []).append(dict(alpha=alpha, **leak))
        # S1 box: phi phi -> N N -> L H L H gives the R-even portal lambda_ind ~ y^2 y_nu^2/(16 pi^2) x [1, ln(M_P^2/M_N^2)]
        ynu2 = MNU * MN * Hi / VEW ** 2
        lam0 = s['y'] ** 2 * ynu2 / (16 * math.pi ** 2)
        Lg = math.log(MPh ** 2 / MN ** 2)
        # S1 at the rolling-slice self-consistent Gamma (y ~ Gamma^{1/2})
        shift = rr['Gamma_shift']
        res.setdefault('s1_induced_portal', []).append(dict(alpha=alpha, ynu2=ynu2, lam_ind_nolog=lam0, lam_ind_logMP=lam0 * Lg,
                                                           window=[s['lambda_IH']['lam_min'], s['lambda_IH']['lam_max_q1']],
                                                           lam_ind_logMP_at_rolling_Gamma=lam0 * Lg * shift))
        print('  S1 induced portal', res['s1_induced_portal'][-1], flush=True)
        # Pauli blocking of phi -> N N by thermal N near the end of reheating: factor tanh(m_I/(4T)), T ~ T_RH a^{-3/8}
        TRH = s['T_RH_GeV'] / Hi
        blk = {('%.1f' % dN): math.tanh(b / (4 * TRH * math.exp(0.375 * dN))) for dN in (0.0, 1.0, 2.0, 3.0)}
        res.setdefault('pauli', []).append(dict(alpha=alpha, T_RH_over_Hi=TRH, block_factor_vs_efolds_before_RH=blk))
        print('  Pauli', res['pauli'][-1], flush=True)
        # self-consistent values on the rolling slice (K fixed by the calibration; Gamma -> Gamma x shift)
        sh = shift
        s2y = s2['y_axial'] * math.sqrt(sh); f2 = MN / s2y
        mu5 = s2['axial_nonperturbative']['mu5_max'] * (s2['f_over_Hi'] / f2)
        res.setdefault('rolling_selfconsistent', []).append(dict(
            alpha=alpha, Gamma_shift=sh, Gamma_over_m=s['Gamma_over_m'] * sh, T_RH_GeV=s['T_RH_GeV'] * math.sqrt(sh),
            y_S1=s['y'] * math.sqrt(sh), y_S1_over_ynat=s['y_over_ynat'] * math.sqrt(sh), y_axial=s2y, f_GeV=f2 * Hi,
            mu5_max_over_4pif=mu5 / (4 * math.pi * f2), z_RH=s['z_RH'] / math.sqrt(sh), N_matter_era=s['N_matter_era'] + (2 / 3) * math.log(1 / sh),
            K=s['K'], H_i=Hi))
        print('  rolling self-consistent', res['rolling_selfconsistent'][-1], flush=True)
        print('  A_ISigma leak', leak, flush=True)
        res['sq_rolling'].append(dict(alpha=alpha, x0_static=x0_static, x0_rolling=rr['x0_rolling_over_He2'], R_vec_factor=fac, R_vec=Rv,
                                      SQ=SQ, R_vec_with_feff_0p4=Rv2, SQ_with_feff_0p4=SQ2, SQ_track=s['SQ_e1e-6']['SQ']))
        print('  SQ rolling', res['sq_rolling'][-1], flush=True)
    for r_ in res['p_feedback']:
        print('  p-feedback', {k: (('%.3g' % v) if isinstance(v, float) else v) for k, v in r_.items()})
    json.dump(res, open(RES / 'sk_channels.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
