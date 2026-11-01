"""fx_chain items 1-2 (Claude, own code): heat budget at e = 1e-6, the relaxed hidden plasma before the phi onset, the
onset condensate density and the M1 homogeneous-condensate (tilted) crossover along the history.

Reuse, stated:
  * q_state / q_solve: copied from tracks/pt_relax/skeptic/code/equilibria.py (our own quantum Hartree solver with a
    circular condensate, Gauss-Legendre quadrature; checked there against pt_relax's independent solver to 4 digits).
  * R values, decay clocks and the release (N_X, E_X): read from tracks/kin_source/results/source_interface.json.
  * M1 homogeneous condensate e(S, F) = nbar^2/2F^2 + U(S, F): the conventions of tracks/m1_vacuum (U in g^2 v^4,
    n/m_X^3 = 100 Omega F^2, a = 0.1, b = 3, c = 0.2); re-implemented here, and its bulk/spinodal numbers reproduced.
Conformal H_i units for the pre-onset plasma (psi canonical, r^2 = chi1^2 + chi2^2, V = lX r^4/4 + g^2 phi^2 r^2/2).
Labels: [D] derived, [C] computed, [E] estimate, [A] argued.
"""
import json, math, os, sys
import numpy as np
from scipy.optimize import brentq, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'results')
TR = '/home/user/TestingGRUT/research/reality_framework/grut_formation_review_20261025/tracks'
SRC = json.load(open(TR + '/kin_source/results/source_interface.json'))

LX, G2, LP, EPST = 1e-3, 1e-2, 3e-2, 0.002
FP = SRC['meta']['fixed_point']
H = FP['H_GeV']; V = FP['v_GeV']; MX = math.sqrt(G2) * V
KAPPA = LP * V * V / H ** 2                     # lambda_phi v^2 / H_i^2 (tachyonic coefficient of a^2), untilted
NX_HEAVY, EX_HEAVY = abs(FP['N_X']), FP['E_X']
NQ_PILOT, RHO_PILOT = 8112.665383318339, 105481.24746711778
A_REF = SRC['release_population']['onset']['a_ref_from_charge']

# ------------------------------------------------------------------ quantum Hartree (copied from pt_relax skeptic)
XG, WG = np.polynomial.legendre.leggauss(400)


def _grid(T, Mmax):
    kmax = 80 * T + 20 * Mmax
    return 0.5 * kmax * (XG + 1), 0.5 * kmax * WG


def bose(x):
    return np.where(x > 700, 0.0, 1.0 / np.expm1(np.minimum(x, 700)))


def q_state(T, rc2, phi=True, it=500):
    D, P = 0.0, 0.0
    for _ in range(it):
        mu = math.sqrt(LX * (rc2 + 2 * D) + G2 * P); M2 = 2 * LX * (rc2 + D) + G2 * P
        Mp2 = G2 * (rc2 + D) + 3 * LP * P
        k, wk = _grid(T, math.sqrt(max(M2, Mp2)))
        w = np.sqrt(k * k + M2); fp = bose((w - mu) / T); fm = bose((w + mu) / T)
        Dn = np.sum(wk * k * k * (fp + fm) / w) / (2 * math.pi ** 2)
        if phi:
            W = np.sqrt(k * k + Mp2); Pn = np.sum(wk * k * k * bose(W / T) / W) / (2 * math.pi ** 2)
        else:
            Pn = 0.0
        if abs(Dn - D) < 1e-13 * (1 + D) and abs(Pn - P) < 1e-13 * (1 + P):
            D, P = Dn, Pn; break
        D, P = 0.5 * (D + Dn), 0.5 * (P + Pn)
    mu = math.sqrt(LX * (rc2 + 2 * D) + G2 * P); M2 = 2 * LX * (rc2 + D) + G2 * P; Mp2 = G2 * (rc2 + D) + 3 * LP * P
    k, wk = _grid(T, math.sqrt(max(M2, Mp2))); c = wk * k * k / (2 * math.pi ** 2)
    w = np.sqrt(k * k + M2); fp = bose((w - mu) / T); fm = bose((w + mu) / T)
    npl, nmi = np.sum(c * fp), np.sum(c * fm)
    d = D / 2
    pot = LX / 4 * (rc2 ** 2 + 8 * rc2 * d + 8 * d * d) + G2 / 2 * P * (rc2 + D) + 3 * LP / 4 * P * P
    E = 0.5 * mu * mu * rc2 + pot + np.sum(c * w * (fp + fm)) - 0.5 * M2 * D
    Ephi, Nphi = 0.0, 0.0
    if phi:
        W = np.sqrt(k * k + Mp2); Ephi = np.sum(c * W * bose(W / T)) - 0.5 * Mp2 * P; Nphi = np.sum(c * bose(W / T))
        E += Ephi
    return dict(T=T, rc2=rc2, mu=mu, M_X=math.sqrt(M2), M_phi=math.sqrt(Mp2) if phi else None, D=D, P=P,
                n_plus=npl, n_minus=nmi, n_phi=Nphi, charge=mu * rc2 + npl - nmi, energy=E, Ephi_frac=Ephi / E)


def q_solve(nQ, rho, phi=True):
    def rc2_for(T):
        f = lambda x: q_state(T, math.exp(x), phi)['charge'] / nQ - 1
        return math.exp(brentq(f, math.log(1e-2), math.log(1e7), xtol=1e-14))
    g = lambda lt: q_state(math.exp(lt), rc2_for(math.exp(lt)), phi)['energy'] / rho - 1
    Thi = 150.0
    while True:                                   # largest T at which a condensate still exists and E > rho
        try:
            if g(math.log(Thi)) > 0:
                break
        except ValueError:
            pass
        Thi *= 0.9
        if Thi < 3.5:
            raise RuntimeError('no bracket')
    lt = brentq(g, math.log(3.0), math.log(Thi), xtol=1e-13)
    st = q_state(math.exp(lt), rc2_for(math.exp(lt)), phi)
    st.update(Xbar_per_net=st['n_minus'] / nQ, r2_total=st['rc2'] + st['D'], phi_included=phi,
              condensate_charge_fraction=st['mu'] * st['rc2'] / nQ, nQ=nQ, rho=rho,
              phi_per_net=st['n_phi'] / nQ)
    return st


# ------------------------------------------------------------------ M1 homogeneous condensate (m1_vacuum conventions)
A_, B_, C_ = 0.1, 3.0, 0.2            # lX/g^2, lP/g^2, eps/g^2 in units g^2 v^4


def U_phi(S):
    return B_ / 4 * (S * S - 1) ** 2 - C_ * (S - 1) ** 3


def dU_phi(S):
    return B_ * S * (S * S - 1) - 3 * C_ * (S - 1) ** 2


def S_of_F2(F2, tphi=0.0, tlin=0.0, tilt=True):
    """phi equation at fixed F: dU_phi(S) + (F^2 + tphi) S - tlin = 0 (one root in (0, 1) for the tilted model;
    tphi = thermal phi mass^2/m_X^2, tlin = the linear thermal push 3 c g^2 <dphi^2>/m_X^2)."""
    c = C_ if tilt else 0.0
    f = lambda S: B_ * S * (S * S - 1) - 3 * c * (S - 1) ** 2 + (F2 + tphi) * S - tlin
    if not tilt and F2 + tphi >= B_ and tlin == 0:
        return 0.0
    lo, hi = (1e-14, 1.0) if tilt or tlin > 0 else (1e-9, 1.0)
    if f(lo) * f(hi) > 0:
        return 0.0
    return brentq(f, lo, hi, xtol=1e-15)


def hom_state(F, tphi=0.0, tX=0.0, tlin=0.0, tilt=True):
    """Homogeneous rotating condensate at field F (r = vF): S from the phi equation, Omega^2 = S^2 + aF^2 + tX,
    nbar = Omega F^2 (n = 100 nbar m_X^3), energy density e = nbar^2/2F^2 + U (units g^2 v^4 = 100 m_X^4)."""
    S = S_of_F2(F * F, tphi, tlin, tilt)
    Om2 = S * S + A_ * F * F + tX
    Om = math.sqrt(Om2)
    nb = Om * F * F
    U = (U_phi(S) if tilt else B_ / 4 * (S * S - 1) ** 2) + 0.5 * S * S * F * F + 0.25 * A_ * F ** 4 \
        + 0.5 * tphi * S * S + 0.5 * tX * F * F - tlin * (S - 1)
    e = nb * nb / (2 * F * F) + U
    return dict(F=F, S=S, Omega=Om, n=100 * nb, e=100 * e, P=100 * (nb * Om - e))


def hom_scan(tphi=0.0, tX=0.0, tlin=0.0, tilt=True, Fs=None):
    Fs = np.linspace(0.3, 4.0, 3701) if Fs is None else Fs
    rows = [hom_state(F, tphi, tX, tlin, tilt) for F in Fs]
    n = np.array([r['n'] for r in rows]); Om = np.array([r['Omega'] for r in rows]); P = np.array([r['P'] for r in rows])
    i = int(np.argmin(Om))
    # refine spinodal (min Omega along n) and zero pressure (bulk)
    fOm = lambda F: hom_state(F, tphi, tX, tlin, tilt)['Omega']
    r = minimize_scalar(fOm, bracket=(Fs[max(i - 2, 0)], Fs[i], Fs[min(i + 2, len(Fs) - 1)]), tol=1e-12)
    sp = hom_state(r.x, tphi, tX, tlin, tilt)
    bulk = None
    sg = np.where(np.sign(P[:-1]) != np.sign(P[1:]))[0]
    if len(sg):
        j = sg[-1]
        Fb = brentq(lambda F: hom_state(F, tphi, tX, tlin, tilt)['P'], Fs[j], Fs[j + 1], xtol=1e-14)
        bulk = hom_state(Fb, tphi, tX, tlin, tilt)
    return dict(spinodal=sp, bulk=bulk)


def hom_at_n(nphys, tphi=0.0, tX=0.0, tlin=0.0, tilt=True):
    """Homogeneous state at physical charge density n (m_X^3): solve for F (n(F) is monotone, m1_vacuum)."""
    f = lambda F: hom_state(F, tphi, tX, tlin, tilt)['n'] - nphys
    F = brentq(f, 1e-3, 20.0, xtol=1e-14)
    return hom_state(F, tphi, tX, tlin, tilt)


# ------------------------------------------------------------------ item 1: R values and clocks
def item1_sources():
    vs = SRC['vector_sources']; md = SRC['mediator_debris']['totals']
    R = dict(vec_e1e6_band=vs['e1e-6_band']['per_sign_totals']['R'], vec_e1e6_broad=vs['e1e-6_broad']['per_sign_totals']['R'],
             vec_e1e7_band=vs['e1e-7_band']['per_sign_totals']['R'], vec_e1e7_broad=vs['e1e-7_broad']['per_sign_totals']['R'],
             med_e1e6=md['e1e-6']['R_mediator_total'], med_e1e7=md['e1e-7']['R_mediator_total'],
             Y_L=md['e1e-6']['R_Y'], casc_e1e6=md['e1e-6']['R_cascade'], casc_e1e7=md['e1e-7']['R_cascade'])
    R['tot_e1e6_band'] = R['vec_e1e6_band'] + R['med_e1e6']
    R['tot_e1e6_broad'] = R['vec_e1e6_broad'] + R['med_e1e6']
    R['tot_e1e7_band'] = R['vec_e1e7_band'] + R['med_e1e7']
    R['tot_e1e7_broad'] = R['vec_e1e7_broad'] + R['med_e1e7']
    R['check_1_over_e_band'] = R['vec_e1e7_band'] / R['vec_e1e6_band']
    # clocks: Rayleigh decay law N_A(a) = N_A0 exp[-(a^2 - 1)/(2 a_unit^2)] (kin_source section 2)
    clocks = {}
    for k in ('e1e-6_band', 'e1e-7_band'):
        p = vs[k]['parent']
        au = p.get('a_unit') or p.get('a_unit_decay')
        clocks[k] = dict(parent=p)
    yl = SRC['mediator_debris']['YL_three_body']
    return R, clocks, yl


def decay_end(a_unit, frac=1e-6):
    """a at which the undecayed fraction is frac (Rayleigh law)."""
    return a_unit * math.sqrt(2 * math.log(1 / frac))


def hard_thermalisation(q, T, s_avail):
    """[E] number of 2->2 collisions of a hard X (conformal momentum q) on the relaxed plasma before the onset.
    Contact amplitudes only (X phi -> X phi |M| = 2 g^2 on 1 phi dof; X X -> X X |M| ~ 4 lX on 2 X dof), massless
    targets, sigma = |M|^2/(16 pi s), s = 2 q <E_t> with <E_t> = 2.70 T; the condensate is ignored (conservative)."""
    z3 = 1.2020569
    n1 = z3 / math.pi ** 2 * T ** 3
    s = 2 * q * 2.701 * T
    G = n1 * (2 * G2) ** 2 / (16 * math.pi * s) + 2 * n1 * (4 * LX) ** 2 / (16 * math.pi * s)
    return dict(Gamma_conf=G, collisions_before_onset=G * s_avail, log2_q_over_3T=math.log2(q / (3 * T)))


# ------------------------------------------------------------------ the chain per (state, R)
def chain_point(label, N, E0, R, phi=True):
    st = q_solve(N, E0 * (1 + R), phi)
    # onset: kappa_eff a^2 = g^2 (rc2 + D) + 3 lP P  (untilted: lP; M1 curvature at phi = 0: lP - 6 eps)
    Mphi2 = G2 * (st['rc2'] + st['D']) + 3 * LP * st['P']
    a_on = math.sqrt(Mphi2 / KAPPA)
    a_on_m1 = math.sqrt(Mphi2 / (KAPPA * (LP - 6 * EPST) / LP))
    out = dict(label=label, R=R, phi_in_equilibrium=phi, T_conf=st['T'], Xbar_per_net_massless=st['Xbar_per_net'],
               phi_per_net=st['phi_per_net'], condensate_charge_fraction=st['condensate_charge_fraction'],
               rc2_conf=st['rc2'], D_conf=st['D'], P_conf=st['P'], Ephi_frac=st['Ephi_frac'],
               energy_per_charge_conf=E0 * (1 + R) / N)
    for tag, a in (('onset_untilted', a_on), ('onset_M1_curvature', a_on_m1), ('a_ref', A_REF)):
        mc = MX * a / H                      # conformal m_X
        o = dict(a=a, a_over_aref=a / A_REF, T_h_over_mX=st['T'] / mc, T_h_GeV=st['T'] * H / a,
                 nQ_over_mX3=N / mc ** 3, n_condensate_over_mX3=st['mu'] * st['rc2'] / mc ** 3,
                 rc2_over_mX2=st['rc2'] / mc ** 2, r2_total_over_mX2=(st['rc2'] + st['D']) / mc ** 2,
                 dm_phi2_thermal_over_mX2=(G2 * st['D'] + 3 * LP * st['P']) / mc ** 2,
                 dm_X2_thermal_over_mX2=(2 * LX * st['D'] + G2 * st['P']) / mc ** 2,
                 phi_fluct_over_mX2=st['P'] / mc ** 2,
                 energy_per_charge_over_mX=E0 * (1 + R) / N / mc,
                 V0_untilted_per_charge_over_mX=75.0 / (N / mc ** 3), V0_M1_per_charge_over_mX=95.0 / (N / mc ** 3),
                 rho_X_over_mX4=E0 * (1 + R) / mc ** 4)
        out[tag] = o
    return out


def main():
    res = {}
    R, clocks, yl = item1_sources()
    res['R'] = R
    # timing
    vs = SRC['vector_sources']
    tim = {}
    for k in ('e1e-6_band', 'e1e-7_band'):
        p = vs[k]['parent']
        tim[k] = {kk: p[kk] for kk in p if isinstance(p[kk], (int, float))}
    tim['Y_L'] = {kk: yl[kk] for kk in yl if isinstance(yl[kk], (int, float))}
    res['clocks_raw'] = tim
    Rlist = [('R0', 0.0), ('med_e1e6', R['med_e1e6']), ('med_e1e7_1.14', R['med_e1e7']),
             ('tot_e1e6_band_4.97', R['tot_e1e6_band']), ('tot_e1e6_broad', R['tot_e1e6_broad']),
             ('tot_e1e7_band_41.7', R['tot_e1e7_band']), ('tot_e1e7_broad', R['tot_e1e7_broad'])]
    rows = []
    for state, N, E0 in (('heavy_release', NX_HEAVY, EX_HEAVY), ('pilot_state', NQ_PILOT, RHO_PILOT)):
        for tag, r in Rlist:
            for phi in (True, False):
                if state == 'pilot_state' and not phi and tag not in ('R0', 'tot_e1e6_band_4.97'):
                    continue
                c = chain_point(state + ':' + tag, N, E0, r, phi)
                rows.append(c)
                o = c['onset_untilted']; m = c['onset_M1_curvature']
                print('%-34s phi=%d T=%.2f Xb=%.4f cond=%.3f | onset a=%.3e T/m=%.2f nQ=%.1f ncond=%.1f r2=%.1f | M1 a=%.3e T/m=%.2f nQ=%.1f'
                      % (c['label'], phi, c['T_conf'], c['Xbar_per_net_massless'], c['condensate_charge_fraction'],
                         o['a'], o['T_h_over_mX'], o['nQ_over_mX3'], o['n_condensate_over_mX3'], o['r2_total_over_mX2'],
                         m['a'], m['T_h_over_mX'], m['nQ_over_mX3']), flush=True)
    res['chain'] = rows
    # hard daughters thermalisation estimate at e = 1e-6 (heavy release, R = 4.97, phi eq)
    c497 = [c for c in rows if c['label'] == 'heavy_release:tot_e1e6_band_4.97' and c['phi_in_equilibrium']][0]
    qmean = vs['e1e-6_band']['per_sign_totals']['mean_q_over_H']
    a_inj = vs['e1e-6_band']['parent']['a_mean_decay']
    s_avail = (c497['onset_untilted']['a'] - a_inj) / 1.9857962750634526   # conformal time a = a_m + C s
    res['hard_daughters_e1e6'] = dict(q_mean_over_H=qmean, T_conf=c497['T_conf'], a_injection=a_inj,
                                      s_available=s_avail, **hard_thermalisation(qmean, c497['T_conf'], s_avail))
    # timing: source completion vs onset
    au_v6 = vs['e1e-6_band']['parent']['a_unit']; au_v7 = vs['e1e-7_band']['parent']['a_unit']; au_y = 3153.0  # Y_L a_unit (kin_source NOTES section 4)
    res['timing'] = dict(vector_e1e6_a_end=decay_end(au_v6), vector_e1e7_a_end=decay_end(au_v7),
                         YL_a_end=decay_end(au_y), cascade_vectors='decay with the primordial vectors (same m_A)',
                         onset_relaxed_R497=c497['onset_untilted']['a'], onset_relaxed_R497_M1=c497['onset_M1_curvature']['a'],
                         a_ref=A_REF,
                         ratio_onset_over_vector_end=c497['onset_untilted']['a'] / decay_end(au_v6),
                         ratio_onset_over_YL_end=c497['onset_untilted']['a'] / decay_end(au_y))
    print('timing', res['timing'])
    print('hard', res['hard_daughters_e1e6'])
    # ------------------------------------------------ item 2: M1 homogeneous condensate, T = 0 reproduction
    z = hom_scan()
    z0 = hom_scan(tilt=False)
    res['hom_T0_M1'] = z; res['hom_T0_untilted'] = z0
    print('T=0 M1: spinodal n=%.3f F=%.4f Om=%.5f S=%.4f; bulk n=%.3f S=%.5f Om=%.6f' % (
        z['spinodal']['n'], z['spinodal']['F'], z['spinodal']['Omega'], z['spinodal']['S'], z['bulk']['n'], z['bulk']['S'], z['bulk']['Omega']))
    print('T=0 untilted: spinodal n=%.3f Om=%.5f; bulk n=%.3f Om=%.6f' % (z0['spinodal']['n'], z0['spinodal']['Omega'], z0['bulk']['n'], z0['bulk']['Omega']))
    n_on_D = math.sqrt(LX) * (LP * 100 / G2) ** 1.5           # sqrt(lX) (lP v^2/g^2)^{3/2}, v = 10 m_X
    n_on_D_M1 = math.sqrt(LX) * ((LP - 6 * EPST) * 100 / G2) ** 1.5
    res['n_onset_pure_condensate'] = dict(untilted=n_on_D, M1_curvature=n_on_D_M1,
                                          formula='n = sqrt(lambda_X) r^3, g^2 r^2 = lambda_phi v^2 (or lambda_phi - 6 eps)')
    print('n_onset [D] pure circular condensate: %.3f (untilted), %.3f (M1 curvature at phi=0) m_X^3' % (n_on_D, n_on_D_M1))
    # ------------------------------------------------ item 2: along the history (heavy release, each R), thermal masses
    hist = []
    for c in rows:
        if not c['phi_in_equilibrium'] or not c['label'].startswith('heavy'):
            continue
        o = c['onset_untilted']
        a_on = o['a']
        prof = []
        for af in np.r_[np.geomspace(0.5, 0.95, 10), np.linspace(1.0, 2.0, 11)]:
            a = af * a_on
            mc = MX * a / H
            ncond = c['onset_untilted']['n_condensate_over_mX3'] / af ** 3
            nQ = c['onset_untilted']['nQ_over_mX3'] / af ** 3
            tphi = c['onset_untilted']['dm_phi2_thermal_over_mX2'] / af ** 2
            tX = c['onset_untilted']['dm_X2_thermal_over_mX2'] / af ** 2
            tlin = 3 * C_ * G2 * c['onset_untilted']['phi_fluct_over_mX2'] / af ** 2
            hT = hom_at_n(ncond, tphi, tX, tlin, True)
            h0 = hom_at_n(ncond, 0, 0, 0, True)
            prof.append(dict(a_over_aon=float(af), T_h_over_mX=c['onset_untilted']['T_h_over_mX'] / af, n_cond=ncond, nQ=nQ,
                             S_eq_thermal=hT['S'], Omega_thermal=hT['Omega'], S_eq_T0=h0['S'], Omega_T0=h0['Omega'],
                             tphi=tphi, tX=tX))
        # thermal spinodal / bulk at the onset temperature (thermal masses frozen at their onset values) [E]
        zt = hom_scan(o['dm_phi2_thermal_over_mX2'], o['dm_X2_thermal_over_mX2'], 3 * C_ * G2 * o['phi_fluct_over_mX2'], True)
        # a at which the condensate density crosses 449.3 and 233.2 (T = 0 thresholds)
        cross = {}
        for name, nth in (('bulk_449', z['bulk']['n']), ('spinodal_233', z['spinodal']['n']),
                          ('thermal_spinodal', zt['spinodal']['n']), ('thermal_bulk', zt['bulk']['n'] if zt['bulk'] else None)):
            if nth is None:
                continue
            af = (o['n_condensate_over_mX3'] / nth) ** (1 / 3)
            cross[name] = dict(n=nth, a_over_aon=af, T_h_over_mX=o['T_h_over_mX'] / af,
                               S_T0=hom_at_n(nth)['S'])
        hist.append(dict(label=c['label'], R=c['R'], profile=prof, thermal_hom_at_onset_T=zt, crossings=cross))
        print('%-34s ncond(on)=%.1f  bulk crossing a/aon=%.3f (T/m=%.2f)  spinodal a/aon=%.3f (T/m=%.2f) thermal sp n=%.1f bulk n=%s'
              % (c['label'], o['n_condensate_over_mX3'], cross['bulk_449']['a_over_aon'], cross['bulk_449']['T_h_over_mX'],
                 cross['spinodal_233']['a_over_aon'], cross['spinodal_233']['T_h_over_mX'], zt['spinodal']['n'],
                 ('%.1f' % zt['bulk']['n']) if zt['bulk'] else None), flush=True)
    res['history_M1'] = hist
    # phi_eq on the T = 0 homogeneous branch for reference
    res['hom_T0_profile'] = [hom_at_n(n) for n in (600, 449.31, 400, 300, 233.2, 200, 164.3, 120, 76.4, 50, 24)]
    os.makedirs(OUT, exist_ok=True)
    json.dump(res, open(os.path.join(OUT, 'heat_onset.json'), 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
