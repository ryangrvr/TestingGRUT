"""fx_formation SKEPTIC (Claude, 1 Nov 2026): independent checks of the kinetics-side numbers.  All [C] arithmetic on
[E] models, m_X = 1 units unless marked.  Own code; imports only bessel_vapour.py (also own, skeptic/code).

1. Spinodal: own linearisation of the M1 condensate (fields psi = 10 F, phi = 10 S; V = 100 U(S, F)),
   dispersion (k^2+z+2aF^2)(k^2+z)(k^2+z+D) + 4 mu^2 z (k^2+z+D) - 4 F^2 S^2 (k^2+z) = 0, z = -omega^2,
   S(F) from dU/dS = 0 by own root finding; fit Gamma_peak/delta and k_peak/sqrt(delta).
2. Pilot evaporation window: net Hertz-Knudsen flux with the computed fugacity deficit, integrated over the window
   in ln a, with H from the visible clock (own g* table [R]); largest object radius that fully evaporates.
3. Low-heat capture freeze-out: own two-phase equilibrium path (vapour at mu_c, liquid entropy neglected) and
   relaxation ODE; 1 - eta_f, T_f, T_vis(T_f), free-streaming of the leftover X.
4. Mass per charge with the M1 relaxed p of fx_chain (gas-only 9.10e-4) instead of the Codex-reproduction 2.26e-3,
   and the H_i refit for each branch; e_crit with the e-dependent cascade.
"""
import json, math, sys, os
import numpy as np
from scipy.optimize import brentq, minimize_scalar
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bessel_vapour as bv

MX = 0.0924089637375362
YQ = 5.07941423e-9
MPL = 2.435e18
OMB, NB = bv.OMB, bv.NB
A3, B3, C3 = 8.142857534830613, 125.78683178383355, -2846.608236300704   # track fit3 (checked against family.json rows below)
a_, b_, c_ = 0.1, 3.0, 0.2


# ---------------------------------------------------------------- 1. spinodal
def U_S(S, F):
    return b_ * (S * S - 1) * S - 3 * c_ * (S - 1) ** 2 + S * F * F


def U_SS(S, F):
    return b_ * (3 * S * S - 1) - 6 * c_ * (S - 1) + F * F


def S_of_F(F):
    # unique root in (0, 1) (fx_chain [D]); bracket
    return brentq(lambda S: U_S(S, F), 1e-12, 1.0 - 1e-15, xtol=1e-15)


def hom(F):
    S = S_of_F(F)
    mu = math.sqrt(S * S + a_ * F * F)
    return S, mu, 100 * mu * F * F


def growth_max(F):
    S, mu, n = hom(F)
    D = U_SS(S, F)
    def gam(k):
        k2 = k * k
        # polynomial in z
        p = np.polymul(np.polymul([1, k2 + 2 * a_ * F * F], [1, k2]), [1, k2 + D])
        p = np.polyadd(p, 4 * mu * mu * np.polymul([1, 0], [1, k2 + D]))
        p = np.polysub(p, 4 * F * F * S * S * np.array([1, k2]))
        z = np.roots(p)
        zr = z[np.abs(z.imag) < 1e-10 * (1 + np.abs(z.real))].real
        pos = zr[zr > 0]                        # z = -omega^2 > 0  -> growth sqrt(z)
        return math.sqrt(pos.max()) if pos.size else 0.0
    ks = np.geomspace(1e-7, 1.0, 400)
    g = np.array([gam(k) for k in ks])
    i = int(np.argmax(g))
    if g[i] == 0:
        return 0.0, None, n
    r = minimize_scalar(lambda lk: -gam(math.exp(lk)), bounds=(math.log(ks[max(i - 2, 0)]), math.log(ks[min(i + 2, 399)])),
                        method='bounded', options=dict(xatol=1e-10))
    return -r.fun, math.exp(r.x), n


def spinodal():
    # spinodal = minimum of mu(F) along the homogeneous branch (d mu/d n = 0)
    r = minimize_scalar(lambda F: hom(F)[1], bounds=(1.5, 2.2), method='bounded', options=dict(xatol=1e-12))
    Fs = r.x
    S, mu, ns = hom(Fs)
    rows = []
    for d in (1e-6, 1e-5, 1e-4, 1e-3):
        F = brentq(lambda F: hom(F)[2] - ns * (1 - d), 0.5, Fs, xtol=1e-15)
        g, k, n = growth_max(F)
        rows.append(dict(delta=d, Gamma=g, k=k, c1=g / d, c2=k / math.sqrt(d)))
    # also the thermodynamic sound speed check: c_s^2 = (n/mu) dmu/dn at delta = 1e-3
    F = brentq(lambda F: hom(F)[2] - ns * (1 - 1e-3), 0.5, Fs, xtol=1e-15)
    h = 1e-6
    S1, m1, n1 = hom(F - h); S2, m2, n2 = hom(F + h)
    cs2 = (hom(F)[2] / hom(F)[1]) * (m2 - m1) / (n2 - n1)
    return dict(F_sp=Fs, S_sp=S, mu_sp=mu, n_sp=ns, rows=rows, cs2_at_delta_1em3=cs2)


# ---------------------------------------------------------------- visible clock [R]
_GT = np.array([[1e-5, 3.36], [5e-4, 3.36], [1e-3, 10.75], [0.05, 10.75], [0.1, 14.25], [0.15, 17.25], [0.17, 30.0],
                [0.2, 55.0], [0.3, 61.75], [1.0, 69.0], [1.5, 75.0], [3.0, 79.0], [5.0, 83.0], [10.0, 86.25], [50.0, 86.25],
                [80.0, 90.0], [150.0, 100.0], [300.0, 106.75], [1e6, 106.75]])
_GS = _GT.copy(); _GS[:2, 1] = 3.91


def gst(T, tab=_GT):
    return float(np.interp(math.log(T), np.log(tab[:, 0]), tab[:, 1]))


def Tvis_of(nQ):
    s = nQ * MX ** 3 / YQ
    f = lambda lt: math.log(2 * math.pi ** 2 / 45 * gst(math.exp(lt), _GS) * math.exp(3 * lt)) - math.log(s)
    return math.exp(brentq(f, math.log(1e-14), math.log(1e7)))


def H_of(Tv):
    return math.sqrt(math.pi ** 2 * gst(Tv) / 90) * Tv ** 2 / MPL * math.sqrt(1 + 0.8e-9 / Tv) / MX   # in m_X


def vbar(T):
    return math.sqrt(8 * T / math.pi) / math.sqrt(1 + 8 * T / math.pi)


# ---------------------------------------------------------------- 2. pilot window
def window_evaporation(Yh, label, stick=(1.0, 0.01)):
    Tev, Ts = bv.crossings(Yh)
    Tg = np.geomspace(Tev * 0.99999, Ts * 1.00001, 120)
    lna, rate = [], []
    for T in Tg:
        mu = bv.vapour_only_mu(Yh, T)
        v = bv.vapour(mu, T)
        mc = bv.mu_coex(T)
        vs = bv.vapour(mc, T)
        nQ = v['n']
        Tv = Tvis_of(nQ)
        Hh = H_of(Tv)
        deficit = 1 - math.exp((mu - mc) / T)
        lna.append(-math.log(nQ) / 3)
        rate.append(vbar(T) * vs['nX'] * deficit / (4 * NB) / Hh)    # dR/d ln a at stick = 1
    lna = np.array(lna); rate = np.array(rate)
    dR = float(np.trapezoid(rate, lna))
    out = dict(label=label, Y_h=Yh, T_ev=Tev, T_s=Ts, dlna_window=float(lna[-1] - lna[0]),
               T_vis_at_Tev_GeV=Tvis_of(bv.vapour(bv.vapour_only_mu(Yh, Tev * 0.99999), Tev * 0.99999)['n']))
    for s in stick:
        R = s * dR
        out['R_max_evaporated_stick%g' % s] = R
        out['Q_max_evaporated_stick%g' % s] = NB * 4 * math.pi / 3 * R ** 3
    return out


# ---------------------------------------------------------------- 3. low-heat freeze-out
def eq_two_phase(Yh, T):
    mc = bv.mu_coex(T) if T > 0.05 else OMB          # vapour-pressure shift < 1e-9 below T = 0.05
    v = bv.vapour(mc, T)
    phi = (Yh * v['s'] - v['n']) / (NB - v['n'] + Yh * v['s'])
    nQ = Yh * (1 - phi) * v['s']
    eta = phi * NB / nQ
    return eta, nQ, v


_PATH = {}


def eq_path(Yh, T_hi=0.6, T_lo=2.5e-3, N=700):
    if Yh not in _PATH:
        rows = []
        for T in np.geomspace(T_hi, T_lo, N):
            eta, nQ, v = eq_two_phase(Yh, T)
            try:
                Tv = Tvis_of(nQ)
            except ValueError:
                break
            rows.append((T, eta, nQ, -math.log(nQ) / 3, Tv, H_of(Tv)))
        _PATH[Yh] = rows
    return _PATH[Yh]


def freezeout(Yh, Q, stick):
    rows = eq_path(Yh)
    R = (3 * Q / (4 * math.pi * NB)) ** (1 / 3)
    x = 1 - rows[0][1]
    Tf = Tvf = None
    for (T0, e0, n0, l0, tv0, h0), (T1, e1, n1, l1, tv1, h1) in zip(rows[:-1], rows[1:]):
        G = (e1 * n1 / Q) * math.pi * R * R * vbar(T1) * stick / h1
        xeq = 1 - e1
        dl = l1 - l0
        x = xeq + (x - xeq) * math.exp(-G * dl) if G * dl < 700 else xeq
        if Tf is None and G < 1:
            Tf, Tvf = T1, tv1
    return dict(Q=Q, stick=stick, one_minus_eta=x, T_f=Tf, T_vis_f_GeV=Tvf,
                floor_formula=1 / (Yh * (2.5 + (1 - OMB) / Tf)) if Tf else None)


def lambda_fs(Th_f, Tvis_f):
    """Comoving free-streaming length (kpc) of NR X decoupling at T_h = Th_f (m_X units) when T_vis = Tvis_f (GeV):
    v = v_f a_f/a, integrated to today with radiation + matter (Omega_m h^2 = 0.143, Omega_r h^2 = 4.15e-5 x 1.68/1.68)."""
    vf = math.sqrt(3 * Th_f)                              # rms, NR
    T0 = 2.348e-13
    af = (3.91 / gst(Tvis_f, _GS)) ** (1 / 3) * T0 / Tvis_f
    H0 = 100 / 2.998e5 * 1.0    # h/Mpc in c=1 (Mpc^-1) per h
    om, orad = 0.143, 4.15e-5
    # comoving lambda = int v(a)/(a^2 H(a)) da  [Mpc]; H(a) = 100 km/s/Mpc sqrt(om a^-3 + orad a^-4) (h absorbed)
    lg = np.linspace(math.log(af), 0.0, 4000)
    a = np.exp(lg)
    v = np.minimum(vf * af / a, 1.0)
    Hh = (100 / 2.998e5) * np.sqrt(om / a ** 3 + orad / a ** 4)      # Mpc^-1
    integrand = v / (a * a * Hh) * a                       # da = a dln a
    return float(np.trapezoid(integrand, lg) * 1e3)        # kpc


def main(out):
    res = {}
    # family fit check against family.json rows (track file read only)
    fam = json.load(open(os.path.join(HERE, '..', '..', 'results', 'family.json')))
    chk = []
    for r in fam['rows']:
        if r['Q'] >= 3e5:
            fit = OMB + A3 * r['Q'] ** (-1 / 3) + B3 * r['Q'] ** (-2 / 3) + C3 / r['Q']
            dfit = OMB + 2 / 3 * A3 * r['Q'] ** (-1 / 3) + 1 / 3 * B3 * r['Q'] ** (-2 / 3)
            chk.append(dict(Q=r['Q'], EQ=r['E_over_Qm'], fit_resid=fit - r['E_over_Qm'], Omega=r.get('Omega'), dfit_minus_Omega=(dfit - r['Omega']) if r.get('Omega') else None))
    res['family_fit_check'] = dict(max_abs_resid=max(abs(c['fit_resid']) for c in chk),
                                   max_abs_first_law=max(abs(c['dfit_minus_Omega']) for c in chk if c['dfit_minus_Omega'] is not None),
                                   sigma_from_A=A3 / ((36 * math.pi) ** (1 / 3) * NB ** (-2 / 3)), n_rows=len(chk))
    print('family:', res['family_fit_check'], flush=True)
    SCR = os.environ.get('SKEP_SCR')
    part = os.path.join(SCR, 'partial.json') if SCR else None
    if part and os.path.exists(part):
        P = json.load(open(part)); res.update({k: P[k] for k in ('spinodal', 'slow_passage', 'window')})
        print('loaded spinodal/slow_passage/window from', part, flush=True)
    else:
        P = None
    if P is None:
        res['spinodal'] = sp = spinodal()
        print('spinodal:', json.dumps(sp), flush=True)
    sp = res['spinodal']
    # slow passage with own H at the spinodal crossing; n_Q there = n_c/0.917 (condensate fraction at onset, heavy e1e6)
    out_sp = []
    c1 = np.mean([r['c1'] for r in sp['rows'][:2]]); c2 = np.mean([r['c2'] for r in sp['rows'][:2]])
    for nQ in (sp['n_sp'], sp['n_sp'] / 0.9174):
        Tv = Tvis_of(nQ); Hh = H_of(Tv)
        for A in (10.0, 60.0):
            d = math.sqrt(6 * Hh * A / c1); k = c2 * math.sqrt(d)
            out_sp.append(dict(nQ=nQ, T_vis=Tv, H=Hh, A=A, delta=d, k=k, Gamma_over_H=c1 * d / Hh,
                               Q_cube=sp['n_sp'] * (2 * math.pi / k) ** 3, Q_half=sp['n_sp'] * (math.pi / k) ** 3))
    res['slow_passage'] = out_sp
    print('slow passage:', json.dumps(out_sp), flush=True)
    # 2. evaporation windows
    if P is None:
        res['window'] = [window_evaporation(0.27680265413087224, 'pilot_e1e6'), window_evaporation(0.2020822595873346, 'heavy_e1e6')]
        if part:
            json.dump({k: res[k] for k in ('spinodal', 'slow_passage', 'window')}, open(part, 'w'), default=float)
    print('window:', json.dumps(res['window']), flush=True)
    # 3. freeze-out
    fo = {}
    for lab, Yh in (('heavy_R0', 0.809575848826171), ('heavy_med_e1e6', 0.4839026445487382)):
        rows = []
        for Q in (8.5e12, 1e14, 1e17, 1e20):
            for s in (1.0, 0.01):
                r = freezeout(Yh, Q, s)
                if r['T_f']:
                    r['lambda_fs_kpc'] = lambda_fs(r['T_f'], r['T_vis_f_GeV'])
                    r['v_rms_f'] = math.sqrt(3 * r['T_f'])
                rows.append(r)
                print(lab, json.dumps(r), flush=True)
        fo[lab] = rows
    res['freezeout'] = fo
    # 4. mass per charge and H_i refits
    REF = 0.9373999905262296
    p_m1 = 9.10e-4
    M_ev = 1 + 2 * p_m1
    def refit(M):
        return (M / REF) ** (-2 / 3)
    res['mass'] = dict(evaporated=dict(p_used_track=0.0045855 / 2.028545, p_fx_chain_M1=p_m1, M_track=1 + 2 * 0.0045855 / 2.028545,
                                       M_corrected=M_ev, dM_corrected=M_ev / REF - 1,
                                       test_corrected=0.0667805 + 2.13356 * p_m1, Hi_refit_factor=refit(M_ev)),
                       lowheat=dict(M_range=[0.7822, 0.8017], Hi_refit_factor=[refit(0.8017), refit(0.7822)],
                                    eta1_only=dict(M=OMB, Hi_refit_factor=refit(OMB))),
                       note='H_i refit factor = (M/0.9374)^(-2/3) at fixed D (Y_Q ~ H_i^{3/2}); REVIEW_M1 refit rule')
    # e_crit with the e-dependent cascade: R(e) = Y_L + (R_vec(1e-6) + casc(1e-6)) (1e-6/e)
    YL, Rv, cs = 0.8055701153147813, 4.053801019129129, 0.10716449056156367
    ec = {}
    for lab, Rc in (('heavy', 2.757083481029306), ('pilot', 4.66760758942212)):
        ec[lab] = dict(R_crit=Rc, e_crit_track=1e-6 * Rv / (Rc - 0.9127346451110923), e_crit_cascade_scaled=1e-6 * (Rv + cs) / (Rc - YL))
    res['e_crit'] = ec
    json.dump(res, open(out, 'w'), indent=1, default=float)
    print(json.dumps(res['mass'], indent=1)); print(json.dumps(ec, indent=1))


if __name__ == '__main__':
    sys.argv = [sys.argv[0], sys.argv[1], os.path.join(HERE, '..', '..', 'results', 'eq_history.json')]
    main(sys.argv[1])
