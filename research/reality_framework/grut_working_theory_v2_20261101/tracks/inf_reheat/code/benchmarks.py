#!/usr/bin/env python3
"""inf_reheat items 3-4 (Claude, own code): reheating benchmarks paired with stand-in inflaton families, and the
consequence for M1r's calibration (H_i, r) and branch (S/Q, estimate).

Inputs: ../results/engine_runs.json (rprh.c: N_X, E_X, N_YL, N_YH versus m_I/H_i at slow decay), ../results/fluid.json
(K(b, Gamma) = a^2 H in radiation; D = K^{3/2}), ../results/channels.json (holder leakage Gamma_hid(b)).

Families (stand-ins, model choices; the inflaton-potential track supplies the real ones):
  F0  the engine's phenomenological oscillator: m_I = 3 H_i, released at rest with amplitude sqrt(6) M_P/3 = 0.816 M_P;
      no inflationary stage (r floor assumes H_CMB = H_i; flatness tested at phi_* = 0.816 M_P, the most lenient choice).
  T   alpha-attractor T-model V = V0 tanh^2(phi/sqrt(6 alpha) M_P), exact background to epsilon_H = 1 (the slice:
      H_i := H_end), N_* from the reheating history (Planck 2018 X eq. 47 form, w = 0 during reheating).
Channels (item 1 decides; see NOTES):
  S1  scalar Majorana Yukawa -(y/2) phi Nbar N (reflection phi -> -phi with N -> i gamma5 N, softly broken by M_N; seesaw):
      Gamma = N_f y^2 m beta^3/(16 pi); one-loop flatness y_nat = (0.008 pi^2 V_*/(N_f phi_*^2))^{1/4} (|Delta eta| <= 0.002
      with the quartic counterterm fixed at phi_*); kinematic blocking y Phi_1 > m/2 until a_unb = (2 y Phi_1/m)^{2/3}.
  S2  axial (pseudo-Goldstone) coupling -(M_N/2) Nbar exp(i gamma5 phi/f) N: |m_N| = M_N for every phi, so no one-loop
      potential and no blocking; Gamma = N_f (M_N/f)^2 m beta/(16 pi) (perturbative part only).
  S3  shift-symmetric gauge coupling (g/4) phi G Gtilde: Gamma = g^2 m^3/(8 pi); inflationary gauge production
      xi_A = g phidot_*/(2 H_*) must stay <~ 2.5 [Q, e.g. Barnaby-Peloso 1011.1500 / Planck 2018 X].
For each: Gamma/H_i fixed by the dark-matter calibration H_req(b, Gamma) = H_end (free-X calibration 1.291e11 on (3, 0.5),
Y ~ H_i^{3/2} N_X / D), the coupling it needs, T_RH, K, the matter era, the plasma time [E], leptogenesis [E], S/Q [E].
Usage: benchmarks.py -> ../results/benchmarks.json.   Labels: C on computed inputs; E where stated.
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from scipy.special import k1

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rh_common as rc

ENG = json.load(open(rc.RES / 'engine_runs.json'))
FLU = json.load(open(rc.RES / 'fluid.json'))
CHA = json.load(open(rc.RES / 'channels.json'))

# ------------------------------------------------------------------------------- engine interpolants (slow decay)
_bs = sorted({r['b']: r for r in ENG['bscan'] + ENG.get('bench_points', [])}.values(), key=lambda r: r['b'])   # slow-decay engine runs
B_ARR = np.array([r['b'] for r in _bs])


def _interp(name, b):
    v = np.array([abs(r[name]) for r in _bs])
    return float(np.interp(b, B_ARR, v)), bool(b < B_ARR[0] or b > B_ARR[-1])


def engine_at(b):
    out = {k: _interp(k, b)[0] for k in ('NX', 'EX', 'NYL', 'NYH', 'eps')}
    out['extrapolated'] = _interp('NX', b)[1]
    return out


_fits = sorted((float(k), v['c_K']) for k, v in FLU['fits'].items())
_fb = [k for k, _ in _fits]
_fc = [c for _, c in _fits]


def cK(b):
    return float(np.exp(np.interp(math.log(b), np.log(_fb), np.log(_fc))))


def K_of(b, G):
    """K(b, Gamma): fluid grid (exact) interpolated in (ln b, ln Gamma); power law below Gamma = 1e-6."""
    if G <= 1e-6:
        return cK(b) * G ** (-1 / 3)
    rows = [g for g in FLU['grid'] if g['K'] is not None]
    bs = sorted({g['b'] for g in rows})
    lk = []
    for bb in bs:
        gg = sorted([g for g in rows if g['b'] == bb], key=lambda g: g['Gamma'])
        lk.append(np.interp(math.log(G), [math.log(g['Gamma']) for g in gg], [math.log(g['K']) for g in gg]))
    return float(np.exp(np.interp(math.log(b), np.log(bs), lk)))


def NMD_of(b, G):
    """e-folds from the slice to inflaton-radiation equality (fluid grid; power law (2/3) ln(1/Gamma) + c below 1e-6)."""
    rows = [g for g in FLU['grid'] if g['b'] == min(_fb, key=lambda x: abs(x - b)) and g.get('N_eq')]
    gg = sorted(rows, key=lambda g: g['Gamma'])
    c = np.mean([g['N_eq'] - (2 / 3) * math.log(1 / g['Gamma']) for g in gg if g['Gamma'] <= 1e-6])
    if G <= 1e-6:
        return (2 / 3) * math.log(1 / G) + c
    return float(np.interp(math.log(G), [math.log(g['Gamma']) for g in gg], [g['N_eq'] for g in gg]))


def H_req(b, G, H_cal=rc.H_CAL_FREEX):
    """H_i needed for the dark matter: Y ~ H_i^{3/2} N_X / D at fixed dimensionless history, D = K^{3/2}."""
    e = engine_at(b)
    D = K_of(b, G) ** 1.5
    return H_cal * ((rc.NX_REF_3_05 / e['NX']) * (D / rc.K_REF ** 1.5)) ** (2 / 3)


# ------------------------------------------------------------------------------- T-model background
def tmodel(alpha, Nstar):
    w = math.sqrt(6 * alpha)
    V = lambda p: math.tanh(p / w) ** 2
    dV = lambda p: 2 * math.tanh(p / w) / (w * math.cosh(p / w) ** 2)
    d2V = lambda p: (2 / w ** 2) * (1 / math.cosh(p / w) ** 4 - 2 * math.tanh(p / w) ** 2 / math.cosh(p / w) ** 2)

    def rhs(N, y):
        p, pi = y
        H2 = V(p) / (3 - pi * pi / 2)
        return [pi, -(3 - pi * pi / 2) * pi - dV(p) / H2]

    def end_and_N(p0, turn=False):
        ev = lambda N, y: y[1] ** 2 / 2 - 1
        ev.terminal = True
        s = solve_ivp(rhs, (0, 400), [p0, -dV(p0) / V(p0)], events=ev, rtol=1e-10, atol=1e-13)
        Ne = s.t_events[0][0]; pe, pie = s.y_events[0][0]
        out = dict(Ne=Ne, pe=pe, pie=pie)
        if turn:
            # cosmic time (units 1/sqrt(V0)); N-time is singular where V = 0
            He2 = V(pe) / (3 - pie ** 2 / 2)
            def rhs_t(t, y):
                p, pd = y
                Hh = math.sqrt((0.5 * pd * pd + V(p)) / 3)
                return [pd, -3 * Hh * pd - dV(p)]
            ev2 = lambda t, y: y[1]
            ev2.terminal = True; ev2.direction = 1
            s2 = solve_ivp(rhs_t, (0, 20 * math.pi * w), [pe, pie * math.sqrt(He2)], events=ev2, rtol=1e-10, atol=1e-13,
                           max_step=0.02 * w)          # m = sqrt(2)/w in these units: ~14 periods
            out['Phi1'] = abs(s2.y_events[0][0][0]); out['t_turn_Hend'] = s2.t_events[0][0] * math.sqrt(He2)
        return out
    xe = 0.5 * math.asinh(math.sqrt(8) / w)
    xs = 0.5 * math.acosh(math.cosh(2 * xe) + 8 * Nstar / w ** 2)
    ps = brentq(lambda p: end_and_N(p)['Ne'] - Nstar, w * xs * 0.9, w * xs * 1.1)
    E = end_and_N(ps, turn=True)
    e = 0.5 * (dV(ps) / V(ps)) ** 2; eta = d2V(ps) / V(ps)
    V0 = rc.AS * 24 * math.pi ** 2 * e / V(ps)
    Hs = math.sqrt(V0 * V(ps) / 3); He = math.sqrt(V0 * V(E['pe']) / (3 - E['pie'] ** 2 / 2))
    m = math.sqrt(2 * V0 / w ** 2)
    return dict(alpha=alpha, w=w, N_star=Nstar, phi_star=ps, phi_end=E['pe'], Phi1=E['Phi1'], eps_star=e, eta_star=eta,
                ns=1 - 6 * e + 2 * eta, r=16 * e, V_star=V0 * V(ps), H_star=Hs * rc.MP, H_end=He * rc.MP, m=m * rc.MP, b=m / He,
                rho_end=3 * He ** 2)


def Nstar_formula(V_star, rho_end, rho_th):
    """Planck 2018 X eq. (47) form with k_* = 0.05/Mpc, w_int = 0 (M_P = 1 units)."""
    return 66.89 - math.log(0.05 / 2.247e-4) + 0.25 * math.log(V_star ** 2 / rho_end) + (1 / 12) * math.log(rho_th / rho_end) \
        - (1 / 12) * math.log(rc.GSTAR / 106.75)


# ------------------------------------------------------------------------------- channel numbers
def y_nat(V_star, phi_star, Nf=1):
    return (0.008 * math.pi ** 2 * V_star / (Nf * phi_star ** 2)) ** 0.25


def lam_nat(V_star, phi_star):
    return (0.004 * math.pi ** 2 * V_star / phi_star ** 2) ** 0.5


def washout(z_rh, Kw=rc.MNU / 1.08e-12):
    f = lambda z: 0.25 * Kw * z ** 3 * k1(z)
    return quad(f, z_rh, 200.0, limit=200)[0]


def plasma_solve(b, G, MPh, chan, tau_unb, GN_lab, E, alpha_g=1 / 40.0, need_m2=1e-3, T_need=0.25, tmax=1e7, ip=None):
    """[E] Early SM plasma in H_i units (see the module docstring for the model).  Densities in units of rho_I(0) = 3 MPh^2;
    matter era H = 1/(1 + 1.5 tau); the N source switches on at tau_unb (S1 blocking)."""
    rho0 = 3 * MPh ** 2

    # exponential integrator on a log grid (the system is linear; 3000 steps per decade)
    t = np.exp(np.linspace(math.log(1e-3), math.log(tmax), int(3000 * math.log10(tmax / 1e-3))))
    rN = np.zeros_like(t); rSn = np.zeros_like(t); rB = np.zeros_like(t); rSB = np.zeros_like(t)
    EB = None
    for i in range(1, len(t)):
        dt = t[i] - t[i - 1]; tm = 0.5 * (t[i] + t[i - 1]); H = 1 / (1 + 1.5 * tm)
        src = G * H * H if tm >= tau_unb else 0.0
        if ip is not None:
            # instant preheating (Felder-Kofman-Linde hep-ph/9812289) on Greene-Kofman production, blocked regime only:
            # per crossing f1 = y^{5/2} (Phi/m)^{1/2} / (2 pi^3) of rho_I goes into N (first crossing in full, later ones
            # limited by the phase space freed by N decays, Pauli), N of mass ~ (2/pi) y Phi decay to L H at
            # y_nu^2 m_N / (8 pi) [E].  Phi ~ Phi1 a^{-3/2} = Phi1 H (matter era).
            Ph = ip['Phi1H'] * H
            if ip['y'] * Ph > b / 2:
                f1 = ip['y'] ** 2.5 * (Ph / b) ** 0.5 / (2 * math.pi ** 3)
                GB = ip['ynu2'] * (2 / math.pi) * ip['y'] * Ph / (8 * math.pi)
                first = (t[i - 1] < math.pi / (2 * b) <= t[i])
                inj = f1 * H * H * (1.0 if first else 0.0) + f1 * H * H * min(1.0, GB * math.pi / b) * (b / math.pi) * dt
                eB = math.exp(-(3 * H + GB) * dt)
                rB[i] = rB[i - 1] * eB + inj
                rSB[i] = rSB[i - 1] * math.exp(-4 * H * dt) + GB * rB[i] * dt
                if EB is None:
                    EB = ip['y'] * Ph / math.pi
            else:
                rB[i] = rB[i - 1] * math.exp(-3 * H * dt); rSB[i] = rSB[i - 1] * math.exp(-4 * H * dt)
        if chan == 'S3':
            lam = 4 * H; e = math.exp(-lam * dt)
            rSn[i] = rSn[i - 1] * e + src * (1 - e) / lam
            continue
        lam = 4 * H + GN_lab; e = math.exp(-lam * dt)
        rN[i] = rN[i - 1] * e + src * (1 - e) / lam
        e4 = math.exp(-4 * H * dt)
        rSn[i] = rSn[i - 1] * e4 + GN_lab * 0.5 * (rN[i] + rN[i - 1]) * (1 - e4) / (4 * H)
    rS = rSn * rho0
    m2 = 0.018 * rS / E ** 2
    if ip is not None and EB is not None:
        m2 = m2 + 0.018 * rSB * rho0 / EB ** 2
        rS = rS + rSB * rho0
    T = (30 * np.maximum(rS, 0) / (math.pi ** 2 * rc.GSTAR)) ** 0.25
    tth = np.where(T > 0, alpha_g ** -2 * np.sqrt(E / np.maximum(T, 1e-300)) / np.maximum(T, 1e-300), np.inf)
    i1 = np.where(m2 >= need_m2)[0]; i2 = np.where((T >= T_need) & (t >= tth))[0]
    return dict(tau_induced=float(t[i1[0]]) if len(i1) else None, tau_thermal=float(t[i2[0]]) if len(i2) else None,
                m2_at_8=float(np.interp(8.0, t, m2)), T_at_8=float(np.interp(8.0, t, T)), tth_at_8=float(np.interp(8.0, t, tth)),
                rhoSM_over_rhoI_at_8=float(np.interp(8.0, t, rS)) / (rho0 / (1 + 12.0) ** 2))


def SQ_estimate(b, G, Hi, e=1e-6, R_inf=0.0):
    """[E] S/Q of the relaxed plasma for a slow-decay history, by rescaling the decisive-pass budget (portals at M_Y):
    R_vec(1e-6) = 3.608 at K = 1.9858, E_X = 188706 (dq_heat), carrying N_A ~ a_*^3 m_A (a_* at H = m_A, radiation:
    (K/m_A)^{1/2}; matter era: (x0/(3 m_A^2))^{1/3}) and the deposit ~ K^{1/2}; R_Y = 0.3786 and R_casc+Sigma = 0.0184
    (x (1e-6/e)^{1/2}) scaled by the engine's N_YL+N_YH and K^{1/2}; S/Q = closed form of dq_entropy + 0.0055 (increments),
    which returns 4.955 at the reference."""
    eg = engine_at(b)
    K = K_of(b, G)
    mA = 2.4e-4 * (e / 1e-6)
    Kref, Eref = rc.K_REF, 188706.29
    NA_ref = (Kref / mA) ** 1.5 * mA
    if G > mA:
        NA = (K / mA) ** 1.5 * mA
    else:
        kb = sorted(X0)
        x0 = float(np.exp(np.interp(math.log(b), np.log(kb), np.log([X0[k] for k in kb]))))   # rho_I a^3 in the matter era
        NA = (x0 / (3 * mA * mA)) * mA
    Rvec_ref_e = 3.608 * (1e-6 / e) ** 1.007
    Rvec = Rvec_ref_e * (NA / NA_ref) * math.sqrt(K / Kref) * (Eref / eg['EX'])
    ny = (eg['NYL'] + eg['NYH']) / (0.56664 + 0.010304)
    RY = 0.3786 * ny * math.sqrt(K / Kref) * (rc.EX_REF_3_05 / eg['EX'])
    Rc = 0.0184 * (eg['NYH'] / 0.010304) * math.sqrt(K / Kref) * (rc.EX_REF_3_05 / eg['EX']) * (1e-6 / e) ** 0.5
    R = Rvec + RY + Rc + R_inf
    EEc = eg['eps'] ** (-4 / 3)
    W = EEc * (1 + R)
    SQ = 0.190255 * (W - 0.97494) ** 0.75013 - 0.00428 + SQ_OFF
    return dict(K=K, R_vec=Rvec, R_Y=RY, R_casc=Rc, R_inf=R_inf, R_total=R, E_over_Ec=EEc, W=W, SQ=SQ,
                vector_crossing='radiation' if G > mA else 'matter era')


X0 = {}
SQ_OFF = 0.0


def x0_of_b(b):
    """comoving inflaton energy rho_I a^3 deep in the matter era (Gamma = 0), H_i units: from fluid.stage1."""
    import fluid
    s = fluid.stage1(b, 0.0, 0.0, 80.0)
    la, I, P, rr = s.y[:, -1]
    return (0.5 * P * P + 0.5 * b * b * I * I) * math.exp(3 * la)


def e_crit(b, G, Hi, R_inf=0.0):
    f = lambda le: SQ_estimate(b, G, Hi, math.exp(le), R_inf)['SQ'] - 3.4716
    try:
        return math.exp(brentq(f, math.log(1e-9), math.log(1e-2)))
    except ValueError:
        return None


# ------------------------------------------------------------------------------- one benchmark
_FAM = {}


def benchmark(alpha, Nf=1, chan='S1', MN_frac=0.25):
    out = dict(family='T(alpha=%g)' % alpha, channel=chan, Nf=Nf)
    if alpha in _FAM:
        fam, G, Nst, why = _FAM[alpha]
        if G is None:
            out.update(fam=fam, solution=None, reason=why)
            return out
    else:
        fam, G, Nst, why = solve_family(alpha)
        _FAM[alpha] = (fam, G, Nst, why)
        if G is None:
            out.update(fam=fam, solution=None, reason=why)
            return out
    return finish(out, fam, G, Nst, Nf, chan, MN_frac)


def solve_family(alpha):
    Nst = 55.0
    for it in range(8):
        fam = tmodel(alpha, Nst)
        b, He = fam['b'], fam['H_end']
        lo, hi = math.log(1e-16), math.log(min(0.5, 0.15 * b))
        f = lambda lg: math.log(H_req(b, math.exp(lg)) / He)
        if f(hi) > 0:
            return fam, None, Nst, 'H_end below the fastest perturbative calibration'
        if f(lo) < 0:
            return fam, None, Nst, 'needs Gamma < 1e-16 H_i'
        G = math.exp(brentq(f, lo, hi))
        rho_th = 3 * (G * He / rc.MP) ** 2
        Nn = Nstar_formula(fam['V_star'], fam['rho_end'], rho_th)
        if abs(Nn - Nst) < 2e-3:
            break
        Nst = Nn
    return fam, G, Nst, None


def finish(out, fam, G, Nst, Nf, chan, MN_frac):
    b, He = fam['b'], fam['H_end']
    MPh = rc.MP / He
    Phi1H = fam['Phi1'] * MPh
    MN = MN_frac * b                                       # H_i units
    beta = math.sqrt(1 - 4 * MN_frac ** 2)
    Gg = G * He
    o = dict(Gamma_over_Hi=G, Gamma_over_m=G / b, Gamma_GeV=Gg, T_RH_GeV=rc.T_RH(Gg), K=K_of(b, G), D=K_of(b, G) ** 1.5,
             N_matter_era=NMD_of(b, G), H_i=He, r=fam['r'], ns=fam['ns'], N_star=Nst, b=b, m_I_GeV=fam['m'], H_star=fam['H_star'],
             Phi1_over_MP=fam['Phi1'], phi_star=fam['phi_star'])
    if chan == 'S1':
        y = math.sqrt(16 * math.pi * G / (Nf * b * beta ** 3))
        o.update(y=y, y_nat=y_nat(fam['V_star'], fam['phi_star'], Nf), y_over_ynat=None, yPhi1_over_half_m=y * Phi1H / (b / 2))
        o['y_over_ynat'] = y / o['y_nat']; o['flatness_tuning'] = o['y_over_ynat'] ** 4
        o['a_unblock'] = max(1.0, o['yPhi1_over_half_m'] ** (2 / 3)); o['tau_unblock'] = (o['a_unblock'] ** 1.5 - 1) / 1.5
        n1 = (y * Phi1H * b) ** 1.5 / (4 * math.pi ** 3)
        o['GK_energy_fraction_first_crossing'] = n1 * y * Phi1H / (0.5 * b * b * Phi1H ** 2)
    elif chan == 'S2':
        yp = math.sqrt(16 * math.pi * G / (Nf * b * beta))
        o.update(y_axial=yp, f_over_Hi=MN / yp, f_GeV=MN / yp * He, Phi1_over_f=Phi1H / (MN / yp))
    elif chan == 'S3':
        g = math.sqrt(8 * math.pi * G / b ** 3)            # 1/H_i
        xiA = g * math.sqrt(2 * fam['eps_star']) * MPh / 2      # xi_A = g phidot/(2H), phidot = sqrt(2 eps) H M_P
        o.update(g_times_Hi=g, inv_g_GeV=He / g, xi_A=xiA, xi_A_ok=xiA < 2.5)
    # N decay and leptogenesis (S1, S2)
    if chan in ('S1', 'S2'):
        MNg = MN * He
        GN = rc.MNU * MNg ** 2 / (8 * math.pi * rc.VEW ** 2) / He          # H_i units, at rest
        gam = (b / 2) / MN
        o.update(M_N_GeV=MNg, Gamma_N_over_Hi=GN, tau_N_lab=gam / GN, Gamma_N_over_Gamma_I=GN / G)
        TR = o['T_RH_GeV']; z = MNg / TR
        eps_need = rc.NB_OBS / (0.53 * TR / fam['m'] * 1.0)
        eps_DI = 3 / (16 * math.pi) * MNg * rc.MNU / rc.VEW ** 2
        o.update(z_RH=z, washout_integral=washout(z), eps_needed_nonthermal=eps_need, eps_DI=eps_DI,
                 nonthermal_max_over_observed=eps_DI / eps_need, nB_s_nonthermal_at_eps_DI=0.53 * eps_DI * TR / fam['m'])
        # T_RH ~ M_N: the N reach equilibrium; strong-washout thermal leptogenesis (Buchmuller-Di Bari-Plumacher
        # hep-ph/0401240): K = mt/m*, m* = 1.08e-3 eV, kappa_f = 2/(K z_B), z_B = 2 + 4 K^0.13 exp(-2.5/K);
        # n_B/s = (28/79) * (n_N^eq/s at z << 1 = 135 zeta(3)/(4 pi^4 g*)) * eps * kappa_f   [E]
        Kw = rc.MNU / 1.08e-12; zB = 2 + 4 * Kw ** 0.13 * math.exp(-2.5 / Kw); kap = 2 / (Kw * zB)
        yN = 135 * 1.2020569 / (4 * math.pi ** 4 * rc.GSTAR)
        o.update(thermal_lepto=dict(K=Kw, kappa_f=kap, nB_s_per_eps=(28 / 79) * yN * kap, eps_needed=rc.NB_OBS / ((28 / 79) * yN * kap),
                                    eps_needed_over_DI=rc.NB_OBS / ((28 / 79) * yN * kap) / eps_DI))
        if chan == 'S2':
            # [E] axial (derivative) production: Landau-Zener at k_eff = h k + mu5 = 0 with gap M_N, mu5 = phidot/(2f);
            # two crossings per half period give n_k ~ 2 pi M_N^2 / |mu5dot| (helicity suppressed) for k < Phi m/(2f)
            fH = MN / o['y_axial']
            mu5max = Phi1H * b / (2 * fH); mu5dot = b * b * Phi1H / (2 * fH)
            nk = min(1.0, 2 * math.pi * MN ** 2 / mu5dot)
            frac = nk * mu5max ** 4 / (4 * math.pi ** 2) / (0.5 * b * b * Phi1H ** 2)
            o['axial_nonperturbative'] = dict(mu5_max=mu5max, n_k=nk, energy_fraction_per_half_period=frac,
                                              Gamma_eff_initial=frac * b / math.pi, ratio_to_perturbative=frac * b / math.pi / G)
        tau_unb = o.get('tau_unblock', 0.0)
        ip = dict(y=o['y'], Phi1H=Phi1H, ynu2=rc.MNU * MNg / rc.VEW ** 2) if chan == 'S1' else None
        o['plasma'] = plasma_solve(b, G, MPh, chan, tau_unb, GN / gam, b / 4, ip=ip)
        if chan == 'S1':
            o['plasma_no_instant_preheating'] = plasma_solve(b, G, MPh, chan, tau_unb, GN / gam, b / 4)
    else:
        o['plasma'] = plasma_solve(b, G, MPh, 'S3', 0.0, 0.0, b / 2, alpha_g=1 / 25.0)
    # Higgs portal window
    mM = fam['m'] / rc.MP
    o['lambda_IH'] = dict(lam_min=(10 / 3) * mM ** 2, lam_max_q1=4 * mM ** 2 / fam['Phi1'] ** 2, lam_max_q350=4 * 350 * mM ** 2 / fam['Phi1'] ** 2,
                          lam_nat=lam_nat(fam['V_star'], fam['phi_star']))
    # holder leakage to the hidden sector at this b
    hr = CHA['hidden']['rows']
    Gh = float(np.exp(np.interp(b, [r['b'] for r in hr], [math.log(max(r['Gamma_hidden'], 1e-300)) for r in hr])))
    eg = engine_at(b)
    Rinf = 3 * MPh ** 2 * o['K'] ** 2 * (Gh / G) / eg['EX']
    o['hidden_leak'] = dict(Gamma_hid=Gh, BR=Gh / G, R_inf=Rinf)
    o['engine_at_b'] = eg
    o['SQ_e1e-6'] = SQ_estimate(b, G, He, 1e-6, Rinf)
    o['e_crit_bulk'] = e_crit(b, G, He, Rinf)
    out.update(fam=fam, solution=o)
    return out


def F0_rows(ys=(1e-4, 2.6e-4, 1e-3, 1e-2, 0.1), Nf=1):
    """engine oscillator (3, Gamma(y)): the H_i the dark matter needs and the resulting r floor (no potential closure)."""
    rows = []
    b = 3.0; beta3 = (1 - 4 * 0.25 ** 2) ** 1.5
    for y in ys:
        G = Nf * y * y * b * beta3 / (16 * math.pi)
        Hq = H_req(b, G)
        MPh = rc.MP / Hq
        rows.append(dict(y=y, Gamma_over_Hi=G, Gamma_over_m=G / b, K=K_of(b, G), D=K_of(b, G) ** 1.5, N_matter_era=NMD_of(b, G),
                         H_req=Hq, r_floor=rc.r_of_H(Hq), T_RH_GeV=rc.T_RH(G * Hq), yPhi_over_half_m=y * math.sqrt(6) / 3 * MPh / 1.5,
                         y_nat_at_0p816MP=y_nat(3 * (Hq / rc.MP) ** 2, math.sqrt(6) / 3, Nf), SQ=SQ_estimate(b, G, Hq)['SQ']))
    # Gamma needed for r = 0.036 at H_CMB = H_i
    Hmax = rc.H_of_r(rc.R_BOUND)
    Gmin = math.exp(brentq(lambda lg: math.log(H_req(b, math.exp(lg)) / Hmax), math.log(1e-16), math.log(0.5)))
    return dict(rows=rows, H_max_r036=Hmax, Gamma_min_for_r036=Gmin, y_min_for_r036=math.sqrt(16 * math.pi * Gmin / (Nf * b * beta3)))


def main():
    global SQ_OFF
    for b in (1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.153, 11.427, 33.185, 101.687):
        X0[b] = x0_of_b(b)
    # calibrate the additive constant to the decisive-pass 4.955 at (3, 0.5): budget R = 4.006 at E/E_c = 15.6036
    W = 15.6036 * (1 + 4.006)
    SQ_OFF = 4.955 - (0.190255 * (W - 0.97494) ** 0.75013 - 0.00428)
    res = dict(method=__doc__, SQ_offset=SQ_OFF, x0=X0)
    # reference check of the S/Q rescaling at (3, 0.5) with the engine's own (3, 0.5) numbers
    res['F0'] = F0_rows()
    print('F0: Gamma_min(r<0.036) = %.3g H_i, y_min = %.3g' % (res['F0']['Gamma_min_for_r036'], res['F0']['y_min_for_r036']))
    for r_ in res['F0']['rows']:
        print('  y=%.2g G/m=%.2e K=%.1f H_req=%.3e r=%.2e T_RH=%.2e yPhi/(m/2)=%.2e y_nat=%.2e SQ=%.2f' % (
            r_['y'], r_['Gamma_over_m'], r_['K'], r_['H_req'], r_['r_floor'], r_['T_RH_GeV'], r_['yPhi_over_half_m'], r_['y_nat_at_0p816MP'], r_['SQ']))
    bms = []
    for alpha in (1e-4, 1e-3, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0):
        for chan in ('S1', 'S2', 'S3'):
            bm = benchmark(alpha, 1, chan)
            bms.append(bm)
            s = bm['solution']
            if s is None:
                print('%-14s %s  no solution: %s (b=%.3g H_end=%.3g)' % (bm['family'], chan, bm['reason'], bm['fam']['b'], bm['fam']['H_end']), flush=True)
                continue
            extra = ''
            if chan == 'S1':
                extra = 'y=%.3g y/ynat=%.3g yPhi/(m/2)=%.3g tau_unb=%.3g' % (s['y'], s['y_over_ynat'], s['yPhi1_over_half_m'], s['tau_unblock'])
            elif chan == 'S2':
                extra = 'y_ax=%.3g f=%.3g GeV Phi1/f=%.3g' % (s['y_axial'], s['f_GeV'], s['Phi1_over_f'])
            else:
                extra = 'g=%.3g/H_i xi_A=%.3g' % (s['g_times_Hi'], s['xi_A'])
            print('%-14s %s b=%.3f H_i=%.3e r=%.2e ns=%.4f N*=%.1f G/Hi=%.3e T_RH=%.2e K=%.1f NMD=%.1f SQ=%.2f ecrit=%s leak R=%.2g pl=%s/%s %s' % (
                bm['family'], chan, s['b'], s['H_i'], s['r'], s['ns'], s['N_star'], s['Gamma_over_Hi'], s['T_RH_GeV'], s['K'], s['N_matter_era'],
                s['SQ_e1e-6']['SQ'], ('%.2e' % s['e_crit_bulk']) if s['e_crit_bulk'] else None, s['hidden_leak']['R_inf'],
                s['plasma']['tau_induced'], s['plasma']['tau_thermal'], extra), flush=True)
    res['benchmarks'] = bms
    json.dump(res, open(rc.RES / 'benchmarks.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
