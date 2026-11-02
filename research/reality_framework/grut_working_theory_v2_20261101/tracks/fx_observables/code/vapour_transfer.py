"""fx_observables item 2: the free-X gas coupled to the objects, and the linear matter power it leaves.  [C]/[E]

Own code (Claude, 1 Nov 2026).  Two gas pictures:

 (A) Codex's gas (f8186fa/94f660d): 6.595 % of the charge free, v_rms = 0.0684 c at T_vis = 2 MeV, cold objects that
     neither absorb nor emit (sticking 0); momentum coupling to Q = 2e5 cores through the elastic sigma = 332/m_X^2
     (tracks/m1_transport SKEPTIC.md); T_gas ~ a^-2.
 (B) Equilibrium vapour (tracks/fx_formation: exchange rates >> H): objects of charge Q coexist with an X vapour at
     mu = Omega(Q) m_X (Gibbs-Thomson); with the hidden entropy per charge S/Q conserved and carried by the vapour,
     n_v (5/2 + Delta/T_h) = (S/Q) n_Q,   n_v = (m T_h / 2 pi)^{3/2} e^{-Delta/T_h},  Delta = (1 - Omega(Q)) m_X,
     so T_h falls only logarithmically (latent heat) and a vapour of (S/Q)/(5/2 + Delta/T_h) ~ 1-5 % of the charge at
     v ~ 0.1 c survives until vapour capture freezes out: Gamma_c = n_obj s pi R^2 vbar = H (relaxation ODE as in
     fx_formation).  After that T_v ~ a^-2; momentum decoupling when n_obj sigma_m vbar = H (sigma_m = pi R^2, or
     332/m_X^2 at Q = 2e5).  Own reconstruction of the history (independent of fx_formation's code; compared in the
     json with its table rows).

 Linear perturbations (sub-horizon Newtonian, smooth radiation, standard Planck-like background; baryons smooth until
 z = 1100 then pressureless; identical in the CDM reference):
   tight coupling (a < a_m): the dark matter (objects + gas) is one fluid with sound speed
       c_eff^2 = gamma n_v T / rho_DM,  gamma = (5/2 + x)/(3/2 + x - x/(5/2 + x)) (two-phase, x = Delta/T) in
       chemical equilibrium, 5/3 after chemical freeze-out;
       delta'' + H delta' + c_eff^2 k^2 delta = 4 pi G a^2 delta_rho ;
   after a_m: cold objects (fraction 1 - f_v of the DM mass) plus a collisionless Maxwellian gas (fraction f_v) with
       1D comoving momentum dispersion sigma_p = a_m sqrt(T(a_m)/m), solved with the exact Newtonian Volterra
       (Gilbert) equation  delta_g(tau) = [delta_m + a_m delta'_m S] e^{-k^2 sigma_p^2 S^2/2}
                                          + int dtau' a(tau') (3/2) H'^2 Omega_m' delta_src(tau') S e^{-k^2 sigma_p^2 S^2/2},
       S(tau, tau') = int_{tau'}^{tau} dtau''/a.  A fluid treatment of the gas (c_s^2 = T/m ~ a^-2) is run as a
       control.
 Initial condition: delta = ln(0.6 k tau) with d delta/d ln a = 1 at k tau = 10 (Hu & Sugiyama 1996 sub-horizon RD
 form [R]); the same in M1 and CDM, and only the ratio P_M1/P_CDM at z = 0 is reported.
 Outputs results/vapour_transfer.json.
"""
import json, math, os, sys, time
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')

# ---------------- background (Planck-2018-like) ----------------
h = 0.674
Omh2_c, Omh2_b = 0.1200, 0.0224
T0_GeV = 2.7255 * 8.617333e-14
Og_h2 = 2.4728e-5
Or_h2 = Og_h2 * (1 + 0.2271 * 3.044)
Om = (Omh2_c + Omh2_b) / h**2
Ob = Omh2_b / h**2
Oc = Omh2_c / h**2
Or = Or_h2 / h**2
OL = 1 - Om - Or
H0_Mpc = h / 2997.92458              # 1/Mpc
H0_GeV = 2.1332e-42 * h
MX = 0.0924089637375362
YQ = 5.07941423e-9
S0_GeV3 = 2891.2 * (1.973269804e-14)**3     # entropy density today [GeV^3]
NQ0 = YQ * S0_GeV3 / MX**3                  # n_Q a^3 in m_X^3
A_EQ = Or / Om
A_REC = 1 / 1101.0
NB = 449.31336788


# radiation with g*(T), g_s(T) steps [R: standard SM values], smoothly interpolated in ln T; entropy conservation
_lgT = np.array([-14, -4.0, -3.4, -3.0, -2.6, -2.2, -1.3, -0.7, -0.5, 0.0, 1.0, 2.0, 3.0])
_gs = np.array([3.909, 3.909, 3.909, 5.6, 10.75, 10.75, 10.75, 17.25, 47.5, 61.75, 86.25, 106.75, 106.75])
_gr = np.array([3.363, 3.363, 3.363, 5.0, 10.75, 10.75, 10.75, 17.25, 47.5, 61.75, 86.25, 106.75, 106.75])
_T = np.logspace(-14, 3, 3000)
_a_of_T = (3.909 / np.interp(np.log10(_T), _lgT, _gs))**(1 / 3) * T0_GeV / _T
_lna_tab = np.log(_a_of_T)[::-1]
_lnT_tab = np.log(_T)[::-1]
MPL = 1.22089e19


def T_of_a(a):
    return math.exp(np.interp(math.log(a), _lna_tab, _lnT_tab))


def rho_r_GeV4(a):
    T = T_of_a(a)
    return np.interp(math.log10(T), _lgT, _gr) * math.pi**2 / 30 * T**4


RHOC0 = 3 * H0_GeV**2 * MPL**2 / (8 * math.pi)
Or_eff0 = rho_r_GeV4(1.0) / RHOC0


def E2(a):
    return rho_r_GeV4(a) / RHOC0 + Om / a**3 + OL


def H_Mpc(a):
    return H0_Mpc * math.sqrt(E2(a))


def H_mX(a):
    return H0_GeV * math.sqrt(E2(a)) / MX


def T_vis_GeV(a):
    return T_of_a(a)




def Om_a(a):
    return (Om / a**3) / E2(a)


# ---------------- object family (from results/objects.json + fit) ----------------
OBJ = json.load(open(os.path.join(RES, 'objects.json')))


def EQ_fit(Q):
    return 0.775865 + 8.1429 * Q**(-1 / 3) + 125.79 * Q**(-2 / 3) - 2847.0 / Q


def Omega_fit(Q):
    return 0.775865 + (2 / 3) * 8.1429 * Q**(-1 / 3) + (1 / 3) * 125.79 * Q**(-2 / 3)


def obj_props(Q):
    rows = {r['Q']: r for r in OBJ['rows']}
    if Q in rows:
        r = rows[Q]
        return dict(Omega=r['Omega'], EQ=r['E_over_Qm'], R=r['r90'])
    return dict(Omega=Omega_fit(Q), EQ=EQ_fit(Q), R=(3 * Q / (4 * math.pi * NB))**(1 / 3))


# ---------------- (B) equilibrium vapour history ----------------
def vapour_T(nQ, SQ, Delta):
    """solve n_v (5/2 + Delta/T) = SQ nQ for T (m_X units); returns T, n_v."""
    target = SQ * nQ
    f = lambda lnT: math.log((math.exp(lnT) / (2 * math.pi))**1.5) - Delta / math.exp(lnT) + math.log(2.5 + Delta / math.exp(lnT)) - math.log(target)
    lo, hi = math.log(1e-5), math.log(0.3)
    if f(hi) < 0:       # vapour cannot hold the entropy: (near) vapour-only, out of this model's range
        return None, None
    lnT = brentq(f, lo, hi, xtol=1e-13)
    T = math.exp(lnT)
    return T, (T / (2 * math.pi))**1.5 * math.exp(-Delta / T)


def history_B(Q, SQ, stick, sigma_m=None, a_lo=1e-14, a_hi=1.0, n=8000):
    p = obj_props(Q)
    Delta = 1 - p['Omega']
    R = p['R']
    sig_c = stick * math.pi * R * R
    sig_m = sigma_m if sigma_m is not None else math.pi * R * R
    lna = np.linspace(math.log(a_lo), math.log(a_hi), n)
    Nv = None; Tv = None; af = None; am = None; Tf = None; valid = False
    rows = []
    frozen = False
    for i, la in enumerate(lna):
        a = math.exp(la)
        nQ = NQ0 / a**3
        Teq, nveq = vapour_T(nQ, SQ, Delta)
        if Teq is None or (Nv is None and Teq > 0.05):
            continue          # start the NR two-phase model only once T_h < 0.05 m_X (equilibrium assumed before)
        Neq = min(nveq / nQ, 1.0)
        if Nv is None:
            Nv, Tv = Neq, Teq
        eta = 1 - Nv
        nobj = eta * nQ / Q
        Hh = H_mX(a)
        if not frozen:
            Tv = Teq
            vb = math.sqrt(8 * Tv / math.pi)
            Gc = nobj * sig_c * vb / Hh
            dl = lna[i] - lna[i - 1] if i > 0 else 0.0
            Nv = Neq + (Nv - Neq) * math.exp(-Gc * dl) if Gc * dl < 700 else Neq
            if Gc < 1.0:
                frozen = True; af = a; Tf = Tv
                valid = (i > 0) and (Teq < 0.05)
        else:
            Tv = Tf * (af / a)**2
            vb = math.sqrt(8 * Tv / math.pi)
            Gc = nobj * sig_c * vb / Hh
        Gm = nobj * sig_m * math.sqrt(8 * Tv / math.pi) / Hh
        if am is None and frozen and Gm < 1.0:
            am = a
        rows.append((a, Teq, Neq, Nv, Tv, Gc, Gm))
        if am is not None and a > 10 * am:
            break
    rows = np.array(rows)
    if af is None or am is None:
        return None
    x = Delta / np.maximum(rows[:, 4], 1e-30)
    out = dict(Q=Q, S_over_Q=SQ, sticking=stick, Delta=Delta, R_mX=R, sigma_c=sig_c, sigma_m=sig_m,
               freeze_out_inside_NR_model=bool(valid), T_eq_at_a_lo=float(rows[0, 1]),
               a_f=af, Tvis_f_GeV=T_vis_GeV(af), T_h_f=Tf, one_minus_eta_f=float(rows[-1, 3]),
               a_m=am, Tvis_m_GeV=T_vis_GeV(am), T_h_m=Tf * (af / am)**2,
               v_rms_at_m=math.sqrt(3 * Tf * (af / am)**2), sigma_p_today=am * math.sqrt(Tf * (af / am)**2),
               EQ=p['EQ'])
    out['sample'] = [dict(a=float(r[0]), Tvis_GeV=T_vis_GeV(r[0]), T_eq=float(r[1]), Nv_eq=float(r[2]), Nv=float(r[3]),
                          T_v=float(r[4]), Gc_over_H=float(r[5]), Gm_over_H=float(r[6])) for r in rows[::max(1, len(rows) // 25)]]
    # c_eff^2(a) table for the perturbation code
    a_arr = rows[:, 0]
    Nv_arr = rows[:, 3]; Tv_arr = rows[:, 4]
    gam = np.where(a_arr < af, (2.5 + x) / (1.5 + x - x / (2.5 + x)), 5 / 3)
    rho = (1 - Nv_arr) * p['EQ'] + Nv_arr * (1 + 1.5 * Tv_arr)
    out['_ceff2'] = (a_arr, gam * Nv_arr * Tv_arr / rho)
    out['f_v_mass'] = float(Nv_arr[-1] / ((1 - Nv_arr[-1]) * p['EQ'] + Nv_arr[-1]))
    # equivalent thermal-WDM mass of the gas alone (match 1D velocity dispersion today; 1 keV <-> 7.8e-8 c) [E]
    out['m_WDM_equiv_keV'] = (7.8e-8 / out['sigma_p_today'])**0.75
    # comoving sound horizon of the coupled fluid up to a_m
    tau_int = lambda a: math.sqrt(np.interp(a, a_arr, out['_ceff2'][1])) / (a * a * H_Mpc(a))
    aa = np.geomspace(a_arr[0], am, 4000)
    out['r_s_Mpc'] = float(np.trapezoid([tau_int(z) for z in aa], aa))
    return out


def history_A():
    """Codex's gas with elastic coupling only (sticking 0) to Q = 2e5 cores."""
    f = 0.065946
    a2 = T0_GeV * (3.909 / 10.75)**(1 / 3) / 2e-3
    v2 = 0.0684
    T2 = v2**2 / 3
    p = obj_props(2e5)
    sig_m = 332.0
    lna = np.linspace(math.log(a2), 0.0, 4000)
    am = None
    rows = []
    for la in lna:
        a = math.exp(la)
        nQ = NQ0 / a**3
        nobj = (1 - f) * nQ / 2e5
        T = T2 * (a2 / a)**2
        Gm = nobj * sig_m * math.sqrt(8 * T / math.pi) / H_mX(a)
        rows.append((a, T, Gm))
        if am is None and Gm < 1:
            am = a
    rows = np.array(rows)
    rho = (1 - f) * p['EQ'] + f
    out = dict(case='A_codex_gas', f_charge=f, a_start=a2, Gm_over_H_start=float(rows[0, 2]), a_m=am, Tvis_m_GeV=T_vis_GeV(am),
               T_h_m=T2 * (a2 / am)**2, sigma_p_today=am * math.sqrt(T2 * (a2 / am)**2), f_v_mass=f / rho)
    out['_ceff2'] = (rows[:, 0], 5 / 3 * f * rows[:, 1] / rho)
    out['m_WDM_equiv_keV'] = (7.8e-8 / out['sigma_p_today'])**0.75
    aa = rows[rows[:, 0] <= am, 0]
    out['r_s_Mpc'] = float(np.trapezoid(np.sqrt(out['_ceff2'][1][:len(aa)]) / (aa * aa * np.array([H_Mpc(z) for z in aa])), aa))
    return out


# ---------------- perturbations ----------------
def a_of_ktau(k, ktau=10.0):
    # tau(a) = int da / (a^2 H)
    f = lambda la: math.log(k * quad(lambda z: 1 / (z * z * H_Mpc(z)), 1e-14, math.exp(la), limit=200)[0]) - math.log(ktau)
    return math.exp(brentq(f, math.log(1e-13), math.log(0.5)))


def tau_of(a):
    return quad(lambda z: 1 / (z * z * H_Mpc(z)), 1e-14, a, limit=200)[0]


def growth(k, hist, gas_mode='gilbert', dN=0.004, reference=False):
    """delta_m(z=0).  Before a_m: DOP853 on the single DM fluid (pressure from hist unless reference); after a_m:
    Heun predictor-corrector marching for objects/baryons with the gas from the Gilbert integral (or a fluid).
    reference=True runs the identical pipeline with c_eff = 0 and f_v = 0 (pure CDM), so integrator errors cancel."""
    a_i = a_of_ktau(k, 10.0)
    d0, s0 = math.log(6.0), 1.0
    fb = Ob / Om
    fc = Oc / Om
    a_m = hist['a_m']
    fv = 0.0 if reference else hist['f_v_mass']
    sigp = 0.0 if reference else hist['sigma_p_today']

    def ceff2(a):
        if reference or a >= a_m:
            return 0.0
        aa, cc = hist['_ceff2']
        if a <= aa[0]:
            return float(cc[0])
        return float(np.interp(math.log(a), np.log(aa), cc))

    def rhs(N, y):
        a = math.exp(N)
        Hc = a * H_Mpc(a)
        dd, ud, db, ub = y
        bar_on = a >= A_REC
        src = fc * dd + (fb * db if bar_on else 0.0)
        acc = 1.5 * Hc * Hc * Om_a(a)
        dud = a / Hc * (acc * src - ceff2(a) * k * k * dd)
        dub = a / Hc * acc * src if bar_on else 0.0
        return [ud / (a * Hc), dud, ub / (a * Hc) if bar_on else 0.0, dub]

    def run_ode(N0, N1, y):
        segs = [N0]
        if N0 < math.log(A_REC) < N1:
            segs.append(math.log(A_REC))
        segs.append(N1)
        for s0_, s1_ in zip(segs[:-1], segs[1:]):
            sol = solve_ivp(rhs, (s0_, s1_), y, method='DOP853', rtol=1e-9, atol=1e-13)
            y = list(sol.y[:, -1])
        return y

    Ni = math.log(a_i)
    Hci = a_i * H_Mpc(a_i)
    y0 = [d0, a_i * Hci * s0, 0.0, 0.0]
    Nm = math.log(a_m)
    if Ni >= Nm:
        Nm = Ni
        ym = list(y0)
    else:
        ym = run_ode(Ni, Nm, list(y0))
    Ngrid = np.linspace(Nm, 0.0, max(int(round((0.0 - Nm) / dN)), 50) + 1)
    dNs = Ngrid[1] - Ngrid[0]
    a = np.exp(Ngrid)
    Hc = a * np.array([H_Mpc(z) for z in a])
    Oma = np.array([Om_a(z) for z in a])
    C = np.concatenate([[0.0], np.cumsum(0.5 * (1 / (a[1:] * Hc[1:]) + 1 / (a[:-1] * Hc[:-1])) * dNs)])
    w_src = a * Hc * 1.5 * Oma
    A = a / Hc * 1.5 * Hc * Hc * Oma          # du/dN = A * src
    B = 1 / (a * Hc)                           # d delta/dN = B * u
    n = len(Ngrid)
    do = np.zeros(n); dg = np.zeros(n); db = np.zeros(n); src = np.zeros(n)
    uo, ub = ym[1], ym[3]
    do[0] = ym[0]; dg[0] = ym[0]; db[0] = ym[2]
    dprime0 = ym[1] / a[0]
    fo = fc * (1 - fv); fg = fc * fv
    bar = a >= A_REC
    src[0] = fo * do[0] + fg * dg[0] + (fb * db[0] if bar[0] else 0.0)
    ug = ym[1]
    for i in range(1, n):
        # gas at i (explicit: the kernel vanishes at tau' = tau)
        if gas_mode == 'gilbert' or reference:
            S = C[i] - C[:i + 1]
            damp = np.exp(-0.5 * (k * sigp * S)**2)
            free = (ym[0] + a[0] * dprime0 * S[0]) * damp[0]
            dg[i] = free + np.trapezoid(w_src[:i + 1] * src[:i + 1] * S * damp, dx=dNs)
        # predictor
        dop = do[i - 1] + dNs * B[i - 1] * uo
        uop = uo + dNs * A[i - 1] * src[i - 1]
        if bar[i - 1] or bar[i]:
            dbp = db[i - 1] + dNs * B[i - 1] * ub
            ubp = ub + dNs * A[i - 1] * src[i - 1]
        else:
            dbp, ubp = 0.0, 0.0
        if gas_mode == 'fluid' and not reference:
            c2a = sigp**2 / a[i - 1]**2; c2b = sigp**2 / a[i]**2
            dgp = dg[i - 1] + dNs * B[i - 1] * ug
            ugp = ug + dNs * (A[i - 1] * src[i - 1] - a[i - 1] / Hc[i - 1] * c2a * k * k * dg[i - 1])
        else:
            dgp = dg[i]
        srcp = fo * dop + fg * dgp + (fb * dbp if bar[i] else 0.0)
        # corrector
        do[i] = do[i - 1] + 0.5 * dNs * (B[i - 1] * uo + B[i] * uop)
        uo = uo + 0.5 * dNs * (A[i - 1] * src[i - 1] + A[i] * srcp)
        if bar[i - 1] or bar[i]:
            db[i] = db[i - 1] + 0.5 * dNs * (B[i - 1] * ub + B[i] * ubp)
            ub = ub + 0.5 * dNs * (A[i - 1] * src[i - 1] + A[i] * srcp)
        if gas_mode == 'fluid' and not reference:
            dg[i] = dg[i - 1] + 0.5 * dNs * (B[i - 1] * ug + B[i] * ugp)
            ug = ug + 0.5 * dNs * (A[i - 1] * src[i - 1] - a[i - 1] / Hc[i - 1] * c2a * k * k * dg[i - 1]
                                   + A[i] * srcp - a[i] / Hc[i] * c2b * k * k * dgp)
        src[i] = fo * do[i] + fg * dg[i] + (fb * db[i] if bar[i] else 0.0)
    return (Oc * (1 - fv) * do[-1] + Oc * fv * dg[-1] + Ob * db[-1]) / Om


def transfer(hist, ks, gas_mode='gilbert'):
    res = []
    for k in ks:
        d1 = growth(k, hist, gas_mode)
        d0 = growth(k, hist, gas_mode, reference=True)
        res.append(dict(k=k, P_ratio=(d1 / d0)**2))
    return res


def main():
    """usage: vapour_transfer.py A | Qlist (comma-separated) | merge"""
    t0 = time.time()
    ks = [1.0, 3.0, 10.0, 20.0, 50.0, 100.0, 300.0, 1000.0, 3000.0]
    mode = sys.argv[1] if len(sys.argv) > 1 else 'A'
    part = os.path.join(os.environ.get('SP', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad'), 'wt2', 'fx_observables', 'parts')
    os.makedirs(part, exist_ok=True)
    if mode == 'A':
        hA = history_A()
        o = {k: v for k, v in hA.items() if not k.startswith('_')}
        o['transfer'] = transfer(hA, [100.0, 1e3, 1e4, 3e4])
        print('A:', o, '%.0fs' % (time.time() - t0), flush=True)
        # CDM control of the marching integrator
        fake = dict(a_m=1e-6, f_v_mass=0.5, sigma_p_today=0.0, _ceff2=(np.array([1e-12, 1e-6]), np.array([0.0, 0.0])))
        o['control_cold_gas_ratio_k3_k100'] = [growth(k, fake) / growth(k, fake, reference=True) for k in (3.0, 100.0)]
        json.dump(o, open(os.path.join(part, 'A.json'), 'w'), indent=1, default=float)
        return
    if mode == 'merge':
        out = dict(scope=__doc__.strip().splitlines()[0], background=dict(h=h, Om=Om, Ob=Ob, Or=Or, a_eq=A_EQ, NQ0_mX3=NQ0))
        out['A'] = json.load(open(os.path.join(part, 'A.json')))
        B = []
        for fn in sorted(os.listdir(part)):
            if fn.startswith('B_'):
                B += json.load(open(os.path.join(part, fn)))
        B.sort(key=lambda r: (r['Q'], r['S_over_Q'], -r['sticking']))
        out['B'] = B
        with open(os.path.join(RES, 'vapour_transfer.json'), 'w') as f:
            json.dump(out, f, indent=1, default=float)
        print('merged', len(B))
        return
    Qs = [float(x) for x in mode.split(',')]
    scen = []
    for Q in Qs:
        for SQ in (1.235, 2.07, 3.4):
            for st in (1.0, 1e-2, 1e-4):
                hB = history_B(Q, SQ, st, sigma_m=332.0 if Q == 2e5 else None)
                if hB is None:
                    scen.append(dict(Q=Q, S_over_Q=SQ, sticking=st, note='no freeze-out in range')); continue
                row = {k: v for k, v in hB.items() if not k.startswith('_')}
                print('B Q=%.0e S/Q=%.2f s=%g: Tvis_f=%.3g GeV T_f=%.4f 1-eta=%.4f f_v=%.4f Tvis_m=%.3g sig_p=%.2e mWDM=%.3g keV r_s=%.3g Mpc' % (
                    Q, SQ, st, row['Tvis_f_GeV'], row['T_h_f'], row['one_minus_eta_f'], row['f_v_mass'], row['Tvis_m_GeV'],
                    row['sigma_p_today'], row['m_WDM_equiv_keV'], row['r_s_Mpc']), flush=True)
                if SQ in (1.235, 2.07):
                    row['transfer'] = transfer(hB, ks)
                    print('  P_M1/P_CDM:', ' '.join('%g:%.4f' % (r['k'], r['P_ratio']) for r in row['transfer']), '%.0fs' % (time.time() - t0), flush=True)
                    if Q in (1e12, 1e14) and SQ == 1.235 and st == 1.0:
                        row['fluid_gas_control'] = transfer(hB, [3.0, 10.0, 50.0, 300.0], 'fluid')
                scen.append(row)
    json.dump(scen, open(os.path.join(part, 'B_%s.json' % mode.replace(',', '_')), 'w'), indent=1, default=float)
    print('done %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
