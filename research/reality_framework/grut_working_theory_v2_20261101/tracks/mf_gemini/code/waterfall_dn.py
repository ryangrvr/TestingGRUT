"""mf_gemini: the mild waterfall of the proposal as written, linear modes + delta N.

Code units M_P = H0 = 1, H0 = sqrt(V0/3) (the proposal's H_i).  Proposal: v^2 = 8/3, lam = 27/16,
m_chi = 3, g^2 sigma_c^2 = lam v^2.  V_SR near the band |sigma| < sigma_c is a model choice; two shapes:
  'lin'  constant slope p through the band (sigma can leave the band on the far side);
  'quad' convex 3 Delta sigma^2/2 near sigma_c (slope p = Delta sigma_c; sigma decelerates, stays in band).
Delta = p/sigma_c is the fraction of sigma_c crossed per e-fold at sigma_c.

Steps:
 1. sigma-only background (chi = 0; at linear order chi modes decouple from the metric).
 2. Linear chi modes (Bunch-Davies at k/aH = 100):  u'' + 3h u' + (k^2/a^2 + g^2 sigma^2 - lam v^2) u = 0.
 3. chi_rms from super-horizon modes; match at chi_rms = f_m chi_trap, chi_trap = (2 p sigma_c/3)^{1/2}
    (below it the back-reaction g^2 chi^2 sigma on sigma is < f_m^2 of the slope).
 4. Nonlinear background from the match with chi = c chi_rms (full equations, then the WKB-averaged sigma
    once g chi > 30 h); N_end(c) at eps_H = 1.
 5. Linear delta N: P_zeta(k) = (dN/dchi_m)^2 P_chi(k, t_m).  Nonlinear delta N: zeta = N(chi_L) - <N>,
    chi_L ~ Gaussian(0, chi_rms).
 6. PBH masses for the peak scale (instant radiation after the end), beta (Gaussian, nonlinear tail),
    f_PBH by radiation-era scaling.  [E] for 6 (threshold, gamma = 0.2, no non-Gaussian collapse model).

Run: python3 -B waterfall_dn.py part 0 7; python3 -B waterfall_dn.py part 7 14; python3 -B waterfall_dn.py combine
     (-> results/waterfall_dn.json; ~5 min per part, one core; per-case files in the scratch dir SCR)
     python3 -B waterfall_dn.py regress  (de Sitter mode regression)
"""
import json
import os
import sys
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import erf, erfc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mfg_common import Model, slow_roll_ic, MP_GEV, eps_star, two_stage

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
GEV_G = 1.78266192e-24   # g per GeV
MSUN_G = 1.98892e33


def sigma_background(m, s_start, n_after=60.0):
    y0 = slow_roll_ic(m, s_start)
    V0c = 0.25 * m.lam * m.v2**2

    def rhs(t, y):
        N, s, ds = y
        h = np.sqrt((0.5 * ds**2 + m.VSR(s) + V0c) / 3.0)
        return [h, ds, -3 * h * ds - m.dVSR(s)]

    def evc(t, y):
        return y[1] - m.sigc
    state = {}

    def evN(t, y):
        if 'Nc' not in state and y[1] <= m.sigc:
            state['Nc'] = y[0]
        return y[0] - (state.get('Nc', 1e9) + n_after)
    evN.terminal = True
    sol = solve_ivp(rhs, (0, 4000.0), y0[:3], method='DOP853', rtol=1e-11, atol=1e-16,
                    dense_output=True, events=[evc, evN])
    tc = sol.t_events[0][0]
    Nc = sol.sol(tc)[0]
    return sol, rhs, tc, Nc


def integrate_mode(m, rhs_bg, y_bg0, t0, k, t_out):
    N0, s0, ds0 = y_bg0
    a0 = np.exp(N0)
    h0 = rhs_bg(t0, y_bg0)[0]
    m2 = m.g2 * s0**2 - m.lam * m.v2
    w = np.sqrt((k / a0)**2 + m2)
    C = a0**1.5 * np.sqrt(2 * w)     # u_true = u_scaled/C ; BD/WKB: u = 1/(a^{3/2} sqrt(2 w_phys))
    du = (-1j * w - h0)
    y0 = [N0, s0, ds0, 1.0, 0.0, du.real, du.imag]

    def rhs(t, y):
        N, s, ds, ur, ui, dur, dui = y
        bg = rhs_bg(t, y[:3])
        h = bg[0]
        mm = (k * np.exp(-N))**2 + m.g2 * s**2 - m.lam * m.v2
        return [bg[0], bg[1], bg[2], dur, dui, -3 * h * dur - mm * ur, -3 * h * dui - mm * ui]
    sol = solve_ivp(rhs, (t0, t_out[-1]), y0, method='DOP853', rtol=1e-9, atol=1e-12, t_eval=t_out)
    return (sol.y[3] + 1j * sol.y[4]) / C


def run_case(sigc, vsr, H_GeV=1.3e11, f_m=1e-2, Nb=12.0, Na=60.0, dlnk=0.25, label=''):
    t_start = time.time()
    hub = H_GeV / MP_GEV
    m = Model(sigc, vsr, hub=hub)
    p_c = m.dVSR(sigc) / 3.0
    # start Npre + 5 e-folds before sigma_c so that modes exiting from Npre e-folds before sigma_c
    # (chi heavy there) are resolved; slow-roll IC
    D = p_c / sigc
    Npre = min(max(4.0, 0.3 / D), 30.0)
    if vsr['kind'] == 'quad':
        s_start = sigc * min(np.exp(D * (Npre + 5.0)), 3.0)
    else:
        s_start = min(sigc + p_c * (Npre + 5.0), 3.0 * sigc)
    sol, rhs_bg, tc, Nc = sigma_background(m, s_start, n_after=Na)
    tt = np.linspace(0, sol.t[-1], 200001)
    NN = sol.sol(tt)[0]
    N0 = NN[0]
    Nouts = np.arange(N0 + 0.3, NN[-1] - 0.5, 0.05)
    touts = np.interp(Nouts, NN, tt)
    ybs = sol.sol(touts)
    hs = np.array([rhs_bg(0, ybs[:, i])[0] for i in range(len(touts))])
    ah = np.exp(ybs[0]) * hs
    Nexit = np.arange(N0 + 0.5, min(Nc + 40.0, NN[-1] - 2.0), dlnk)
    ks, P = [], []
    for Ne in Nexit:
        te = np.interp(Ne, NN, tt)
        k = np.exp(Ne) * rhs_bg(0, sol.sol(te))[0]
        ts = max(np.interp(Ne - np.log(100.0), NN, tt), 0.0)   # WKB start (adiabatic: chi heavy or k/aH >= 100)
        sel = touts > ts + 1e-9
        Pk = np.zeros(len(touts))
        u = integrate_mode(m, rhs_bg, sol.sol(ts), ts, k, touts[sel])
        Pk[sel] = k**3 * np.abs(u)**2 / (2 * np.pi**2)
        ks.append(k)
        P.append(Pk)
    ks = np.array(ks)
    P = np.array(P)
    mask = ks[:, None] < ah[None, :]
    chirms = np.sqrt((P * mask).sum(axis=0) * dlnk) * hub        # in M_P
    sig_t = ybs[1]
    chi_trap = np.sqrt(2.0 * p_c * sigc / 3.0)
    target = f_m * chi_trap
    inband = np.abs(sig_t) < sigc
    ok_idx = np.where((chirms >= target) & inband & (Nouts > Nc))[0]
    base = dict(label=label, H_GeV=H_GeV, sigc=sigc, vsr=vsr, g_phys=float(m.g * hub), p_c=float(p_c),
                Delta=float(p_c / sigc), chi_trap=float(chi_trap), N_sigma_c=float(Nc))
    if ok_idx.size == 0:
        base.update(ok=False, reason='chi_rms does not reach f_m chi_trap while sigma is inside the band '
                    '(sigma leaves |sigma| < sigma_c first: no waterfall in this crossing)',
                    chirms_max_over_trap=float(chirms.max() / chi_trap), runtime_s=time.time() - t_start)
        return base
    im = int(ok_idx[0])
    tm, Nm = touts[im], Nouts[im]
    gr = (np.log(chirms[im + 2]) - np.log(chirms[im - 2])) / (touts[im + 2] - touts[im - 2])
    ybg = sol.sol(tm)
    cache = {}

    def endN(c):
        if c in cache:
            return cache[c]
        y0 = np.array([ybg[0], ybg[1], ybg[2], c * chirms[im], c * chirms[im] * gr])
        r2 = two_stage(m, y0)
        cache[c] = (r2['N_end'], r2)
        return cache[c]
    N1, r1 = endN(1.0)
    if not np.isfinite(N1):
        base.update(ok=False, reason='classical trajectory from chi_rms leaves the band', runtime_s=time.time() - t_start)
        return base
    cs = np.exp(np.linspace(np.log(1e-4), np.log(6.0), 23))
    Ncur = np.array([endN(c)[0] for c in cs])
    good = np.isfinite(Ncur)
    dN = (endN(1.02)[0] - endN(0.98)[0]) / (np.log(1.02) - np.log(0.98))
    s_eff = -1.0 / dN
    # Gaussian average over chi_L = |g| chi_rms; outside the c grid: asymptotic log law (small c)
    gq, wq = np.polynomial.hermite_e.hermegauss(240)
    wq = wq[gq > 0] / wq[gq > 0].sum()
    gq = gq[gq > 0]
    lnc = np.log(cs[good])
    Ng = Ncur[good]
    Nint = np.interp(np.log(gq), lnc, Ng)
    lo = np.log(gq) < lnc[0]
    Nint[lo] = Ng[0] - (np.log(gq[lo]) - lnc[0]) / s_eff
    Nmean = float((wq * Nint).sum())
    varN = float((wq * (Nint - Nmean)**2).sum())

    def tail(x):
        tgt = Nmean + x
        if tgt >= Ng[0]:
            cstar = cs[good][0] * np.exp(-(tgt - Ng[0]) * s_eff)
        else:
            cstar = float(np.exp(np.interp(-tgt, -Ng, lnc)))
        return float(erf(cstar / np.sqrt(2.0)))
    dNdchi = dN / chirms[im]
    exited = ks < ah[im]
    Pzeta = dNdchi**2 * P[exited, im] * hub**2
    Nexit_rel = Nexit[exited] - Nc
    Nbefore = N1 - Nexit[exited]
    ipk = int(np.argmax(Pzeta))
    hend = r1['h_end']
    Hend_GeV = hend * H_GeV
    gamma = 0.2

    def mpbh(nb):
        return gamma * 4 * np.pi * MP_GEV**2 / Hend_GeV * np.exp(2 * nb) * GEV_G

    def fpbh(beta, Mg):
        MH = Mg / gamma / GEV_G
        Hf = 4 * np.pi * MP_GEV**2 / MH
        gf = 106.75
        Tf = (90.0 / (np.pi**2 * gf))**0.25 * np.sqrt(Hf * MP_GEV)
        Teq = 0.7965e-9
        return beta * (Tf / Teq) * (gf / 3.909)**(1.0 / 3.0) * (0.315 / 0.264)
    Mpk = mpbh(Nbefore[ipk])
    Ppk = float(Pzeta[ipk])
    zc = [0.5, 0.67, 1.0]
    tails = {f'{x}': tail(x) for x in zc}
    bg = {f'{x}': float(0.5 * erfc(x / np.sqrt(2 * Ppk))) for x in zc}
    base.update(
        ok=True, N_match_after_sigc=float(Nm - Nc), chi_rms_match=float(chirms[im]), growth_rate_match=float(gr),
        N_end_after_sigc=float(N1 - Nc), N_mean_after_sigc=float(Nmean - Nc), varN=varN, s_eff=float(s_eff),
        H_end_over_H0=float(hend), end_state=r1.get('end_state'), N_switch_after_sigc=float(r1.get('N_switch', np.nan) - Nc),
        Pzeta_peak=Ppk, Pzeta_peak_Nexit_after_sigc=float(Nexit_rel[ipk]), Pzeta_peak_Nbefore_end=float(Nbefore[ipk]),
        sum_Pzeta_dlnk=float(Pzeta.sum() * dlnk), M_PBH_peak_g=float(Mpk),
        beta_tail_nonlinear=tails, beta_gaussian_peak=bg,
        f_PBH_tail_zc067=float(fpbh(tails['0.67'], Mpk)), f_PBH_gauss_zc067=float(fpbh(bg['0.67'], Mpk)),
        frac_band_exit_in_cgrid=float(1 - good.mean()),
        spectrum=dict(N_exit_after_sigc=[float(x) for x in Nexit_rel[::2]], Pzeta=[float(x) for x in Pzeta[::2]]),
        Ncurve=dict(c=[float(x) for x in cs], N_after_sigc=[float(x - Nc) for x in Ncur]),
        runtime_s=time.time() - t_start)
    return base


def regression():
    """massless and m^2 = 2 H^2 (conformal) modes in exact de Sitter vs the analytic spectra."""
    out = {}
    for m2 in [0.0, 2.0, -4.5]:
        def rhs_bg(t, y):
            return [1.0, 0.0, 0.0]
        class M:
            pass
        mm = M()
        mm.g2 = 0.0
        mm.lam = -m2       # g^2 sigma^2 - lam v^2 = m2 with v2 = 1
        mm.v2 = 1.0
        N0 = 3.0
        k = np.exp(N0 + np.log(100.0))
        tout = np.array([N0 + np.log(100.0) + 3.0])
        u = integrate_mode(mm, rhs_bg, np.array([N0, 0.0, 0.0]), N0, k, tout)
        P = k**3 * np.abs(u[0])**2 / (2 * np.pi**2)
        x = np.exp(-3.0)                       # k/aH at read-out
        if m2 == 0.0:
            exact = (1 / (2 * np.pi))**2 * (1 + x**2)
        elif m2 == 2.0:
            exact = (1 / (2 * np.pi))**2 * x**2      # conformally coupled: P = (H/2pi)^2 (k/aH)^2
        else:
            exact = None
        out[f'm2={m2}'] = dict(P=float(P), exact=None if exact is None else float(exact))
    return out


def todo_list():
    H = 1.3e11
    p_star = np.sqrt(2 * eps_star(H))
    p_c = p_star * np.exp(0.015 * 36.0)
    todo = []
    for D in [0.01, 0.03, 0.1, 0.3, 1.0, 3.0]:
        todo.append((p_c / D, dict(kind='quad', eta=D), f'quad_D{D:g}', H))
    for D in [0.03, 0.06, 0.1, 0.2]:
        todo.append((p_c / D, dict(kind='lin', p=p_c), f'lin_D{D:g}', H))
    for D in [0.1, 1.0]:
        todo.append((p_star / D, dict(kind='quad', eta=D), f'quad_D{D:g}_pstar', H))
    for Hx in [1.25e11, 1.55e11]:
        pcx = np.sqrt(2 * eps_star(Hx)) * np.exp(0.015 * 36.0)
        todo.append((pcx / 0.3, dict(kind='quad', eta=0.3), f'quad_D0.3_H{Hx:.3g}', Hx))
    return todo


SCR = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/mf_gemini/wdn'


def run_part(i0, i1):
    os.makedirs(SCR, exist_ok=True)
    for sigc, vsr, lab, H in todo_list()[i0:i1]:
        fn = os.path.join(SCR, lab + '.json')
        if os.path.exists(fn):
            continue
        r = run_case(sigc, vsr, H_GeV=H, label=lab)
        keys = ['ok', 'Delta', 'N_end_after_sigc', 'N_mean_after_sigc', 'varN', 's_eff', 'Pzeta_peak',
                'Pzeta_peak_Nexit_after_sigc', 'Pzeta_peak_Nbefore_end', 'M_PBH_peak_g', 'beta_tail_nonlinear',
                'f_PBH_tail_zc067', 'H_end_over_H0', 'runtime_s', 'reason']
        print(lab, {k: r.get(k) for k in keys}, flush=True)
        with open(fn, 'w') as f:
            json.dump(r, f)


def combine():
    H = 1.3e11
    p_star = np.sqrt(2 * eps_star(H))
    out = dict(model='proposal as written: v^2 = 8/3 M_P^2, lam = 27/16 (H0/M_P)^2, m_chi = 3 H0',
               H_GeV=H, eps_star=float(eps_star(H)), p_star=float(p_star),
               p_c=float(p_star * np.exp(0.54)), cases=[])
    out['regression_de_Sitter_modes'] = regression()
    for sigc, vsr, lab, Hx in todo_list():
        with open(os.path.join(SCR, lab + '.json')) as f:
            out['cases'].append(json.load(f))
    with open(os.path.join(RES, 'waterfall_dn.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'one':
        H = 1.3e11
        pc = np.sqrt(2 * eps_star(H)) * np.exp(0.54)
        D = float(sys.argv[2])
        r = run_case(pc / D, dict(kind=sys.argv[3], eta=D, p=pc), H_GeV=H, label='one')
        print({k: v for k, v in r.items() if k not in ('spectrum', 'Ncurve')})
        print(r.get('spectrum'))
        print(r.get('Ncurve'))
    elif len(sys.argv) > 1 and sys.argv[1] == 'part':
        run_part(int(sys.argv[2]), int(sys.argv[3]))
    elif len(sys.argv) > 1 and sys.argv[1] == 'regress':
        print(regression())
    elif len(sys.argv) > 1 and sys.argv[1] == 'combine':
        combine()
    else:
        run_part(0, 100)
        combine()
