"""qm_micro item 3 (Claude, 2 Nov 2026): coarsening of the Q-matter object population from the spinodal to capture
freeze-out, with the computed sticking.  Own code.  [C] for the numerics of the stated model, [E] for the model.

Model (m_X = 1; all object quantities per unit comoving net charge):
  * background: the two-phase equilibrium history at fixed S/Q (qm_vapour.history: T_h(a), n_Q(a), mu_c(T), vapour
    n_X, n_Xbar, (e+P), the sticking-weighted one-sided fluxes Phi_X(mu_c), Phi_Xbar(-mu_c), H(a) from the visible clock);
  * objects: Lagrangian cohorts (Q_c, N_c), bulk radius R = (3Q/4 pi n_B)^{1/3}, Gibbs-Thomson
    Omega_c = mu_c(T) + 2 sigma/(n_B R) (sigma = 98.8);
  * exchange with the vapour (ripening, global condensation/evaporation, and capture), by detailed balance with the
    sticking-limited fluxes of item 1:
       dQ_c/dt = 4 pi R^2 { Phi_X [e^{d/T} - e^{(Omega_c - mu_c)/T}] - Phi_Xbar [e^{-d/T} - e^{-(Omega_c - mu_c)/T}] },
    d = mu_v - mu_c the vapour's chemical-potential offset, fixed IMPLICITLY each step by charge conservation
       sum_c N_c Q_c' + [n_X e^{d/T} - n_Xbar e^{-d/T}]/n_Q = 1        (exchange rates are 1e5-1e11 H);
    cohorts that shrink below Q_min dissolve into the vapour;
  * mergers ("hidden-sector kinetic encounters"): Brownian coagulation with the Fuchs-type harmonic kernel
       K = K_c K_fm/(K_c + K_fm),  K_c = 4 pi (R_i+R_j)(D_i+D_j),  D = T/zeta,  zeta = (4 pi/3) R^2 (e+P)_vap vbar (Epstein),
       K_fm = pi (R_i+R_j)^2 sqrt(v_i^2+v_j^2),  v = sqrt(8 T/(pi M)),  M = Omega_b Q,   merger probability 1;
  * freeze-out: the capture rate per vapour X, Gamma_c = sum_c n_c 4 pi R_c^2 Phi_X / n_X, against H; after
    Gamma_c = H the run continues on the equilibrium T(a) (variant a) and, separately, with the decoupled vapour cooling
    as T_v ~ a^-2 and capture only (variant b).
Initial state at the spinodal T_sp(S/Q): lognormal in Q (median Q_0, sigma_lnQ = 0.7) holding eta_eq(T_sp); Q_0 =
4e11 or 1.2e13 (fx_formation's first lumps), or the drag-limited bicontinuous estimate (option 'hydro').
"""
import json, math, os, sys, time
import numpy as np
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

NB = 449.31336788454803
OMB = 0.7758647074595397
SIGMA = 98.8
QMIN = 56709.0
MX = 0.0924089637375362
FA, FB = 8.142857534830613, 125.78683178383355     # E/Q fit of fx_formation (family.json fit3)


def EoverQ(Q):
    return OMB + FA * Q ** (-1 / 3) + FB * Q ** (-2 / 3)


class History:
    def __init__(self, path):
        d = json.load(open(path))
        rows = [r for r in d['rows'] if r['phase'] == 'two-phase']
        rows.sort(key=lambda r: -r['T'])
        self.SQ = d['S_over_Q']; self.YQ = d['Y_Q']
        self.lna = -np.log(np.array([r['nQ'] for r in rows])) / 3
        self.lna -= self.lna[0]
        self.keys = ['T', 'eta', 'nQ', 'mu_c', 'n_X', 'n_Xbar', 'eP_vap', 'Phi_X', 'Phi_Xbar', 'Phi_X_bare', 'H', 'T_vis_GeV', 'phi_vol', 's_vap']
        self.tab = {k: np.array([r[k] for r in rows]) for k in self.keys}
        self.logk = {'T', 'nQ', 'n_X', 'n_Xbar', 'eP_vap', 'Phi_X', 'Phi_Xbar', 'Phi_X_bare', 'H', 'T_vis_GeV', 'phi_vol', 's_vap'}

    stick_mode = 'computed'          # 'computed' | ('scale', x) | 'geometric' (s = 1, as fx_formation assumed)

    def at(self, la):
        out = {}
        for k in self.keys:
            y = self.tab[k]
            if k in self.logk:
                out[k] = float(np.exp(np.interp(la, self.lna, np.log(np.maximum(y, 1e-300)))))
            else:
                out[k] = float(np.interp(la, self.lna, y))
        out['trap'] = (self.stick_mode == 'trap')
        if self.stick_mode == 'geometric':
            f = out['Phi_X_bare'] / out['Phi_X']
            out['Phi_X'] *= f; out['Phi_Xbar'] *= f          # Xbar: same factor (approximate; Xbar is negligible below T ~ 1)
        elif isinstance(self.stick_mode, tuple):
            out['Phi_X'] *= self.stick_mode[1]; out['Phi_Xbar'] *= self.stick_mode[1]
        return out

    def la_of_T(self, T):
        return float(np.interp(-math.log(T), -np.log(self.tab['T']), self.lna))


def radius(Q):
    return (3 * Q / (4 * math.pi * NB)) ** (1 / 3)


_EOS = None


def trap_factor(R, T):
    """[E] probability that the phonon left by a captured X decays (Beliaev, spontaneous only: a lower bound) before it
    can re-emit an X: 1 - exp(-L/mfp), L = 2 x (4R/3) (two mean chords: by reciprocity the phonon lies in the escape
    cone, escaping with probability ~s per wall hit), mfp(T) from results/eos_phonon.json (w = Delta + 1.5 T)."""
    global _EOS
    if _EOS is None:
        d = json.load(open(os.path.join(HERE, '..', 'results', 'eos_phonon.json')))['rows']
        _EOS = (np.log([r['T'] for r in d]), np.log([r['mfp'] for r in d]))
    lt = math.log(T)
    lm = np.interp(lt, *_EOS) if _EOS[0][0] <= lt <= _EOS[0][-1] else (
        _EOS[1][0] + (lt - _EOS[0][0]) * (_EOS[1][1] - _EOS[1][0]) / (_EOS[0][1] - _EOS[0][0]) if lt < _EOS[0][0]
        else _EOS[1][-1] + (lt - _EOS[0][-1]) * (_EOS[1][-1] - _EOS[1][-2]) / (_EOS[0][-1] - _EOS[0][-2]))
    return -np.expm1(-(8.0 / 3.0) * R / math.exp(lm))


def exchange_rate(Q, d, s):
    """dQ/dt per cohort (array Q) at vapour offset d, background s."""
    R = radius(Q)
    g = 2 * SIGMA / (NB * R) / s['T']
    T = s['T']
    return 4 * math.pi * R * R * (s['Phi_X'] * (np.expm1(d / T) - np.expm1(g)) - s['Phi_Xbar'] * (np.expm1(-d / T) - np.expm1(-g)))


RMIN = (3 * QMIN / (4 * math.pi * NB)) ** (1 / 3)


def radial_velocity(R, d, s):
    g = 2 * SIGMA / (NB * np.maximum(R, RMIN)) / s['T']
    T = s['T']
    v = (s['Phi_X'] * (np.expm1(d / T) - np.expm1(g)) - s['Phi_Xbar'] * (np.expm1(-d / T) - np.expm1(-g))) / NB
    return v * trap_factor(R, T) if s.get('trap') else v


def refine(lq, N, nnew=240):
    """Re-resolve a cohort list whose charge-carrying part has thinned out: smooth CDF of the number through the
    cohort midpoints, re-sampled on nnew cohorts; the charge is conserved exactly by a final rescale of N."""
    o = np.argsort(lq); lq, N = lq[o], N[o]
    if lq.size < 3:
        return lq, N
    C = np.cumsum(N) - 0.5 * N
    sp = np.median(np.diff(lq)) if lq.size > 2 else 0.05
    x = np.linspace(lq[0] - 0.5 * sp, lq[-1] + 0.5 * sp, nnew + 1)
    Cx = np.interp(x, np.r_[lq[0] - 0.5 * sp, lq, lq[-1] + 0.5 * sp], np.r_[0.0, C, N.sum()])
    Nn = np.diff(Cx); xn = 0.5 * (x[1:] + x[:-1])
    keep = Nn > 0
    xn, Nn = xn[keep], Nn[keep]
    Nn *= (N * np.exp(lq)).sum() / (Nn * np.exp(xn)).sum()
    return xn, Nn


def coag_kernel(Q, s, epstein=1.0):
    T = s['T']
    R = radius(Q)
    M = OMB * Q
    vbar = math.sqrt(8 * T / math.pi) / math.sqrt(1 + 8 * T / math.pi)
    zeta = epstein * 4 * math.pi / 3 * R * R * s['eP_vap'] * vbar
    D = T / zeta
    v2 = 8 * T / (math.pi * M)
    Rs = R[:, None] + R[None, :]
    Kc = 4 * math.pi * Rs * (D[:, None] + D[None, :])
    Kfm = math.pi * Rs ** 2 * np.sqrt(v2[:, None] + v2[None, :])
    return Kc * Kfm / (Kc + Kfm)


def merge_close(lq, N, nmax, grid=0.01):
    """Vectorised re-binning when more than nmax cohorts exist: cohorts are aggregated on a uniform ln Q grid (spacing
    chosen to leave <= nmax bins), each bin becoming one cohort at its charge-weighted mean charge (number and charge
    conserved exactly)."""
    if lq.size <= nmax:
        return lq, N
    span = lq.max() - lq.min()
    h = max(span / (nmax - 1), grid)
    b = np.floor((lq - lq.min()) / h).astype(np.int64)
    ub, inv = np.unique(b, return_inverse=True)
    Nn = np.zeros(ub.size); Wn = np.zeros(ub.size)
    Q = np.exp(lq)
    np.add.at(Nn, inv, N); np.add.at(Wn, inv, N * Q)
    return np.log(Wn / Nn), Nn


def stats(lq, N):
    Q = np.exp(lq)
    w = N * Q
    eta = w.sum()
    if eta <= 0:
        return dict(eta=0.0)
    meanQ = float((w * Q).sum() / eta)                     # charge-weighted <Q>
    Qeff = float((eta / (N * Q ** (2 / 3)).sum()) ** 3)      # capture-weighted: n_obj <R^2> = (eta/Qeff) Qeff^{2/3}
    mu = float((w * lq).sum() / eta)
    sd = float(math.sqrt(max((w * (lq - mu) ** 2).sum() / eta, 0.0)))
    EQ = float((w * EoverQ(Q)).sum() / eta)
    return dict(eta=float(eta), Q_charge_weighted=meanQ, Q_capture_eff=Qeff, Q_number_mean=float(eta / N.sum()),
                lnQ_charge_mean=mu, sigma_lnQ_charge=sd, n_obj_per_charge=float(N.sum()), EoverQ_charge_weighted=EQ)


def run(H, Q0, sig0=0.7, epstein=1.0, coag=True, nmax=300, T_start=None, la_stop=None, verbose=True, label='', every=25):
    t0 = time.time()
    la = H.la_of_T(T_start) if T_start else 0.0
    s = H.at(la)
    lq = np.linspace(math.log(Q0) - 4 * sig0, math.log(Q0) + 4 * sig0, 161)
    N = np.exp(-0.5 * ((lq - math.log(Q0)) / sig0) ** 2)
    N *= s['eta'] / (N * np.exp(lq)).sum()
    la_end = H.lna[-1] if la_stop is None else la_stop
    dla = 1e-3
    out_rows = []
    frozen = None
    nstep = 0
    nref = 0
    coag_on = coag
    while la < la_end - 1e-12:
        dla = min(dla * 1.25, 0.02, la_end - la)
        while True:
            s1 = H.at(la + dla)
            sm = H.at(la + 0.5 * dla)
            dt = dla / sm['H']
            Q = np.exp(lq)
            T = s1['T']
            R0 = radius(Q)
            def newQ(d, nsub=4):
                # integrate the RADIUS: dR/dt = (dQ/dt)/(4 pi R^2 n_B) is bounded and smooth (non-stiff), RK2 midpoint
                Rn = R0.copy(); h = dt / nsub
                for _ in range(nsub):
                    Rm = np.maximum(Rn + 0.5 * h * radial_velocity(Rn, d, s1), RMIN)
                    Rn = Rn + h * radial_velocity(Rm, d, s1)
                    Rn = np.maximum(Rn, 0.5 * RMIN)
                Qn = 4 * math.pi / 3 * NB * Rn ** 3
                return np.where(Qn < QMIN, 0.0, Qn)
            def g(d):
                eta_n = (N * newQ(d)).sum()
                vfrac = max(1 - eta_n * s1["nQ"] / NB, 0.0)        # vapour volume fraction (liquid excluded); g stays monotonic
                fv = vfrac * (s1['n_X'] * math.exp(d / T) - s1['n_Xbar'] * math.exp(-d / T)) / s1['nQ']
                return eta_n + fv - 1.0
            try:
                d = brentq(g, -30 * T, 30 * T, xtol=1e-18 * max(T, 1e-3), rtol=1e-15, maxiter=400)
            except ValueError:
                dla *= 0.5
                if dla < 1e-9:
                    print('no bracket: g(-30T)=%r g(+30T)=%r g(0)=%r la=%.6f T=%.4g ncoh=%d eta=%.6f' % (
                        g(-30 * T), g(30 * T), g(0.0), la, T, lq.size, (N * np.exp(lq)).sum()), flush=True)
                    raise
                continue
            Qn = newQ(d)
            keep = Qn > 0
            wq = N * Q
            dissolved = wq[~keep].sum() / max(wq.sum(), 1e-300)
            mshift = (wq[keep] * np.abs(np.log(Qn[keep] / Q[keep]))).sum() / max(wq[keep].sum(), 1e-300) if keep.any() else 0.0
            if (mshift > 0.02 or dissolved > 0.02) and dla > 1e-8:
                dla *= 0.5
                continue
            break
        lq_new, N_new = np.log(Qn[keep]), N[keep]
        # coagulation (explicit, switched off once its rate stays below 1e-4 H)
        coag_rate = 0.0
        if coag_on and lq_new.size > 0:
            Qc = np.exp(lq_new)
            K = coag_kernel(Qc, sm, epstein)
            rate_i = (K * N_new[None, :]).sum(1) * s1['nQ']
            coag_rate = float((rate_i * N_new * Qc).sum() / (N_new * Qc).sum() / sm['H'])
            if coag_rate < 1e-4 and s1['T'] < T_start * 0.9 if T_start else coag_rate < 1e-4:
                coag_on = False
            nsub = max(1, int(math.ceil(rate_i.max() * dt / 0.05)))
            for isub in range(nsub):
                if isub:
                    Qc = np.exp(lq_new)
                    K = coag_kernel(Qc, sm, epstein)
                n = Qc.size
                iu, ju = np.triu_indices(n)
                fac = np.where(iu == ju, 0.5, 1.0)
                dN = fac * K[iu, ju] * N_new[iu] * N_new[ju] * s1['nQ'] * dt / nsub
                loss = np.zeros(n)
                np.add.at(loss, iu, dN); np.add.at(loss, ju, dN)
                scale = np.minimum(1.0, 0.5 * N_new / np.maximum(loss, 1e-300))
                dN = dN * np.minimum(scale[iu], scale[ju])
                loss = np.zeros(n); np.add.at(loss, iu, dN); np.add.at(loss, ju, dN)
                N_new = N_new - loss
                sel = dN > 1e-12 * dN.max() if dN.size else dN > 0
                lq_new = np.r_[lq_new, np.log(Qc[iu][sel] + Qc[ju][sel])]; N_new = np.r_[N_new, dN[sel]]
                lq_new, N_new = merge_close(lq_new, N_new, nmax)
        lq, N = lq_new, N_new
        # re-resolve when fewer than 80 cohorts carry 99.9% of the charge
        if lq.size:
            wq = N * np.exp(lq); o = np.argsort(-wq); cum = np.cumsum(wq[o]) / wq.sum()
            if np.searchsorted(cum, 0.999) + 1 < 80:
                lq, N = refine(lq, N)
                nref = nref + 1
        la += dla
        nstep += 1
        st = stats(lq, N)
        fv = 1 - st.get('eta', 0.0)
        Qc = np.exp(lq)
        tf = trap_factor(radius(Qc), s1['T']) if (s1.get('trap') and lq.size) else 1.0
        Gc = (N * s1['nQ'] * 4 * math.pi * radius(Qc) ** 2 * tf).sum() * s1['Phi_X'] / s1['n_X'] if lq.size else 0.0
        row = dict(lna=la, T=s1['T'], T_vis_GeV=s1['T_vis_GeV'], H=s1['H'], d_over_T=d / T, f_v=fv, f_v_eq=1 - s1['eta'],
                   Gamma_capture_over_H=Gc / s1['H'], coag_rate_over_H=coag_rate, phi_vol=st.get('eta', 0) * s1['nQ'] / NB,
                   R_c_from_d=(2 * SIGMA / (NB * d)) if d > 0 else None, ncoh=int(lq.size), **st)
        if frozen is None and Gc / s1['H'] < 1.0 and s1['T'] < 0.1:
            frozen = dict(row)
            frozen['lq'] = lq.tolist(); frozen['N'] = N.tolist()
        if nstep % every == 0 or la >= la_end - 1e-12:
            out_rows.append(row)
            if verbose and nstep % 200 == 0:
                print('  %s step %d lna=%.3f T=%.4g Tvis=%.3g GeV eta=%.4f (eq %.4f) <Q>=%.3e Qeff=%.3e sig=%.2f Gc/H=%.2e coag/H=%.2e ncoh=%d t=%.0fs' % (
                    label, nstep, la, s1['T'], s1['T_vis_GeV'], st.get('eta', 0), s1['eta'], st.get('Q_charge_weighted', 0),
                    st.get('Q_capture_eff', 0), st.get('sigma_lnQ_charge', 0), Gc / s1['H'], coag_rate, lq.size, time.time() - t0), flush=True)
    # variant (b): decoupled vapour from the freeze-out point, T_v ~ a^-2, capture only, distribution frozen
    vb = None
    if frozen is not None:
        lqf, Nf = np.array(frozen['lq']), np.array(frozen['N'])
        la0, Tf = frozen['lna'], frozen['T']
        sf = H.at(la0)
        sflux = sf['Phi_X'] / sf['Phi_X_bare']                       # thermal-flux sticking at T_f
        Nv = frozen['f_v']
        stick = json.load(open(os.path.join(HERE, '..', 'results', 'vapour_gt.json')))['kinetic']
        Tk = np.array([r['T'] for r in stick]); sk = np.array([r['s_flux_avg'] for r in stick])
        x = np.linspace(la0, H.lna[-1], 4000)
        for a0, a1 in zip(x[:-1], x[1:]):
            s = H.at(0.5 * (a0 + a1))
            Tv = Tf * math.exp(-2 * (0.5 * (a0 + a1) - la0))
            sv = float(np.interp(math.log(Tv), np.log(Tk), sk)) if Tv > Tk[0] else sk[0] * math.sqrt(Tv / Tk[0])
            if H.stick_mode == 'geometric':
                sv = 1.0
            elif isinstance(H.stick_mode, tuple):
                sv *= H.stick_mode[1]
            vb_ = math.sqrt(8 * Tv / math.pi)
            tf = trap_factor(radius(np.exp(lqf)), Tv) if H.stick_mode == 'trap' else 1.0
            Gc = (Nf * s['nQ'] * math.pi * radius(np.exp(lqf)) ** 2 * tf).sum() * sv * vb_
            Nv *= math.exp(-Gc / s['H'] * (a1 - a0))
        if H.stick_mode == 'trap':
            wq = Nf * np.exp(lqf)
            sflux *= float((wq * trap_factor(radius(np.exp(lqf)), Tf)).sum() / wq.sum())
        vb = dict(f_v_final=Nv, f_v_at_freeze=frozen['f_v'], sticking_at_freeze=sflux)
    res = dict(label=label, S_over_Q=H.SQ, Q0=Q0, sigma0=sig0, epstein=epstein, coag=coag, T_start=T_start, steps=nstep, refines=nref,
               rows=out_rows, freeze=({k: v for k, v in frozen.items() if k not in ('lq', 'N')} if frozen else None),
               final=out_rows[-1] if out_rows else None, variant_b=vb, runtime_s=time.time() - t0,
               state=dict(lna=la, lq=lq.tolist(), N=N.tolist()))
    if frozen:
        lqf, Nf = np.array(frozen['lq']), np.array(frozen['N'])
        h, e = np.histogram(lqf / math.log(10), bins=40, weights=Nf * np.exp(lqf))
        res['freeze_charge_histogram_log10Q'] = dict(edges=e.tolist(), charge_fraction=h.tolist())
    return res


def wagner_test():
    """Pure ripening at fixed background (T = 0.5, no global drive): R_c^2 vs t against Wagner's R_c^2 = K_W t/2."""
    class Fixed:
        pass
    import qm_vapour as qv
    T = 0.5
    mc = 0.775905
    FX, dFX = qv.flux(mc, T, 'X', True, True); FB, dFB = qv.flux(-mc, T, 'Xbar', True, True)
    import bessel_vapour_copy as bv
    v = bv.vapour(mc, T)
    KW = 2 * SIGMA * (dFX + dFB) / NB ** 2
    nQ = 1.0                       # charge density chosen so the vapour holds ~4% of the charge
    s = dict(T=T, Phi_X=FX, Phi_Xbar=FB, n_X=v['nX'], n_Xbar=v['nXb'], nQ=nQ)
    Q0 = 1e12
    lq = np.linspace(math.log(Q0) - 2.8, math.log(Q0) + 2.8, 201)
    N = np.exp(-0.5 * ((lq - math.log(Q0)) / 0.7) ** 2)
    eta0 = 1 - v['n'] / nQ
    N *= eta0 / (N * np.exp(lq)).sum()
    t = 0.0
    dt = 0.002 * radius(Q0) ** 2 / KW
    rows = []
    for it in range(20000):
        Q = np.exp(lq)
        def g(d):
            Qn = Q + exchange_rate(Q, d, s) * dt
            Qn = np.where(Qn < QMIN, 0.0, Qn)
            return (N * Qn).sum() + (s['n_X'] * math.exp(d / T) - s['n_Xbar'] * math.exp(-d / T)) / nQ - 1
        d = brentq(g, -1e-3, 1e-3, xtol=1e-20, rtol=1e-15, maxiter=400)
        Qn = Q + exchange_rate(Q, d, s) * dt
        keep = Qn > QMIN
        lq, N = np.log(Qn[keep]), N[keep]
        wq = N * np.exp(lq); o = np.argsort(-wq); cum = np.cumsum(wq[o]) / wq.sum()
        if np.searchsorted(cum, 0.999) + 1 < 80:
            lq, N = refine(lq, N)
        t += dt
        dt *= 1.003
        if it % 500 == 0:
            Rc = 2 * SIGMA / (NB * (d - 0.0)) if d > 0 else float('nan')
            Rbar = float((N * radius(np.exp(lq))).sum() / N.sum())
            Qq = np.exp(lq); Rq = radius(Qq)
            rows.append(dict(t=t, Rc=Rc, Rbar_number=Rbar, Rc2_over_KWt=Rc ** 2 / (KW * t), Rbar2_over_KWt=Rbar ** 2 / (KW * t), n=int(lq.size),
                             Rmax_over_Rc=float(Rq.max() / Rc), Rbar_over_Rc=float(Rbar / Rc), Q_charge_weighted=float((N * Qq * Qq).sum() / (N * Qq).sum())))
    return dict(K_W=KW, rows=rows, wagner_prediction=dict(Rc2_over_KWt=0.5))


def T_spinodal(SQ):
    d = json.load(open(os.path.join(HERE, '..', '..', 'dq_entropy', 'results', 'summary.json')))['benchmark']['SQ_of_R_table']
    sq = np.array([r['S_per_Q_final_T0'] for r in d]); Ton = np.array([r['T_h_over_mX_onset'] for r in d])
    if SQ < sq[0]:
        Ton_v = Ton[0] + (SQ - sq[0]) * (Ton[1] - Ton[0]) / (sq[1] - sq[0])   # linear extrapolation below the R = 0 row
    else:
        Ton_v = float(np.interp(SQ, sq, Ton))
    # spinodal crossing at a/a_on = 0.84-0.87 (fx_chain): T_sp = T_on / 0.855
    return Ton_v / 0.855, Ton_v


def hydro_Q0(H, T_sp):
    """Drag-limited coarsening of the bicontinuous liquid network (volume fraction ~0.5 at the spinodal) until the
    volume fraction falls to 0.3: L^2 = 6 sigma t/((e+P)_vap vbar), Q0 = n_B (pi/6) L^3  [E]."""
    la0 = H.la_of_T(T_sp)
    s0 = H.at(la0)
    phi0 = s0['phi_vol']
    if phi0 <= 0.3:
        return None, dict(phi0=phi0)
    la1 = la0 + math.log(phi0 / 0.3) / 3
    x = np.linspace(la0, la1, 200)
    tb = sum((b - a) / H.at(0.5 * (a + b))['H'] for a, b in zip(x[:-1], x[1:]))
    vbar = math.sqrt(8 * s0['T'] / math.pi) / math.sqrt(1 + 8 * s0['T'] / math.pi)
    L = math.sqrt(6 * SIGMA * tb / (s0['eP_vap'] * vbar))
    return NB * math.pi / 6 * L ** 3, dict(phi0=phi0, dlna=la1 - la0, t=tb, L=L, T_end=H.at(la1)['T'])


def main(out, which):
    SCR = os.environ.get('QM_SCR', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/qm_micro')
    if which == 'wagner':
        r = wagner_test()
        json.dump(r, open(out, 'w'), indent=1, default=float)
        for row in r['rows']:
            print(row)
        return
    SQ = float(which)
    H = History(os.path.join(SCR, 'hist_%s.json' % which))
    Tsp, Ton = T_spinodal(SQ)
    Tstart = min(Tsp, H.tab['T'][0])
    res = dict(S_over_Q=SQ, T_spinodal=Tsp, T_onset=Ton, T_start=Tstart, runs=[])
    Qh, hinfo = hydro_Q0(H, Tstart)
    res['hydro_stage'] = dict(Q0=Qh, **hinfo)
    print('S/Q=%.2f T_sp=%.3f (onset %.3f) start %.3f; hydro Q0=%s %s' % (SQ, Tsp, Ton, Tstart, '%.3e' % Qh if Qh else None, hinfo), flush=True)
    cases = [('Q0_4e11', 4e11, 1.0, True, 'computed'), ('Q0_1.2e13', 1.2e13, 1.0, True, 'computed'),
             ('Q0_4e11_nocoag', 4e11, 1.0, False, 'computed'), ('Q0_4e11_epstein1.4', 4e11, 1.4, True, 'computed'),
             ('Q0_4e11_stick_x0.1', 4e11, 1.0, True, ('scale', 0.1)), ('Q0_4e11_stick_geometric', 4e11, 1.0, True, 'geometric'),
             ('Q0_4e11_trap_beliaev', 4e11, 1.0, True, 'trap')]
    if Qh:
        cases.append(('Q0_hydro', Qh, 1.0, True, 'computed'))
    only = os.environ.get('QM_CASES')
    if only:
        cases = [c for c in cases if c[0] in only.split(',')]
        prev = os.path.join(os.path.dirname(out), os.path.basename(out))
        if os.path.exists(out):
            old = json.load(open(out))
            res['runs'] = [r for r in old['runs'] if r['case'] not in only.split(',')]
    for lab, Q0, ep, cg, sm in cases:
        H.stick_mode = sm
        r = run(H, Q0, epstein=ep, coag=cg, T_start=Tstart, label='%s %s' % (which, lab), verbose=False)
        r['stick_mode'] = sm if isinstance(sm, str) else 'scale %g' % sm[1]
        r.pop('state', None)
        r['case'] = lab
        res['runs'].append(r)
        f = r['freeze']; fb = r['variant_b']
        print('%s %s: freeze T_f=%.4g (Tvis %.3g GeV) f_v=%.4f <Q>=%.3e Qeff=%.3e | final f_v=%.4f <Q>=%.3e | variant b f_v=%.4f | %.0fs' % (
            which, lab, f['T'], f['T_vis_GeV'], f['f_v'], f['Q_charge_weighted'], f['Q_capture_eff'], r['final']['f_v'],
            r['final']['Q_charge_weighted'], fb['f_v_final'], r['runtime_s']), flush=True)
        json.dump(res, open(out, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
