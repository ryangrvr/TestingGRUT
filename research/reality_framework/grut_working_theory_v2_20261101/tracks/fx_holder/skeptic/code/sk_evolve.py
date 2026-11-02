#!/usr/bin/env python3
"""fx_holder_skeptic (Claude; own code) -- independent evolution through reheating for the trilinear holder.

Different from the track's rptri.c in every layer:
  * integrator: scipy solve_ivp DOP853 (8th order) instead of a hand-written DP5(4);
  * background: PRESCRIBED Einstein-frame inflaton + radiation, H^2 = (J^2/2 + b^2 I^2/2 + rho_r)/3, R = 2 b^2 I^2 - J^2
    (the hidden energy, the inflaton back-reaction and the xi_X|X|^2 F-term are O(eta^2) ~ 3e-15 and dropped) instead of
    the constraint-solved Jordan-frame engine with exact back-reaction;
  * holder written WITHOUT g_I or eta^2:  G phihat(t) = G phihat(0) I(t)/I(0) = (10i - mu) * scale * I(t)/I(0), so the
    trilinear is  u(t) = mu + scale (10 i - mu) I(t)/I_ref  (I_ref = sqrt6/b_ref);  Option 1A u = mu + lock R;
  * own preparations (sk_common.mp_stationary / mp_held, 40-digit Newton), own window means on dt = 0.05 samples;
  * optional phase-blind counterterm cport (|u|^2 - mu^2)(|S|^2 + 2|Y|^2): the operators the track's one-loop running
    induces (phi^2|S|^2, phi|S|^2, phi^2|Y|^2, phi|Y|^2 with the fixed ratios of Tr M^4), cport = 1/(4 pi^2) per unit ln Q.
Usage: sk_evolve.py TAG MODE LAM B GAM SCALE CPORT END [BREF]  (MODE in orig, 1A, tri; SCALE multiplies the holder term
       (lock for 1A); BREF = inflaton mass that fixed g_I; writes ../results/runs/TAG.json)
"""
import json, sys, time
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sk_common as sc

RES = Path(__file__).resolve().parent.parent / 'results' / 'runs'
RES.mkdir(parents=True, exist_ok=True)


def make_rhs(mode, p, b, gam, scale, Iref):
    lx, ly, ls, kap, MM2, h2 = p['lx'], p['ly'], p['ls'], p['kap'], p['M'] ** 2, p['h'] ** 2 / 2
    mu, cp = p['mu'], p['cport']
    w = (10j - mu) * scale / Iref
    lock = sc.LOCK * (scale if mode == '1A' else 1.0)
    b2 = b * b

    def rhs(t, y):
        I, J, rad = y[0], y[1], y[2]
        X = complex(y[4], y[5]); Y = complex(y[6], y[7]); S = complex(y[8], y[9])
        vX = complex(y[10], y[11]); vY = complex(y[12], y[13]); vS = complex(y[14], y[15])
        H = ((0.5 * J * J + 0.5 * b2 * I * I + rad) / 3.0) ** 0.5
        R = 2 * b2 * I * I - J * J
        if mode == 'tri':
            u = mu + w * I
        elif mode == '1A':
            u = mu + lock * R
        else:
            u = mu + 0j
        du = cp * (u.real * u.real + u.imag * u.imag - mu * mu)
        x2 = X.real * X.real + X.imag * X.imag; y2 = Y.real * Y.real + Y.imag * Y.imag; s2 = S.real * S.real + S.imag * S.imag
        Xc = X.conjugate()
        gX = (2 * lx * x2 + sc.XIX * R) * X + 3 * kap * Y * Xc * Xc
        gY = (MM2 + 2 * ly * y2 + 2 * du) * Y + kap * X * X * X + 2 * u.conjugate() * S * Y.conjugate()
        gS = (2 * ls * (s2 - h2) + du) * S + u * Y * Y
        h3 = 3 * H
        aX = -h3 * vX - gX; aY = -h3 * vY - gY; aS = -h3 * vS - gS
        return [J, -(h3 + gam) * J - b2 * I, gam * J * J - 4 * H * rad, H,
                vX.real, vX.imag, vY.real, vY.imag, vS.real, vS.imag,
                aX.real, aX.imag, aY.real, aY.imag, aS.real, aS.imag]
    return rhs


def u_of(mode, p, b, scale, Iref, I, J):
    R = 2 * b * b * I * I - J * J
    if mode == 'tri':
        return p['mu'] + (10j - p['mu']) * scale * I / Iref
    if mode == '1A':
        return p['mu'] + sc.LOCK * scale * R
    return p['mu'] + 0j


def prep(mode, p, b, scale, Iref):
    I0 = np.sqrt(6.0) / b
    u0 = u_of(mode, p, b, scale, Iref, I0, 0.0)
    if mode == 'orig':
        z = sc.mp_held(p['mu'], p, -10.0, [76.4, -2.1, -0.95, 28.3])
        return np.array([complex(c) for c in z]), u0, 0.0
    # continuation in the cubic coefficient from the known lambda-dependent u = 10 i point to u0 (the branch is smooth and
    # unique over |u| in [5.8, 45], sk_point.py s-scan)
    guess = np.array([87.2, -5.75, -28.4j]) if p['ly'] > 0.75 else np.array([90.8, -6.9, -28.6j])
    pp = dict(p); pp['cport'] = 0.0
    z, _ = sc.mp_stationary(10j, pp, -10.0, guess)
    n = 40
    for k in range(1, n + 1):
        f = k / n
        uk = 10j + f * (u0 - 10j)
        pk = dict(p); pk['cport'] = p['cport'] * f
        z, res = sc.mp_stationary(uk, pk, -10.0, np.array([complex(c) for c in z]))
    return np.array([complex(c) for c in z]), u0, float(res)


def window(t, f, a, b_):
    m = (t >= a - 1e-9) & (t <= b_ + 1e-9)
    return float(np.trapezoid(f[m], t[m]) / (t[m][-1] - t[m][0]))


def main():
    tag, mode = sys.argv[1], sys.argv[2]
    lam, b, gam, scale, cport, end = map(float, sys.argv[3:9])
    bref = float(sys.argv[9]) if len(sys.argv) > 9 else b       # inflaton history that fixed g_I (default: the same)
    p = sc.par(lam, cport)
    Iref = np.sqrt(6.0) / bref
    t0 = time.time()
    z0, u0, pres = prep(mode, p, b, scale, Iref)
    y0 = np.zeros(16)
    y0[0] = np.sqrt(6.0) / b; y0[3] = 0.0
    y0[4:10] = np.c_[z0.real, z0.imag].ravel()
    dt = 0.05
    te = np.arange(0.0, end + dt / 2, dt)
    sol = solve_ivp(make_rhs(mode, p, b, gam, scale, Iref), (0.0, end), y0, method='DOP853', rtol=1e-10, atol=1e-12, t_eval=te)
    if not sol.success:
        raise RuntimeError(sol.message)
    t, Y = sol.t, sol.y
    a = np.exp(Y[3]); a3 = a ** 3
    N = []
    for k in range(3):
        z = Y[4 + 2 * k] + 1j * Y[5 + 2 * k]; v = Y[10 + 2 * k] + 1j * Y[11 + 2 * k]
        N.append(-2 * a3 * (np.conj(z) * v).imag)
    N = np.array(N)
    gauss = np.max(np.abs(N[0] + 3 * N[1] + 6 * N[2]) / np.maximum(np.abs(N[0]) + 3 * np.abs(N[1]) + 6 * np.abs(N[2]), 1.0))
    hi = min(end, 1000.0)
    out = dict(tag=tag, mode=mode, lam=lam, b=b, gam=gam, scale=scale, cport=cport, end=end, bref=bref, u0=[u0.real, u0.imag],
               prep_amps=np.abs(z0).tolist(), prep_delta=float(np.angle(z0[2] * np.exp(-6j * np.angle(z0[0])))), prep_residual=pres,
               nfev=int(sol.nfev), seconds=None, max_gauss=float(gauss))
    out['NX_400_1000'] = window(t, N[0], 400, hi)
    out['NX_640_800'] = window(t, N[0], 640, 800)
    if end >= 1000:
        out['NX_800_1000'] = window(t, N[0], 800, 1000)
    starts = np.arange(400.0, hi - 80 + 1e-9, 20.0)
    r80 = np.array([window(t, N[0], s, s + 80) for s in starts])
    out['NX_run80'] = [float(r80.min()), float(r80.max())]
    out['run80'] = r80.tolist()
    # exchange partner
    n6 = np.array([window(t, 6 * N[2], s, s + 80) for s in starts]); n3 = np.array([window(t, 3 * N[1], s, s + 80) for s in starts])
    out['corr_NX_6NS'] = float(np.corrcoef(r80, n6)[0, 1]); out['corr_NX_3NY'] = float(np.corrcoef(r80, n3)[0, 1]) if np.std(n3) > 0 else None
    # late Y excitation proxies: std of N_Y and a free-quantum number with the forced part removed (own derivation:
    # M^2 Y + c conj(Y) = -kap X^3 with c = 2 conj(u) S  =>  Y_p = (-M^2 J + c conj(J))/(M^4 - |c|^2))
    m = t >= 400
    X = Y[4] + 1j * Y[5]; Yc = Y[6] + 1j * Y[7]; S = Y[8] + 1j * Y[9]; vY = Y[12] + 1j * Y[13]
    uu = np.array([u_of(mode, p, b, scale, Iref, I, J) for I, J in zip(Y[0], Y[1])])
    c = 2 * np.conj(uu) * S; Jx = p['kap'] * X ** 3
    Yp = (-p['M'] ** 2 * Jx + c * np.conj(Jx)) / (p['M'] ** 4 - np.abs(c) ** 2)
    mL2 = p['M'] ** 2 - np.abs(c); mH2 = p['M'] ** 2 + np.abs(c)
    # light/heavy directions: V2 = M^2|w|^2 + 2|u S| Re(w^2) with Y - Yp = e^{-i beta} w, beta = -arg(u conj(S))/2 + ... :
    rot = np.exp(0.5j * np.angle(uu * np.conj(S)))
    dy = (Yc - Yp) * rot; dv = vY * rot                      # (forced-part velocity neglected late: X slow)
    NYL = a3 * (dv.imag ** 2 + mL2 * dy.imag ** 2) / np.sqrt(mL2)
    NYH = a3 * (dv.real ** 2 + mH2 * dy.real ** 2) / np.sqrt(mH2)
    out['NYL_400_1000'] = window(t[m], NYL[m], 400, hi); out['NYH_400_1000'] = window(t[m], NYH[m], 400, hi)
    mm = t >= 0.9 * hi
    out['NY_std_late'] = float(N[1][mm].std())
    # yield at the M1 H_i
    rad = Y[2][-1]; aend = a[-1]
    rho_phys = rad * sc.MR_GEV ** 2 * sc.H_GEV ** 2
    T = (30 * rho_phys / (np.pi ** 2 * sc.GSTAR)) ** 0.25
    s_ent = 2 * np.pi ** 2 / 45 * sc.GSTAR * T ** 3
    for w in ('NX_400_1000', 'NX_640_800'):
        Yx = abs(out[w]) * sc.H_GEV ** 3 / (aend ** 3 * s_ent)
        out['Y_over_YQ_' + w[3:]] = float(Yx / sc.YQ)
    out['Rr_end'] = float(rad * aend ** 4)
    out['H_new_GeV'] = float(sc.H_GEV * out['Y_over_YQ_400_1000'] ** (-2 / 3))
    keep = [0.4, 0.8, 1.2, 2, 3, 5, 10, 20, 50, 100]
    out['NX_frac_at'] = {str(k): float(np.interp(k, t, N[0]) / out['NX_400_1000']) for k in keep}
    early = t <= 6
    out['u_min_im_early'] = float(uu[early].imag.min()); out['mL2_min_early'] = float(mL2[early].min())
    out['seconds'] = time.time() - t0
    json.dump(out, open(RES / (tag + '.json'), 'w'), indent=1, default=float)
    print(tag, 'NX(400-1000)=%.6f NX(640-800)=%.6f run80=%s Y/YQ=%.6f H=%.5e NYL=%.4f NYstd=%.4g gauss=%.1e corr6NS=%.6f t=%.0fs' % (
        out['NX_400_1000'], out['NX_640_800'], np.round(out['NX_run80'], 3), out['Y_over_YQ_400_1000'], out['H_new_GeV'], out['NYL_400_1000'],
        out['NY_std_late'], gauss, out['corr_NX_6NS'], out['seconds']))


if __name__ == '__main__':
    main()
