"""qm_micro SKEPTIC (Claude, 2 Nov 2026): independent recomputation of the planar sticking of X on bulk M1 Q-matter.
Own code, written from the field equations, not from scatter1d.py / qm_common.py.

Method (different from qm_micro's Numerov block-tridiagonal solve):
  * bulk state by nested root finding of P(Omega) = max_S P(S; Omega) = 0 (brentq in S, brentq in Omega);
  * planar interface by scipy collocation (solve_bvp) with Omega as a free parameter and a pinning condition
    F(0) = F_B/2, Dirichlet ends at +-L;
  * linear scattering by two-sided shooting with DOP853 (rtol 1e-12) on the complex 3-channel ODE Y'' = M(x) Y,
    the left (bulk) solutions started as exact bulk eigenmodes (decaying into the bulk or outgoing phonon), the right
    solutions as exact vacuum eigenmodes (decaying closed channels, outgoing open channels, one incoming X), matched
    in value and derivative at x_m;
  * outgoing directions from the numerically differentiated bulk dispersion (group velocity), not from a norm formula.
Units m_X = 1; rescaled fields psi = sqrt2 g X/m_X, S = phi/v; U = Uphi + S^2|psi|^2/2 + a|psi|^4/4, a = 0.1, b = 3, c = 0.2.
Linearisation: psi = e^{-i Om t}(F + u + i v), S -> S + s, (u, v, s) = (U, iV, Sx) e^{-iwt} + c.c.:
  U''  = (S^2 + 3aF^2 - Om^2 - w^2 + kp^2) U + 2 Om w V + 2 S F Sx
  V''  = 2 Om w U + (S^2 + aF^2 - Om^2 - w^2 + kp^2) V
  Sx'' = 2 S F U + (Uphi''(S) + F^2 - w^2 + kp^2) Sx          (derived independently in SKEPTIC.md)
"""
import json, math, sys, time
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_bvp, solve_ivp

a, b, c = 0.1, 3.0, 0.2


def Uphi(S):
    return b / 4 * (S * S - 1) ** 2 - c * (S - 1) ** 3


def dUphi(S):
    return b * S * (S * S - 1) - 3 * c * (S - 1) ** 2


def d2Uphi(S):
    return b * (3 * S * S - 1) - 6 * c * (S - 1)


def P_of(Om):
    # stationary S of P(S) = (Om^2 - S^2)^2/(4a) - Uphi(S) in the condensed branch (S small)
    g = lambda S: -S * (Om * Om - S * S) / a - dUphi(S)
    S = brentq(g, 1e-3, 0.3, xtol=1e-16)
    return (Om * Om - S * S) ** 2 / (4 * a) - Uphi(S), S


def bulk():
    Om = brentq(lambda o: P_of(o)[0], 0.70, 0.85, xtol=1e-16)
    S = P_of(Om)[1]
    F = math.sqrt((Om * Om - S * S) / a)
    return dict(Omega=Om, S=S, F=F, n=100 * Om * F * F)


def profile(L=40.0, nodes=3001, tol=1e-10):
    """Planar interface by collocation on two half-lines glued at x = 0 (pin F(0) = F_B/2 there), Omega free,
    exact linear asymptotic (Robin) conditions at x = -L (bulk stable subspace) and x = +L (vacuum)."""
    B = bulk()
    FB, SB, OmB = B['F'], B['S'], B['Omega']

    def f1(F, Fp, S, Sp, Om):
        return np.vstack([Fp, (S * S + a * F * F - Om * Om) * F, Sp, dUphi(S) + F * F * S])

    def f(t, y, p):
        Om = p[0]
        yl, yr = y[:4], y[4:]
        return np.vstack([-L * f1(*yl, Om), L * f1(*yr, Om)])

    def bc(ya, yb, p):
        Om = p[0]
        yl0, yr0, yl1, yr1 = ya[:4], ya[4:], yb[:4], yb[4:]
        # bulk Jacobian at the current Omega's bulk point is close to the exact one; use the exact bulk (Omega_b)
        K = np.array([[SB * SB + 3 * a * FB * FB - Om * Om, 2 * SB * FB], [2 * SB * FB, d2Uphi(SB) + FB * FB]])
        lam, E = np.linalg.eigh(K)
        d = np.array([yl1[0] - FB, yl1[2] - SB]); dp = np.array([yl1[1], yl1[3]])
        cl = [E[:, i] @ (dp - math.sqrt(lam[i]) * d) for i in range(2)]      # decaying into the bulk: d' = +kappa d
        kF = math.sqrt(1 - Om * Om); kS = math.sqrt(d2Uphi(1.0))
        return np.array([*(yl0 - yr0), yl0[0] - 0.5 * FB, cl[0], cl[1], yr1[1] + kF * yr1[0], yr1[3] + kS * (yr1[2] - 1.0)])

    t = np.linspace(0, 1, nodes)
    xl, xr = -L * t, L * t
    pr = lambda x: 0.5 * (1 - np.tanh(x / 1.5))
    dpr = lambda x: -0.5 / 1.5 / np.cosh(x / 1.5) ** 2
    y0 = np.vstack([FB * pr(xl), FB * dpr(xl), 1 - (1 - SB) * pr(xl), -(1 - SB) * dpr(xl),
                    FB * pr(xr), FB * dpr(xr), 1 - (1 - SB) * pr(xr), -(1 - SB) * dpr(xr)])
    sol = solve_bvp(f, bc, t, y0, p=[OmB], tol=tol, max_nodes=500000)
    if not sol.success:
        raise RuntimeError(sol.message)
    Om = sol.p[0]

    class Prof:
        def __init__(self, s):
            self.s = s
        def sol(self, x):
            x = np.asarray(x, float)
            tt = np.abs(x) / L
            yy = self.s.sol(np.clip(tt, 0, 1))
            if np.ndim(x) == 0:
                return yy[:4] if x < 0 else yy[4:]
            out = np.where(x < 0, yy[:4], yy[4:])
            return out
    pf = Prof(sol)
    xx = np.linspace(-L, L, 200001)
    yy = pf.sol(xx)
    sig = 100 * np.trapezoid(yy[1] ** 2 + yy[3] ** 2, xx)
    Veff = Uphi(yy[2]) + 0.5 * yy[2] ** 2 * yy[0] ** 2 + a / 4 * yy[0] ** 4 - 0.5 * Om * Om * yy[0] ** 2
    sig_dir = 100 * np.trapezoid(0.5 * yy[1] ** 2 + 0.5 * yy[3] ** 2 + Veff, xx)
    return dict(sol=pf, Omega=Om, Omega_b=OmB, bulk=B, sigma=sig, sigma_direct=sig_dir, L=L, x_pin=0.0, nodes=sol.x.size,
                rms_res=float(np.max(sol.rms_residuals)))


class Scatter:
    """Linear scattering as a complex linear two-point BVP (scipy collocation, real/imaginary split), with exact
    radiation (Robin) conditions Y' = Z Y built from the asymptotic eigenmodes at both ends.  No shooting, so no
    loss of precision from exponentially growing closed channels."""

    def __init__(self, P, Lb=36.0, Lv=36.0):
        self.P = P
        self.sol = P['sol']   # object with .sol(x) -> (F, F', S, S')
        self.Om = P['Omega_b']
        self.B = P['bulk']
        self.Lb, self.Lv = Lb, Lv
        xs = np.linspace(-Lb, Lv, 14001)
        self.xs = xs
        y = self.sol.sol(xs)
        self.Fs, self.Ss = y[0], y[2]
        from scipy.interpolate import CubicSpline
        self.Fi = CubicSpline(xs, self.Fs); self.Si = CubicSpline(xs, self.Ss)

    def Mat(self, F, S, w, kp):
        Om = self.Om
        M = np.zeros((3, 3))
        M[0, 0] = S * S + 3 * a * F * F - Om * Om - w * w + kp * kp
        M[1, 1] = S * S + a * F * F - Om * Om - w * w + kp * kp
        M[2, 2] = d2Uphi(S) + F * F - w * w + kp * kp
        M[0, 1] = M[1, 0] = 2 * Om * w
        M[0, 2] = M[2, 0] = 2 * S * F
        return M

    def modes(self, F, S, w, kp):
        """eigenmodes of the constant matrix; for open ones the group velocity of e^{+ikx} from the numerically
        differentiated dispersion lambda_i(w): v_g = -2k / (dlambda/dw)."""
        M = self.Mat(F, S, w, kp)
        lam, E = np.linalg.eigh(M)
        dw = 1e-6
        lp, Ep = np.linalg.eigh(self.Mat(F, S, w + dw, kp))
        lm, Em = np.linalg.eigh(self.Mat(F, S, w - dw, kp))
        out = []
        for i in range(3):
            e = E[:, i]
            ip = int(np.argmax(np.abs(Ep.T @ e))); im = int(np.argmax(np.abs(Em.T @ e)))
            dl = (lp[ip] - lm[im]) / (2 * dw)
            if lam[i] > 0:
                out.append(dict(vec=e, kind='closed', kap=math.sqrt(lam[i])))
            else:
                k = math.sqrt(-lam[i])
                out.append(dict(vec=e, kind='open', k=k, vg_plus=-2 * k / dl))
        return out

    def solve(self, w, kp, inc='X', tol=1e-9):
        Lb, Lv = self.Lb, self.Lv
        mL = self.modes(self.B['F'], self.B['S'], w, kp)
        mR = self.modes(0.0, 1.0, w, kp)
        qL = []
        for m in mL:
            qL.append(m['kap'] if m['kind'] == 'closed' else (1j * m['k'] if m['vg_plus'] < 0 else -1j * m['k']))
        EL = np.array([m['vec'] for m in mL]).T
        ZL = (EL * np.array(qL)[None, :]) @ EL.T
        chX = np.array([1.0, -1.0, 0.0]) / math.sqrt(2); chB = np.array([1.0, 1.0, 0.0]) / math.sqrt(2)
        qR, names = [], []
        iinc = None
        for j, m in enumerate(mR):
            v = m['vec']
            nm = 'X' if abs(abs(v @ chX) - 1) < 1e-10 else ('Xbar' if abs(abs(v @ chB) - 1) < 1e-10 else 'phi')
            names.append(nm)
            if m['kind'] == 'closed':
                qR.append(-m['kap'])
            else:
                qR.append(1j * m['k'] if m['vg_plus'] > 0 else -1j * m['k'])
                if nm == inc:
                    iinc = j
        if iinc is None:
            return dict(closed=True)
        ER = np.array([m['vec'] for m in mR]).T
        ZR = (ER * np.array(qR)[None, :]) @ ER.T
        ein = ER[:, iinc]; qin = -qR[iinc]
        g = (qin - qR[iinc]) * ein        # Y' - ZR Y at x = Lv for the incoming wave of unit amplitude there

        Fi, Si = self.Fi, self.Si

        def f(x, y):
            F = Fi(x); S = Si(x)
            Om = self.Om
            m00 = S * S + 3 * a * F * F - Om * Om - w * w + kp * kp
            m11 = S * S + a * F * F - Om * Om - w * w + kp * kp
            m22 = d2Uphi(S) + F * F - w * w + kp * kp
            m01 = 2 * Om * w; m02 = 2 * S * F
            # y = [Yr(3), Yi(3), Ypr(3), Ypi(3)]
            Yr, Yi, Ypr, Ypi = y[0:3], y[3:6], y[6:9], y[9:12]
            def Mv(Y):
                return np.vstack([m00 * Y[0] + m01 * Y[1] + m02 * Y[2], m01 * Y[0] + m11 * Y[1], m02 * Y[0] + m22 * Y[2]])
            return np.vstack([Ypr, Ypi, Mv(Yr), Mv(Yi)])

        def bc(ya, yb):
            Ya = ya[0:3] + 1j * ya[3:6]; Ypa = ya[6:9] + 1j * ya[9:12]
            Yb = yb[0:3] + 1j * yb[3:6]; Ypb = yb[6:9] + 1j * yb[9:12]
            r1 = Ypa - ZL @ Ya
            r2 = Ypb - ZR @ Yb - g
            return np.concatenate([r1.real, r1.imag, r2.real, r2.imag])

        x = np.linspace(-Lb, Lv, 3001)
        y0 = np.zeros((12, x.size))
        sol = solve_bvp(f, bc, x, y0, tol=tol, max_nodes=300000)
        if not sol.success:
            raise RuntimeError(sol.message)
        yb = sol.y[:, -1]; ya = sol.y[:, 0]
        Yb = yb[0:3] + 1j * yb[3:6]
        Ya = ya[0:3] + 1j * ya[3:6]
        bR = ER.T @ (Yb - ein)
        Jin = abs(qin.imag)
        out = {}
        for j, m in enumerate(mR):
            if m['kind'] == 'open':
                out['out_' + names[j]] = abs(bR[j]) ** 2 * abs(qR[j].imag) / Jin
        cL = EL.T @ Ya
        tl = []
        for j, m in enumerate(mL):
            if m['kind'] == 'open':
                tl.append(dict(k=m['k'], frac=abs(cL[j]) ** 2 * m['k'] / Jin))
        trans = sum(t['frac'] for t in tl)
        return dict(reflected=out.get('out_' + inc, 0.0), out=out, trans=trans, trans_list=tl,
                    unitarity=sum(out.values()) + trans - 1, nodes=sol.x.size)


def bulk_dispersion(B, ks):
    """w(k) of the three bulk branches at kp = 0, k_x = k: roots of det(M(w) + k^2) = 0 in w (scan + brentq)."""
    Om, F, S = B['Omega'], B['F'], B['S']
    def det(w, k):
        M = np.array([[S * S + 3 * a * F * F - Om * Om - w * w + k * k, 2 * Om * w, 2 * S * F],
                      [2 * Om * w, S * S + a * F * F - Om * Om - w * w + k * k, 0.0],
                      [2 * S * F, 0.0, d2Uphi(S) + F * F - w * w + k * k]])
        return np.linalg.det(M)
    rows = []
    for k in ks:
        ws = np.linspace(1e-6, 4.0, 4001)
        dv = np.array([det(w, k) for w in ws])
        roots = []
        for i in range(ws.size - 1):
            if dv[i] == 0 or dv[i] * dv[i + 1] < 0:
                roots.append(brentq(det, ws[i], ws[i + 1], args=(k,), xtol=1e-15))
        rows.append(dict(k=float(k), w=roots))
    return rows


def main(out):
    t0 = time.time()
    P = profile()
    B = P['bulk']
    print('bulk', B, 'profile Omega-Omega_b = %.2e sigma = %.5f nodes %d pin %.1e rms %.1e (%.0fs)' % (
        P['Omega'] - P['Omega_b'], P['sigma'], P['nodes'], P['x_pin'], P['rms_res'], time.time() - t0), flush=True)
    sc = Scatter(P)
    res = dict(bulk=B, profile=dict(Omega=P['Omega'], Omega_minus_Omega_b=P['Omega'] - P['Omega_b'], sigma=P['sigma'], nodes=P['nodes'], L=P['L']))
    rows = []
    for T in (0.005, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0):
        E = 1 + T; p = math.sqrt(E * E - 1)
        for th in (0.0, 60.0, 85.0):
            kp = p * math.sin(math.radians(th))
            r = sc.solve(E - sc.Om, kp, 'X')
            rX = dict(T=T, theta=th, s=1 - r['reflected'], unitarity=r['unitarity'], out=r['out'], trans=r['trans_list'])
            if th == 0.0:
                r2 = Scatter(P, Lb=28.0, Lv=30.0).solve(E - sc.Om, kp, 'X', tol=1e-10)
                rX['s_ctrl_domain_xm'] = 1 - r2['reflected']
                q = sc.solve(E + sc.Om, kp, 'Xbar')
                rX['Xbar_destroyed'] = 1 - q['reflected']; rX['Xbar_out'] = q['out']; rX['Xbar_unitarity'] = q['unitarity']
            rows.append(rX)
            print('T=%.3f th=%4.0f s=%.6f uni=%.1e %s' % (T, th, rX['s'], rX['unitarity'],
                  ('ctrl %.6f Xbar destroyed %.6f' % (rX['s_ctrl_domain_xm'], rX['Xbar_destroyed'])) if th == 0 else ''), flush=True)
    res['X'] = rows
    # threshold coefficient
    thr = []
    for T in (1e-5, 1e-4, 1e-3):
        r = sc.solve(1 + T - sc.Om, 0.0, 'X')
        thr.append(dict(T=T, s=1 - r['reflected'], s_over_sqrtT=(1 - r['reflected']) / math.sqrt(T)))
    res['threshold'] = thr
    print('threshold', thr, flush=True)
    ks = np.r_[np.linspace(0.01, 0.2, 5), np.linspace(0.25, 2.0, 36)]
    disp = bulk_dispersion(B, ks)
    # phonon branch convexity
    ph = np.array([[r['k'], r['w'][0]] for r in disp])
    cs_small = ph[0, 1] / ph[0, 0]
    curv = np.gradient(np.gradient(ph[:, 1], ph[:, 0]), ph[:, 0])
    res['dispersion'] = dict(rows=disp, c_s_smallk=cs_small, phonon_curvature=curv.tolist(),
                             gaps_k0=bulk_dispersion(B, [1e-6])[0]['w'])
    print('c_s(k->0) = %.6f gaps %s; phonon w/k at k = 0.4, 0.6, 1.0: %s' % (cs_small, res['dispersion']['gaps_k0'],
          [float(np.interp(k, ph[:, 0], ph[:, 1]) / k) for k in (0.4, 0.6, 1.0)]), flush=True)
    res['runtime_s'] = time.time() - t0
    json.dump(res, open(out, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main(sys.argv[1])
