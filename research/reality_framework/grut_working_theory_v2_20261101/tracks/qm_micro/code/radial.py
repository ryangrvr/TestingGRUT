"""qm_micro item 1d (Claude, 2 Nov 2026): finite-object partial-wave cross-check.  Own code (method independent of
tracks/m1_transport/skeptic/code/sk_xcore.py: second-order finite differences, direct sparse solve with
Riccati-Hankel radiation conditions, instead of RK4 shooting with QR re-orthonormalisation).

Radial form of the planar operator of scatter1d.py: y = r (U, V, Sx), y'' = [M(r; w) + l(l+1)/r^2] y, y(0) = 0.
Objects: M1 Q-balls from fx_formation's solver (qball_m1_copy.py, imported read-only) at Q = 2e5, 1e6, 1e7, 1e8.
(1) Undamped (the physical linear problem): for an X below the Xbar and phi thresholds only the elastic channel is
    open, S_l = e^{2 i delta_l}; sigma_el = (pi/k^2) sum (2l+1) |1 - S_l|^2.  Compared with sk_xcore at Q = 2e5,
    k = 0.0684 (delta_0 = -0.356, delta_1 = -0.027, delta_2 = +0.003, sigma_el = 332/m_X^2).
(2) Interior absorber: a smooth damping -i gamma(r) w on the diagonal, gamma = g0 (1 - r/R_d)^2 for r < R_d (well
    inside the surface), which removes the transmitted phonons as a thermalising interior would.  Then
    sigma_abs = (pi/k^2) sum (2l+1)(1 - |S_l|^2) is compared with the geometric-optics value pi R^2 s_ang(E),
    s_ang = int_0^1 s(E, mu) 2 mu dmu from the planar tables (results/stick_thermal.json).
"""
import json, math, os, sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy.interpolate import CubicSpline
from scipy.special import spherical_jn, spherical_yn, spherical_kn
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import qm_common as qc
FXF = os.path.join(HERE, '..', '..', 'fx_formation', 'code')
sys.path.insert(0, FXF)
import qball_m1_copy as qb          # read-only import of fx_formation's copy of the m1_vacuum solver

a_ = qc.a_
NB = 449.31336788454803


def objects(Qs):
    out = {}
    # the continuation sequence of fx_formation/code/family.py
    seq = [(2e5, 60.), (3e5, 60.), (5e5, 60.), (1e6, 60.), (2e6, 70.), (5e6, 80.), (1e7, 90.), (2e7, 100.), (3e7, 110.),
           (5e7, 120.), (7e7, 130.), (1e8, 140.)]
    sv, u, w, Om, ob = qb.solve_fixed_Q(200000, 60.0, 0.025)
    cur = (sv, u, w, Om)
    for Q, L in seq:
        if Q > max(Qs):
            break
        if Q != 2e5:
            try:
                sv, u, w, Om, ob = qb.solve_fixed_Q(Q, L, 0.025, seed=cur)
            except RuntimeError as e:
                print('continuation stopped at Q=%.3g: %s' % (Q, e), flush=True)
                break
            cur = (sv, u, w, Om)
        if Q in Qs:
            r = sv.r
            out[Q] = dict(r=r, F=u / r, S=1 + w / r, Omega=Om, E_over_Q=ob['E_over_Qm'], r50=ob['r50'], r90=ob['r90'],
                          L=sv.L, R_gibbs=(3 * Q / (4 * math.pi * NB)) ** (1 / 3))
            print('object Q=%.3g Omega=%.6f E/Q=%.6f r90=%.2f R_gibbs=%.2f' % (Q, Om, ob['E_over_Qm'], ob['r90'], out[Q]['R_gibbs']), flush=True)
    return out


class Radial:
    def __init__(self, obj, rmax, h):
        self.h = h
        self.r = h * np.arange(1, int(round(rmax / h)) + 1)
        rr = np.r_[0.0, obj['r'], obj['L']]
        Fs = CubicSpline(rr, np.r_[obj['F'][0] + (obj['F'][0] - obj['F'][1]) * 0, obj['F'], 0.0])
        Ss = CubicSpline(rr, np.r_[obj['S'][0], obj['S'], 1.0])
        inside = self.r < obj['L']
        self.F = np.where(inside, Fs(np.minimum(self.r, obj['L'])), 0.0)
        self.S = np.where(inside, Ss(np.minimum(self.r, obj['L'])), 1.0)
        self.Om = obj['Omega']

    def solve(self, w, l, gamma=None):
        r, F, S, Om, h = self.r, self.F, self.S, self.Om, self.h
        n = r.size
        cen = l * (l + 1) / r ** 2
        damp = np.zeros(n, complex) if gamma is None else -1j * gamma * w
        d0 = S * S + 3 * a_ * F * F - Om * Om - w * w + cen + damp
        d1 = S * S + a_ * F * F - Om * Om - w * w + cen + damp
        d2 = qc.Uphi_pp(S) + F * F - w * w + cen + damp
        # vacuum channels at r_max: B = (1,-1,0)/sqrt2 (X, energy Om+w), A = (1,1,0)/sqrt2 (Xbar), phi
        mB, mA, mF = 1 - (Om + w) ** 2, 1 - (Om - w) ** 2, 6 - w * w
        E = np.array([[1, 1, 0], [-1, 1, 0], [0, 0, math.sqrt(2)]]) / math.sqrt(2)     # columns B, A, phi
        rN, rN1 = r[-1], r[-1] + h
        def ratio_out(m):
            if m < 0:
                q = math.sqrt(-m)
                hp = lambda x: x * (spherical_jn(l, x) + 1j * spherical_yn(l, x))
                return hp(q * rN1) / hp(q * rN), q
            kap = math.sqrt(m)
            kk = lambda x: x * spherical_kn(l, x)
            return kk(kap * rN1) / kk(kap * rN), None
        rhos, qs = zip(*(ratio_out(m) for m in (mB, mA, mF)))
        K = (E * np.array(rhos)[None, :]) @ E.T                 # y_{N+1} = K y_N (outgoing/decaying)
        qB = qs[0]
        hm = lambda x: x * (spherical_jn(l, x) - 1j * spherical_yn(l, x))
        hp = lambda x: x * (spherical_jn(l, x) + 1j * spherical_yn(l, x))
        # incoming X in channel B, amplitude 1:  y_B = hm + S hp
        inc_N, inc_N1 = hm(qB * rN), hm(qB * rN1)
        rhs_end = (inc_N1 - rhos[0] * inc_N) * E[:, 0]
        i2 = 1 / h ** 2
        base = 3 * np.arange(n)
        rows, cols, vals = [], [], []
        # equation j: (y_{j+1} - 2 y_j + y_{j-1})/h^2 - M_j y_j = 0
        diag = [d0, d1, d2]
        for a in range(3):
            rows.append(base + a); cols.append(base + a); vals.append(-2 * i2 - diag[a])
            rows.append(base[1:] + a); cols.append(base[:-1] + a); vals.append(np.full(n - 1, i2, complex))
            rows.append(base[:-1] + a); cols.append(base[1:] + a); vals.append(np.full(n - 1, i2, complex))
        off = [(0, 1, -2 * Om * w * np.ones(n)), (1, 0, -2 * Om * w * np.ones(n)), (0, 2, -2 * S * F), (2, 0, -2 * S * F)]
        for a, b, v in off:
            rows.append(base + a); cols.append(base + b); vals.append(v.astype(complex))
        # last row: y_{N+1} = K y_N + rhs_end  ->  add K/h^2 to the last diagonal block
        for a in range(3):
            for b in range(3):
                rows.append(np.array([base[-1] + a])); cols.append(np.array([base[-1] + b])); vals.append(np.array([K[a, b] * i2]))
        A = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(3 * n, 3 * n))
        rhs = np.zeros(3 * n, complex)
        rhs[3 * (n - 1):] = -rhs_end * i2
        y = splu(A).solve(rhs).reshape(n, 3)
        yc = E.T @ y[-1]
        S_l = (yc[0] - inc_N) / hp(qB * rN)
        other = {}
        if qs[1] is not None:
            other['Xbar'] = abs(yc[1] / hp(qs[1] * rN)) ** 2 * qs[1] / qB
        if qs[2] is not None:
            other['phi'] = abs(yc[2] / hp(qs[2] * rN)) ** 2 * qs[2] / qB
        return S_l, other, qB


def cross_sections(R, w, gamma=None, lmax=None, tol=1e-9):
    out = []
    l = 0
    sig_el = sig_abs = 0.0
    k = None
    while True:
        S, oth, k = R.solve(w, l, gamma)
        el = (math.pi / k ** 2) * (2 * l + 1) * abs(1 - S) ** 2
        ab = (math.pi / k ** 2) * (2 * l + 1) * (1 - abs(S) ** 2 - sum(oth.values()))
        sig_el += el; sig_abs += ab
        out.append(dict(l=l, S=[S.real, S.imag], absS=abs(S), delta=float(np.angle(S) / 2), other=oth))
        l += 1
        if lmax is not None and l > lmax:
            break
        if lmax is None and l > k * (R.Rg + 15) + 10 and (el + abs(ab)) < tol * max(sig_el + abs(sig_abs), 1e-300):
            break
        if l > 400:
            break
    return sig_el, sig_abs, out, k


def main(out):
    t0 = time.time()
    st = json.load(open(os.path.join(HERE, '..', 'results', 'stick_thermal.json')))
    Tk = np.array(st['T_kin']); s_ang = np.array(st['s_X_angle_avg'])
    s_ang_of = lambda T: float(np.exp(np.interp(math.log(T), np.log(Tk), np.log(s_ang))))
    objs = objects([2e5, 1e6, 1e7, 1e8])
    res = dict(objects={('%.0e' % Q): {k: v for k, v in o.items() if k not in ('r', 'F', 'S')} for Q, o in objs.items()},
               absorber_model='gamma = 0.2 * (1 - tanh((r - (R_gibbs - 5))/2))/2 on all three diagonals (-i gamma w)')
    # (1) undamped, Q = 2e5, k = 0.0684: compare with sk_xcore
    o = objs[2e5]
    rows = []
    for rmax, h in ((60.0, 0.01), (90.0, 0.01), (60.0, 0.005)):
        R = Radial(o, rmax, h); R.Rg = o['R_gibbs']
        k0 = 0.0684; E = math.sqrt(1 + k0 * k0); w = E - o['Omega']
        se, sa, pw, k = cross_sections(R, w, None, lmax=8)
        rows.append(dict(rmax=rmax, h=h, k=k, sigma_el=se, sigma_abs=sa, delta=[p['delta'] for p in pw[:4]], absS=[p['absS'] for p in pw[:4]]))
        print('Q=2e5 k=%.4f rmax=%g h=%g: sigma_el=%.2f sigma_abs=%.1e delta0..3=%s' % (k, rmax, h, se, sa, ['%.4f' % p['delta'] for p in pw[:4]]), flush=True)
    res['undamped_Q2e5_k0.0684'] = rows
    res['sk_xcore_reference'] = dict(delta=[-0.356, -0.027, 0.003], sigma_el=332.0)
    scan = []
    R = Radial(o, 60.0, 0.01); R.Rg = o['R_gibbs']
    for k0 in (0.02, 0.04, 0.0684, 0.1, 0.15, 0.2, 0.3, 0.5):
        E = math.sqrt(1 + k0 * k0)
        se, sa, pw, k = cross_sections(R, E - o['Omega'], None, lmax=10)
        scan.append(dict(k=k0, sigma_el=se, sigma_abs=sa))
        print('   undamped Q=2e5 k=%.3f sigma_el=%.1f (|1-|S|^2| max %.1e)' % (k0, se, max(abs(1 - p['absS'] ** 2) for p in pw)), flush=True)
    res['undamped_Q2e5_kscan'] = scan
    # (2) absorbing interior: flat damping gamma0 = 0.2 inside R_d = R - 5 (tanh edge, width 2)
    ab = []
    # model spheres of bulk Q-matter with the planar interface glued at R0 (an approximate background, O(1/R0))
    P = qc.planar_profile(Lb=60, Lv=60, h=0.01)
    for R0 in (100.0, 200.0):
        rr = np.arange(0.01, R0 + 60.0, 0.01)
        objs['glued_R%g' % R0] = dict(r=rr, F=np.interp(rr - R0, P['x'], P['F'], left=P['F_B'], right=0.0),
                                      S=np.interp(rr - R0, P['x'], P['S'], left=P['S_B'], right=1.0), Omega=P['Omega_b'],
                                      L=R0 + 60.0, R_gibbs=R0 + P['x_gibbs'], E_over_Q=None, r50=None, r90=None)
    for key in list(objs):
        o = objs[key]
        Rg = o['R_gibbs']
        rmax = max(60.0, Rg + 50.0)
        R = Radial(o, rmax, 0.01); R.Rg = Rg
        Rd = Rg - 5.0
        gam = 0.2 * 0.5 * (1 - np.tanh((R.r - Rd) / 2.0))
        Tlist = (0.003, 0.02, 0.1, 0.3, 1.0) if not isinstance(key, str) else (0.1, 0.3)
        Q = key
        for T in Tlist:
            E = 1 + T; k0 = math.sqrt(E * E - 1)
            if (k0 * (Rg + 30) > 400):
                continue
            se, sa, pw, k = cross_sections(R, E - o['Omega'], gam)
            geo = math.pi * Rg ** 2 * s_ang_of(T)
            ab.append(dict(Q=Q if not isinstance(Q, str) else None, model=Q if isinstance(Q, str) else 'Q-ball', T=T, k=k, R_gibbs=Rg, R_damp=Rd, sigma_abs=sa, sigma_el=se, geo_piR2_sang=geo,
                           ratio=sa / geo, s_ang=s_ang_of(T), lmax=len(pw) - 1, kR=k * Rg))
            print('   absorber %s T=%.3f kR=%.2f: sigma_abs=%.4g  pi R^2 s_ang=%.4g  ratio %.3f  sigma_el=%.4g (l<=%d)  t=%.0fs' % (
                Q if isinstance(Q, str) else '%.0e' % Q, T, k * Rg, sa, geo, sa / geo, se, len(pw) - 1, time.time() - t0), flush=True)
    res['absorber'] = ab
    res['runtime_s'] = time.time() - t0
    json.dump(res, open(out, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main(sys.argv[1])
