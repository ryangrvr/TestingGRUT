"""Skeptic part B: Floquet chart of phi_k on the frozen pure-quartic X orbit, by a different route from pt_matching.

[D] For H = |pi|^2/2 + lam |psi|^4/4, u = |psi|^2 obeys u'' = 4E - 3 lam u^2 (any charge q).  Hence u = -(2/lam) P(s+s0)
with P the Weierstrass function (g2 = 4 lam E), and the linearised phi_k equation
    phi'' = -(g^2 u - lam_phi M^2 + k^2) phi = [ (2 g^2/lam) P + B ] phi,   B = lam_phi M^2 - k^2
is a Lame equation with n(n+1) = 2 g^2/lam_X = 20, i.e. n = 4: a finite-gap potential with exactly four instability
gaps (some possibly closed) below the semi-infinite unstable region.  k = 0 at mass M and mode k at mass M' with
lam_phi M'^2 - k^2 = lam_phi M^2 share the same multiplier.

Numerics: u(s) from the scalar ODE u'' = 4E - 3 lam u^2 (not the 2-D orbit), with the 2x2 fundamental matrix of phi
for a vector of B values integrated together by scipy DOP853 over one radial period T_r."""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sk_common import *
from scipy.integrate import solve_ivp

t0 = time.time()


def monodromy_traces(u0, up0, E, Tr, Bs, rtol=1e-12):
    nB = len(Bs)
    def rhs(t, y):
        u, up = y[0], y[1]
        F = y[2:].reshape(4, nB)          # rows: phi_a, phi_a', phi_b, phi_b'
        c = G2 * u - Bs
        d = np.empty_like(F)
        d[0] = F[1]; d[1] = -c * F[0]; d[2] = F[3]; d[3] = -c * F[2]
        return np.r_[up, 4 * E - 3 * LX * u * u, d.ravel()]
    F0 = np.zeros((4, nB)); F0[0] = 1; F0[3] = 1
    sol = solve_ivp(rhs, (0, Tr), np.r_[u0, up0, F0.ravel()], method='DOP853', rtol=rtol, atol=1e-12)
    y = sol.y[:, -1]
    F = y[2:].reshape(4, nB)
    tr = F[0] + F[3]; det = F[0] * F[3] - F[2] * F[1]
    return tr, det, abs(y[0] - u0) / u0


def chart(y4, info, bmin, bmax, nb, chunk=400):
    E = energy4(y4); u0 = y4[0] ** 2 + y4[1] ** 2; up0 = 2 * (y4[0] * y4[2] + y4[1] * y4[3])
    umin = info['rmin'] ** 2; Tr = info['T_r']
    bs = np.linspace(bmin, bmax, nb)
    tr = np.empty(nb); det = np.empty(nb); clos = 0
    for i in range(0, nb, chunk):
        Bs = bs[i:i + chunk] * G2 * umin
        a, d, c = monodromy_traces(u0, up0, E, Tr, Bs)
        tr[i:i + chunk] = a; det[i:i + chunk] = d; clos = max(clos, c)
    nu = np.where(np.abs(tr) > 2, np.arccosh(np.maximum(np.abs(tr), 2) / 2) / Tr, 0.0)
    return bs, tr, det, nu, clos


def gaps(bs, tr):
    un = np.abs(tr) > 2
    out = []; i = 0
    while i < len(bs):
        if un[i]:
            j = i
            while j + 1 < len(bs) and un[j + 1]: j += 1
            out.append([float(bs[i]), float(bs[j]), int(np.sign(tr[i]))])
            i = j + 1
        else:
            i += 1
    return out


res = {}
y4 = tau20_state(); q = charge(y4)
lateA = 31.494634874218747   # pericentre of the late (A-matched) orbit, part A
for name, yy in [('tau20', y4), ('late', np.array([lateA, 0.0, 0.0, q / lateA]))]:
    info = orbit_info(yy, nsamp=1024)
    umin, umax = info['rmin'] ** 2, info['rmax'] ** 2
    bmax = 1.15 * umax / umin
    bs, tr, det, nu, clos = chart(yy, info, -60.0, bmax, 12001)
    gp = gaps(bs, tr)
    # refine every gap edge by bisection on |tr| - 2
    E = energy4(yy); u0 = yy[0] ** 2 + yy[1] ** 2; up0 = 2 * (yy[0] * yy[2] + yy[1] * yy[3])
    def trb(b): return monodromy_traces(u0, up0, E, info['T_r'], np.array([b * G2 * umin]))[0][0]
    def edge(a, b):
        fa = abs(trb(a)) - 2
        for _ in range(60):
            m = 0.5 * (a + b); fm = abs(trb(m)) - 2
            if np.sign(fm) == np.sign(fa): a, fa = m, fm
            else: b = m
        return 0.5 * (a + b)
    db = bs[1] - bs[0]
    ref = []
    for lo, hi, sg in gp:
        L = edge(lo - db, lo) if lo > bs[0] else lo
        H = edge(hi, hi + db) if hi < bs[-1] else None
        ref.append(dict(b_lo=L, b_hi=H, sign=sg,
                        M_over_Mperi_lo=float(np.sqrt(L)) if L > 0 else None,
                        M_over_Mperi_hi=(float(np.sqrt(H)) if (H is not None and H > 0) else None),
                        k_at_M0_lo=float(np.sqrt(-L * G2 * umin)) if L < 0 else None,
                        k_at_M0_hi=(float(np.sqrt(-H * G2 * umin)) if (H is not None and H < 0) else None),
                        max_nu_over_omega_r=float(nu[(bs >= lo) & (bs <= hi)].max() / info['omega_r'])))
    # exponents at the track's quoted points (k = 0)
    pts = {}
    for m in [0.75, 0.87, 1.0, 1.75, 3.0]:
        t_ = trb(m * m)
        pts[f'{m}'] = float(np.arccosh(abs(t_) / 2) / info['T_r'] / info['omega_r']) if abs(t_) > 2 else 0.0
    res[name] = dict(gaps=ref, n_finite_gaps=sum(1 for g in ref if g['b_hi'] is not None), g2_umin=G2 * umin,
                     nu_over_omega_r_at_M_over_Mperi=pts, max_det_err=float(np.max(np.abs(det - 1))), u_closure=clos,
                     omega_r=info['omega_r'], T_r=info['T_r'],
                     table=dict(b=bs[::20], trace=tr[::20], nu_over_omega_r=nu[::20] / info['omega_r']))
    print(name, 'g2 umin', G2 * umin, 'det err', res[name]['max_det_err'], 'closure', clos, time.time() - t0)
    for g in ref: print('   gap', g)
    print('   nu/omega_r at k=0:', pts)
save_json('skB_floquet.json', res)
print('elapsed', time.time() - t0)
