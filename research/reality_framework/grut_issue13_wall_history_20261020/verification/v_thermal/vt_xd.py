"""X and D curvatures along phi_min(T), analytic second derivatives (Daleckii-Krein divided differences),
plus a coarse global 3-field search.  Independent of the module (which used finite differences).
Writes results/xd_v.json."""
import json
import os
import numpy as np
from scipy.optimize import brentq, minimize
import vt_core as V

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
P = V.BENCH
MU = V.MU0
PI = np.pi


def Fp_parts(m2, T, mu, kind):
    """first derivative w.r.t. m^2 of the single-dof function used on eigenvalues of a given matrix.
    kind 'cwJ': CW + thermal; 'daisy': -(T/12pi) Re(m2)^{3/2}; 'daisy_sub': +(T/12pi) Re(m2)^{3/2}"""
    m2 = np.atleast_1d(np.asarray(m2, float))
    if kind == 'cwJ':
        out = V.cwp(m2, mu)
        if T > 0:
            out = out + V.vthp(m2, T)
        return out
    if kind == 'daisy':
        return -T / (8 * PI) * V.re_pow12(m2)
    if kind == 'daisy_sub':
        return T / (8 * PI) * V.re_pow12(m2)


def divdiff(a, b, T, mu, kind):
    if abs(a - b) > 1e-7:
        return float((Fp_parts(a, T, mu, kind) - Fp_parts(b, T, mu, kind))[0] / (a - b))
    h = 1e-6
    return float((Fp_parts(a + h, T, mu, kind) - Fp_parts(a - h, T, mu, kind))[0] / (2 * h))


def pieces(scheme, T):
    Pi = V.debye(P)
    Pd = np.array([Pi['phi'], Pi['X'], Pi['X'], Pi['D']]) * T**2
    zero = np.zeros(4)
    if scheme == 'oneloop':
        return [(zero, 'cwJ')]
    if scheme == 'AE':
        return [(zero, 'cwJ'), (Pd, 'daisy'), (zero, 'daisy_sub')] if T > 0 else [(zero, 'cwJ')]
    if scheme == 'Parwani':
        return [(Pd, 'cwJ')]


def diag_masses(phi):
    return np.array([P['mp'] + 3 * P['lp'] * phi**2, P['mx'] + P['c'] * phi**2,
                     P['mx'] + P['c'] * phi**2, P['md'] + P['d'] * phi**2])


def curv_dir(phi, T, scheme, direction):
    """d^2 V / d q1^2 (direction 'X') or d^2V/dD^2 ('D') at (phi, 0, 0, 0)."""
    m = diag_masses(phi)
    if direction == 'X':
        tree = P['mx'] + P['c'] * phi**2
        j, coup = 1, 2 * P['c'] * phi
        Mpp = np.array([2 * P['c'], 6 * P['lx'], 2 * P['lx'], 2 * P['al']])
    else:
        tree = P['md'] + P['d'] * phi**2
        j, coup = 3, 2 * P['d'] * phi
        Mpp = np.array([2 * P['d'], 2 * P['al'], 2 * P['al'], 6 * P['ld']])
    tot = tree
    for shift, kind in pieces(scheme, T):
        mm = m + shift
        Fp = Fp_parts(mm, T, MU, kind)
        tot += float(np.sum(Fp * Mpp))
        tot += 2 * coup**2 * divdiff(mm[0], mm[j], T, MU, kind)
    return tot


def dVds_scheme(phi, T, scheme):
    if T == 0:
        # T = 0: CW only
        s = phi**2
        out = P['mp'] / 2 + P['lp'] * s / 2
        for n, m2, dm in V.masses_phi(np.array([phi]), P):
            out = out + n * dm * V.cwp(m2, MU)
        return float(out[0])
    return float(V.dVds(np.array([phi]), T, P, MU, scheme)[0])


def phi_min(T, scheme):
    g0 = dVds_scheme(1e-6, T, scheme)
    if g0 >= 0:
        return 0.0
    return brentq(lambda x: dVds_scheme(x, T, scheme), 1e-6, 1.5, xtol=1e-15)


# ------------------------------------------------ full 4-field V (eigenvalues) for global search
def hessian(phi, q1, D):
    lp, lx, c, d, ld, al = P['lp'], P['lx'], P['c'], P['d'], P['ld'], P['al']
    n = np.broadcast(phi, q1, D).shape
    H = np.zeros(n + (4, 4))
    H[..., 0, 0] = P['mp'] + 3 * lp * phi**2 + c * q1**2 + d * D**2
    H[..., 1, 1] = P['mx'] + c * phi**2 + 3 * lx * q1**2 + al * D**2
    H[..., 2, 2] = P['mx'] + c * phi**2 + lx * q1**2 + al * D**2
    H[..., 3, 3] = P['md'] + d * phi**2 + 3 * ld * D**2 + al * q1**2
    H[..., 0, 1] = H[..., 1, 0] = 2 * c * phi * q1
    H[..., 0, 3] = H[..., 3, 0] = 2 * d * phi * D
    H[..., 1, 3] = H[..., 3, 1] = 2 * al * q1 * D
    return H


def V4(phi, q1, D, T, scheme):
    lp, lx, c, d, ld, al = P['lp'], P['lx'], P['c'], P['d'], P['ld'], P['al']
    R2 = q1**2
    V0 = (P['mp'] * phi**2 / 2 + P['mx'] * R2 / 2 + P['md'] * D**2 / 2 + lp * phi**4 / 4 + lx * R2**2 / 4
          + ld * D**4 / 4 + c * phi**2 * R2 / 2 + d * phi**2 * D**2 / 2 + al * R2 * D**2 / 2)
    H = hessian(phi, q1, D)
    Pi = V.debye(P)
    Pd = np.diag([Pi['phi'], Pi['X'], Pi['X'], Pi['D']]) * T**2
    ev = np.linalg.eigvalsh(H)
    out = V0 + 0.0
    if scheme in ('oneloop', 'AE'):
        out = out + np.sum(V.cw(ev, MU), axis=-1)
        if T > 0:
            out = out + np.sum(V.vth(ev.reshape(-1), T).reshape(ev.shape), axis=-1)
    if scheme == 'AE' and T > 0:
        evP = np.linalg.eigvalsh(H + Pd)
        out = out - T / (12 * PI) * (np.sum(V.re_pow32(evP), axis=-1) - np.sum(V.re_pow32(ev), axis=-1))
    if scheme == 'Parwani':
        evP = np.linalg.eigvalsh(H + Pd)
        out = out + np.sum(V.cw(evP, MU), axis=-1)
        if T > 0:
            out = out + np.sum(V.vth(evP.reshape(-1), T).reshape(evP.shape), axis=-1)
    return out


def main():
    res = {}
    Tc_AE = 1.5634661935  # placeholder overwritten below with own value
    tr = json.load(open(os.path.join(OUT, 'transition_v.json')))
    Tcs = {s: tr['fixed']['mu0'][s]['Tc'] for s in ['oneloop', 'AE', 'Parwani']}
    for scheme in ['oneloop', 'AE', 'Parwani']:
        Tc = Tcs[scheme]
        Ts = sorted(set(list(np.linspace(0, 3, 31)) + [Tc * (1 - 1e-3), Tc, Tc * (1 + 1e-3)]))
        rows = []
        for T in Ts:
            pm = phi_min(T, scheme) if T > 0 else phi_min(0, scheme)
            mX = curv_dir(pm, T, scheme, 'X')
            mD = curv_dir(pm, T, scheme, 'D')
            Pi = V.debye(P)
            lead_X = P['c'] * pm**2 + Pi['X'] * T**2
            rows.append(dict(T=T, phi_min=pm, m2X_eff=mX, m2D_eff=mD, m2X_leading=lead_X))
        res[scheme] = dict(rows=rows, min_m2X=min(r['m2X_eff'] for r in rows),
                           T_at_min_m2X=min(rows, key=lambda r: r['m2X_eff'])['T'],
                           min_m2D=min(r['m2D_eff'] for r in rows),
                           T_at_min_m2D=min(rows, key=lambda r: r['m2D_eff'])['T'])
        r = min(rows, key=lambda r: r['m2X_eff'])
        res[scheme]['ratio_min_to_leading'] = r['m2X_eff'] / r['m2X_leading']
        print(scheme, res[scheme]['min_m2X'], res[scheme]['T_at_min_m2X'], res[scheme]['ratio_min_to_leading'],
              res[scheme]['min_m2D'], res[scheme]['T_at_min_m2D'], flush=True)
    # finite-difference cross-check of analytic curvature using the 4x4 eigenvalue V (AE)
    chk = []
    for T, ph in [(1.0, 0.77), (1.4, 0.45), (Tcs['AE'], 0.0), (2.0, 0.0), (1.2, 0.707)]:
        h = 1e-2
        fX = (V4(ph, h, 0.0, T, 'AE') - 2 * V4(ph, 0.0, 0.0, T, 'AE') + V4(ph, -h, 0.0, T, 'AE')) / h**2
        fD = (V4(ph, 0.0, h, T, 'AE') - 2 * V4(ph, 0.0, 0.0, T, 'AE') + V4(ph, 0.0, -h, T, 'AE')) / h**2
        chk.append(dict(T=T, phi=ph, fdX=float(fX), anX=curv_dir(ph, T, 'AE', 'X'),
                        fdD=float(fD), anD=curv_dir(ph, T, 'AE', 'D')))
    res['fd_check_AE'] = chk
    print(chk)
    # coarse global search
    glob = []
    g = np.linspace(0, 1.6, 33)
    PH, Q, DD = np.meshgrid(g, g, g, indexing='ij')
    for scheme in ['AE', 'Parwani']:
        for T in [0.5, 1.0, 1.4, Tcs[scheme], 1.7, 2.5]:
            Vg = V4(PH, Q, DD, T, scheme)
            i = np.unravel_index(np.argmin(Vg), Vg.shape)
            # off-axis minimum
            mask = (Q > 0) | (DD > 0)
            Voff = np.where(mask, Vg, np.inf)
            j = np.unravel_index(np.argmin(Voff), Vg.shape)
            # Nelder-Mead from off-axis starts
            nm = []
            for st in [(0.5, 0.5, 0.5), (0.1, 1.0, 0.1), (1.0, 0.2, 1.0), (0.0, 0.8, 0.0), (0.3, 0.0, 1.2)]:
                r = minimize(lambda x: float(V4(x[0], x[1], x[2], T, scheme)), st, method='Nelder-Mead',
                             options=dict(xatol=1e-7, fatol=1e-15, maxiter=20000, maxfev=40000))
                nm.append([abs(r.x[0]), abs(r.x[1]), abs(r.x[2])])
            glob.append(dict(scheme=scheme, T=T, grid_min=[float(PH[i]), float(Q[i]), float(DD[i])],
                             off_axis_gap=float(Voff[j] - Vg[i]),
                             nm_max_absR=max(x[1] for x in nm), nm_max_absD=max(x[2] for x in nm),
                             nm_phis=[x[0] for x in nm]))
            print(glob[-1], flush=True)
    res['global'] = glob
    with open(os.path.join(OUT, 'xd_v.json'), 'w') as f:
        json.dump(res, f, indent=1, default=float)


if __name__ == '__main__':
    main()
