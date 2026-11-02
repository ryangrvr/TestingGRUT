"""qm_micro item 1a (Claude, 2 Nov 2026): linear (Bogoliubov) scattering of X and Xbar on the planar M1 Q-matter
interface, exact at linear order.  Own code.

Rotating frame psi = e^{-i Omega t} (F(x) + u + i v), S -> S(x) + s, Omega = Omega_b, plane waves
e^{i k_perp y}; with (u, v, s) = (U, iV, Sx) e^{-i w t} + c.c. the linearised field equations become the REAL
SYMMETRIC system  Y'' = M(x; w, k_perp) Y,  Y = (U, V, Sx):
    M = [[S^2+3aF^2-Om^2-w^2+kp^2, 2 Om w, 2 S F],
         [2 Om w, S^2+aF^2-Om^2-w^2+kp^2, 0],
         [2 S F, 0, Uphi''(S)+F^2-w^2+kp^2]].
(Rotated to (U-V, U+V)/sqrt2 this is the canonical operator of bf8ab3b / sk_xcore.py with nu = w.)
Lab-frame content: psi = (U-V) e^{-i(Om+w)t} + (U+V)^* e^{-i(Om-w)t}: in vacuum the channel (1,-1,0)/sqrt2 is an X
of energy Om + w, (1,1,0)/sqrt2 an Xbar of energy w - Om (open for w > Om + sqrt(1+kp^2)), (0,0,1) a phi of energy w.
Every quantum in every channel carries rotating-frame energy w, and the time-averaged rotating-frame energy flux is
2 w Im(Y^dag Y'), so the Wronskian J = Im(Y^dag Y') is the quantum-number flux in every channel (positive-norm
modes; the norm sign -e^T dM/dw e is checked for every channel).  An incoming X of energy E has w = E - Omega_b; an
incoming Xbar of energy E has w = E + Omega_b.
Discretisation: Numerov in W = (1 - h^2 M/12) Y:  W_{j+1} + W_{j-1} = C_j W_j, C_j = 2 + h^2 M (1 - h^2 M/12)^{-1}
(real symmetric, so Im(W_j^dag W_{j+1}) is conserved exactly).  Radiation conditions from the exact discrete modes
of the asymptotic matrices (z + 1/z = c_i): decaying or outgoing (group velocity away from the interface) on both
sides, one incoming channel on the vacuum side.  Block-tridiagonal sparse solve.
"""
import json, math, os, sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import qm_common as qc

a_ = qc.a_


class Interface:
    def __init__(self, Lb=60.0, Lv=60.0, h=0.01, hs=None, profile=None):
        P = qc.planar_profile(Lb=Lb, Lv=Lv, h=h) if profile is None else profile
        self.P = P
        self.Om = P['Omega_b']
        x, F, S = P['x'], P['F'], P['S']
        if hs is not None and abs(hs - h) > 1e-15:
            # resample on a finer/coarser grid by cubic interpolation of the profile
            from scipy.interpolate import CubicSpline
            xs = np.arange(x[0], x[-1] + 0.5 * hs, hs)
            F = CubicSpline(x, F)(xs); S = CubicSpline(x, S)(xs); x = xs; h = hs
        self.x, self.F, self.S, self.h = x, F, S, h
        self.N = x.size

    def Mfield(self, w, kp):
        F, S, Om = self.F, self.S, self.Om
        n = self.N
        M = np.zeros((n, 3, 3))
        M[:, 0, 0] = S * S + 3 * a_ * F * F - Om * Om - w * w + kp * kp
        M[:, 1, 1] = S * S + a_ * F * F - Om * Om - w * w + kp * kp
        M[:, 2, 2] = qc.Uphi_pp(S) + F * F - w * w + kp * kp
        M[:, 0, 1] = M[:, 1, 0] = 2 * Om * w
        M[:, 0, 2] = M[:, 2, 0] = 2 * S * F
        return M

    def Mw(self, w):
        Om = self.Om
        return np.array([[-2 * w, 2 * Om, 0.0], [2 * Om, -2 * w, 0.0], [0.0, 0.0, -2 * w]])


def modes(Mc, h, Mw, side, scheme='numerov'):
    """Discrete modes of the constant matrix Mc.  side = 'L' (bulk, x -> -inf) or 'R' (vacuum, x -> +inf).
    Returns eigvecs (columns), allowed z (decaying or outgoing on that side), incoming z (R only), info."""
    lam, E = np.linalg.eigh(Mc)
    c = 2 + h * h * lam / (1 - h * h * lam / 12) if scheme == 'numerov' else 2 + h * h * lam
    zs, zin, info = [], [], []
    for i in range(3):
        e = E[:, i]
        normsign = -float(e @ Mw @ e)
        if c[i] > 2 or c[i] < -2 or lam[i] >= 0:
            r = np.roots([1, -c[i], 1])
            r = np.sort(np.abs(r.real)) * np.sign(r.real[np.argsort(np.abs(r.real))])
            # decaying away from the interface: L side needs |z| > 1, R side |z| < 1
            z = r[1] if side == 'L' else r[0]
            zs.append(complex(z)); zin.append(None)
            info.append(dict(eig=float(lam[i]), open=False, kappa=float(math.log(abs(r[1])) / h), normsign=normsign))
        else:
            th = math.acos(c[i] / 2)
            zp, zm = complex(math.cos(th), math.sin(th)), complex(math.cos(th), -math.sin(th))
            # positive norm: energy flux ~ J = |amp|^2 Im z ; outgoing on L needs J < 0, on R J > 0
            pos = normsign > 0
            if side == 'L':
                z = zm if pos else zp
                zi = zp if pos else zm
            else:
                z = zp if pos else zm
                zi = zm if pos else zp
            zs.append(z); zin.append(zi)
            info.append(dict(eig=float(lam[i]), open=True, k=float(th / h), k_cont=float(math.sqrt(-lam[i])), normsign=normsign))
    return E, np.array(zs), zin, info


def solve(I, w, kp, inc, scheme='numerov'):
    """inc: 'X' (incoming channel (1,-1,0)/sqrt2) or 'Xbar' ((1,1,0)/sqrt2).  Returns flux fractions.
    scheme 'numerov' (continuum to O(h^4)) or 'fd2' (the exact frequency-domain image of a second-order lattice)."""
    h = I.h
    M = I.Mfield(w, kp)
    n = I.N
    if scheme == 'numerov':
        T = np.eye(3)[None] - h * h / 12 * M
        Tinv = np.linalg.inv(T)
        C = 2 * np.eye(3)[None] + h * h * np.einsum('nij,njk->nik', M, Tinv)
    else:
        C = 2 * np.eye(3)[None] + h * h * M
    Mw = I.Mw(w)
    EL, zL, _, infoL = modes(M[0], h, Mw, 'L', scheme)
    ER, zR, zinR, infoR = modes(M[-1], h, Mw, 'R', scheme)
    KL = (EL * (1 / zL)[None, :]) @ EL.T          # W_{-1} = KL W_0
    KR = (ER * zR[None, :]) @ ER.T                # outgoing/decaying part: W_{N+1} = KR W_N
    # incoming channel vector in the vacuum basis
    ein = np.array([1.0, -1.0, 0.0]) / math.sqrt(2) if inc == 'X' else np.array([1.0, 1.0, 0.0]) / math.sqrt(2)
    # identify it among ER columns
    ic = int(np.argmax(np.abs(ER.T @ ein)))
    if abs(abs(ER[:, ic] @ ein) - 1) > 1e-10:
        raise RuntimeError('vacuum eigenvectors are not the channel vectors (degenerate?)')
    ein = ER[:, ic]
    if not infoR[ic]['open']:
        return dict(w=w, kp=kp, inc=inc, closed=True)
    z_in = zinR[ic]; z_out = zR[ic]
    # assemble block tridiagonal: rows j: -W_{j-1} + C_j W_j - W_{j+1} = rhs
    Cb = C.astype(complex)
    Cb[0] = Cb[0] - KL
    Cb[-1] = Cb[-1] - KR
    rows, cols, vals = [], [], []
    base = 3 * np.arange(n)
    for a in range(3):
        for b in range(3):
            rows.append(base + a); cols.append(base + b); vals.append(Cb[:, a, b].astype(complex))
    for a in range(3):
        rows.append(base[1:] + a); cols.append(base[:-1] + a); vals.append(-np.ones(n - 1, complex))
        rows.append(base[:-1] + a); cols.append(base[1:] + a); vals.append(-np.ones(n - 1, complex))
    A = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(3 * n, 3 * n))
    rhs = np.zeros(3 * n, complex)
    Ain = 1.0
    rhs[3 * (n - 1):] = Ain * (z_in - z_out) * ein
    Wsol = splu(A).solve(rhs).reshape(n, 3)
    # fluxes (discrete Wronskian, units of 1/h dropped consistently)
    J = np.imag(np.einsum('ni,ni->n', Wsol[:-1].conj(), Wsol[1:]))
    out = {}
    inc_flux = abs(Ain) ** 2 * abs(z_in.imag)
    b = ER.T @ (Wsol[-1] - Ain * ein)
    names_R = []
    for i in range(3):
        v = ER[:, i]
        nm = 'X' if abs(v @ np.array([1, -1, 0]) / math.sqrt(2)) > 0.99 else ('Xbar' if abs(v @ np.array([1, 1, 0]) / math.sqrt(2)) > 0.99 else 'phi')
        names_R.append(nm)
        if infoR[i]['open']:
            out['out_' + nm] = abs(b[i]) ** 2 * abs(zR[i].imag) / inc_flux
    cL = EL.T @ Wsol[0]
    trans = []
    for i in range(3):
        if infoL[i]['open']:
            trans.append(dict(k=infoL[i]['k'], k_cont=infoL[i]['k_cont'], frac=abs(cL[i]) ** 2 * abs(zL[i].imag) / inc_flux,
                              normsign=infoL[i]['normsign'], vec=[float(t) for t in EL[:, i]]))
    trans.sort(key=lambda d: d['k'])
    tot = sum(v for k, v in out.items()) + sum(t['frac'] for t in trans)
    res = dict(w=w, kp=kp, inc=inc, closed=False, inc_channel_k=infoR[ic]['k'],
               reflected=out.get('out_' + inc, 0.0), out=out, trans=trans, trans_total=sum(t['frac'] for t in trans),
               unitarity=tot - 1.0, J_spread=float((J.max() - J.min()) / inc_flux),
               bulk_channels=infoL, vac_channels=[dict(name=nm, **inf) for nm, inf in zip(names_R, infoR)])
    return res


def bulk_branches(Om, SB, FB, w, kp=0.0):
    """Bulk propagating channels at (w, kp): k_x, norm sign, group velocity sign; continuum operator."""
    D = qc.Uphi_pp(SB) + FB * FB
    M = np.array([[SB * SB + 3 * a_ * FB * FB - Om * Om - w * w + kp * kp, 2 * Om * w, 2 * SB * FB],
                  [2 * Om * w, SB * SB + a_ * FB * FB - Om * Om - w * w + kp * kp, 0.0],
                  [2 * SB * FB, 0.0, D - w * w + kp * kp]])
    lam, E = np.linalg.eigh(M)
    Mw = np.array([[-2 * w, 2 * Om, 0.0], [2 * Om, -2 * w, 0.0], [0.0, 0.0, -2 * w]])
    out = []
    for i in range(3):
        if lam[i] < 0:
            e = E[:, i]
            q = math.sqrt(-lam[i])
            out.append(dict(kx=q, norm=-float(e @ Mw @ e), vg=-2 * q / float(e @ Mw @ e)))
    return out


def scan(out):
    t0 = time.time()
    I = Interface(h=0.01)
    Om = I.Om
    P = I.P
    res = dict(profile=dict(Omega=Om, Omega_minus_Omega_b=P['Omega'] - P['Omega_b'], sigma_planar=P['sigma_equip'],
                            sigma_direct=P['sigma_direct'], x_gibbs=P['x_gibbs'], S_B=P['S_B'], F_B=P['F_B'], n_B=P['n_B'],
                            grid=dict(Lb=60, Lv=60, h=I.h)))
    # bulk propagating channels along the X and Xbar lines (kp = 0)
    res['bulk_channels'] = [dict(w=float(w), channels=bulk_branches(Om, P['S_B'], P['F_B'], float(w))) for w in np.linspace(0.05, 9.0, 180)]
    # X sticking at normal incidence over a wide energy range
    Ts = np.unique(np.r_[np.geomspace(1e-4, 0.02, 12), np.linspace(0.02, 1.0, 50), np.geomspace(1.0, 8.0, 25)])
    rowsX = []
    for T in Ts:
        E = 1 + T
        r = solve(I, E - Om, 0.0, 'X')
        rowsX.append(dict(T=float(T), E=E, w=E - Om, s=1 - r['reflected'], reflected=r['reflected'], out=r['out'],
                          trans=[dict(k=t['k'], frac=t['frac']) for t in r['trans']], unitarity=r['unitarity'],
                          J_spread=r['J_spread']))
    print('X normal incidence: %d energies, max |unitarity| %.1e, t=%.0fs' % (len(rowsX), max(abs(r['unitarity']) for r in rowsX), time.time() - t0), flush=True)
    for r in rowsX[::6]:
        print('   T=%.4g s=%.6f out=%s trans=%s' % (r['T'], r['s'], {k: round(v, 6) for k, v in r['out'].items()}, [(round(t['k'], 4), round(t['frac'], 6)) for t in r['trans']]), flush=True)
    res['X_normal'] = rowsX
    rowsB = []
    for T in Ts:
        E = 1 + T
        r = solve(I, E + Om, 0.0, 'Xbar')
        rowsB.append(dict(T=float(T), E=E, w=E + Om, reflected=r['reflected'], destroyed=1 - r['reflected'], out=r['out'],
                          trans=[dict(k=t['k'], frac=t['frac']) for t in r['trans']], unitarity=r['unitarity'],
                          J_spread=r['J_spread']))
    print('Xbar normal incidence: max |unitarity| %.1e t=%.0fs' % (max(abs(r['unitarity']) for r in rowsB), time.time() - t0), flush=True)
    for r in rowsB[::6]:
        print('   T=%.4g destroyed=%.6f out=%s trans=%s' % (r['T'], r['destroyed'], {k: round(v, 6) for k, v in r['out'].items()}, [(round(t['k'], 4), round(t['frac'], 6)) for t in r['trans']]), flush=True)
    res['Xbar_normal'] = rowsB
    # oblique incidence grid
    Tg = np.unique(np.r_[np.geomspace(1e-3, 0.02, 6), np.linspace(0.02, 1.0, 25), np.geomspace(1.0, 8.0, 10)])
    thg = np.radians(np.r_[np.linspace(0, 85, 18), 88, 89.5])
    grid = {}
    for name in ('X', 'Xbar'):
        tab = np.zeros((Tg.size, thg.size)); uni = 0.0
        for i, T in enumerate(Tg):
            E = 1 + T; p = math.sqrt(E * E - 1)
            for j, th in enumerate(thg):
                kp = p * math.sin(th)
                r = solve(I, (E - Om) if name == 'X' else (E + Om), kp, name)
                tab[i, j] = 1 - r['reflected']
                uni = max(uni, abs(r['unitarity']))
        grid[name] = dict(T=Tg.tolist(), theta_deg=np.degrees(thg).tolist(), one_minus_R=tab.tolist(), max_unitarity=uni)
        print('%s oblique grid done, max |unitarity| %.1e, t=%.0fs' % (name, uni, time.time() - t0), flush=True)
    res['oblique'] = grid
    # resolution control
    I2 = Interface(h=0.01, hs=0.005)
    ctrl = []
    for T in (0.02, 0.1, 0.3, 1.0, 4.0):
        r1 = solve(I, 1 + T - Om, 0.0, 'X'); r2 = solve(I2, 1 + T - Om, 0.0, 'X')
        q1 = solve(I, 1 + T + Om, 0.0, 'Xbar'); q2 = solve(I2, 1 + T + Om, 0.0, 'Xbar')
        ctrl.append(dict(T=T, s_h01=1 - r1['reflected'], s_h005=1 - r2['reflected'], dbar_h01=1 - q1['reflected'], dbar_h005=1 - q2['reflected']))
    res['resolution_control'] = ctrl
    I3 = Interface(Lb=90.0, Lv=90.0, h=0.01)
    dom = []
    for T in (0.02, 0.3, 1.0):
        r3 = solve(I3, 1 + T - Om, 0.0, 'X'); q3 = solve(I3, 1 + T + Om, 0.0, 'Xbar')
        dom.append(dict(T=T, s_L90=1 - r3['reflected'], dbar_L90=1 - q3['reflected']))
    res['domain_control'] = dom
    print('controls', ctrl, dom)
    res['runtime_s'] = time.time() - t0
    json.dump(res, open(out, 'w'), indent=1, default=float)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'test':
        I = Interface()
        for T in (0.02, 0.1, 0.3, 1.0):
            t0 = time.time()
            r = solve(I, 1 + T - I.Om, 0.0, 'X')
            print('X T=%.3f s=%.8f out=%s trans=%s uni=%.2e Jspread=%.1e (%.2fs)' % (T, 1 - r['reflected'], r['out'], [(t['k'], t['frac'], t['normsign']) for t in r['trans']], r['unitarity'], r['J_spread'], time.time() - t0))
            q = solve(I, 1 + T + I.Om, 0.0, 'Xbar')
            print('   Xbar destroyed=%.8f out=%s trans=%s uni=%.2e' % (1 - q['reflected'], q['out'], [(round(t['k'], 4), t['frac'], t['normsign']) for t in q['trans']], q['unitarity']))
            print('   bulk channels X:', bulk_branches(I.Om, I.P['S_B'], I.P['F_B'], 1 + T - I.Om))
    else:
        scan(sys.argv[1])
