"""Track m1_vacuum, item 2: independent spherically symmetric Q-ball solver for GRUT-M1.
Independent code (Claude, 1 Nov 2026): Numerov (4th-order) discretisation of
    u'' = (S^2 + a F^2 - Omega^2) u,          u = rho F
    w'' = rho [ U_phi'(S) + F^2 S ],           w = rho (S - 1)
with Dirichlet u = w = 0 at rho = 0 and rho = L, solved by damped Newton on the sparse Jacobian.
Two modes: fixed Omega (unknowns u, w) and fixed Qbar = g^2 Q (extra unknown Omega, extra equation
Qbar = 4 pi Omega int u^2 drho).  Units: rho = m_X r, m_X = g v, a = 0.1, b = 3, c = 0.2 (see potential_m1.py).
Energy E/(Q m_X) = E_b/Qbar with E_b = 4 pi int rho^2 [ (F'^2 + S'^2)/2 + Omega^2 F^2/2 + U ] drho.
Not Codex's code: Codex used scipy.solve_bvp on (F, F', S, S', j) with Omega as a parameter.
"""
import json, math, sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.interpolate import CubicSpline

a, b, c = 0.1, 3.0, 0.2
G2 = 0.01
V_GEV = 0.924089637375362
MX_GEV = 0.1 * V_GEV


def Uphi(S):
    return b / 4 * (S * S - 1)**2 - c * (S - 1)**3


def Uphi_p(S):
    return b * (S * S - 1) * S - 3 * c * (S - 1)**2


def Uphi_pp(S):
    return b * (3 * S * S - 1) - 6 * c * (S - 1)


class Solver:
    def __init__(self, L, h):
        self.N = int(round(L / h))
        self.h = h
        self.L = self.N * h
        self.r = np.arange(1, self.N) * h          # interior nodes 1..N-1
        self.n = self.N - 1

    # ---- right-hand sides and partial derivatives at interior nodes
    def G(self, u, w, Om):
        r = self.r
        F = u / r
        S = 1 + w / r
        Gu = (S * S + a * F * F - Om * Om) * u
        Gw = r * (Uphi_p(S) + F * F * S)
        dGu_du = S * S + 3 * a * F * F - Om * Om
        dGu_dw = 2 * S * F
        dGu_dOm = -2 * Om * u
        dGw_dw = Uphi_pp(S) + F * F
        dGw_du = 2 * F * S
        return Gu, Gw, dGu_du, dGu_dw, dGu_dOm, dGw_du, dGw_dw

    def residual_and_jac(self, u, w, Om, Qbar=None):
        n, h = self.n, self.h
        Gu, Gw, Au, Aw, AOm, Bu, Bw = self.G(u, w, Om)
        # pad with boundary values (u_0 = w_0 = 0, G_0 = 0; u_N = w_N = 0, G_N = 0)
        up = np.r_[0.0, u, 0.0]
        wp = np.r_[0.0, w, 0.0]
        Gup = np.r_[0.0, Gu, 0.0]
        Gwp = np.r_[0.0, Gw, 0.0]
        k = h * h / 12
        Ru = up[2:] - 2 * up[1:-1] + up[:-2] - k * (Gup[2:] + 10 * Gup[1:-1] + Gup[:-2])
        Rw = wp[2:] - 2 * wp[1:-1] + wp[:-2] - k * (Gwp[2:] + 10 * Gwp[1:-1] + Gwp[:-2])
        # Jacobian blocks: d Ru_i / d u_j  = delta_{j,i±1} (1 - k Au_j) + delta_{j,i} (-2 - 10 k Au_i)
        idx = np.arange(n)
        rows, cols, vals = [], [], []
        def add(rr, cc, vv):
            rows.append(rr); cols.append(cc); vals.append(vv)
        # Ru wrt u
        add(idx, idx, -2 - 10 * k * Au)
        add(idx[:-1], idx[1:], 1 - k * Au[1:])
        add(idx[1:], idx[:-1], 1 - k * Au[:-1])
        # Ru wrt w (offset n in column)
        add(idx, idx + n, -10 * k * Aw)
        add(idx[:-1], idx[1:] + n, -k * Aw[1:])
        add(idx[1:], idx[:-1] + n, -k * Aw[:-1])
        # Rw wrt w
        add(idx + n, idx + n, -2 - 10 * k * Bw)
        add(idx[:-1] + n, idx[1:] + n, 1 - k * Bw[1:])
        add(idx[1:] + n, idx[:-1] + n, 1 - k * Bw[:-1])
        # Rw wrt u
        add(idx + n, idx, -10 * k * Bu)
        add(idx[:-1] + n, idx[1:], -k * Bu[1:])
        add(idx[1:] + n, idx[:-1], -k * Bu[:-1])
        if Qbar is None:
            J = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(2 * n, 2 * n))
            return np.r_[Ru, Rw], J
        # fixed-charge mode: extra unknown Omega (index 2n), extra equation charge
        AOmp = np.r_[0.0, AOm, 0.0]
        dRu_dOm = -k * (AOmp[2:] + 10 * AOmp[1:-1] + AOmp[:-2])
        add(idx, np.full(n, 2 * n), dRu_dOm)
        C = 4 * math.pi * Om * h * np.sum(u * u) - Qbar
        add(np.full(n, 2 * n), idx, 8 * math.pi * Om * h * u)
        add(np.array([2 * n]), np.array([2 * n]), np.array([4 * math.pi * h * np.sum(u * u)]))
        J = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(2 * n + 1, 2 * n + 1))
        return np.r_[Ru, Rw, C], J

    def solve(self, u, w, Om, Qbar=None, tol=1e-10, maxit=60, verbose=False):
        # tol 1e-10 on the max Numerov residual: the round-off floor is ~1e-13 * max|u|, i.e. ~5e-12 at Q = 1e7;
        # Newton converges quadratically so the solution error is far below the discretisation error.
        n = self.n
        x = np.r_[u, w] if Qbar is None else np.r_[u, w, Om]
        for it in range(maxit):
            if Qbar is None:
                R, J = self.residual_and_jac(x[:n], x[n:2 * n], Om)
            else:
                R, J = self.residual_and_jac(x[:n], x[n:2 * n], x[2 * n], Qbar)
            nr = np.max(np.abs(R))
            if verbose:
                print('  it', it, 'res', nr, 'Om', x[2 * n] if Qbar is not None else Om, flush=True)
            if nr < tol:
                break
            dx = spsolve(J, -R)
            # damping: keep the step bounded and the fields in a sane range
            lam = 1.0
            for _ in range(30):
                xn = x + lam * dx
                if Qbar is None:
                    Rn, _ = self.residual_and_jac(xn[:n], xn[n:2 * n], Om)
                else:
                    if not (0 < xn[2 * n] < 1.5):
                        lam *= 0.5; continue
                    Rn, _ = self.residual_and_jac(xn[:n], xn[n:2 * n], xn[2 * n], Qbar)
                if np.all(np.isfinite(Rn)) and np.max(np.abs(Rn)) < nr * (1 - 1e-4 * lam) or lam < 1e-3:
                    break
                lam *= 0.5
            x = xn
        else:
            raise RuntimeError('Newton did not converge, residual %g' % nr)
        u, w = x[:n], x[n:2 * n]
        if Qbar is not None:
            Om = float(x[2 * n])
        return u, w, Om, nr, it

    # ---- observables
    def observables(self, u, w, Om):
        h, r = self.h, self.r
        F = u / r
        S = 1 + w / r
        # derivatives of u and w by 4th-order central differences on the padded arrays
        up = np.r_[0.0, u, 0.0]
        wp = np.r_[0.0, w, 0.0]
        def d1(y):
            # 4th-order central, 2nd-order near ends (ends are zero anyway at L; near 0 use y(0)=0 symmetry: u odd)
            yy = np.r_[-y[2], -y[1], y, 0.0, 0.0]   # odd extension at 0 (u(-r) = -u(r)), zero beyond L
            return (-yy[4:] + 8 * yy[3:-1] - 8 * yy[1:-3] + yy[:-4]) / (12 * h)
        du = d1(up)[1:-1]
        dw = d1(wp)[1:-1]
        Fp = (du - u / r) / r
        Sp = (dw - w / r) / r
        U = Uphi(S) + 0.5 * S * S * F * F + a / 4 * F**4
        grad = 0.5 * (Fp * Fp + Sp * Sp)
        kin = 0.5 * Om * Om * F * F
        I = lambda z: 4 * math.pi * h * np.sum(r * r * z)      # trapezoid; endpoints are zero
        Egrad, Ekin, EU = I(grad), I(kin), I(U)
        E = Egrad + Ekin + EU
        Qbar = Om * I(F * F)
        # central values by even extrapolation F(0) = (4F(h) - F(2h))/3 + O(h^4)
        F0 = (4 * F[0] - F[1]) / 3
        S0 = (4 * S[0] - S[1]) / 3
        virial = (Egrad + 3 * (EU - Ekin)) / E
        E_virial = 2 * Om * Qbar - 2 * EU          # from virial identity Egrad = 3(Ekin - EU) and Ekin = Om Qbar / 2
        # radii containing 50/90/99% of charge
        cum = np.cumsum(r * r * F * F) * h * 4 * math.pi * Om / Qbar
        rad = {p: float(np.interp(p, cum, r)) for p in (0.5, 0.9, 0.99)}
        return dict(Omega=Om, Qbar=Qbar, Q=Qbar / G2, E_b=E, E_over_Qm=E / Qbar, E_virial_over_Qm=E_virial / Qbar,
                    mass_GeV=E * V_GEV / 0.1, Egrad=Egrad, Ekin=Ekin, EU=EU, virial_residual=virial,
                    F0=F0, S0=S0, r50=rad[0.5], r90=rad[0.9], r99=rad[0.99], L=self.L, h=self.h)

    def guess(self, Fc, Sc, R, d):
        r = self.r
        prof = 1 / (1 + np.exp((r - R) / d))
        F = Fc * prof
        S = 1 - (1 - Sc) * prof
        return r * F, r * (S - 1)

    def interp_from(self, other, u, w):
        """Map a solution from another grid (cubic spline in u, w, zero outside)."""
        cu = CubicSpline(np.r_[0.0, other.r, other.L], np.r_[0.0, u, 0.0])
        cw = CubicSpline(np.r_[0.0, other.r, other.L], np.r_[0.0, w, 0.0])
        rr = self.r
        un = np.where(rr < other.L, cu(np.minimum(rr, other.L)), 0.0)
        wn = np.where(rr < other.L, cw(np.minimum(rr, other.L)), 0.0)
        return un, wn


def solve_fixed_Q(Q, L=60.0, h=0.025, seed=None, Om0=0.88, verbose=False):
    sv = Solver(L, h)
    if seed is None:
        R = (3 * G2 * Q / (4 * math.pi * Om0 * 2.6**2))**(1 / 3)
        u, w = sv.guess(2.6, 0.13, R, 0.5)
        Om = Om0
    else:
        s0, u0, w0, Om = seed
        u, w = sv.interp_from(s0, u0, w0)
    u, w, Om, res, it = sv.solve(u, w, Om, Qbar=G2 * Q, verbose=verbose)
    ob = sv.observables(u, w, Om)
    ob.update(newton_residual=res, newton_iterations=it)
    return sv, u, w, Om, ob


def main():
    out = {}
    t0 = time.time()
    # ---------------- primary object Q = 200000 at three resolutions (Richardson check) and two domains
    runs = {}
    seed = None
    for (L, h) in [(60.0, 0.05), (60.0, 0.025), (60.0, 0.0125), (90.0, 0.025)]:
        sv, u, w, Om, ob = solve_fixed_Q(200000, L, h, seed=seed)
        seed = (sv, u, w, Om)
        runs['L%g_h%g' % (L, h)] = ob
        print('Q=2e5 L=%g h=%g  Om=%.10f E/Qm=%.10f (virial form %.10f) S0=%.8f F0=%.8f vir=%.1e res=%.1e it=%d  t=%.1fs' % (
            L, h, Om, ob['E_over_Qm'], ob['E_virial_over_Qm'], ob['S0'], ob['F0'], ob['virial_residual'], ob['newton_residual'], ob['newton_iterations'], time.time() - t0), flush=True)
    out['primary_resolution_runs'] = runs
    # Richardson (4th order) on Omega and E/Qm from h = 0.05, 0.025, 0.0125 at L = 60
    o1, o2, o3 = runs['L60_h0.05'], runs['L60_h0.025'], runs['L60_h0.0125']
    def rich(k):
        return (16 * o3[k] - o2[k]) / 15, (o3[k] - o2[k]) / (o2[k] - o1[k])
    prim = {}
    for k in ('Omega', 'E_over_Qm', 'S0', 'F0', 'mass_GeV'):
        val, ratio = rich(k)
        prim[k] = val; prim[k + '_convergence_ratio'] = ratio
    prim['E_over_Qm_virial_form_h0.0125'] = o3['E_virial_over_Qm']
    prim['codex'] = dict(Omega=0.8817204919524214, E_over_Qm=0.9373999905262296, S0=0.1257717809745518,
                         F0=2.6309605904405573, mass_GeV=17324.832346421026)
    prim['diff_vs_codex'] = {k: prim[k] - prim['codex'][k] for k in prim['codex']}
    out['primary'] = prim
    print('PRIMARY Richardson:', json.dumps(prim, indent=1), flush=True)

    # ---------------- profile comparison with Codex's PROFILE_M1.csv (rho step 0.05, L = 160)
    sv, u, w, Om = seed  # L=90 h=0.025
    cod = np.loadtxt(sys.argv[1], delimiter=',', skiprows=1)
    rr = cod[:, 0]
    cu = CubicSpline(np.r_[0.0, sv.r, sv.L], np.r_[0.0, u, 0.0])
    cw = CubicSpline(np.r_[0.0, sv.r, sv.L], np.r_[0.0, w, 0.0])
    m = (rr > 0) & (rr < 60)
    Fm = cu(rr[m]) / rr[m]
    Sm = 1 + cw(rr[m]) / rr[m]
    out['profile_vs_codex'] = dict(max_abs_dF=float(np.max(np.abs(Fm - cod[m, 1]))), max_abs_dS=float(np.max(np.abs(Sm - cod[m, 3]))),
                                   rho_range=[float(rr[m][0]), float(rr[m][-1])])
    print('profile diff vs Codex: dF max %.2e dS max %.2e' % (out['profile_vs_codex']['max_abs_dF'], out['profile_vs_codex']['max_abs_dS']), flush=True)

    # ---------------- neighbours: first law dE/dQ = Omega m_X, dOmega/dQ (h = 0.025, L = 60)
    base = runs['L60_h0.025']
    s60 = Solver(60.0, 0.025)
    sv0, u0, w0, Om0, _ = solve_fixed_Q(200000, 60.0, 0.025, seed=seed)
    nb = {}
    for dQ in (500, 1000, 2000):
        res = {}
        for sgn in (-1, 1):
            _, _, _, _, ob = solve_fixed_Q(200000 + sgn * dQ, 60.0, 0.025, seed=(sv0, u0, w0, Om0))
            res[sgn] = ob
        dE_dQ = (res[1]['E_b'] - res[-1]['E_b']) / (2 * dQ) / G2     # E_b/Qbar units -> dE/dQ in m_X
        dOm_dQ = (res[1]['Omega'] - res[-1]['Omega']) / (2 * dQ)
        nb[str(dQ)] = dict(dE_dQ_over_mX=dE_dQ, Omega=Om0, first_law_rel=dE_dQ / Om0 - 1, dOmega_dQ=dOm_dQ,
                           dE_dQ_GeV=dE_dQ * MX_GEV)
        print('dQ=%d dE/dQ=%.10f Om=%.10f rel=%.2e dOm/dQ=%.4e' % (dQ, dE_dQ, Om0, dE_dQ / Om0 - 1, dOm_dQ), flush=True)
    nb['codex'] = dict(neighbor_dE_dQ_GeV=0.08147890035738237, Omega_mX_GeV=0.08147887696747388,
                       relative_difference=2.870671436117078e-07, neighbor_dOmega_dQ=-2.0500726844691065e-07)
    out['first_law'] = nb

    # ---------------- E(Q)/Q for Q = 1e5 ... 1e7 (fixed-Q continuation), approach to the bulk value
    Om_bulk = 0.7758647074595397
    scan = []
    cur = (sv0, u0, w0, Om0)
    for Q in (150000, 100000, 80000, 70000):
        svq, uq, wq, Omq, ob = solve_fixed_Q(Q, 60.0, 0.025, seed=cur)
        cur = (svq, uq, wq, Omq)
        ob['surface_coefficient'] = (ob['E_over_Qm'] - Om_bulk) * Q**(1 / 3)
        scan.append(ob)
        print('Q=%d Om=%.8f E/Qm=%.8f S0=%.5f r90=%.3f (E/Q-Ob)Q^1/3=%.4f' % (Q, Omq, ob['E_over_Qm'], ob['S0'], ob['r90'], ob['surface_coefficient']), flush=True)
    cur = (sv0, u0, w0, Om0)
    for Q, L in ((300000, 60.0), (500000, 60.0), (1000000, 60.0), (2000000, 70.0), (5000000, 80.0), (10000000, 90.0), (30000000, 110.0)):
        svq, uq, wq, Omq, ob = solve_fixed_Q(Q, L, 0.025, seed=cur)
        cur = (svq, uq, wq, Omq)
        ob['surface_coefficient'] = (ob['E_over_Qm'] - Om_bulk) * Q**(1 / 3)
        scan.append(ob)
        print('Q=%d Om=%.8f E/Qm=%.8f S0=%.5f r90=%.3f (E/Q-Ob)Q^1/3=%.4f  t=%.0fs' % (Q, Omq, ob['E_over_Qm'], ob['S0'], ob['r90'], ob['surface_coefficient'], time.time() - t0), flush=True)
    scan.append(dict(Q=200000, **{k: base[k] for k in ('Omega', 'E_over_Qm', 'S0', 'F0', 'r90', 'mass_GeV')}))
    scan.sort(key=lambda d: d['Q'])
    out['EQ_scan'] = scan
    out['codex_Q1e6'] = dict(Omega=0.8342966820499844, E_over_Qm=0.8670265120251068, S0=0.12122433502210861, mass_GeV=80121.0215092106)

    # ---------------- fixed-Omega branch toward Omega -> 1: Q(Omega), E(Omega); Q_min and E/Q = 1 crossing
    branch = []
    svb = Solver(100.0, 0.05)
    ub, wb = svb.interp_from(sv0, u0, w0)
    Om = Om0
    ub, wb, _, _, _ = svb.solve(ub, wb, Om)
    for Om in list(np.arange(0.882, 0.9801, 0.002)) + list(np.arange(0.981, 0.9951, 0.001)):
        try:
            ub, wb, _, res, it = svb.solve(ub, wb, float(Om))
        except RuntimeError as e:
            print('branch stopped at Omega', Om, e); break
        ob = svb.observables(ub, wb, float(Om))
        branch.append(dict(Omega=float(Om), Q=ob['Q'], E_over_Qm=ob['E_over_Qm'], S0=ob['S0'], F0=ob['F0'], r90=ob['r90'], r99=ob['r99'],
                           E_minus_OmQ_over_mX=(ob['E_b'] - Om * ob['Qbar']) / G2))
    out['branch_fixed_Omega'] = branch
    Qs = np.array([d['Q'] for d in branch]); Oms = np.array([d['Omega'] for d in branch]); EQ = np.array([d['E_over_Qm'] for d in branch])
    i = int(np.argmin(Qs))
    # parabolic refinement of Q_min
    if 0 < i < len(Qs) - 1:
        pa = np.polyfit(Oms[i - 1:i + 2], Qs[i - 1:i + 2], 2)
        Om_min = -pa[1] / (2 * pa[0]); Q_min = float(np.polyval(pa, Om_min))
    else:
        Om_min, Q_min = float(Oms[i]), float(Qs[i])
    # E/Q = 1 crossing on the stable side (Omega < Om_min)
    cross = None
    for j in range(len(branch) - 1):
        if (EQ[j] - 1) * (EQ[j + 1] - 1) < 0 and Oms[j] < Om_min:
            t = (1 - EQ[j]) / (EQ[j + 1] - EQ[j])
            cross = dict(Omega=float(Oms[j] + t * (Oms[j + 1] - Oms[j])), Q=float(Qs[j] + t * (Qs[j + 1] - Qs[j])))
    out['Q_min'] = dict(Omega=float(Om_min), Q=float(Q_min), grid_point=dict(Omega=float(Oms[i]), Q=float(Qs[i])))
    out['E_over_Q_equals_mX'] = cross
    out['dQ_dOmega_sign_at_reference'] = 'negative (stable branch) if Omega_ref < Omega(Q_min)'
    print('Q_min', out['Q_min'], 'E/Q=1 crossing', cross, flush=True)
    for d in branch[::5]:
        print('  Om=%.3f Q=%.4e E/Qm=%.6f S0=%.4f r90=%.2f' % (d['Omega'], d['Q'], d['E_over_Qm'], d['S0'], d['r90']))
    out['runtime_s'] = time.time() - t0
    json.dump(out, open(sys.argv[2], 'w'), indent=1, default=float)
    # save the primary profile (h = 0.025, L = 60) for the Hessian script
    np.savez(sys.argv[3], r=sv0.r, u=u0, w=w0, Om=Om0, h=sv0.h, L=sv0.L)


if __name__ == '__main__':
    main()
