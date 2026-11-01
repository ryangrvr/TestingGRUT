"""m1_vacuum SKEPTIC, part 2: independent Q-ball solver on a DIFFERENT discretisation from both Codex
(scipy solve_bvp collocation on (F,F',S,S',j)) and the track (Numerov on u = rho F, w = rho(S-1), fixed-Q Newton).
Independent code (Claude skeptic, 1 Nov 2026).

Method: cell-centred grid rho_i = (i+1/2)h, conservative second-order discretisation of (rho^2 F')' = rho^2 G_F,
(rho^2 S')' = rho^2 G_S; regularity at the origin is automatic (zero flux through rho = 0); Dirichlet at rho = L
through ghost cells. Fixed-Omega damped Newton (sparse LU). Charge Qbar = 4 pi Omega h sum rho_i^2 F_i^2 (midpoint
rule); the Omega that gives Qbar = g^2 Q is found by a secant iteration. Second-order errors: two-level Richardson
over h, h/2, h/4 (error ratios must be 4). Energy E_b = 4 pi int rho^2 [F'^2/2 + S'^2/2 + Omega^2 F^2/2 + U].
Also: DOP853 shooting from (F0, S0, Omega) as a divergence-radius diagnostic of the central values, the fixed-Omega
branch Q(Omega), E(Omega) for Q_min and the E/Q = m_X crossing, and the first law dE/dQbar = Omega from the branch.
Units as in the track: rho = m_X r, a = 0.1, b = 3, c = 0.2, Qbar = g^2 Q, E/(Q m_X) = E_b/Qbar, mass = E_b v/g.
"""
import json, math, sys, time
import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import splu
from scipy.integrate import solve_ivp

a, b, c = 0.1, 3.0, 0.2
G2 = 0.01
V_GEV = 0.924089637375362


def Uphi(S):   return b / 4 * (S * S - 1)**2 - c * (S - 1)**3
def Uphi_p(S): return b * (S * S - 1) * S - 3 * c * (S - 1)**2
def Uphi_pp(S): return b * (3 * S * S - 1) - 6 * c * (S - 1)


class CC:
    """cell-centred conservative FD solver at fixed Omega"""
    def __init__(self, L, h):
        self.N = N = int(round(L / h)); self.h = h; self.L = N * h
        self.r = (np.arange(N) + 0.5) * h
        self.r2 = self.r**2
        Af = (np.arange(N + 1) * h)**2            # face areas rho_{j}^2, j = 0..N (face N at L)
        main = -(Af[1:] + Af[:-1]) / h**2
        main[-1] = -(2 * Af[-1] + Af[-2]) / h**2  # ghost: X_N = -X_{N-1} (Dirichlet 0 at face N)
        off = Af[1:-1] / h**2
        self.Dm = sps.diags([off, main, off], [-1, 0, 1], format='csr')
        self.srcS = np.zeros(N); self.srcS[-1] = 2 * Af[-1] / h**2       # S ghost = 2 - S_{N-1}
        self.Af = Af

    def residual(self, F, S, Om):
        """unknowns are F and w = S - 1 (so the stiff Dirichlet stencil never acts on O(1) values; the ghost for w
        is w_N = -w_{N-1}, i.e. the same matrix as for F and no source term)"""
        GF = (S * S + a * F * F - Om * Om) * F
        GS = Uphi_p(S) + F * F * S
        return np.r_[self.Dm @ F - self.r2 * GF, self.Dm @ (S - 1.0) - self.r2 * GS]

    def jacobian(self, F, S, Om):
        N = self.N
        dFF = sps.diags(-self.r2 * (S * S + 3 * a * F * F - Om * Om))
        dFS = sps.diags(-self.r2 * 2 * S * F)
        dSS = sps.diags(-self.r2 * (Uphi_pp(S) + F * F))
        dSF = sps.diags(-self.r2 * 2 * F * S)
        return sps.bmat([[self.Dm + dFF, dFS], [dSF, self.Dm + dSS]], format='csc')

    def solve(self, F, S, Om, tol=1e-10, maxit=60):
        """damped Newton (2-norm line search). tol is on the max-norm of the residual, whose round-off floor is
        ~ rho^2 * 1e-16 * |fields| ~ 1e-12 here; Newton is quadratic so the field error is far below that."""
        N = self.N
        x = np.r_[F, S]
        for it in range(maxit):
            R = self.residual(x[:N], x[N:], Om)
            nr = np.max(np.abs(R)); n2 = np.linalg.norm(R)
            if nr < tol:
                return x[:N], x[N:], nr, it
            dx = splu(self.jacobian(x[:N], x[N:], Om)).solve(-R)
            lam = 1.0
            while lam > 1e-4:
                xn = x + lam * dx
                Rn = self.residual(xn[:N], xn[N:], Om)
                if np.all(np.isfinite(Rn)) and np.linalg.norm(Rn) < (1 - 1e-4 * lam) * n2:
                    break
                lam *= 0.5
            if lam <= 1e-4 and nr < 1e-6:
                return x[:N], x[N:], nr, it          # round-off floor reached (never triggered after the w = S-1 reformulation)
            x = xn
        raise RuntimeError('Newton failed, residual %g' % nr)

    def observables(self, F, S, Om):
        h, r2, N = self.h, self.r2, self.N
        Fg = np.r_[F, -F[-1]]; Sg = np.r_[S, 2 - S[-1]]           # ghosts at L (w = S-1 odd about the face)
        Fe = np.r_[F[0], Fg]; Se = np.r_[S[0], Sg]                 # even ghost at origin
        dF = (Fe[2:] - Fe[:-2]) / (2 * h); dS = (Se[2:] - Se[:-2]) / (2 * h)
        U = Uphi(S) + 0.5 * S * S * F * F + a / 4 * F**4
        I = lambda z: 4 * math.pi * h * np.sum(r2 * z)
        Egrad = I(0.5 * (dF * dF + dS * dS))
        # face-based gradient energy (alternative second-order form)
        dFf = (Fg[1:] - Fg[:-1]) / h; dSf = (Sg[1:] - Sg[:-1]) / h
        Egrad_face = 4 * math.pi * h * np.sum(self.Af[1:] * 0.5 * (dFf**2 + dSf**2))
        Ekin = I(0.5 * Om * Om * F * F); EU = I(U)
        Qbar = Om * I(F * F)
        E = Egrad + Ekin + EU
        F0 = (9 * F[0] - F[1]) / 8; S0 = (9 * S[0] - S[1]) / 8     # even extrapolation to rho = 0
        cum = np.cumsum(r2 * F * F) * 4 * math.pi * h * Om / Qbar
        return dict(Omega=Om, Qbar=Qbar, Q=Qbar / G2, E_b=E, E_b_face=Egrad_face + Ekin + EU, E_over_Qm=E / Qbar,
                    mass_GeV=E * V_GEV / 0.1, Egrad=Egrad, Ekin=Ekin, EU=EU, F0=F0, S0=S0,
                    virial=(Egrad + 3 * (EU - Ekin)) / E, r90=float(np.interp(0.9, cum, self.r)), r99=float(np.interp(0.99, cum, self.r)),
                    L=self.L, h=self.h)

    def guess(self, Fc=2.63, Sc=0.126, R=5.0, d=1.0):
        """tanh guess; at fixed Omega the Newton basin is narrow (R = 5.9, d = 0.6 diverges; R = 5, d = 1 converges in 6 steps)"""
        p = 1 / (1 + np.exp((self.r - R) / d))
        return Fc * p, 1 - (1 - Sc) * p

    def interp(self, other, F, S):
        return np.interp(self.r, other.r, F, right=0.0), np.interp(self.r, other.r, S, right=1.0)


def fixed_Q(sv, Qbar, Om0, F, S, tol=1e-10):
    """Omega such that the discrete charge equals Qbar. Phase 1: walk along the fixed-Omega branch in steps of 0.004
    (dQ/dOmega < 0 on the thin-wall side) until the target is bracketed; phase 2: Illinois (safeguarded secant) in
    the bracket. Returns the solution and observables."""
    def q_at(Om, F, S):
        F, S, _, _ = sv.solve(F, S, Om); return sv.observables(F, S, Om)['Qbar'] - Qbar, F, S
    Om1 = Om0; q1, F, S = q_at(Om1, F, S)
    step = 0.004 if q1 > 0 else -0.004          # too much charge -> raise Omega
    Om2 = Om1 + step; q2, F, S = q_at(Om2, F, S)
    while q1 * q2 > 0:
        Om1, q1 = Om2, q2
        Om2 = Om2 + step; q2, F, S = q_at(Om2, F, S)
    # Illinois
    side = 0
    for it in range(60):
        Om3 = Om2 - q2 * (Om2 - Om1) / (q2 - q1)
        q3, F, S = q_at(Om3, F, S)
        if abs(q3) < tol:
            Om2, q2 = Om3, q3; break
        if q3 * q2 < 0:
            Om1, q1 = Om2, q2; side = 0
        else:
            q1 *= 0.5; side += 1
        Om2, q2 = Om3, q3
    ob = sv.observables(F, S, Om2); ob['secant_iterations'] = it; ob['charge_residual'] = q2
    return F, S, Om2, ob


def shoot(F0, S0, Om, rmax=40.0):
    """DOP853 outward integration from the series at rho0; returns the radius where the solution leaves the vacuum"""
    r0 = 1e-3
    GF0 = (S0 * S0 + a * F0 * F0 - Om * Om) * F0; GS0 = Uphi_p(S0) + F0 * F0 * S0
    y0 = [F0 + GF0 * r0**2 / 6, GF0 * r0 / 3, S0 + GS0 * r0**2 / 6, GS0 * r0 / 3]
    def rhs(r, y):
        F, Fp, S, Sp = y
        return [Fp, -2 * Fp / r + (S * S + a * F * F - Om * Om) * F, Sp, -2 * Sp / r + Uphi_p(S) + F * F * S]
    # divergence = the fields leave the range of the true solution (|F| <= 2.64, |S - 1| <= 0.88)
    def leave(r, y): return 3.0 - max(abs(y[0]), abs(y[2] - 1))
    leave.terminal = True
    sol = solve_ivp(rhs, (r0, rmax), y0, method='DOP853', rtol=1e-13, atol=1e-18, events=[leave], dense_output=True)
    if len(sol.t_events[0]):
        rdiv = float(sol.t_events[0][0])
    else:
        rdiv = float(sol.t[-1])          # integration ended (failure = blow-up) or reached rmax without leaving
    return rdiv, sol


def main():
    out = {}
    t0 = time.time()
    QBAR = 2000.0
    # ---------- primary: h = 0.04, 0.02, 0.01 at L = 60, then L = 80 at h = 0.02
    grids = [(60.0, 0.04), (60.0, 0.02), (60.0, 0.01), (80.0, 0.02)]
    runs = {}
    prev = None
    for (L, h) in grids:
        sv = CC(L, h)
        if prev is None:
            F, S = sv.guess(); Om = 0.8817
        else:
            F, S = sv.interp(prev[0], prev[1], prev[2]); Om = prev[3]
        F, S, Om, ob = fixed_Q(sv, QBAR, Om, F, S)
        prev = (sv, F, S, Om)
        runs['L%g_h%g' % (L, h)] = ob
        print('L=%g h=%g Om=%.12f E/Qm=%.12f (face %.12f) F0=%.10f S0=%.10f mass=%.8f vir=%.1e  t=%.1fs' % (
            L, h, Om, ob['E_over_Qm'], ob['E_b_face'] / QBAR, ob['F0'], ob['S0'], ob['mass_GeV'], ob['virial'], time.time() - t0), flush=True)
    o1, o2, o3 = runs['L60_h0.04'], runs['L60_h0.02'], runs['L60_h0.01']
    prim = {}
    for k in ('Omega', 'E_over_Qm', 'F0', 'S0', 'mass_GeV'):
        r1 = (4 * o2[k] - o1[k]) / 3; r2 = (4 * o3[k] - o2[k]) / 3
        prim[k] = (16 * r2 - r1) / 15
        prim[k + '_level1'] = r2
        prim[k + '_error_ratio_2nd_order_expect_4'] = (o2[k] - o1[k]) / (o3[k] - o2[k])
    # face-form energy, extrapolated the same way
    e1, e2, e3 = o1['E_b_face'] / QBAR, o2['E_b_face'] / QBAR, o3['E_b_face'] / QBAR
    prim['E_over_Qm_face_form'] = (16 * (4 * e3 - e2) / 3 - (4 * e2 - e1) / 3) / 15
    prim['L80_vs_L60_Omega_change_at_h0.02'] = runs['L80_h0.02']['Omega'] - o2['Omega']
    prim['codex'] = dict(Omega=0.8817204919524214, E_over_Qm=0.9373999905262296, S0=0.1257717809745518,
                         F0=2.6309605904405573, mass_GeV=17324.832346421026)
    prim['track'] = dict(Omega=0.8817204919524271, E_over_Qm=0.9373999905262296, S0=0.1257717809745547,
                         F0=2.630960590440569, mass_GeV=17324.832346421023)
    prim['mine_minus_codex'] = {k: prim[k] - prim['codex'][k] for k in prim['codex']}
    prim['mine_minus_codex_relative'] = {k: (prim[k] - prim['codex'][k]) / prim['codex'][k] for k in prim['codex']}
    out['primary_runs'] = runs; out['primary'] = prim
    print('PRIMARY (2-level Richardson):', json.dumps({k: prim[k] for k in ('Omega', 'E_over_Qm', 'F0', 'S0', 'mass_GeV', 'mine_minus_codex', 'mine_minus_codex_relative')}, indent=1), flush=True)
    print('error ratios:', {k: prim[k + '_error_ratio_2nd_order_expect_4'] for k in ('Omega', 'E_over_Qm', 'F0', 'S0')}, flush=True)

    # ---------- shooting diagnostic of the central values (DOP853)
    sh = {}
    cases = {'codex': (2.6309605904405573, 0.1257717809745518, 0.8817204919524214),
             'mine': (prim['F0'], prim['S0'], prim['Omega']),
             'codex_S0_plus_1e-9': (2.6309605904405573, 0.1257717809745518 + 1e-9, 0.8817204919524214),
             'codex_S0_plus_1e-6': (2.6309605904405573, 0.1257717809745518 + 1e-6, 0.8817204919524214),
             'codex_F0_plus_1e-6': (2.6309605904405573 + 1e-6, 0.1257717809745518, 0.8817204919524214)}
    for name, (F0, S0, Om) in cases.items():
        rdiv, sol = shoot(F0, S0, Om)
        sh[name] = dict(F0=F0, S0=S0, Omega=Om, divergence_radius=rdiv, implied_precision_S_channel=math.exp(-math.sqrt(6) * rdiv))
        print('shoot %-22s rho_div = %.3f  (S-channel growth e^{sqrt6 rho}: implied |delta| ~ %.1e)' % (name, rdiv, math.exp(-math.sqrt(6) * rdiv)), flush=True)
    out['shooting'] = sh

    # ---------- E/Q at other charges: Codex's Q = 1e6 control and a few track rows (h = 0.02 and 0.01, level-1 Richardson)
    ctrl = {}
    for Q, L in ((1000000, 60.0), (100000, 60.0), (500000, 60.0), (70000, 60.0)):
        vals = {}
        prevq = None
        for h in (0.02, 0.01):
            sv = CC(L, h)
            if prevq is None:
                F, S = sv.interp(prev[0], prev[1], prev[2]); Om = prev[3]
            else:
                F, S = sv.interp(prevq[0], prevq[1], prevq[2]); Om = prevq[3]
            F, S, Om, ob = fixed_Q(sv, G2 * Q, Om, F, S)
            prevq = (sv, F, S, Om); vals[h] = ob
        ctrl[str(Q)] = {k: (4 * vals[0.01][k] - vals[0.02][k]) / 3 for k in ('Omega', 'E_over_Qm', 'S0', 'F0', 'mass_GeV', 'r90')}
        print('Q=%d  Om=%.9f E/Qm=%.9f S0=%.6f mass=%.4f GeV  t=%.0fs' % (Q, ctrl[str(Q)]['Omega'], ctrl[str(Q)]['E_over_Qm'], ctrl[str(Q)]['S0'], ctrl[str(Q)]['mass_GeV'], time.time() - t0), flush=True)
    ctrl['codex_Q1e6'] = dict(Omega=0.8342966820499844, E_over_Qm=0.8670265120251068, S0=0.12122433502210861, mass_GeV=80121.0215092106)
    ctrl['track'] = {'100000': dict(Omega=0.916103, E_over_Qm=0.978960), '500000': dict(Omega=0.850935, E_over_Qm=0.892765),
                     '70000': dict(Omega=0.941722, E_over_Qm=1.001065)}
    out['controls'] = ctrl

    # ---------- fixed-Omega branch: Q(Omega), E(Omega); first law; Q_min; E/Q = 1 crossing   (h = 0.02, L = 100)
    svb = CC(100.0, 0.02)
    F, S = svb.interp(prev[0], prev[1], prev[2])
    F, S, _, _ = svb.solve(F, S, prev[3])
    branch = []
    Om_ref = 0.8817204919524214
    Oms = [Om_ref - 1e-3, Om_ref - 5e-4, Om_ref, Om_ref + 5e-4, Om_ref + 1e-3] + list(np.arange(0.884, 0.9801, 0.004)) + list(np.arange(0.982, 0.9941, 0.002))
    for Om in Oms:
        try:
            F, S, nr, it = svb.solve(F, S, float(Om))
        except RuntimeError as e:
            print('branch stopped at', Om, e); break
        ob = svb.observables(F, S, float(Om))
        branch.append(dict(Omega=float(Om), Qbar=ob['Qbar'], Q=ob['Q'], E_b=ob['E_b'], E_b_face=ob['E_b_face'], E_over_Qm=ob['E_over_Qm'], S0=ob['S0'], F0=ob['F0'], r99=ob['r99']))
    out['branch'] = branch
    Om_arr = np.array([d['Omega'] for d in branch]); Q_arr = np.array([d['Qbar'] for d in branch])
    E_arr = np.array([d['E_b_face'] for d in branch])      # the face-form energy is the discrete functional whose
    E_cc = np.array([d['E_b'] for d in branch])            # stationarity gives the discrete equations (up to the last cell)
    # first law from the branch around the reference: dE/dQbar vs Omega (central differences on the 5 close points)
    i = 2
    fl = dict(dE_dQ_central=(E_arr[i + 1] - E_arr[i - 1]) / (Q_arr[i + 1] - Q_arr[i - 1]), Omega=Om_arr[i],
              dE_dQ_wide=(E_arr[i + 2] - E_arr[i - 2]) / (Q_arr[i + 2] - Q_arr[i - 2]),
              dE_dQ_central_cellcentred_energy=(E_cc[i + 1] - E_cc[i - 1]) / (Q_arr[i + 1] - Q_arr[i - 1]),
              dQ_dOmega=(Q_arr[i + 1] - Q_arr[i - 1]) / (Om_arr[i + 1] - Om_arr[i - 1]) / G2,
              Qbar_at_reference_Omega=Q_arr[i])
    fl['first_law_rel_central'] = fl['dE_dQ_central'] / fl['Omega'] - 1
    fl['first_law_rel_wide'] = fl['dE_dQ_wide'] / fl['Omega'] - 1
    fl['first_law_rel_cellcentred_energy'] = fl['dE_dQ_central_cellcentred_energy'] / fl['Omega'] - 1
    fl['dOmega_dQ'] = 1 / fl['dQ_dOmega']
    fl['note'] = 'the discrete first law holds to O(dOmega^2) only for the variationally consistent (face-form) energy; the cell-centred central-difference energy violates it at O(h^2) ~ 3e-5'
    out['first_law_from_branch'] = fl
    print('first law from branch (face energy): dE/dQ/Om-1 = %.2e (central, dOm=5e-4), %.2e (wide, 1e-3); cell-centred energy: %.2e; dOmega/dQ = %.5e ; Qbar at Om_ref = %.6f' % (
        fl['first_law_rel_central'], fl['first_law_rel_wide'], fl['first_law_rel_cellcentred_energy'], fl['dOmega_dQ'], Q_arr[i]), flush=True)
    # Q_min and E/Q = 1 crossing
    j = int(np.argmin(Q_arr))
    pa = np.polyfit(Om_arr[j - 1:j + 2], Q_arr[j - 1:j + 2], 2)
    Om_min = -pa[1] / (2 * pa[0]); Q_min = np.polyval(pa, Om_min) / G2
    EQ = E_arr / Q_arr
    cross = None
    for jj in range(len(branch) - 1):
        if (EQ[jj] - 1) * (EQ[jj + 1] - 1) < 0 and Om_arr[jj] < Om_min:
            t = (1 - EQ[jj]) / (EQ[jj + 1] - EQ[jj])
            cross = dict(Omega=float(Om_arr[jj] + t * (Om_arr[jj + 1] - Om_arr[jj])), Q=float((Q_arr[jj] + t * (Q_arr[jj + 1] - Q_arr[jj])) / G2))
    out['Q_min'] = dict(Omega=float(Om_min), Q=float(Q_min), track=dict(Omega=0.9728, Q=56710))
    out['E_over_Q_equals_mX'] = dict(mine=cross, track=dict(Omega=0.9402, Q=71263))
    print('Q_min = %.0f at Omega = %.4f ; E/Q = 1 at Q = %.0f, Omega = %.4f' % (Q_min, Om_min, cross['Q'], cross['Omega']), flush=True)
    for d in branch[::4]:
        print('  Om=%.4f Q=%.5e E/Qm=%.6f S0=%.4f r99=%.1f' % (d['Omega'], d['Q'], d['E_over_Qm'], d['S0'], d['r99']))
    out['runtime_s'] = time.time() - t0
    json.dump(out, open(sys.argv[1], 'w'), indent=1, default=float)
    np.savez(sys.argv[2], r=prev[0].r, F=prev[1], S=prev[2], Om=prev[3], h=prev[0].h, L=prev[0].L)


if __name__ == '__main__':
    main()
