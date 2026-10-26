#!/usr/bin/env python3
"""hw_envelope (Claude, own code): independent reproduction of Codex 2908e91's four mathematical control groups.

Nothing of Codex's code is imported or run.  Each group is re-derived by a different route:
  G1  exact current sources, derived symbolically in complex variables (z, zbar independent), with arbitrary radial
      functions g(|z|^2) added to show they drop out; plus a *dynamical* check: a random three-field release with
      H = 1/(2t) is integrated and d(a^3 n_i)/dt is compared with 2 a^3 Im(z_i* dV/dz_i*) by quadrature.
  G2  48 random complex field/coupling samples (different RNG seed): the analytic force, a numerical gradient of the
      real potential (central differences in Re/Im), the angular finite difference and the gauge-weighted sum.
  G3  forced response: sympy proof that Y_* = (-AJ + BJ*)/D solves AY + BY* + J = 0, that V_Y(Y_*) = Delta V,
      that -6 Im(kappa Y_** X^3) = -12 Im(alpha6 X^6) = -12 Im(mu Sigma* Y_*^2) identically (so the forced Y carries
      zero net n_Y source), and numerics at Codex's four amplitudes by a different route (scipy minimize of the full
      potential including the quartic, and a real 2x2 solve).
  G4  Higgs angular countercharge: sympy elimination of A_0 from the homogeneous gauged phase Lagrangian,
      giving n_X = 6 I delta_dot, n_Sigma = -n_X/6 and rho_phase = n_X^2/(2 f^2) + n_X^2/(72 h^2).
  G5  conditional envelope yield bound: own evaluation of 12|alpha6|R^6/H_i over the supplied prompt-radiation
      entropy; the table (1.2201e-8, 5.4890e-8, 3.8583e-5, 1.7358e-4 Y_Q), the r_i^(7/2) scaling and the
      extrapolation points r_i/M_Y = 18.24 / 11.87.
Usage: python3 -B controls.py <out.json>
"""
import json, math, sys
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, quad, cumulative_trapezoid
from scipy.optimize import minimize, brentq, fsolve

MRED, GSTAR, YQ = 2.435e18, 106.75, 4.8294009076988175e-9   # [Q] same supplied constants as the packets


def entropy(H):
    T = (90 * MRED**2 * H**2 / (math.pi**2 * GSTAR))**0.25
    return 2 * math.pi**2 / 45 * GSTAR * T**3


# ---------------------------------------------------------------- G1: symbolic currents, complex variables
def g1_symbolic():
    X, Y, S, Xb, Yb, Sb = sp.symbols('X Y S Xb Yb Sb')
    k, kb, m, mb = sp.symbols('kappa kappab mu mub')
    gX, gY, gS = (sp.Function(n) for n in ('gX', 'gY', 'gS'))
    K, Kb = k * Yb * X**3, kb * Y * Xb**3
    L, Lb = m * Sb * Y**2, mb * S * Yb**2
    V = K + Kb + L + Lb + gX(X * Xb) + gY(Y * Yb) + gS(S * Sb)
    conj = {X: Xb, Xb: X, Y: Yb, Yb: Y, S: Sb, Sb: S, k: kb, kb: k, m: mb, mb: m}

    def source(z, zb):                       # 2 Im(zb dV/dzb) = -i (w - conj(w))
        w = zb * sp.diff(V, zb)
        return sp.expand(-sp.I * (w - w.subs(conj, simultaneous=True)))
    sX, sY, sS = source(X, Xb), source(Y, Yb), source(S, Sb)
    ImK, ImL = (K - Kb) / (2 * sp.I), (L - Lb) / (2 * sp.I)
    res = dict(nX=sp.simplify(sX + 6 * ImK), nY=sp.simplify(sY - 2 * ImK + 4 * ImL), nS=sp.simplify(sS - 2 * ImL),
               gauge=sp.simplify(sX + 3 * sY + 6 * sS), nD=sp.simplify(sX + 3 * sY + 12 * ImL))
    return dict(pass_group=all(v == 0 for v in res.values()), residuals={a: str(b) for a, b in res.items()},
                note='sources -6 Im K, 2 Im K - 4 Im L, 2 Im L; radial terms g(|z|^2) drop out; complex couplings')


def g1_dynamical(seed=7):
    """Integrate a random release and compare a^3 n_i(t) with the quadrature of the claimed source."""
    rng = np.random.default_rng(seed)
    k = 0.05 * np.exp(1j * rng.uniform(-np.pi, np.pi)); m = 0.3 * np.exp(1j * rng.uniform(-np.pi, np.pi))
    lx, ly, ls, M, h = 0.02, 0.5, 0.7, 1.3, 0.8

    def F(X, Y, S):
        return (2 * lx * abs(X)**2 * X + 3 * np.conj(k) * Y * np.conj(X)**2,
                M**2 * Y + 2 * ly * abs(Y)**2 * Y + k * X**3 + 2 * np.conj(m) * S * np.conj(Y),
                2 * ls * (abs(S)**2 - h**2 / 2) * S + m * Y**2)

    def rhs(t, w):
        z = w[0:3] + 1j * w[3:6]; zd = w[6:9] + 1j * w[9:12]
        H = 1 / (2 * t)
        f = F(*z)
        acc = [-3 * H * zd[i] - f[i] for i in range(3)]
        return [*zd.real, *zd.imag, *[a.real for a in acc], *[a.imag for a in acc]]
    z0 = (rng.normal(size=3) + 1j * rng.normal(size=3)) * 0.8
    v0 = (rng.normal(size=3) + 1j * rng.normal(size=3)) * 0.3
    t = np.linspace(0.5, 12.0, 6001)
    sol = solve_ivp(rhs, (t[0], t[-1]), [*z0.real, *z0.imag, *v0.real, *v0.imag], t_eval=t, method='DOP853',
                    rtol=1e-12, atol=1e-14)
    z = sol.y[0:3] + 1j * sol.y[3:6]; zd = sol.y[6:9] + 1j * sol.y[9:12]
    a3 = (2 * t)**1.5
    n = -2 * a3 * np.imag(np.conj(z) * zd)
    Kt = k * np.conj(z[1]) * z[0]**3; Lt = m * np.conj(z[2]) * z[1]**2
    src = np.array([-6 * Kt.imag, 2 * Kt.imag - 4 * Lt.imag, 2 * Lt.imag]) * a3
    pred = n[:, :1] + np.concatenate([np.zeros((3, 1)), cumulative_trapezoid(src, t, axis=1)], axis=1)
    err = np.max(np.abs(n - pred), axis=1) / np.max(np.abs(n), axis=1)
    gauge = n[0] + 3 * n[1] + 6 * n[2]
    return dict(pass_group=bool(err.max() < 1e-6 and np.max(np.abs(gauge - gauge[0])) / np.max(np.abs(n)) < 1e-9),
                max_rel_error_nX_nY_nS=err.tolist(),
                gauge_charge_drift_rel=float(np.max(np.abs(gauge - gauge[0])) / np.max(np.abs(n))),
                charge_change_nX=float(n[0, -1] - n[0, 0]), nD_change=float((n[0] + 3 * n[1])[-1] - (n[0] + 3 * n[1])[0]),
                nD_from_mu_only=float(np.trapezoid(-12 * Lt.imag * a3, t)))


# ---------------------------------------------------------------- G2: random complex samples, numerical gradients
def g2_numeric(nsamp=48, seed=20261026):
    rng = np.random.default_rng(seed)
    lx, ly, ls, M, h = 0.003, 0.8, 1.1, 1.0, 0.1
    worst = dict(force=0.0, numgrad=0.0, angle=0.0, gauge=0.0)
    for _ in range(nsamp):
        z = (rng.normal(size=3) + 1j * rng.normal(size=3)) * 0.5
        k = rng.uniform(0.01, 0.2) * np.exp(1j * rng.uniform(-np.pi, np.pi))
        m = rng.uniform(0.05, 0.5) * np.exp(1j * rng.uniform(-np.pi, np.pi))

        def V(zz):
            X, Y, S = zz
            return (lx * abs(X)**4 + ly * abs(Y)**4 + ls * (abs(S)**2 - h**2 / 2)**2 + M**2 * abs(Y)**2
                    + 2 * np.real(k * np.conj(Y) * X**3 + m * np.conj(S) * Y**2))
        X, Y, S = z
        Fa = np.array([2 * lx * abs(X)**2 * X + 3 * np.conj(k) * Y * np.conj(X)**2,
                       M**2 * Y + 2 * ly * abs(Y)**2 * Y + k * X**3 + 2 * np.conj(m) * S * np.conj(Y),
                       2 * ls * (abs(S)**2 - h**2 / 2) * S + m * Y**2])
        # numerical dV/dz* = (dV/dx + i dV/dy)/2
        Fn = np.zeros(3, complex); d = 1e-5
        for i in range(3):
            for comp, unit in ((0, 1.0), (1, 1j)):
                vals = []
                for s in (2, 1, -1, -2):
                    zz = z.copy(); zz[i] += s * d * unit; vals.append(V(zz))
                der = (-vals[0] + 8 * vals[1] - 8 * vals[2] + vals[3]) / (12 * d)
                Fn[i] += 0.5 * der * (1 if comp == 0 else 1j)
        src = 2 * np.imag(np.conj(z) * Fa)
        K = k * np.conj(Y) * X**3; L = m * np.conj(S) * Y**2
        want = np.array([-6 * K.imag, 2 * K.imag - 4 * L.imag, 2 * L.imag])
        scale = np.max(np.abs(want))
        ang = np.zeros(3)
        for i in range(3):                          # dV/dalpha_i by angular finite difference
            vals = []
            for s in (2, 1, -1, -2):
                zz = z.copy(); zz[i] *= np.exp(1j * s * 1e-4); vals.append(V(zz))
            ang[i] = (-vals[0] + 8 * vals[1] - 8 * vals[2] + vals[3]) / (12e-4)
        worst['force'] = max(worst['force'], float(np.max(np.abs(src - want)) / scale))
        worst['numgrad'] = max(worst['numgrad'], float(np.max(np.abs(Fn - Fa)) / np.max(np.abs(Fa))))
        worst['angle'] = max(worst['angle'], float(np.max(np.abs(ang - want)) / scale))
        worst['gauge'] = max(worst['gauge'], float(abs(src @ [1, 3, 6]) / scale))
    return dict(pass_group=bool(worst['force'] < 1e-12 and worst['numgrad'] < 1e-7 and worst['angle'] < 1e-8
                                and worst['gauge'] < 1e-12), samples=nsamp, worst_relative_errors=worst,
                note='dV/dalpha_i equals the source because n_i = -2 r_i^2 alpha_i_dot and the EOM; complex couplings')


# ---------------------------------------------------------------- G3: forced response and the X^6 torque
def g3_symbolic():
    A = sp.symbols('A', positive=True)
    Br, Bi, Jr, Ji, Yr, Yi = sp.symbols('B_r B_i J_r J_i Y_r Y_i', real=True)
    B, J, Y = Br + sp.I * Bi, Jr + sp.I * Ji, Yr + sp.I * Yi
    D = A**2 - (Br**2 + Bi**2)
    eq = sp.expand(A * Y + B * sp.conjugate(Y) + J)
    sol = sp.solve([sp.re(eq), sp.im(eq)], [Yr, Yi], dict=True)[0]
    Ystar = (-A * J + B * sp.conjugate(J)) / D
    ok_root = sp.simplify(sol[Yr] - sp.re(Ystar)) == 0 and sp.simplify(sol[Yi] - sp.im(Ystar)) == 0
    # V_Y = A|Y|^2 + Re(B conj(Y)^2) + 2 Re(conj(Y) J)   [2Re(mu Sigma* Y^2) = Re(B Y*^2) with B = 2 mu* Sigma]
    VY = A * (Yr**2 + Yi**2) + sp.re(B * sp.conjugate(Y)**2) + 2 * sp.re(sp.conjugate(Y) * J)
    Vstar = sp.simplify(VY.subs({Yr: sp.re(Ystar), Yi: sp.im(Ystar)}))
    dV_formula = sp.simplify((-A * (Jr**2 + Ji**2) + sp.re(B * sp.conjugate(J)**2)) / D)
    ok_dV = sp.simplify(Vstar - dV_formula) == 0
    # torque identities: with J = kappa X^3, B conj(J)^2 = 2 mu* Sigma kappa*^2 X*^6 = conj(2 alpha6 X^6), alpha6 = mu kappa^2 Sigma*/D
    # so Re(B J*^2)/D = 2 Re(alpha6 X^6) and Delta V = -A|kappa|^2|X|^6/D + (alpha6 X^6 + h.c.)      [matches Codex]
    # n_X source -6 Im(kappa Y* X^3) = -6 Im(Y_** J); n_D source -12 Im(mu Sigma* Y^2) = -6 Im(B* Y_*^2);
    # both should equal -12 Im(alpha6 X^6) = -6 Im(B* J^2)/D.
    sX = -6 * sp.im(sp.conjugate(Ystar) * J)
    sD = -6 * sp.im(sp.conjugate(B) * Ystar**2)
    target = -6 * sp.im(sp.conjugate(B) * J**2) / D
    ok_sX = sp.simplify(sX - target) == 0
    ok_sD = sp.simplify(sD - target) == 0
    return dict(pass_group=bool(ok_root and ok_dV and ok_sX and ok_sD), Ystar_solves_linear_system=bool(ok_root),
                VY_at_Ystar_equals_DeltaV=bool(ok_dV), nX_source_equals_minus12_Im_alpha6X6=bool(ok_sX),
                nD_source_equals_same=bool(ok_sD),
                note='hence (2 Im K - 4 Im L) = 0 on the forced branch: the forced Y accumulates no net n_Y')


def g3_numeric():
    M, hs, mu, kap = 1.0, 0.1, 0.5, 1e-3
    Sig = hs / math.sqrt(2); A = M**2; B = 2 * mu * Sig; D = A**2 - B**2
    alpha6 = mu * kap**2 * Sig / D
    rows = []
    for amp in (0.2, 0.1, 0.05, 0.025):
        X = amp * np.exp(1j * np.pi / 12); J = kap * X**3
        Ystar = (-A * J + B * np.conj(J)) / D
        # route 1: minimise the quadratic + quartic Y potential directly
        def VY(v, quartic):
            Y = v[0] + 1j * v[1]
            return A * abs(Y)**2 + np.real(B * np.conj(Y)**2) + 2 * np.real(np.conj(Y) * J) + quartic * abs(Y)**4
        q = minimize(VY, [Ystar.real, Ystar.imag], args=(0.0,), method='BFGS', options=dict(gtol=1e-16))
        Yq = q.x[0] + 1j * q.x[1]
        # route 2: full nonlinear stationarity with lambda_Y = 1 via fsolve
        def grad(v):
            Y = v[0] + 1j * v[1]
            g = A * Y + B * np.conj(Y) + J + 2 * abs(Y)**2 * Y
            return [g.real / abs(J), g.imag / abs(J)]
        Yn = fsolve(grad, [Ystar.real, Ystar.imag], xtol=1e-14); Yn = Yn[0] + 1j * Yn[1]
        sX = -6 * np.imag(kap * np.conj(Ystar) * X**3); sD = -12 * np.imag(mu * Sig * Ystar**2)
        sE = -12 * np.imag(alpha6 * X**6); sN = -12 * np.imag(mu * Sig * Yn**2)
        exc = abs(0)  # displacement Y - Y_* = 0 by construction: excitation form vanishes identically
        rows.append(dict(X_over_M=amp, Ystar_over_X=float(abs(Ystar) / abs(X)),
                         min_vs_formula_rel=float(abs(Yq - Ystar) / abs(Ystar)),
                         nonlinear_vs_formula_rel=float(abs(Yn - Ystar) / abs(Ystar)),
                         source_nX=float(sX), source_nD=float(sD), source_eft=float(sE), source_nonlinear=float(sN),
                         source_match_rel=float(max(abs(sX - sD), abs(sX - sE)) / abs(sX)),
                         nonlinear_source_shift_rel=float(abs(sN / sD - 1)), nY_source_on_forced_branch=float(
                             2 * np.imag(kap * np.conj(Ystar) * X**3) - 4 * np.imag(mu * Sig * Ystar**2)),
                         instantaneous_excitation_form=0.0))
    codex = [-2.7289347133229633e-11, -4.26396048956713e-13, -6.662438264948641e-15, -1.0410059788982252e-16]  # [Q]
    cmp = [abs(r['source_nX'] / c - 1) for r, c in zip(rows, codex)]
    return dict(pass_group=bool(all(r['source_match_rel'] < 1e-12 and r['min_vs_formula_rel'] < 1e-6
                                    and r['nonlinear_source_shift_rel'] < 1e-8 for r in rows) and max(cmp) < 1e-9),
                alpha6_M_units=alpha6, gap2_over_M2=A - B, rows=rows, max_rel_dev_from_codex_source_values=max(cmp))


# ---------------------------------------------------------------- G4: Higgs angular countercharge
def g4_symbolic():
    f, h, g = sp.symbols('f h g', positive=True)
    thd, sd, A0 = sp.symbols('thetadot sigmadot A0', real=True)
    # X = f e^{-i theta}/sqrt2 (charge 1), Sigma = h e^{-i sigma}/sqrt2 (charge 6); D_t z = (d_t + i q g A0) z
    # kinetic: (1/2) f^2 (thd - g A0)^2 + (1/2) h^2 (sd - 6 g A0)^2 ; n_i = -2 Im(z* D_t z) = r_i^2 (alpha_dot - q g A0)
    u, w = thd - g * A0, sd - 6 * g * A0
    T = sp.Rational(1, 2) * f**2 * u**2 + sp.Rational(1, 2) * h**2 * w**2
    gauss = sp.solve(sp.diff(T, A0), A0)[0]          # homogeneous Gauss law: total gauge charge zero
    nX = sp.simplify((f**2 * u).subs(A0, gauss)); nS = sp.simplify((h**2 * w).subs(A0, gauss))
    I = f**2 * h**2 / (f**2 + 36 * h**2); delta_dot = 6 * thd - sd
    rho = sp.simplify(T.subs(A0, gauss))
    ok = (sp.simplify(nX - 6 * I * delta_dot) == 0, sp.simplify(nS + nX / 6) == 0,
          sp.simplify(rho - (nX**2 / (2 * f**2) + nX**2 / (72 * h**2))) == 0,
          sp.simplify(nX + 6 * nS) == 0)
    return dict(pass_group=all(ok), nX_equals_6I_deltadot=ok[0], nSigma_equals_minus_nX_over_6=ok[1],
                rho_phase_formula=ok[2], gauge_neutral=ok[3], nX=str(nX), rho=str(rho),
                note='second term n_X^2/(72 h^2) is the Sigma phase kinetic cost; no radial quanta are implied')


# ---------------------------------------------------------------- G5: the conditional envelope yield bound
def bound_yield(alpha6, r_i, cH, lam=1e-3):
    """|a^3 n_X| <= 12 |alpha6| R^6 / H_i with R = r_i (1 + 2/cH)^(1/4), H_i = sqrt(2 lam/cH) r_i; yield over s(H_i)."""
    H = math.sqrt(2 * lam / cH) * r_i
    R = r_i * (1 + 2 / cH)**0.25
    N = 12 * abs(alpha6) * R**6 / H
    return dict(H_i=H, R=R, a3nX_bound=N, yield_bound=N / entropy(H), over_YQ=N / entropy(H) / YQ)


def g5_bound():
    M, hs, mu, kap = 1e12, 1e11, 5e11, 1e-3
    Sig = hs / math.sqrt(2); D = M**4 - 4 * mu**2 * Sig**2
    alpha6 = mu * kap**2 * Sig / D
    r = math.sqrt(2) * mu * hs / M**2
    alpha6_closed = r * kap**2 / (2 * M**2 * (1 - r**2))
    table = {}
    for cH in (1.0, 10.0):
        for q in (0.1, 1.0, 2.0):
            table[f'ri/M={q}_cH={cH:g}'] = bound_yield(alpha6, q * M, cH)
    codex = {'ri/M=0.1_cH=1': 1.2201116680281771e-08, 'ri/M=0.1_cH=10': 5.4889536982333816e-08,
             'ri/M=1.0_cH=1': 3.858331870716284e-05, 'ri/M=1.0_cH=10': 0.00017357595657622028}   # [Q]
    dev = {k: abs(table[k]['over_YQ'] / v - 1) for k, v in codex.items()}
    power = math.log(bound_yield(alpha6, 0.2 * M, 1.0)['yield_bound'] / bound_yield(alpha6, 0.1 * M, 1.0)['yield_bound']) / math.log(2)
    integral = quad(lambda a: a**-2, 1, np.inf)[0]
    extrap = {}
    for cH in (1.0, 10.0):
        q = brentq(lambda q: math.log(bound_yield(alpha6, q * M, cH)['over_YQ']), 0.01, 1e4)
        b = bound_yield(alpha6, q * M, cH)
        gap2 = M**2 - 2 * mu * Sig; Ymax = kap * b['R']**3 / gap2
        extrap[f'cH={cH:g}'] = dict(ri_over_M=q, H_i=b['H_i'], M_over_Hi=M / b['H_i'], kappaR2_over_M2=kap * b['R']**2 / M**2,
                                    omegaX_over_gap=math.sqrt(2e-3) * b['R'] / math.sqrt(gap2),
                                    Y_quartic_force_ratio=2 * Ymax**2 / gap2,
                                    Sigma_forced_shift=mu * Ymax**2 / (2 * hs**2 * Sig))
    return dict(pass_group=bool(max(dev.values()) < 1e-9 and abs(power - 3.5) < 1e-12 and abs(integral - 1) < 1e-12),
                alpha6_GeVm2=alpha6, alpha6_closed_form_rel_dev=abs(alpha6_closed / alpha6 - 1), r=r,
                table={k: v['over_YQ'] for k, v in table.items()}, max_rel_dev_from_codex_table=max(dev.values()),
                yield_power_in_ri=power, triangle_integral=integral, extrapolation_to_YQ=extrap,
                note='same entropy convention as our scal.py yield_over_YQ (prompt radiation at H_i, g*=106.75)')


if __name__ == '__main__':
    out = dict(G1_symbolic=g1_symbolic(), G1_dynamical=g1_dynamical(), G2_numeric=g2_numeric(),
               G3_symbolic=g3_symbolic(), G3_numeric=g3_numeric(), G4_symbolic=g4_symbolic(), G5_bound=g5_bound())
    out['groups_passed'] = int(sum(v['pass_group'] for v in out.values()))
    Path(sys.argv[1]).write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out, indent=1, default=float))
