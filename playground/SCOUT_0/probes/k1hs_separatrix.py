# K1-HS: does the early-GR-restoring mu>1 K1 branch admit a globally regular, EFT-stable continuation to a = 1?
# Exact K1 ODE (Linder 2020, PRIMARY-TEXT): 2(R-1) aB' + [2(R-1) + aB](aB + 2 aM) = 0,  R' = -aM R, ' = d/dln a.
# EFT stability (VERIFIED vs published Horndeski eqs; Peirone et al. PRD 97 043519 (2018) App. Eq. A18, alpha_T = 0):
#   alpha = alpha_K + 3 aB^2/2 > 0 (alpha_K free), N = alpha c_s^2 = D/2 - (2 - aB) H'/H - rho_m/(H^2 M_*^2) > 0.
# Background LambdaCDM (Om0 = 0.3). Two independent constructions:
#   (I)  high-order analytic series about a = 0 (t = a^n) + forward integration, a0 -> 0 and order -> high;
#   (II) backward shooting from a = 1 (backward integration is ODE-attracting toward the early branch).
# Requires numpy, scipy, sympy, mpmath.
import numpy as np, sympy as sp, mpmath as mp
from scipy.integrate import solve_ivp, quad

Om0 = 0.3
Om = lambda a: Om0*a**-3/(Om0*a**-3 + 1 - Om0)
mp.mp.dps = 40

# alpha_M histories: (label, n, aM(t) as sympy in t = a^n, R(t) exact as sympy)
t = sp.symbols('t', positive=True)
cS = sp.symbols('c', positive=True)
rr = (1 - Om0)/Om0
HIST = {
    "OmegaDE-tracking": (3, cS*t/(Om0 + (1 - Om0)*t), (1 + rr*t)**(-cS/(3*(1 - Om0)))),
    "a^1":              (1, cS*t, sp.exp(-cS*t)),
    "a^1.5":            (sp.Rational(3, 2), cS*t, sp.exp(-cS*t/sp.Rational(3, 2))),
    "a^2/(1+a^2)":      (2, cS*t/(1 + t), (1 + t)**(-cS/2)),
}

def series_coeffs(n, aM, R, cval, order):
    """Analytic (strong-branch) series aB = sum_{k>=1} b_k t^k solving the exact ODE order by order in t = a^n."""
    aMv, Rv = aM.subs(cS, cval), R.subs(cS, cval)
    bs = sp.symbols('b1:%d' % (order + 1))
    aB = sum(bs[k]*t**(k + 1) for k in range(order))
    expr = sp.expand(sp.series(2*(Rv - 1)*n*t*sp.diff(aB, t) + (2*(Rv - 1) + aB)*(aB + 2*aMv), t, 0, order + 2).removeO())
    sol = {}
    c2 = sp.Poly(expr, t).coeff_monomial(t**2)
    roots = sp.solve(c2, bs[0]); b1 = max(roots, key=lambda r: float(r))     # strong branch = larger root
    sol[bs[0]] = b1
    for k in range(1, order):
        ck = sp.Poly(sp.expand(expr.subs(sol)), t).coeff_monomial(t**(k + 2))
        sol[bs[k]] = sp.solve(ck, bs[k])[0]
    return [float(sol[b]) for b in bs], [float(r) for r in roots]

def make_rhs(n, aM, R, cval):
    fM = sp.lambdify(t, aM.subs(cS, cval)); fR = sp.lambdify(t, R.subs(cS, cval))
    def aMa(a): return float(fM(a**float(n)))
    def Ra(a): return float(fR(a**float(n)))
    def rhs(x, y):
        a = np.exp(x); B = y[0]; d = Ra(a) - 1; M = aMa(a)
        return [-(2*d + B)*(B + 2*M)/(2*d)]
    return rhs, aMa, Ra

def diagnostics(sol_fun, rhs, aMa, Ra, xs):
    out = []
    for x in xs:
        a = np.exp(x); B = sol_fun(x); M = aMa(a); R = Ra(a)
        Bp = rhs(x, [B])[0]; A = B + 2*M; D = (2 - B)*A + 2*Bp
        mu = R*((2 + 2*M)*A + 2*Bp)/D; Sg = R*((2 + M)*A + 2*Bp)/D
        N = 0.5*D + (2 - B)*1.5*Om(a) - 3*Om(a)*R
        out.append((a, M, B, R, mu, Sg, N))
    return np.array(out)

for label, (n, aM, R) in HIST.items():
    for cval in ((0.1, 0.3) if label == "OmegaDE-tracking" else (0.05, 0.2)):
        rhs, aMa, Ra = make_rhs(n, aM, R, cval)
        print(f"\n===== {label}, c = {cval} =====")
        # (I) series + forward integration, convergence in a0 and order
        coeffs, roots = series_coeffs(n, aM, R, cval, 5)
        print(f"  leading roots b1/c candidates: {[round(r/cval, 4) for r in roots]}; strong-branch series coeffs b1..b5 = {[f'{v:.4g}' for v in coeffs]}")
        res = {}
        for a0 in (1e-1, 3e-2, 1e-2):
            for order in (1, 3, 5):
                tt = a0**float(n); B0 = sum(coeffs[k]*tt**(k + 1) for k in range(order))
                s = solve_ivp(rhs, (np.log(a0), 0.0), [B0], rtol=1e-12, atol=1e-16, method='DOP853', dense_output=True)
                if s.status != 0:
                    res[(a0, order)] = f"blow-up at a={np.exp(s.t[-1]):.3f}"
                else:
                    d = diagnostics(lambda x: s.sol(x)[0], rhs, aMa, Ra, np.linspace(np.log(a0), 0, 300))
                    res[(a0, order)] = f"aB(1)={d[-1,2]:+.5f} mu-1={d[-1,4]-1:+.4f} minN={d[:,6].min():+.2e}"
        for k, v in res.items(): print(f"  (I) forward from a0={k[0]:.0e}, series order {k[1]}: {v}")
        # (II) backward shooting from a = 1: family labelled by aB(1)
        a_end = 1e-4
        fam = []
        for B1 in np.linspace(-1.0, 1.5, 51):
            s = solve_ivp(rhs, (0.0, np.log(a_end)), [B1], rtol=1e-11, atol=1e-15, method='DOP853', dense_output=True)
            if s.status != 0 or not np.all(np.isfinite(s.y)):
                fam.append((B1, None)); continue
            xs = np.linspace(np.log(a_end), 0.0, 400)
            d = diagnostics(lambda x: s.sol(x)[0], rhs, aMa, Ra, xs)
            ratio0 = d[0, 2]/d[0, 1]
            fam.append((B1, (ratio0, d[-1, 4] - 1, d[-1, 5] - 1, d[:, 6].min(), np.max(np.abs(2*d[:, 5] - d[:, 4] - 1)), d[:, 3].min(), np.max(np.abs(d[:, 4])))))
        b_strong, b_weak = max(roots)/cval, min(roots)/cval
        print(f"  (II) backward shooting, aB(1) grid -> early ratio aB/aM at a=1e-4 (strong={b_strong:+.3f}, weak={b_weak:+.3f}):")
        viable_strong = []
        for B1, info in fam[::5]:
            if info is None: print(f"     aB(1)={B1:+.2f}: backward integration singular (not GR-restoring)"); continue
            r0, m1, s1, mN, k1, Rmin, mumax = info
            print(f"     aB(1)={B1:+.2f}: early ratio {r0:+.3f}  mu(1)-1={m1:+.4f} Sig(1)-1={s1:+.4f}  minN={mN:+.2e}  K1res={k1:.0e}  minR={Rmin:.3f}")
        for B1, info in fam:
            if info is None: continue
            r0, m1, s1, mN, k1, Rmin, mumax = info
            if abs(r0 - b_strong) < 0.05 and m1 > 0 and mN > 0 and Rmin > 0 and np.isfinite(mumax):
                viable_strong.append((B1, m1, s1, mN))
        if viable_strong:
            lo, hi = viable_strong[0], viable_strong[-1]
            print(f"  => GR-restoring (strong-branch origin), EFT-stable (N>0, R>0, alpha_K free), mu(1)>1 members: {len(viable_strong)} grid points,"
                  f" aB(1) in [{lo[0]:+.2f}, {hi[0]:+.2f}], mu(1)-1 in [{min(v[1] for v in viable_strong):+.4f}, {max(v[1] for v in viable_strong):+.4f}]")
        else:
            print("  => no grid member with strong-branch origin, mu(1)>1 and N>0 throughout")
