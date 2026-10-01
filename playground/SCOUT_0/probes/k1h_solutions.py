# K1-H steps 2-5: does the exact K1 constraint admit early-GR-restoring, stable, non-GR luminal Horndeski histories?
# Inputs: Linder (2020) mu, Sigma (PRIMARY-TEXT as supplied by the external audit); R' = -alpha_M R.
# Stability (RECONSTRUCTED, not primary-text-verified here): Bellini-Sawicki form with alpha_T = 0,
#   Dkin = alpha_K + (3/2) alpha_B^2 > 0 (alpha_K free -> satisfiable), and
#   N := Dkin*c_s^2 = (1/2)[(2-alpha_B)A + 2 alpha_B'] - (2-alpha_B) H'/H - rho_m/(H^2 M_*^2),
#   sign-checked: alpha_B = alpha_M = 0 gives N = phidot^2/(H^2 M^2) (quintessence c_s^2 = 1);
#   the bracket equals Linder's mu/Sigma denominator D. Background: LambdaCDM, Om0 = 0.3; rho_m/(H^2 M_*^2) = 3 Om(a) R.
# Requires numpy, scipy, sympy.
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp

print("=== early-time asymptotics: alpha_M = c a^n, alpha_B = b a^n, delta = R - 1 = -(c/n) a^n ===")
b, c, n = sp.symbols('b c n', real=True)
eq = sp.expand(-2*b*c - (-(-2*c/n + b)*(b + 2*c)))     # 2 delta aB' + (2 delta + aB)(aB + 2 aM) = 0 at O(a^{2n})
roots = sp.solve(sp.Eq(eq, 0), b)
print("  leading-order condition:", sp.factor(eq), "= 0  ->  b =", [sp.simplify(r) for r in roots])
for nv in (3, 1):
    print(f"  n={nv}: b/c =", [float(sp.N(r.subs({n: nv, c: 1}))) for r in roots],
          " (No Slip has b/c = -2, Only Run b/c = 0: both are distinct)")
Nlead = sp.expand((b + 2*c) + n*b - sp.Rational(3, 2)*b + 3*c/n)  # leading O(a^n) part of N in matter era
for nv in (3,):
    for r in roots:
        val = sp.simplify(Nlead.subs(b, r).subs(n, nv))
        print(f"  n={nv}, b={sp.simplify(r.subs(n, nv))}: leading N/a^n = {val}  -> stable iff c > 0" if True else "")

Om0 = 0.3
Om = lambda a: Om0*a**-3/(Om0*a**-3 + 1 - Om0)
OL = lambda a: 1 - Om(a)
dlnH = lambda a: -1.5*Om(a)                           # H'/H for LambdaCDM

def run(cM, branch, a0=1e-2, pert=0.0):
    aM = lambda a: cM*OL(a)/(1 - Om0)                   # alpha_M proportional to Omega_DE(a)/Omega_Lambda (n = 3 early)
    rts = [float(sp.N(r.subs({n: 3, c: 1}))) for r in roots]
    bc = rts[branch]
    x0 = np.log(a0)
    delta0 = np.expm1(-cM*np.log(1 + (1-Om0)/Om0*a0**3)/(3*(1-Om0)))   # exact integral of alpha_M for this choice
    aB0 = bc*aM(a0)*(1 + pert)
    def rhs(x, y):
        a = np.exp(x); d, B = y; M = aM(a)
        return [-M*(1 + d), -(2*d + B)*(B + 2*M)/(2*d)]
    sol = solve_ivp(rhs, (x0, 0.0), [delta0, aB0], rtol=1e-10, atol=1e-14, dense_output=True)
    xs = np.linspace(x0, 0.0, 400); out = []
    for x in xs:
        a = np.exp(x); d, B = sol.sol(x); M = aM(a); R = 1 + d
        Bp = rhs(x, [d, B])[1]; A = B + 2*M; D = (2 - B)*A + 2*Bp
        mu = R*((2 + 2*M)*A + 2*Bp)/D; Sg = R*((2 + M)*A + 2*Bp)/D
        N = 0.5*D - (2 - B)*dlnH(a) - 3*Om(a)*R
        out.append((a, M, B, R, mu, Sg, N))
    return np.array(out), sol.status

for cM in (0.1, 0.3):
    for br in (0, 1):
        arr, st = run(cM, br)
        a, M, B, R, mu, Sg, N = arr.T
        print(f"\n  c_M={cM}, branch b/c={[float(sp.N(r.subs({n: 3, c: 1}))) for r in roots][br]:+.3f}: solver status {st}")
        print(f"    max|2Sigma-mu-1| = {np.max(np.abs(2*Sg - mu - 1)):.1e}   min N (needs > 0) = {N.min():+.3e}   "
              f"z=0: alpha_M={M[-1]:.3f} alpha_B={B[-1]:+.3f} R={R[-1]:.4f} mu-1={mu[-1]-1:+.4f} Sigma-1={Sg[-1]-1:+.4f}")
        for ai in (0.1, 0.5, 1.0):
            k = np.argmin(abs(a - ai))
            print(f"      a={a[k]:.2f}: alpha_B/alpha_M={B[k]/M[k]:+.3f}  mu-1={mu[k]-1:+.5f}  Sigma-1={Sg[k]-1:+.5f}  N={N[k]:+.4f}")

print("\n=== attractor test: perturb the initial alpha_B by +-5% (branch with b/c > 0, c_M = 0.3) ===")
base, _ = run(0.3, 0)
for p in (0.05, -0.05):
    arr, st = run(0.3, 0, pert=p)
    print(f"  pert {p:+.2f}: status {st}; alpha_B(z=0) = {arr[-1,2]:+.4f} vs {base[-1,2]:+.4f}; max|2Sigma-mu-1| = {np.max(np.abs(2*arr[:,5]-arr[:,4]-1)):.1e}")
print("  (K1 holds identically along ANY solution of the constraint ODE: a one-parameter family per alpha_M history)")
