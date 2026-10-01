# K1-H robustness: leading-order branch signs for general n, and a scan of alpha_M histories alpha_M = c * a^n
# (pure power law, plus the Omega_DE-tracking shape) on a LambdaCDM background; same reconstructed stability N.
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
b, c, n = sp.symbols('b c n', positive=False)
nn = sp.symbols('n', positive=True)
roots = [c*(1 - sp.sqrt(4*nn + 1))/nn, c*(1 + sp.sqrt(4*nn + 1))/nn]
print("=== leading order (alpha_M = c a^n, c > 0): sign of mu-1 and of N per branch ===")
for k, bb in enumerate(roots):
    mu1 = sp.simplify(-c/nn + (2*c + bb)**2/(2*(bb + 2*c + nn*bb)))
    N = sp.simplify(bb*(nn - sp.Rational(1, 2)) + c*(2 + 3/nn))
    vals = [(nv, float(mu1.subs({nn: nv, c: 1})), float(N.subs({nn: nv, c: 1}))) for nv in (0.5, 1, 1.5, 2, 3, 4, 6)]
    print(f"  branch {k} (b/c = {'(1-sqrt(4n+1))/n' if k == 0 else '(1+sqrt(4n+1))/n'}):",
          "  ".join(f"n={v[0]}: (mu-1)/(c a^n)={v[1]:+.3f}, N/(c a^n)={v[2]:+.3f}" for v in vals))
print("  c < 0 flips both signs of N at leading order -> early-time gradient instability for alpha_M < 0")

Om0 = 0.3
Om = lambda a: Om0*a**-3/(Om0*a**-3 + 1 - Om0)
def run(aM, aMint, bc, a0=1e-2, pert=0.0):
    def rhs(x, y):
        a = np.exp(x); d, B = y; M = aM(a)
        return [-M*(1 + d), -(2*d + B)*(B + 2*M)/(2*d)]
    y0 = [np.expm1(-aMint(a0)), bc*aM(a0)*(1 + pert)]
    sol = solve_ivp(rhs, (np.log(a0), 0.0), y0, rtol=1e-10, atol=1e-14, dense_output=True)
    if sol.status != 0: return None
    res = []
    for x in np.linspace(np.log(a0), 0, 300):
        a = np.exp(x); d, B = sol.sol(x); M = aM(a); R = 1 + d; Bp = rhs(x, [d, B])[1]
        A = B + 2*M; D = (2 - B)*A + 2*Bp
        mu = R*((2 + 2*M)*A + 2*Bp)/D; Sg = R*((2 + M)*A + 2*Bp)/D
        N = 0.5*D + (2 - B)*1.5*Om(a) - 3*Om(a)*R
        res.append((mu, Sg, N))
    return np.array(res)

print("\n=== numerical scan (start a0 = 0.01 on the leading-order branch, +-5% perturbations) ===")
for nv in (1.0, 1.5, 3.0):
    for cv in (0.05, 0.2, -0.05):
        aM = lambda a, c_=cv, n_=nv: c_*a**n_
        aMi = lambda a, c_=cv, n_=nv: c_*a**n_/n_
        for k in (0, 1):
            bc = float(roots[k].subs({nn: nv, c: 1}))
            outs = [run(aM, aMi, bc, pert=p) for p in (0.0, 0.05, -0.05)]
            tags = []
            for o in outs:
                if o is None: tags.append("BLOW-UP"); continue
                mu, Sg, N = o.T
                tags.append(f"mu-1={mu[-1]-1:+.4f} Sig-1={Sg[-1]-1:+.4f} K1res={np.max(abs(2*Sg-mu-1)):.0e} minN={N.min():+.1e}")
            print(f"  n={nv} c={cv:+.2f} branch {k} (b/c={bc:+.3f}): [0] {tags[0]} | [+5%] {tags[1]} | [-5%] {tags[2]}")
