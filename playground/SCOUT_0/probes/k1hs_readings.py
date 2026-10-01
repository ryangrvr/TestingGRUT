# K1-HS follow-up: (1) exact identities on K1 from Linder's mu, Sigma; (2) stability under two readings of the denominator:
#   reading BS: alpha c_s^2 = N = D/2 - (2-aB)H'/H - rho_m/(H^2 M_*^2)   (published Horndeski eq., verified)
#   reading L : alpha c_s^2 = D/2  (the identification implied by mu = R(1 + A^2/D) having the standard QSA form
#               mu = R(1 + A^2/(2 alpha c_s^2)); OWED: check against Linder's text which reading his Eqs. assume).
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
aB, aM, aBp, R = sp.symbols('alpha_B alpha_M alpha_Bp R')
A = aB + 2*aM; D = (2 - aB)*A + 2*aBp
mu = R*((2 + 2*aM)*A + 2*aBp)/D; Sg = R*((2 + aM)*A + 2*aBp)/D
print("mu - R(1 + A^2/D) =", sp.simplify(mu - R*(1 + A**2/D)), ";  Sigma - R(1 + (aM+aB)A/D) =", sp.simplify(Sg - R*(1 + (aM + aB)*A/D)))
aBp_K1 = sp.solve(2*(R - 1)*aBp + (2*(R - 1) + aB)*A, aBp)[0]
print("on K1: D =", sp.factor(D.subs(aBp, aBp_K1)), "   mu - 1 =", sp.factor(sp.simplify(mu.subs(aBp, aBp_K1) - 1)),
      "   Sigma - 1 =", sp.factor(sp.simplify(Sg.subs(aBp, aBp_K1) - 1)))
print("=> on K1 (R<1 for alpha_M>0 histories): sign(mu-1) = sign(alpha_M/alpha_B); reading L gives alpha c_s^2 = alpha_B A R/(2(1-R)),")
print("   i.e. reading-L gradient stability <=> alpha_B (alpha_B + 2 alpha_M) > 0  (strong side alpha_B>0: always; weak side needs alpha_B < -2 alpha_M)")

Om0 = 0.3
Om = lambda a: Om0*a**-3/(Om0*a**-3 + 1 - Om0)
def hist(kind, c):
    if kind == "a^1":   return (lambda a: c*a), (lambda a: np.exp(-c*a)), 1.0, c
    if kind == "a^1.5": return (lambda a: c*a**1.5), (lambda a: np.exp(-c*a**1.5/1.5)), 1.5, c
    if kind == "a^2/(1+a^2)": return (lambda a: c*a*a/(1 + a*a)), (lambda a: (1 + a*a)**(-c/2)), 2.0, c
    rr = (1 - Om0)/Om0
    return (lambda a: c*a**3/(Om0 + (1 - Om0)*a**3)), (lambda a: (1 + rr*a**3)**(-c/(3*(1 - Om0)))), 3.0, c/Om0
def traj(kind, c, B1, a_end=1e-4):
    aMa, Ra, n, c_eff = hist(kind, c)
    def rhs(x, y):
        a = np.exp(x); d = Ra(a) - 1; M = aMa(a); B = y[0]
        return [-(2*d + B)*(B + 2*M)/(2*d)]
    s = solve_ivp(rhs, (0.0, np.log(a_end)), [B1], rtol=1e-11, atol=1e-15, method='DOP853', dense_output=True)
    if s.status != 0: return None
    rows = []
    for x in np.linspace(np.log(a_end), 0, 500):
        a = np.exp(x); B = s.sol(x)[0]; M = aMa(a); Rv = Ra(a); Bp = rhs(x, [B])[0]
        Av = B + 2*M; Dv = (2 - B)*Av + 2*Bp
        muv = Rv*(1 + Av**2/Dv); Sv = Rv*(1 + (M + B)*Av/Dv)
        N = 0.5*Dv + (2 - B)*1.5*Om(a) - 3*Om(a)*Rv
        rows.append((a, M, B, Rv, muv, Sv, N, 0.5*Dv))
    r = np.array(rows)
    n_ = n; strong = (1 + np.sqrt(4*n_ + 1))/n_
    return r, r[0, 2]/r[0, 1], strong
print("\n=== backward-shooting family, strong-branch origin, alpha_B(1) > 0 (mu(1) > 1 by the K1 identity) ===")
for kind, cs in (("OmegaDE", (0.1, 0.3)), ("a^1", (0.05, 0.2)), ("a^1.5", (0.05, 0.2)), ("a^2/(1+a^2)", (0.05, 0.2))):
    for c in cs:
        good_BS, good_L, n_tot = [], [], 0
        for B1 in np.linspace(0.05, 1.5, 30):
            out = traj(kind, c, B1)
            if out is None: continue
            r, ratio0, strong = out; n_tot += 1
            if abs(ratio0 - strong) > 0.02*abs(strong) or r[:, 2].min() <= 0: continue
            m1 = r[-1, 4] - 1; okBS = r[:, 6].min() > 0; okL = r[:, 7].min() > 0
            if okBS: good_BS.append((B1, m1))
            if okL: good_L.append((B1, m1))
        def rng(g): return f"{len(g)}/{n_tot}, mu(1)-1 in [{min(v[1] for v in g):+.4f}, {max(v[1] for v in g):+.4f}]" if g else f"0/{n_tot}"
        print(f"  {kind:12} c={c}: strong-origin, alpha_B>0 throughout, mu>1, EFT-stable members — reading BS: {rng(good_BS)};  reading L: {rng(good_L)}")

print("\n=== K1-H weak-branch solutions (forward, exact weak start) under the two readings ===")
for kind, c in (("OmegaDE", 0.1), ("OmegaDE", 0.3), ("a^1", 0.2), ("a^1.5", 0.2)):
    aMa, Ra, n, c_eff = hist(kind, c)
    weak = (1 - np.sqrt(4*n + 1))/n
    def rhs(x, y):
        a = np.exp(x); d = Ra(a) - 1; M = aMa(a); B = y[0]
        return [-(2*d + B)*(B + 2*M)/(2*d)]
    a0 = 1e-2
    s = solve_ivp(rhs, (np.log(a0), 0), [weak*aMa(a0)], rtol=1e-11, atol=1e-15, dense_output=True)
    Ns, Ls, mus = [], [], []
    for x in np.linspace(np.log(a0), 0, 400):
        a = np.exp(x); B = s.sol(x)[0]; M = aMa(a); Rv = Ra(a); Bp = rhs(x, [B])[0]; Av = B + 2*M; Dv = (2 - B)*Av + 2*Bp
        Ns.append(0.5*Dv + (2 - B)*1.5*Om(a) - 3*Om(a)*Rv); Ls.append(0.5*Dv); mus.append(Rv*(1 + Av**2/Dv))
    print(f"  {kind:8} c={c}: weak branch mu(1)-1={mus[-1]-1:+.4f}  min N (BS) = {min(Ns):+.2e}   min D/2 (L) = {min(Ls):+.2e}"
          f"  -> {'stable' if min(Ns) > 0 else 'UNSTABLE'} (BS) / {'stable' if min(Ls) > 0 else 'UNSTABLE'} (L)")
