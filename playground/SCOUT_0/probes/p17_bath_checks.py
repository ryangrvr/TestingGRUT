# P-17 numerical/symbolic checks (requires numpy, scipy, sympy). See P17_RESULT.md.
# C0  Drude: friction kernel from J(w)/w is exponential (p=1 edge) -- REDISCOVERED-KNOWN; quantum symmetrized
#     noise kernel of the same bath is NOT CM at low T (Matsubara weights change sign) -> keep objects separate.
# C1  Pathwise identity: full autonomous Hamiltonian (system + finite harmonic bath) vs GLE driven by the
#     free-bath force F computed from the bath initial microstate (deterministic microstate control).
# C3  O(beta) retained-mean discriminator S-D = -12 beta int g(t-s) mu(s) v(s) ds (covariance-only), across
#     bath regimes; short-time order; convergence to the C-B value -12 beta T a t^2 in the white + overdamped limit.
import numpy as np, sympy as sp
from scipy.linalg import expm
from scipy.integrate import solve_ivp, quad

print("=== C0 Drude: J(w) = g0 w k^2/(w^2+k^2)  ->  gamma(t) = (2/pi) int_0^inf J(w)/w cos(wt) dw ===")
w, t, kap, g0 = sp.symbols('omega t kappa gamma0', positive=True)
gam = sp.simplify(2/sp.pi*sp.integrate(g0*kap**2/(w**2+kap**2)*sp.cos(w*t), (w, 0, sp.oo)))
print("  gamma(t) =", gam, " (exponential = P-06b p=1 edge; classical noise C_F = T gamma by FDT)")
# quantum symmetrized force correlation (hbar = 1): C_sym(t) = (1/pi) int_0^inf J(w) coth(w/2T) cos(wt) dw.
# Exact Matsubara expansion (coth(x/2T) = 2T/x + sum_n 4T x/(x^2+nu_n^2), nu_n = 2 pi n T, partial fractions):
#   C_sym(t) = (g0 k^2/2) cot(k/2T) e^{-k t} + sum_{n>=1} [2 T g0 k^2 nu_n/(nu_n^2 - k^2)] e^{-nu_n t}.
# The exponentials are distinct, so the Bernstein measure is unique: any negative weight => NOT CM.
def matsubara(T, k=1.0, g=1.0, nmax=200000):
    n = np.arange(1, nmax+1); nu = 2*np.pi*n*T
    return g*k**2/2/np.tan(k/(2*T)), nu, 2*T*g*k**2*nu/(nu**2 - k**2)
def Csym(tv, T, k=1.0, g=1.0):
    c0, nu, cn = matsubara(T, k, g)
    return c0*np.exp(-k*tv) + np.sum(cn*np.exp(-nu*tv))
for T in (5.0, 0.3, 0.05):
    c0, nu, cn = matsubara(T)
    neg = [int(i+1) for i in np.where(cn < 0)[0][:3]] + (["c0"] if c0 < 0 else [])
    vals = [Csym(tv, T) for tv in (0.5, 2.0, 5.0, 10.0)]
    print(f"  T/kappa={T:4.2f}: negative Bernstein weights: {neg if neg else 'none'};  C_sym(t=0.5,2,5,10) =",
          ", ".join(f"{v:+.3e}" for v in vals), "->", "NOT CM" if neg else "CM (all weights > 0)")
print("  classical limit (T >> kappa): C_sym -> T gamma(t), CM;  low T (nu_1 = 2 pi T < kappa): Matsubara weights",
      "turn negative -> quantum noise kernel NOT CM though gamma(t) is exactly exponential")

print("\n=== C1 pathwise identity: autonomous Hamiltonian bath  ==  GLE with free-bath force F (one microstate) ===")
rng = np.random.default_rng(17)
Nb, kq, M, beta, a, T = 40, 2.3, 1.0, 1.0, 1.0, 0.7
def drude_modes(Nb, g=1.0, k=2.0, wmax=12.0):
    wj = (np.arange(Nb)+0.5)*wmax/Nb; dw = wmax/Nb
    J = g*wj*k**2/(wj**2+k**2); cj = np.sqrt(2/np.pi*J*wj*dw)      # gamma(t) = sum c^2/w^2 cos(w t)
    return wj, cj
def twoexp_modes(Nb, wmax=12.0):   # CM, non-Drude: J/w = 0.5*L(k=1) + 0.5*L(k=4)
    wj = (np.arange(Nb)+0.5)*wmax/Nb; dw = wmax/Nb
    Jw = 0.5*1/(wj**2+1) + 0.5*16/(wj**2+16); cj = np.sqrt(2/np.pi*Jw*wj**2*dw)
    return wj, cj
def osc_modes(Nb, wmax=12.0):      # non-CM: resonant bath, J/w peaked at w=3
    wj = (np.arange(Nb)+0.5)*wmax/Nb; dw = wmax/Nb
    Jw = 1/((wj-3)**2+0.25) + 1/((wj+3)**2+0.25); cj = np.sqrt(2/np.pi*0.3*Jw*wj**2*dw)
    return wj, cj
def run_pathwise(modes, label, tmax=6.0, nt=6000):
    wj, cj = modes
    xi = rng.normal(0, np.sqrt(T)/wj); eta = rng.normal(0, np.sqrt(T), size=Nb)   # ONE bath microstate
    def rhs(_, y):
        qv, pv, xs, ps = y[0], y[1], y[2:2+Nb], y[2+Nb:]
        dV = kq*qv + 4*beta*qv**3
        return np.concatenate(([pv/M, -dV + np.sum(cj*(xs - cj*qv/wj**2))], ps, -wj**2*(xs - cj*qv/wj**2)))
    y0 = np.concatenate(([a, 0.0], cj*a/wj**2 + xi, eta))
    ts = np.linspace(0, tmax, nt+1); dt = ts[1]
    sol = solve_ivp(rhs, (0, tmax), y0, t_eval=ts, rtol=1e-11, atol=1e-12, method='DOP853')
    qH = sol.y[0]
    F = (cj[None, :]*(xi[None, :]*np.cos(np.outer(ts, wj)) + eta[None, :]/wj*np.sin(np.outer(ts, wj)))).sum(1)
    gk = (cj**2/wj**2*np.cos(np.outer(ts, wj))).sum(1)                      # gamma(t)
    # GLE: M q'' = -V'(q) - int_0^t gamma(t-s) q'(s) ds + F(t), trapezoid memory + Heun step
    q, v = np.zeros(nt+1), np.zeros(nt+1); q[0] = a
    def acc(n, vn):   # acceleration at step n given v history up to n-1 and current vn
        hist = v[:n+1].copy(); hist[n] = vn
        if n == 0: mem = 0.0
        else:
            wts = np.ones(n+1); wts[0] = wts[-1] = 0.5
            mem = dt*np.dot(wts*gk[n::-1][:n+1], hist)
        return (-(kq*q[n] + 4*beta*q[n]**3) - mem + F[n])/M
    a0 = acc(0, v[0])
    for n in range(nt):
        qp, vp = q[n] + dt*v[n], v[n] + dt*a0
        q[n+1] = qp; an = acc(n+1, vp)
        q[n+1] = q[n] + dt*(v[n]+vp)/2; vc = v[n] + dt*(a0+an)/2
        q[n+1] = q[n] + dt*(v[n]+vc)/2; v[n+1] = vc; a0 = acc(n+1, vc)
    err = np.max(np.abs(qH - q))
    print(f"  {label:<26} max|q_H - q_GLE| on [0,{tmax}] = {err:.2e}  (|q| <= {np.max(np.abs(qH)):.2f}; GLE dt={dt:.0e}, O(dt^2))")
    return err
run_pathwise(drude_modes(Nb), "Drude (CM edge)")
run_pathwise(twoexp_modes(Nb), "two-exponential (CM)")
run_pathwise(osc_modes(Nb), "resonant (non-CM)")
print("  deterministic microstate: q_H(t) is a fixed trajectory; the 'exogenous' match is the Dirac law on that F path")

print("\n=== C3 O(beta) discriminator S-D(t) = -12 beta int_0^t g(t-s) mu(s) v(s) ds  (covariance-only) ===")
def disc(A, b, Dm, Sig0, x0, tv, ns=2000):
    # mu(s) = e^{As} x0 ; Sigma(s) via exact Van Loan-free quadrature on a grid ; g(tau) = e1^T e^{A tau} b
    s = np.linspace(0, tv, ns+1); ds = s[1]
    E = expm(A*ds); n = A.shape[0]
    Sig = Sig0.copy(); mu = x0.copy(); vs, mus = [Sig[0, 0]], [mu[0]]
    # exact one-step covariance increment: Q_ds = int_0^ds e^{Au} D e^{A^T u} du (Van Loan)
    big = np.block([[-A, Dm], [np.zeros((n, n)), A.T]])*ds
    EB = expm(big); Qd = EB[n:, n:].T @ EB[:n, n:]
    for _ in range(ns):
        Sig = E @ Sig @ E.T + Qd; mu = E @ mu; vs.append(Sig[0, 0]); mus.append(mu[0])
    gs = np.array([(expm(A*(tv-si)) @ b)[0] for si in s])
    integrand = gs*np.array(mus)*np.array(vs)
    return -12*beta*np.trapezoid(integrand, s)
k1, Tq, aq = 2.3, 1.0, 1.0
def overdamped_white():
    return np.array([[-k1]]), np.array([1.0]), np.array([[2*Tq]]), np.zeros((1, 1)), np.array([aq])
def inertial_white(Mm, g=1.0):
    A = np.array([[0, 1/Mm], [-k1, -g/Mm]]); return A, np.array([0, 1.0]), np.diag([0, 2*g*Tq]), np.zeros((2, 2)), np.array([aq, 0])
def inertial_drude(Mm, kk, g=1.0):
    # z = -int g kk e^{-kk(t-s)} qdot ds + F ; F stationary OU, Var = g kk T at t=0 (shifted-Gibbs preparation)
    A = np.array([[0, 1/Mm, 0], [-k1, 0, 1], [0, -g*kk/Mm, -kk]])
    return A, np.array([0, 1.0, 0]), np.diag([0, 0, 2*g*kk**2*Tq]), np.diag([0, 0, g*kk*Tq]), np.array([aq, 0, 0])
CB = lambda tv: disc(*overdamped_white(), tv)
print("  overdamped white (= C-B, 1 site, K11=2.3): S-D / (-12 beta T a t^2) at t=1e-3,1e-2:",
      ", ".join(f"{CB(tv)/(-12*beta*Tq*aq*tv**2):.4f}" for tv in (1e-3, 1e-2)), " (record: 1 as t->0)")
for Mm in (1.0,):
    tv = 1e-2; r = disc(*inertial_white(Mm), tv)/(-2*beta*aq*Tq*tv**5/(5*Mm**3))
    print(f"  inertial white  M={Mm}: leading order t^5, ratio to -2 beta a g0 T t^5/(5M^3) at t=1e-2: {r:.4f}")
for kk in (2.0,):
    tv = 1e-2; r = disc(*inertial_drude(1.0, kk), tv)/(-beta*aq*Tq*kk*tv**6/10)
    print(f"  inertial Drude  M=1, kappa={kk}: leading order t^6, ratio to -beta a T gamma(0) t^6/(10 M^3) at t=1e-2: {r:.4f}")
print("  limit to C-B at fixed t (Ohmic white, M -> 0; then Drude kappa -> inf): S-D(t) / C-B(t)")
for tv in (0.05, 0.2, 1.0):
    row = [disc(*inertial_white(Mm), tv)/CB(tv) for Mm in (1.0, 0.1, 0.01, 0.001)]
    rowd = disc(*inertial_drude(0.001, 1e4), tv)/CB(tv)
    print(f"    t={tv:4}: M=1,0.1,0.01,0.001 ->", ", ".join(f"{x:.4f}" for x in row), f" | Drude M=1e-3,kappa=1e4: {rowd:.4f}")
