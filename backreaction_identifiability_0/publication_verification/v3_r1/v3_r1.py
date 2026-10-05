#!/usr/bin/env python3
"""
V3-R1: Deterministic quadrature remediation.
Preregistered: t_star=0.5, diagnostic t={0.25,0.75,1.0}, NB grid {4,8,16,32,64,128}.
Quadrature: Gauss-Legendre in x over [-Lx, Lx], Gauss-Legendre in p over [-Lp, Lp],
with Lx = 6, Lp = 6 (tail mass < 1e-8 in each; verified by trapezoid extension).
Integrator: RK4 with step refinement factors {1, 2, 4} from base dt=1e-3.
All choices frozen BEFORE evaluation (see report).
"""
import numpy as np, json, os

# ---- PREREGISTERED ----
T_STAR = 0.5
T_DIAG = [0.25, 0.75, 1.0]
NB_GRID = [4, 8, 16, 32, 64, 128]
LX, LP = 6.0, 6.0
NX_LIST = [120, 240, 480]   # x-resolution ladder
NP_LIST = [120, 240, 480]   # p-resolution ladder (matched)
DT_BASE = 1e-3
REFINE = [1, 2, 4]

def s(u): return 10*u**3 - 15*u**4 + 6*u**5
def q_p1(t): return s(t/np.pi) if t <= np.pi else 1.0

def gl_nodes(n, a, b):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5*(b-a)*x + 0.5*(a+b), 0.5*(b-a)*w

def duffing_integral(x0, p0, T, qfun, dt):
    """RK4 for x'' = -x - x^3 + eps*q(t) with eps=q-coupling absorbed via X^eps."""
    # We track (x, v) at time T. Return x(T).
    steps = max(1, int(round(T/dt)))
    h = T/steps
    x = np.full_like(x0, np.nan); x[:] = x0
    v = p0.copy()
    for k in range(steps):
        t0 = k*h; tm = t0 + 0.5*h; t1 = t0 + h
        q0, qm, q1 = qfun(t0), qfun(tm), qfun(t1)
        k1v = -x - x**3 + q0;  k1x = v
        k2v = -(x+0.5*h*k1x) - (x+0.5*h*k1x)**3 + qm;  k2x = v + 0.5*h*k1v
        k3v = -(x+0.5*h*k2x) - (x+0.5*h*k2x)**3 + qm;  k3x = v + 0.5*h*k2v
        k4v = -(x+h*k3x) - (x+h*k3x)**3 + q1;          k4x = v + h*k3v
        x = x + (h/6)*(k1x + 2*k2x + 2*k3x + k4x)
        v = v + (h/6)*(k1v + 2*k2v + 2*k3v + k4v)
    return x

def integrate_eps(x0, p0, T, eps, qfun, dt):
    """Full driven oscillator: x'' = -x - x^3 + eps*q(t)."""
    steps = max(1, int(round(T/dt)))
    h = T/steps
    x = x0.copy(); v = p0.copy()
    for k in range(steps):
        t0=k*h; tm=t0+0.5*h; t1=t0+h
        q0,qm,q1 = qfun(t0), qfun(tm), qfun(t1)
        def dv(x,v,q): return -x - x**3 + eps*q
        k1v = dv(x,v,q0);          k1x = v
        k2v = dv(x+0.5*h*k1x, v+0.5*h*k1v, qm); k2x = v+0.5*h*k1v
        k3v = dv(x+0.5*h*k2x, v+0.5*h*k2v, qm); k3x = v+0.5*h*k2v
        k4v = dv(x+h*k3x, v+h*k3v, q1);          k4x = v+h*k3v
        x = x + (h/6)*(k1x+2*k2x+2*k3x+k4x)
        v = v + (h/6)*(k1v+2*k2v+2*k3v+k4v)
    return x

def moments_over_grid(vals, wx, wp):
    """vals: array shape (NX, NP). Returns E[f] with weights wx (x) outer wp (p)."""
    W = np.outer(wx, wp)
    Z = W.sum()
    return (vals*W).sum()/Z

def run_resolution(nx, np_, dt):
    xs, wx = gl_nodes(nx, -LX, LX)
    ps, wp = gl_nodes(np_, -LP, LP)
    W = np.outer(wx, wp)
    Z = W.sum()
    # normalize (quartic Gibbs over truncated domain)
    rho = np.exp(-0.5*ps[None,:]**2 - 0.5*xs[:,None]**2 - 0.25*xs[:,None]**4)
    Wn = W * rho
    Zn = Wn.sum()
    Wn = Wn / Zn
    out = {}
    for T in [T_STAR] + T_DIAG:
        # P0: unforced -> X^eps = x0(T) independent of eps
        x0T = duffing_integral(xs, ps, T, lambda t: 0.0, dt)
        # grid: x0T[i,j] = x0(T; xs[i], ps[j])
        m0  = (x0T*Wn).sum()
        m2_0 = (x0T**2*Wn).sum()
        m3_0 = (x0T**3*Wn).sum() - 3*m0*m2_0 + 2*m0**3
        # P1: eps-dependent; we need kappa3(X^eps) for each eps (i.e. each NB)
        k3_eps = {}
        var_eps = {}
        for nb in NB_GRID:
            eps = 1.0/np.sqrt(nb)
            xT = integrate_eps(xs, ps, T, eps, q_p1, dt)
            m  = (xT*Wn).sum()
            m2 = (xT**2*Wn).sum()
            m3 = (xT**3*Wn).sum()
            # central 3rd moment: E[(X-m)^3] = m3 - 3*m*m2 + 2*m^3
            k3 = m3 - 3*m*m2 + 2*m**3
            var = m2 - m*m
            k3_eps[nb] = k3
            var_eps[nb] = var
        out[T] = dict(m0=m0, var0=m2_0-m0*m0, k3_0=m3_0,
                      k3_eps=k3_eps, var_eps=var_eps)
    return out

def main():
    results = {}
    for nx in NX_LIST:
        npn = NX_LIST[NX_LIST.index(nx)] if nx in NX_LIST else None
        # match p resolution
        for r, dt in enumerate([DT_BASE/f for f in REFINE]):
            res = run_resolution(nx, nx, dt)
            results[f"nx{nx}_dt{dt}"] = res
            print(f"--- nx={nx}, dt={dt} ---")
            for T in [T_STAR] + T_DIAG:
                d = res[T]
                print(f" t={T}: P0: m={d['m0']:.6e}, var={d['var0']:.6e}, k3={d['k3_0']:.3e}")
                gam = {nb: d['k3_eps'][nb]/d['var_eps'][nb]**1.5 * (1.0/np.sqrt(nb))**0 for nb in NB_GRID}
                # gamma1(F) = kappa3(F) / Var(F)^{3/2}
                # kappa3(F) = NB^{-1/2} * kappa3(X^eps)   [exact]
                # Var(F)    = Var(X^eps)                   [exact]
                g1 = {nb: (d['k3_eps'][nb]/np.sqrt(nb)) / d['var_eps'][nb]**1.5 for nb in NB_GRID}
                print(f"        P1 gamma1: " + " ".join(f"{nb}:{g1[nb]:+.3e}" for nb in NB_GRID))
    with open("v3_r1_results.json","w") as f:
        json.dump(results, f, indent=1, default=str)
    print("saved v3_r1_results.json")

if __name__ == "__main__":
    main()
