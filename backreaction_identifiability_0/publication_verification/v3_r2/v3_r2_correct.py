#!/usr/bin/env python3
"""
V3-R2: Correct tensor quadrature with true 2-D evolution.
Every (x_i, p_j) pair evolved independently. P0 stationarity firewall first.
All scientific choices frozen from V3-R1 (no retuning).
"""
import numpy as np, json, time

T_STAR = 0.5
T_DIAG = [0.25, 0.75, 1.0]
NB_GRID = [4, 8, 16, 32, 64, 128]
L = 6.0
NX_LIST = [120, 240, 480]
DT_LIST = [1e-3, 5e-4, 2.5e-4]

# Reference values (high-precision, from independent 1-D quadrature)
REF_M2 = 0.467919916974
REF_M4 = 0.532080083026
REF_M6 = 0.871679667895

def s(u): return 10*u**3 - 15*u**4 + 6*u**5
def q_p1(t): return s(t/np.pi) if t <= np.pi else 1.0

def integrate_tensor(X0, P0, T, eps, qfun, dt):
    """Evolve full 2-D arrays. x'' = -x - x^3 + eps*q(t). Returns x(T)."""
    assert X0.ndim == 2 and X0.shape == P0.shape
    steps = max(1, int(round(T / dt)))
    h = T / steps
    x = X0.copy()
    v = P0.copy()
    for k in range(steps):
        t0 = k * h
        q0 = qfun(t0)
        # RK4 for coupled system (x, v)
        k1v = -x - x**3 + eps * q0
        k1x = v
        x2 = x + 0.5 * h * k1x
        v2 = v + 0.5 * h * k1v
        qm = qfun(t0 + 0.5 * h)
        k2v = -x2 - x2**3 + eps * qm
        k2x = v2
        x3 = x + 0.5 * h * k2x
        v3 = v + 0.5 * h * k2v
        k3v = -x3 - x3**3 + eps * qm
        k3x = v3
        x4 = x + h * k3x
        v4 = v + h * k3v
        q1 = qfun(t0 + h)
        k4v = -x4 - x4**3 + eps * q1
        k4x = v4
        x = x + (h / 6.0) * (k1x + 2*k2x + 2*k3x + k4x)
        v = v + (h / 6.0) * (k1v + 2*k2v + 2*k3v + k4v)
    return x

def moments(f, Wn):
    m1 = (f * Wn).sum()
    m2 = (f**2 * Wn).sum()
    m3 = (f**3 * Wn).sum()
    m4 = (f**4 * Wn).sum()
    mu = m1
    var = m2 - mu * mu
    k3 = m3 - 3*mu*m2 + 2*mu**3
    return m1, var, k3, m4

def run_grid(nx, dt):
    xg, wg = np.polynomial.legendre.leggauss(nx)
    xs = L * xg
    wx = L * wg
    # Correct 2-D tensor: every (xs[i], xs[j]) pair (using same nodes for x and p)
    X0 = np.broadcast_to(xs[:, None], (nx, nx)).copy()
    P0 = np.broadcast_to(xs[None, :], (nx, nx)).copy()
    assert X0.shape == (nx, nx) and P0.shape == (nx, nx)
    rho = np.exp(-0.5 * P0**2 - 0.5 * X0**2 - 0.25 * X0**4)
    W = wx[:, None] * wx[None, :]
    Wn_raw = W * rho
    Wn = Wn_raw / Wn_raw.sum()
    assert Wn.shape == (nx, nx)
    norm_raw = Wn_raw.sum()

    # --- P0 STATIONARITY FIREWALL (before any P1) ---
    times_p0 = [0.0, 0.25, 0.5, 0.75, 1.0]
    p0 = {}
    for T in times_p0:
        xT = integrate_tensor(X0, P0, T, 0.0, lambda t: 0.0, dt)
        assert xT.shape == Wn.shape, f"xT shape {xT.shape} != Wn shape {Wn.shape}"
        m1, var, k3, m4 = moments(xT, Wn)
        p0[T] = dict(mean=m1, var=var, k3=k3, m4=m4)
    var_ref = p0[0.0]['var']
    var_drift = max(abs(p0[t]['var'] - var_ref) for t in times_p0)
    k3_max = max(abs(p0[t]['k3']) for t in times_p0)
    mean_max = max(abs(p0[t]['mean']) for t in times_p0)
    # Gibbs identity at t=0
    gibbs_res = p0[0.0]['var'] + p0[0.0]['m4'] - 1.0

    # --- INDEPENDENT 1-D CHECK (separate pipeline) ---
    xm, wxm = np.polynomial.legendre.leggauss(4000)
    xm = L * xm
    wxm = L * wxm
    rx = np.exp(-0.5 * xm**2 - 0.25 * xm**4)
    wxn = wxm * rx
    wxn /= wxn.sum()
    ref_1d_m2 = (xm**2 * wxn).sum()
    ref_1d_m4 = (xm**4 * wxn).sum()
    ref_1d_m6 = (xm**6 * wxn).sum()

    # --- P1 ---
    p1 = {}
    for T in [T_STAR] + T_DIAG:
        row = {}
        for nb in NB_GRID:
            eps = 1.0 / np.sqrt(nb)
            xT = integrate_tensor(X0, P0, T, eps, q_p1, dt)
            assert xT.shape == Wn.shape
            m1, var, k3, m4 = moments(xT, Wn)
            k3F = k3 / np.sqrt(nb)
            g1F = k3F / var**1.5
            row[nb] = dict(mean=m1, var=var, k3=k3, k3F=k3F, g1F=g1F, eps=eps)
        p1[T] = row

    # --- ODD-EPSILON CHECK ---
    odd = {}
    for nb in [4, 16]:
        eps_p = 1.0 / np.sqrt(nb)
        eps_m = -eps_p
        xTp = integrate_tensor(X0, P0, T_STAR, eps_p, q_p1, dt)
        xTm = integrate_tensor(X0, P0, T_STAR, eps_m, q_p1, dt)
        _, varp, k3p, _ = moments(xTp, Wn)
        _, varm, k3m, _ = moments(xTm, Wn)
        odd[nb] = dict(k3_plus=k3p, k3_minus=k3m,
                       ratio=k3m / k3p if k3p != 0 else None,
                       sum=k3p + k3m)

    return dict(nx=nx, dt=dt, norm_raw=norm_raw,
                p0=p0, var_drift=var_drift, k3_max=k3_max, mean_max=mean_max,
                gibbs_res=gibbs_res,
                ref_1d=dict(m2=ref_1d_m2, m4=ref_1d_m4, m6=ref_1d_m6),
                p1=p1, odd=odd)

def main():
    t0 = time.time()
    results = {}
    for nx in NX_LIST:
        for dt in DT_LIST:
            key = f"nx{nx}_dt{dt}"
            print(f"--- {key} ---", flush=True)
            r = run_grid(nx, dt)
            results[key] = r
            v0 = r['p0'][0.0]['var']
            print(f"  P0 var: t=0:{v0:.8f} t=0.25:{r['p0'][0.25]['var']:.8f} "
                  f"t=0.5:{r['p0'][0.5]['var']:.8f} t=0.75:{r['p0'][0.75]['var']:.8f} "
                  f"t=1.0:{r['p0'][1.0]['var']:.8f}", flush=True)
            print(f"  P0 var drift: {r['var_drift']:.3e}  k3 max: {r['k3_max']:.3e}  "
                  f"mean max: {r['mean_max']:.3e}", flush=True)
            print(f"  Gibbs identity res: {r['gibbs_res']:.3e}", flush=True)
            print(f"  1D check: m2={r['ref_1d']['m2']:.10f} m4={r['ref_1d']['m4']:.10f} "
                  f"m6={r['ref_1d']['m6']:.10f}", flush=True)
            row = r['p1'][T_STAR]
            print(f"  P1 gamma1 at t*={T_STAR}: " +
                  " ".join(f"{nb}:{row[nb]['g1F']:.4e}" for nb in NB_GRID), flush=True)
            print(f"  Odd-eps: " + " ".join(
                f"nb{nb}:sum={r['odd'][nb]['sum']:.2e}" for nb in [4, 16]), flush=True)
    with open("v3_r2_results.json", "w") as f:
        json.dump(results, f, indent=1, default=str)
    print(f"\nTotal time: {time.time()-t0:.0f}s", flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
