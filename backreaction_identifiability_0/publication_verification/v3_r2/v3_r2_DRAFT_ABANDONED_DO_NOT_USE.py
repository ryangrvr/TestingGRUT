#!/usr/bin/env python3
"""
V3-R2: CORRECT TENSOR QUADRATURE. Full 2-D (nx, np_) grids; P0 firewall first.
Preregistered choices frozen from V3-R1.
"""
import numpy as np, json, time

T_STAR = 0.5
T_DIAG = [0.25, 0.75, 1.0]
NB_GRID = [4, 8, 16, 32, 64, 128]
LX = 6.0
NX_LIST = [120, 240, 480]
DT_LIST = [1e-3, 5e-4, 2.5e-4]

def s(u): return 10*u**3 - 15*u**4 + 6*u**5
def q_p1(t): return s(t/np.pi) if t <= np.pi else 1.0

def integrate_tensor(X0, P0, T, eps, qfun, dt):
    steps = max(1, int(round(T/dt)))
    h = T/steps
    x = X0.copy(); v = P0.copy()
    for k in range(steps):
        t0=k*h; tm=t0+0.5*h; t1=t0+h
        q0,qm,q1 = qfun(t0), qfun(tm), qfun(t1)
        k1v = -x - x**3 + eps*q0;  k1x = v
        x2 = x + 0.5*h*k1x; v2 = v + 0.5*h*k1v
        k2v = -x2 - x2**3 + eps*qm;  k2x = v2
        x3 = x + 0.5*h*k2x; v3 = v + 0.5*h*k2v
        k3v = -x3 - x3**3 + eps*qm;  k3x = v3
        x4 = x + h*k3x; v4 = v + h*k3v
        k4v = -x4 - x4**3 + eps*q1;  k4x = v4
        x = x + (h/6)*(k1x + 2*k2x + 2*k3x + k4x)
        v = v + (h/6)*(k1v + 2*k2v + 2*k3v + k4v)
    return x

def run_resolution(nx, dt):
    xs, wx = np.polynomial.legendre.leggauss(nx)
    xs = 0.5*LX*xs + 0.5*LX; wx = 0.5*LX*wx
    X0 = np.broadcast_to(xs[:, None], (nx, nx)).copy()
    P0 = np.broadcast_to(xs[None, :], (nx, nx)).copy()
    W = wx[:, None] * wx[None, :]
    rho = np.exp(-0.5*P0**2 - 0.5*X0**2 - 0.25*X0**4)
    Wn = W * rho
    Z_raw = Wn.sum()
    Wn /= Wn.sum()
    assert X0.shape == (nx, nx) and P0.shape == (nx, nx) and Wn.shape == (nx, nx)
    out = {"nx": nx, "dt": dt, "Z_raw": Z_raw}
    times = [0.0] + [T_STAR] + T_DIAG
    # P0 firewall
    p0 = {}
    for T in times:
        xT = integrate_tensor(X0, P0, T, 0.0, lambda t: 0.0, dt)
        assert xT.shape == Wn.shape
        m = (xT*Wn).sum(); m2 = (xT**2*Wn).sum(); m3 = (xT**3*Wn).sum()
        var = m2 - m*m; k3 = m3 - 3*m*m2 + 2*m*m*m
        p0[T] = dict(mean=m, var=var, k3=k3)
    out["p0"] = p0
    # Gibbs identity E[x^2]+E[x^4]=1
    x2 = (X0**2*Wn).sum(); x4 = (X0**4*Wn).sum()
    out["gibbs_identity_residual"] = x2 + x4 - 1.0
    # P1
    p1 = {}
    for T in [T_STAR] + T_DIAG:
        dT = {}
        for nb in NB_GRID:
            eps = 1.0/np.sqrt(nb)
            xT = integrate_tensor(X0, P0, T, eps, q_p1, dt)
            assert xT.shape == Wn.shape
            m = (xT*Wn).sum(); m2 = (xT**2*Wn).sum(); m3 = (xT**3*Wn).sum()
            var = m2 - m*m; k3 = m3 - 3*m*m2 + 2*m*m*m
            dT[nb] = dict(mean=m, var=var, k3=k3)
        p1[T] = dT
    out["p1"] = p1
    return out

def main():
    results = {}
    t0wall = time.time()
    for nx in NX_LIST:
        for dt in DT_LIST:
            key = f"nx{nx}_dt{dt}"
            print(f"--- {key} ---", flush=True)
            r = run_resolution(nx, dt)
            results[key] = r
            # P0 firewall report
            v0 = r["p0"][0.0]["var"]
            vars_by_t = {T: r["p0"][T]["var"] for T in r["p0"]}
            k3max = max(abs(r["p0"][T]["k3"]) for T in r["p0"])
            print(f"nx={nx}: P0 variance by time (should be constant ~{v0:.6f}): "
                  + " ".join(f"t={T}:{v:.6f}" for T, v in vars_by_t.items()), flush=True)
            print(f"nx={nx}: Gibbs identity residual: {r['gibbs_identity_residual']:.3e}", flush=True)
            print(f"nx={nx}: P0 k3 max abs: {k3max:.3e}", flush=True)
    with open("v3_r2_results.json", "w") as f:
        json.dump(results, f, indent=1, default=str)
    print(f"DONE in {time.time()-t0wall:.0f}s", flush=True)

if __name__ == "__main__":
    main()
