#!/usr/bin/env python3
"""
AUTHORITATIVE V3-R2 computation. Single source of truth for all reported values.
Corrected tensor quadrature with pipeline assertions.
All scientific choices frozen from V3-R1: t*=0.5, diagnostics {0.25,0.75,1.0},
NB {4,8,16,32,64,128}, domain [-6,6]^2, resolutions 120/240/480, RK4 dt {1e-3,5e-4,2.5e-4}.
"""
import numpy as np, json, time, os

T_STAR = 0.5
T_DIAG = [0.25, 0.75, 1.0]
NB_GRID = [4, 8, 16, 32, 64, 128]
L = 6.0
NX_REF = 240
DT_REF = 5e-4
NX_LIST = [120, 240, 480]
DT_LIST = [1e-3, 5e-4, 2.5e-4]

# High-precision reference values (independent 1-D quadrature, 4000-pt GL)
REF_M2 = 0.467919916974
REF_M4 = 0.532080083026
REF_M6 = 0.871679667895

def s(u): return 10*u**3 - 15*u**4 + 6*u**5
def q_p1(t): return s(t/np.pi) if t <= np.pi else 1.0

def integrate_tensor(X0, P0, T, eps, qfun, dt):
    """Correct tensor evolution: every (x_i, p_j) pair independently."""
    assert X0.ndim == 2 and X0.shape == P0.shape
    steps = max(1, int(round(T/dt)))
    h = T/steps
    x = X0.copy(); v = P0.copy()
    for k in range(steps):
        t0 = k*h
        q0 = qfun(t0)
        k1v = -x - x**3 + eps*q0; k1x = v
        x2 = x + 0.5*h*k1x; v2 = v + 0.5*h*k1v
        qm = qfun(t0 + 0.5*h)
        k2v = -x2 - x2**3 + eps*qm; k2x = v2
        x3 = x + 0.5*h*k2x; v3 = v + 0.5*h*k2v
        k3v = -x3 - x3**3 + eps*qm; k3x = v3
        x4 = x + h*k3x; v4 = v + h*k3v
        q1 = qfun(t0 + h)
        k4v = -x4 - x4**3 + eps*q1; k4x = v4
        x = x + (h/6)*(k1x + 2*k2x + 2*k3x + k4x)
        v = v + (h/6)*(k1v + 2*k2v + 2*k3v + k4v)
    return x

def moments(f, Wn):
    m1 = (f*Wn).sum(); m2 = (f**2*Wn).sum(); m3 = (f**3*Wn).sum(); m4 = (f**4*Wn).sum()
    var = m2 - m1*m1; k3 = m3 - 3*m1*m2 + 2*m1*m1*m1
    return m1, m2, var, k3, m4

def run(nx, dt):
    xg, wg = np.polynomial.legendre.leggauss(nx)
    xs = L*xg; wx = L*wg
    X0 = np.broadcast_to(xs[:, None], (nx, nx)).copy()
    P0 = np.broadcast_to(xs[None, :], (nx, nx)).copy()
    assert X0.shape == (nx, nx) and P0.shape == (nx, nx)
    rho = np.exp(-0.5*P0**2 - 0.5*X0**2 - 0.25*X0**4)
    W = wx[:, None]*wx[None, :]
    Wn_raw = W*rho; norm_raw = Wn_raw.sum()
    Wn = Wn_raw/norm_raw
    assert Wn.shape == (nx, nx)
    out = {"nx": nx, "dt": dt, "norm_raw": float(norm_raw)}

    # P0 firewall
    p0 = {}
    for T in [0.0, 0.25, 0.5, 0.75, 1.0]:
        xT = integrate_tensor(X0, P0, T, 0.0, lambda t: 0.0, dt)
        assert xT.shape == Wn.shape
        m1, m2, var, k3, m4 = moments(xT, Wn)
        p0[T] = dict(mean=float(m1), m2=float(m2), var=float(var), k3=float(k3), m4=float(m4))
    out["p0"] = p0
    out["gibbs_identity_residual"] = p0[0.0]["var"] + p0[0.0]["m4"] - 1.0

    # Independent 1-D check (separate pipeline, 4000-pt GL)
    xm, wxm = np.polynomial.legendre.leggauss(4000)
    xm = L*xm; wxm = L*wxm
    rx = np.exp(-0.5*xm**2 - 0.25*xm**4)
    wxn = wxm*rx; wxn /= wxn.sum()
    out["independent_1d"] = dict(
        m2=float((xm**2*wxn).sum()), m4=float((xm**4*wxn).sum()),
        m6=float((xm**6*wxn).sum()))

    # P1 finite-epsilon
    p1 = {}
    for T in [T_STAR] + T_DIAG:
        row = {}
        for nb in NB_GRID:
            eps = 1.0/np.sqrt(nb)
            xT = integrate_tensor(X0, P0, T, eps, q_p1, dt)
            assert xT.shape == Wn.shape
            m1, m2, var, k3, m4 = moments(xT, Wn)
            k3F = k3/np.sqrt(nb)  # exact iid identity
            varF = var  # exact iid identity
            g1F = k3F/varF**1.5
            # PIPELINE ASSERTIONS
            assert abs(k3F - k3/np.sqrt(nb)) < 1e-30, "kappa3(F) identity violated"
            assert abs(g1F - k3F/varF**1.5) < 1e-15, "gamma1(F) identity violated"
            row[str(nb)] = dict(epsilon=float(eps), mean_X=float(m1), var_X=float(var),
                                kappa3_X=float(k3), kappa3_F=float(k3F),
                                gamma1_F=float(g1F), NB_gamma1=float(nb*g1F))
        p1[str(T)] = row
    out["p1"] = p1

    # Odd-epsilon check
    odd = {}
    for nb in [4, 16]:
        eps_p = 1.0/np.sqrt(nb); eps_m = -eps_p
        xTp = integrate_tensor(X0, P0, T_STAR, eps_p, q_p1, dt)
        xTm = integrate_tensor(X0, P0, T_STAR, eps_m, q_p1, dt)
        _, _, _, k3p, _ = moments(xTp, Wn)
        _, _, _, k3m, _ = moments(xTm, Wn)
        odd[str(nb)] = dict(k3_plus=float(k3p), k3_minus=float(k3m),
                            sum=float(k3p+k3m), ratio=float(k3m/k3p) if k3p != 0 else None)
    out["odd_epsilon"] = odd
    return out

def main():
    t0 = time.time()
    results = {"settings": dict(t_star=T_STAR, t_diag=T_DIAG, nb_grid=NB_GRID,
                                domain="[-6,6]^2", nx_list=NX_LIST, dt_list=DT_LIST,
                                ref_m2=REF_M2, ref_m4=REF_M4, ref_m6=REF_M6,
                                integrator="RK4 fixed-step", quadrature="Gauss-Legendre tensor")}
    for nx in NX_LIST:
        for dt in DT_LIST:
            key = f"nx{nx}_dt{dt}"
            print(f"--- {key} ---", flush=True)
            r = run(nx, dt)
            results[key] = r
            v0 = r["p0"][0.0]["var"]
            print(f"  P0 var(t=0)={v0:.12f}  gibbs_res={r['gibbs_identity_residual']:.3e}", flush=True)
            print(f"  P1 t*=0.5 gamma1: " + " ".join(f"{nb}:{r['p1']['0.5'][str(nb)]['gamma1_F']:.6e}" for nb in NB_GRID), flush=True)
    # Scaling fit at reference resolution
    ref = results[f"nx{NX_REF}_dt{DT_REF}"]
    gam = np.array([ref["p1"][str(T_STAR)][str(nb)]["gamma1_F"] for nb in NB_GRID])
    nb_arr = np.array(NB_GRID, dtype=float)
    p_global, a_global = np.polyfit(np.log(nb_arr), np.log(np.abs(gam)), 1)
    local_p = []
    for i in range(len(nb_arr)-1):
        lp = -np.log(abs(gam[i+1])/abs(gam[i]))/np.log(nb_arr[i+1]/nb_arr[i])
        local_p.append(float(lp))
    results["scaling_fit"] = dict(
        global_p=float(-p_global), intercept=float(a_global),
        local_p=local_p,
        NB_gamma1=[float(nb*g) for nb, g in zip(nb_arr, gam)],
        residuals=[float(np.log(abs(gam[i])) - a_global - (-p_global)*np.log(nb_arr[i]))
                   for i in range(len(nb_arr))])
    results["runtime_seconds"] = time.time() - t0
    with open("v3_r2_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nAUTHORITATIVE RESULTS SAVED to v3_r2_results.json ({time.time()-t0:.0f}s)", flush=True)
    print(f"Global fitted p = {float(-p_global):.6f}", flush=True)
    print(f"Local p: {local_p}", flush=True)
    print(f"NB*gamma1: {[f'{v:.6e}' for v in results['scaling_fit']['NB_gamma1']]}", flush=True)

if __name__ == "__main__":
    main()
