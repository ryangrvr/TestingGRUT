#!/usr/bin/env python3
"""v_production independent check, part (5): spectator energy fraction and number after instant reheating.

Inflation/reheating are outside-action preparation data [P]. Written from scratch by the verifier.
Units: H_e = 1, a_e = 1. Fields (phi, q1, q2, D), X = (q1 + i q2)/sqrt2, r = sqrt(q1^2+q2^2) (|X|^2 = r^2/2).
Quartic part (H >> v, m_D0):  U4 = lphi/4 phi^4 + gX^2 phi^2 r^2/2 + lX r^4/4 + hD^2 phi^2 D^2/2 + lD D^4/4.
Starobinsky-Yokoyama equilibrium P ∝ exp(-8 pi^2 U/(3 H^4)).
"""
import json, math, os, time
import numpy as np
from scipy import special

lphi = 0.01; gX = 0.1; lX = 0.001; hD = 0.01; lD = 0.001
g2 = gX**2; h2 = hD**2
MP = 2.435323e18; gstar = 106.75; Ymax = 4.373e-12
K8 = 8 * math.pi**2 / 3
rng = np.random.default_rng(20261021)
here = os.path.dirname(os.path.abspath(__file__))
out = {}


def U4(p, r, d):
    return lphi / 4 * p**4 + g2 * p**2 * r**2 / 2 + lX / 4 * r**4 + h2 * p**2 * d**2 / 2 + lD / 4 * d**4


def gradU(p, r, d):
    return (lphi * p**3 + g2 * p * r**2 + h2 * p * d**2,
            g2 * p**2 * r + lX * r**3,
            h2 * p**2 * d + lD * d**3)


# ---------- (a) equilibrium moments by tensor Gauss-Legendre quadrature (phi>=0, r>=0, D>=0 with measure r)
def gl(n, a, b):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w


def moments(n):
    p, wp = gl(n, 0, 7.5); r, wr = gl(n, 0, 14.0); d, wd = gl(n, 0, 12.0)
    P, R, Dd = np.meshgrid(p, r, d, indexing="ij")
    W = wp[:, None, None] * wr[None, :, None] * wd[None, None, :] * R * np.exp(-K8 * U4(P, R, Dd))
    Z = W.sum()
    gu = gradU(P, R, Dd)
    res = dict(phi2=(W * P**2).sum() / Z, q2sum=(W * R**2).sum() / Z, D2=(W * Dd**2).sum() / Z,
               U=(W * U4(P, R, Dd)).sum() / Z,
               phi_dot_gradU=(W * (P * gu[0] + R * gu[1] + Dd * gu[2])).sum() / Z,
               gradU2=(W * (gu[0]**2 + gu[1]**2 + gu[2]**2)).sum() / Z)
    res["sum_chi2"] = res["phi2"] + res["q2sum"] + res["D2"]
    return res


m120 = moments(120); m180 = moments(180)
out["equilibrium_moments_H_units"] = m180
out["equilibrium_moments_n120"] = m120
out["exact_identity_U_equals_4x3_over_32pi2"] = 4 * 3 / (32 * math.pi**2)
# single-field check 0.13177/sqrt(lambda)
out["single_field_phi2_coef"] = special.gamma(0.75) / special.gamma(0.25) * math.sqrt(3 / (2 * math.pi**2))

S = m180["sum_chi2"]; Ubar = m180["U"]
c_eps_frozen = S / 6 + Ubar / 3
c_eps_SR = S / 6 + Ubar / 3 - m180["phi_dot_gradU"] / 9 + m180["gradU2"] / 54
out["c_eps_minimal_no_velocity"] = c_eps_frozen
out["c_eps_with_slow_roll_velocity"] = c_eps_SR
out["frozen_over_potential_term"] = (S / 6) / (Ubar / 3)
out["critic_c_over_3_with_c0.07"] = 0.07 / 3
out["ratio_to_critic"] = c_eps_frozen / (0.07 / 3)

# ---------- (b) sanity check of the conformal argument in cosmic time (single quartic field, RD, H = 1/(2t))
from scipy.integrate import solve_ivp
def phys_check(lam, phi0, phidot0):
    te = 0.5
    def rhs(t, y):
        return [y[1], -1.5 / t * y[1] - lam * y[0]**3]
    tmax = 4.0e4
    sol = solve_ivp(rhs, [te, tmax], [phi0, phidot0], method="DOP853", rtol=1e-11, atol=1e-14, dense_output=True)
    ts = np.linspace(0.8 * tmax, tmax, 200001)
    ph, pd = sol.sol(ts)
    a4 = (ts / te)**2
    rho_a4 = np.mean((0.5 * pd**2 + lam / 4 * ph**4) * a4)
    Ec = 0.5 * (phi0 + phidot0)**2 + lam / 4 * phi0**4
    return dict(lam=lam, phi0=phi0, phidot0=phidot0, rho_a4_late=rho_a4, Ec_conformal=Ec,
                Ec_without_frozen_term=lam / 4 * phi0**4, ratio=rho_a4 / Ec)
out["cosmic_time_check"] = [phys_check(0.01, 3.0, 0.0), phys_check(0.001, 6.0, 0.0),
                            phys_check(0.01, 3.0, -0.01 * 27 / 3)]

# ---------- (c) sample P_eq exactly by rejection, then conformal dynamics chi'' + dV4 = 0
def sample(n):
    cp = K8 * lphi / 4; cd = K8 * lD / 4
    sig2 = 3 / (4 * math.pi**2 * lX)   # s = r^2 half-Gaussian: exp(-K8 lX s^2/4) = exp(-s^2/(2 sig2))
    outp = []; outr = []; outd = []
    while sum(len(x) for x in outp) < n:
        m = 4 * n
        p = (rng.gamma(0.25, size=m) / cp)**0.25 * rng.choice([-1, 1], size=m)
        d = (rng.gamma(0.25, size=m) / cd)**0.25 * rng.choice([-1, 1], size=m)
        s = np.abs(rng.normal(0, math.sqrt(sig2), size=m)); r = np.sqrt(s) * rng.choice([-1, 1], size=m)
        acc = rng.random(m) < np.exp(-K8 * (g2 * p**2 * r**2 / 2 + h2 * p**2 * d**2 / 2))
        outp.append(p[acc]); outr.append(r[acc]); outd.append(d[acc])
    return np.concatenate(outp)[:n], np.concatenate(outr)[:n], np.concatenate(outd)[:n]


Nsamp = 2000
p0, r0, d0 = sample(Nsamp)
big = sample(400000)
out["MC_moments_400k"] = dict(phi2=float(np.mean(big[0]**2)), q2sum=float(np.mean(big[1]**2)), D2=float(np.mean(big[2]**2)),
                              U=float(np.mean(U4(*big))))
del big


def evolve(p0, r0, d0, slow_roll, dt=0.02, T1=1000.0, T2=3000.0):
    x = np.stack([p0, r0, d0]).astype(float)
    v = x.copy()                       # chi' = a' phi0 = H phi0 (a_e = 1, H_e = 1)
    if slow_roll:
        g = np.stack(gradU(*x)); v = v - g / 3.0
    def acc(x):
        return -np.stack(gradU(*x))
    E0 = 0.5 * (v**2).sum(0) + U4(*x)
    a = acc(x)
    nsteps = int(round(T2 / dt)); n1 = int(round(T1 / dt))
    v2sum = np.zeros_like(x); cross = np.zeros_like(x); cnt = 0
    sgn = np.sign(x)
    for i in range(nsteps):
        v += 0.5 * dt * a
        x += dt * v
        a = acc(x)
        v += 0.5 * dt * a
        if i >= n1:
            v2sum += v**2; cnt += 1
            s = np.sign(x)
            cross += (s != sgn) & (s != 0)
            sgn = np.where(s != 0, s, sgn)
        elif i == n1 - 1:
            sgn = np.sign(x)
    E1 = 0.5 * (v**2).sum(0) + U4(*x)
    v2 = v2sum / cnt
    omega = math.pi * np.maximum(cross, 1) / (T2 - T1)
    J_i = v2 / omega
    return dict(E0=E0, E1=E1, J_i=J_i, omega=omega, crossings=cross)


t0 = time.time()
res = evolve(p0, r0, d0, slow_roll=True)
Jtot = res["J_i"].sum(0)
out["dyn"] = dict(Nsamp=Nsamp, dt=0.02, window=[1000, 3000], runtime_s=time.time() - t0,
                  max_rel_energy_drift=float(np.max(np.abs(res["E1"] / res["E0"] - 1))),
                  mean_Ec=float(np.mean(res["E0"])),
                  c_eps_sampled=float(np.mean(res["E0"]) / 3.0),
                  J_mean=float(np.mean(Jtot)), J_sem=float(np.std(Jtot) / math.sqrt(Nsamp)),
                  shares=[float(np.mean(res["J_i"][k]) / np.mean(Jtot)) for k in range(3)],
                  min_crossings=[int(res["crossings"][k].min()) for k in range(3)],
                  mean_quantum_energy_over_H=float(np.mean(res["E0"]) / np.mean(Jtot)))
# single-field calibration of the estimator: quartic field alone, J = (4/3) E / (0.8472 sqrt(lam) A)
cal = evolve(np.array([3.0]), np.array([0.0]), np.array([0.0]), slow_roll=False)
lam = lphi; E = cal["E0"][0]; A = (4 * E / lam)**0.25
out["calibration_single_phi"] = dict(J_est=float(cal["J_i"][0][0]), J_exact=float((4 / 3) * E / (0.84721 * math.sqrt(lam) * A)))

# ---------- (d) yields and thresholds
s_coef = 2 * math.pi**2 * gstar / 45
Ta_coef = (90 / (math.pi**2 * gstar))**0.25          # T a = Ta_coef sqrt(M_P H)
YJ_coef = out["dyn"]["J_mean"] / (s_coef * Ta_coef**3)
def floor_coef(ceps):
    xi = (ceps * gstar / 4)**0.25
    return (4 * special.zeta(3) / math.pi**2) * xi**3 / s_coef, xi
Yf, xiE = floor_coef(c_eps_frozen)
Yf_SR, _ = floor_coef(out["dyn"]["c_eps_sampled"])
thr = lambda coef: MP * (Ymax / coef)**(2 / 3)
out["yields"] = dict(xi_E_coef=xiE, Y_floor_coef=Yf, Y_J_coef=YJ_coef,
                     J_over_neq=YJ_coef / Yf,
                     H_threshold_floor_GeV=thr(Yf), H_threshold_J_GeV=thr(YJ_coef),
                     H_threshold_floor_SRvelocity_GeV=thr(Yf_SR),
                     Y_floor_at_1e13=Yf * (1e13 / MP)**1.5, Y_J_at_1e13=YJ_coef * (1e13 / MP)**1.5,
                     critic_Y_at_1e13=1e-9, critic_threshold=1e13 * (Ymax / 1e-9)**(2 / 3))
json.dump(out, open(os.path.join(here, "vp_spectator.json"), "w"), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
