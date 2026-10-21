"""Part 4: Z2 sign memory of a charge-restored (phi ~ 0) patch embedded in a +v background.

Toy model (stated as a toy): spherically symmetric 3D field S = phi/v obeying the unchanged phi equation
    S_tt = (1/r)(r S)_rr - (S^2 - 1 + a x(r,t)) S,
with the charged-field background x = F^2 = 2|X|^2/v^2 FROZEN in space and prescribed in time
(x_in(t) w(r), w a smooth ball of radius R).  This isolates the question "tachyonic growth vs boundary-driven
relaxation"; the X dynamics, non-spherical modes and expansion are NOT included.
Units: length/time ell = 1/(sqrt(lam_phi) v) (= 1/m_X at the benchmark, 1/(100 GeV)).
Run: python3 z2_memory.py  (about 3-5 minutes, 1 core)
"""
import os
import numpy as np
from common import A, LAM_PHI, V_GEV, GEV_TO_FM, save_json

a = A
DR = 0.05
DT = 0.02
HERE = os.path.dirname(os.path.abspath(__file__))
PROFILE = os.path.normpath(os.path.join(HERE, "../../../../grut_candidate_1_20261018/code/results_profile.csv"))


def ball(r, R, delta=0.6):
    return 0.5 * (1 - np.tanh((r - R) / delta))


def static_S(r, xprof, iters=15000):
    """positive static solution of -lap S + (S^2-1+a x)S = 0, S'(0)=0, S(rmax)=1 (gradient flow)."""
    S = 1 - 0.99 * (xprof > 1.0 / a)
    u = r * S
    dt = 0.4 * DR * DR
    for _ in range(iters):
        S = np.where(r > 0, u / np.where(r > 0, r, 1), 0)
        S[0] = S[1]
        lap = np.zeros_like(u)
        lap[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / DR**2
        u[1:-1] += dt * (lap[1:-1] - (S[1:-1] ** 2 - 1 + a * xprof[1:-1]) * u[1:-1])
        u[0] = 0.0
        u[-1] = r[-1]
    S = u / np.where(r > 0, r, 1); S[0] = S[1]
    return S


def evolve(r, S0, xin_of_t, wprof, t_end, seed=None, rec_every=0.25, Rin=None, seed_vel=None):
    """leapfrog in u=rS; returns times, S(0,t), interior volume-mean S(t) (r<Rin), and final S."""
    u = r * S0
    if seed is not None:
        u = u + r * seed
    def acc(u, t):
        S = u / np.where(r > 0, r, 1)
        x = xin_of_t(t) * wprof
        out = np.zeros_like(u)
        out[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / DR**2 - (S[1:-1] ** 2 - 1 + a * x[1:-1]) * u[1:-1]
        return out
    nsteps = int(round(t_end / DT))
    u_old = u - 0.5 * DT * DT * acc(u, 0.0)
    if seed_vel is not None:
        u_old = u_old - DT * r * seed_vel
    u_old[0] = 0.0; u_old[-1] = r[-1]
    ts, c, cin = [], [], []
    if Rin is None:
        Rin = 0.5 * r[-1]
    vin = (r < Rin) * r * r
    every = max(1, int(round(rec_every / DT)))
    for k in range(nsteps):
        t = k * DT
        un = 2 * u - u_old + DT * DT * acc(u, t)
        un[0] = 0.0; un[-1] = r[-1]
        u_old, u = u, un
        if k % every == 0:
            ts.append(t + DT); c.append(u[1] / r[1])
            Sn = u / np.where(r > 0, r, 1); Sn[0] = Sn[1]
            cin.append(np.sum(vin * Sn) / np.sum(vin))
    S = u / np.where(r > 0, r, 1); S[0] = S[1]
    return np.array(ts), np.array(c), np.array(cin), S


def run_case(R, x_before, x_after, eps, tauQ=0.0, t_after=60.0, delta=0.6, noise=None, theta=None):
    m_t = np.sqrt(max(1e-12, 1 - a * x_after))
    t_end = tauQ + R + t_after
    rmax = R + 0.6 * t_end + 10
    r = np.arange(0, rmax + DR / 2, DR)
    w = ball(r, R, delta)
    S0 = static_S(r, x_before * w)
    if tauQ <= 0:
        xin = lambda t: x_after
    else:
        xin = lambda t: x_before + (x_after - x_before) * min(1.0, t / tauQ)
    seed = -eps * w if noise is None else noise(r) * eps * w
    t_after = max(t_after, 30.0 + np.log(1.0 / max(eps, 1e-300)) / m_t)
    seed_vel = None
    if theta is not None:      # oscillating seed: S = -eps cos(w_in t + theta) inside, w_in = sqrt(a x_before - 1)
        w_in = np.sqrt(max(0.0, a * x_before - 1))
        seed = -eps * np.cos(theta) * w
        seed_vel = eps * w_in * np.sin(theta) * w
    ts, c, cin, Sf = evolve(r, S0, xin, w, t_end, seed=seed, Rin=0.5 * R, seed_vel=seed_vel)
    late = c[ts > t_end - 10]
    j = np.argmax(np.abs(c) > 0.5) if np.any(np.abs(c) > 0.5) else None
    dec_sign = None if j is None else int(np.sign(c[j]))
    dec_time = None if j is None else float(ts[j])
    ji = np.argmax(np.abs(cin) > 0.5) if np.any(np.abs(cin) > 0.5) else None
    dec_sign_in = None if ji is None else int(np.sign(cin[ji]))
    neg = cin < -0.5
    neg_duration = float(neg.sum() * (ts[1] - ts[0])) if neg.any() else 0.0
    return dict(R=R, x_before=x_before, x_after=x_after, eps=eps, tauQ=tauQ, m_tach=m_t,
                S_center_initial=float(S0[0]), center_decision_sign=dec_sign, center_decision_time=dec_time,
                inner_half_volume_decision_sign=dec_sign_in,
                interior_time_below_minus_half=neg_duration,
                S_center_final_mean=float(np.mean(late)), S_center_min=float(c.min()),
                interior_mean_min=float(cin.min()), final_center_negative=bool(np.mean(late) < -0.3),
                flipped=bool(dec_sign == -1))


def first_flip_R(x_after, eps, tauQ=0.0, Rs=None, x_before=3.4):
    rows = []
    for R in Rs:
        rows.append(run_case(R, x_before, x_after, eps, tauQ=tauQ))
    flipR = [row["R"] for row in rows if row["flipped"]]
    return (min(flipR) if flipR else None), rows


if __name__ == "__main__":
    out = dict(status="computed (toy: frozen X background, spherical symmetry, flat space)")
    # -- Dirichlet estimate R_* = pi/m_tach (sharp ball) as a function of the patch density (x<1) --
    dir_rows = []
    for x in [0.0, 0.25, 0.5, 0.75, 0.9, 0.99]:
        mt = np.sqrt(1 - a * x)
        R_bar = np.pi / mt
        dir_rows.append(dict(x_after=x, m_tach_over_unit=mt, m_tach_GeV=mt * np.sqrt(LAM_PHI) * V_GEV,
                             R_star_bar=R_bar, R_star_fm=R_bar / (np.sqrt(LAM_PHI) * V_GEV) * GEV_TO_FM))
    out["dirichlet_R_star"] = dir_rows
    # -- sudden quench scans (coherent negative seed = worst case for memory) --
    Rs = [1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 28]
    scans = []
    for x_after in [0.0, 0.5, 0.8]:
        for eps in [1e-2, 1e-4, 1e-6]:
            Rf, rows = first_flip_R(x_after, eps, Rs=Rs)
            mt = np.sqrt(1 - a * x_after)
            m_in = np.sqrt(a * 3.4 - 1)
            R_heur = 2.0 + np.log(0.5 / eps) / m_in + np.log(1 / eps) / mt
            scans.append(dict(x_after=x_after, eps=eps, R_first_flip=Rf, pi_over_mtach=np.pi / mt,
                              ln_inv_eps_over_mtach=np.log(1 / eps) / mt, R_heuristic=R_heur,
                              R_heuristic_formula="2 + ln(0.5/eps)/m_in + ln(1/eps)/m_tach, m_in=sqrt(a x_before-1)",
                              rows=rows))
            print("sudden", x_after, eps, Rf, np.pi / mt, np.log(1 / eps) / mt, flush=True)
    out["sudden_quench_scans"] = scans
    # -- finite quench time (charge removed linearly over tauQ) --
    slow = []
    thetas = np.linspace(0, 2 * np.pi, 8, endpoint=False)
    for tauQ in [0.0, 10.0, 100.0]:
        for R in [4.0, 8.0, 16.0, 32.0]:
            rows = [run_case(R, 3.4, 0.0, 1e-4, tauQ=tauQ, theta=th, t_after=30.0) for th in thetas]
            neg = [rw["center_decision_sign"] == -1 for rw in rows]
            slow.append(dict(tauQ=tauQ, R=R, x_after=0.0, eps=1e-4, n_phases=len(thetas),
                             fraction_negative_center_decision=float(np.mean(neg)),
                             fraction_negative_inner_half=float(np.mean([rw["inner_half_volume_decision_sign"] == -1 for rw in rows])),
                             rows=rows))
            print("ramp", tauQ, R, np.mean(neg), flush=True)
    out["ramp_quench_scans"] = slow
    # -- control: charged interior (x>1) with large negative kick: oscillation about 0, never -v --
    ctrl = []
    for kick in [0.3, 0.6]:
        rr = run_case(6.0, 3.4, 3.4, kick, t_after=150.0)
        ctrl.append(rr)
    out["charged_interior_control"] = ctrl
    # -- archived soliton S,F profile, charge removed (sudden), coherent negative seed --
    arch = []
    try:
        d = np.loadtxt(PROFILE, delimiter=",", skiprows=1)
        rp, Fp, Sp = d[:, 0], d[:, 1], d[:, 3]
        for eps in [0.0, 1e-2, 0.05, 0.1, 0.3]:
            t_end = 80.0
            rmax = 60.0
            r = np.arange(0, rmax + DR / 2, DR)
            S0 = np.interp(r, rp, Sp, right=1.0)
            x0 = np.interp(r, rp, Fp**2, right=0.0)
            w = x0 / x0[0]
            for xa_frac in [0.0]:
                ts, c, cin, Sf = evolve(r, S0, lambda t: 0.0, w, t_end, seed=-eps * (r < 3.3), Rin=1.65)
                late = c[ts > t_end - 10]
                j = np.argmax(np.abs(c) > 0.5) if np.any(np.abs(c) > 0.5) else None
                arch.append(dict(eps_negative_seed_inside_r3p3=eps, x_after=0.0, S_center_initial=float(S0[0]),
                                 center_decision_sign=None if j is None else int(np.sign(c[j])),
                                 center_decision_time=None if j is None else float(ts[j]),
                                 S_center_min=float(c.min()), interior_mean_min=float(cin.min()),
                                 S_center_final_mean=float(late.mean()),
                                 flipped=bool(j is not None and c[j] < 0)))
                print("archived", eps, arch[-1], flush=True)
    except OSError as e:
        arch.append(dict(error=str(e)))
    out["archived_soliton_charge_removed"] = arch
    # -- closed bubble fate: a flipped interior is bounded by +v and collapses --
    big = run_case(16.0, 3.4, 0.0, 1e-2, t_after=150.0)
    out["closed_bubble_long_run"] = big
    save_json("z2_memory.json", out)
    print("done")
