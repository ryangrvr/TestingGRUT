"""Item 1: exact U(1)_X, Noether charge conservation in FRW, absence of anomalies, inflationary dilution and
the chronology of charge supply. Writes results/charge_conservation.json.

(a) Algebraic check: U depends on X only through |X|^2, so q1 dU/dq2 - q2 dU/dq1 = 0 identically (derived; checked
    at random field points, including the symmetric counterterms allowed by EFT rule 2).
(b) Lattice check: a KDK leapfrog of the full four-field action on a periodic 3D lattice in a radiation-dominated
    FRW background conserves the comoving charge Q = sum (q2 pi1 - q1 pi2) dx^3 to roundoff (derived exactly for
    a U(1)-invariant force and drift).
(c) Dilution: charge present N e-folds before the end of inflation is diluted by e^{-3N}; the charge floor
    rho >= m* |n_Q| bounds N for a required present yield (inflation and reheating are OUTSIDE-ACTION preparation data).
"""
import math
import time

import numpy as np

from common import (GSTAR, HD, LAM_D, LAM_PHI, LAM_X, GX, MD0, MP, MSTAR, M_X, V, Y_B, S_SM_COEF, dump)

T0 = time.time()
rng = np.random.default_rng(20261021)
out = dict(status="derived + computed; inflation/reheating inputs are outside-action preparation data")

# ---------------------------------------------------------------------------------- (a) algebraic Noether check


def dU(q, alpha=0.0, mX0sq=0.0):
    """Gradient of U(h=phi-v, q1, q2, D) including U(1)_X-symmetric counterterms (m_X0^2 |X|^2, alpha |X|^2 D^2)."""
    phi, q1, q2, D = q
    X2 = 0.5 * (q1 * q1 + q2 * q2)
    dphi = LAM_PHI * (phi * phi - V * V) * phi + 2 * GX ** 2 * phi * X2 + HD ** 2 * phi * D * D
    common = GX ** 2 * phi * phi + 2 * LAM_X * X2 + alpha * D * D + mX0sq
    dq1 = common * q1
    dq2 = common * q2
    dD = (MD0 ** 2 + HD ** 2 * phi * phi) * D + LAM_D * D ** 3 + 2 * alpha * X2 * D
    return np.array([dphi, dq1, dq2, dD])


res = []
for _ in range(2000):
    q = rng.normal(0, 800, 4)
    g = dU(q, alpha=rng.normal(0, 1e-3), mX0sq=rng.normal(0, 1e2))
    scale = abs(q[1] * g[2]) + abs(q[2] * g[1]) + 1e-300
    res.append(abs(q[1] * g[2] - q[2] * g[1]) / scale)
out["noether_identity_max_rel"] = float(max(res))
out["noether_statement"] = (
    "Under X -> e^{i a} X every term of U depends on |X|^2 only, so with j^mu = i(X d^mu X* - X* d^mu X) the field "
    "equation Box X = dU/dX* gives nabla_mu j^mu = i(X dU/dX - X* dU/dX*) = i(|X|^2 U' - |X|^2 U') = 0 on shell. "
    "In FRW, d(a^3 n_Q)/dt = -(spatial divergence), whose volume average vanishes: comoving net charge is constant.")
out["anomaly_statement"] = (
    "U(1)_X is global and acts only on the gauge-singlet scalar X. No fermion carries X charge, so there is no "
    "ABJ (gauge) anomaly, no mixed U(1)_X-gravitational anomaly (both need chiral fermions) and no instanton/"
    "sphaleron violation (no gauge field couples to the current). EFT rule 2 admits only U(1)_X-symmetric "
    "counterterms, so running cannot generate violation. Gravity-only production is C-symmetric (hidden_heating). "
    "Hence Candidate 1's dynamics cannot generate a net hidden yield. Quantum-gravity violation of global symmetries "
    "(black holes, wormholes, Planck-suppressed charged operators) is NOT in the action; adding it would be an action "
    "extension (GRUT_CANDIDATE_1.md P4 claims no quantum-gravity protection).")

# ---------------------------------------------------------------------------------- (b) FRW lattice check
Nl = 16
dx = 0.01           # GeV^-1 (= 1/m_X)
dt = 0.002          # GeV^-1
t_start, t_end = 0.5, 5.0   # radiation era a = (t/t_start)^(1/2), H(t_start) = 1 GeV


def a_of(t):
    return math.sqrt(t / t_start)


def lap(f):
    return (np.roll(f, 1, 0) + np.roll(f, -1, 0) + np.roll(f, 1, 1) + np.roll(f, -1, 1) + np.roll(f, 1, 2)
            + np.roll(f, -1, 2) - 6 * f) / dx ** 2


fields = np.zeros((4, Nl, Nl, Nl))
fields[0] = V + rng.normal(0, 50, (Nl, Nl, Nl))
amp = 300.0 * (1 + 0.3 * rng.normal(size=(Nl, Nl, Nl)))
theta = rng.uniform(0, 2 * math.pi, (Nl, Nl, Nl))
fields[1] = amp * np.cos(theta)
fields[2] = amp * np.sin(theta)
fields[3] = rng.normal(0, 30, (Nl, Nl, Nl))
a0 = a_of(t_start)
mom = np.zeros_like(fields)   # pi = a^3 dq/dt
omega = 60.0
mom[1] = a0 ** 3 * (-omega * fields[2]) + rng.normal(0, 1e4, (Nl, Nl, Nl))
mom[2] = a0 ** 3 * (omega * fields[1]) + rng.normal(0, 1e4, (Nl, Nl, Nl))
mom[0] = rng.normal(0, 1e4, (Nl, Nl, Nl))


def force(f, a):
    g = dU(f)
    return np.array([a * lap(f[i]) - a ** 3 * g[i] for i in range(4)])


def charge(f, p):
    return float(np.sum(f[2] * p[1] - f[1] * p[2]) * dx ** 3)


Q0 = charge(fields, mom)
t = t_start
qs = [(t, Q0, Q0 / (Nl * dx * a0) ** 3)]
nsteps = int(round((t_end - t_start) / dt))
F = force(fields, a_of(t))
for i in range(nsteps):
    mom += 0.5 * dt * F
    fields += dt * mom / a_of(t + 0.5 * dt) ** 3
    t += dt
    F = force(fields, a_of(t))
    mom += 0.5 * dt * F
    if (i + 1) % 250 == 0:
        a = a_of(t)
        qs.append((t, charge(fields, mom), charge(fields, mom) / (Nl * dx * a) ** 3))
qs = np.array(qs)
out["lattice_FRW_check"] = dict(
    grid=Nl, dx_GeVinv=dx, dt_GeVinv=dt, t_range_GeVinv=[t_start, t_end], a_final=a_of(t_end),
    Q_comoving_initial=Q0, max_rel_change_comoving_Q=float(np.max(np.abs(qs[:, 1] / Q0 - 1))),
    physical_density_ratio_final_over_initial=float(qs[-1, 2] / qs[0, 2]), expected_a_minus3=a_of(t_end) ** -3,
    max_abs_field_final_GeV=float(np.max(np.abs(fields))), all_finite=bool(np.all(np.isfinite(fields))),
    note="leapfrog in cosmic time with pi = a^3 dq/dt; each kick and drift conserves Q exactly for a U(1)-invariant "
         "potential and lattice Laplacian (summation by parts); residual is floating-point roundoff")

# ---------------------------------------------------------------------------------- (c) inflationary dilution
rows = []
for H in (1e2, 1e3, 1e6, 1e9, 1e13):
    rho_inf = 3 * H * H * MP * MP
    T_rh = (rho_inf / (math.pi ** 2 / 30 * GSTAR)) ** 0.25        # instantaneous reheating (most permissive)
    s_rh = S_SM_COEF * T_rh ** 3
    for r in (1.0, 0.0897, 1e-3):
        nQ_end = r * Y_B * s_rh
        Nmax = math.log(rho_inf / (MSTAR * nQ_end)) / 3.0
        rows.append(dict(H_inf_GeV=H, r=r, T_RH_instant_GeV=T_rh, nQ_end_GeV3=nQ_end, N_max_efolds=Nmax))
out["inflation_dilution"] = dict(
    label="OUTSIDE-ACTION preparation data (H_inf, reheating); charge floor rho_Q >= m* |n_Q| (quoted/derived)",
    condition="charge present N e-folds before the end of inflation needs n_Q,i = n_Q,end e^{3N}; requiring "
              "m* n_Q,i < 3 H^2 M_P^2 (subdominant to the inflaton) gives N < N_max. Instantaneous reheating is "
              "the most permissive case; any later entropy production lowers N_max further.",
    rows=rows)
out["chronology"] = (
    "The net hidden yield Y_Q = r Y_B is an initial datum supplied by OUTSIDE-ACTION physics at a time T_Q that the "
    "action does not fix: either after the last ~N_max e-folds of inflation (N_max ~ 10-20 for r ~ 1, see table; "
    "observable scales leave the horizon ~50-60 e-folds before the end), or during/after reheating. T_Q (visible "
    "temperature at supply), the spatial form of the supplied charge (homogeneous condensate, plasma, localized) and "
    "any later entropy dilution are free inputs. Y_Q is constant only between T_Q and today if no entropy is injected.")
out["runtime_s"] = time.time() - T0
dump("charge_conservation.json", out)
print({k: out[k] for k in ("noether_identity_max_rel",)}, out["lattice_FRW_check"]["max_rel_change_comoving_Q"],
      out["lattice_FRW_check"]["physical_density_ratio_final_over_initial"], out["lattice_FRW_check"]["expected_a_minus3"])
for r in rows:
    print(r)
print("runtime", out["runtime_s"])
