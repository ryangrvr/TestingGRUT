"""Item 3a of track higgs_reheat: the post-inflation background R(t), H(t), a(t), rho_r(t) for (m_I, Gamma_I) = (3, 0.5) H_i.

Own code.  Two independent integrations:
  (1) our engine rpreheat.c (copied unchanged from grut_formation_review_20261025/tracks/rp_reheating/code; mode 0:
      oscillating quadratic inflaton released at rest with H = H_i, Gamma dI/dt decay into radiation, exact Jordan-frame
      trace equation for R).  Hidden fields are set to a static spectator point (X = Y = 0, Sigma at its vev) because
      their contribution to R is eta2 V ~ 1e-10 H_i^2 (eta2 = (H_i/M_P)^2 ~ 3e-15) [D]; this is stated, not hidden.
  (2) an own scipy integration of the same inflaton + radiation system (R = 2 b^2 I^2 - I'^2, units M_P = H_i = 1),
      used by modes.py and patch.py through dense output.
Writes ../results/background.json (downsampled table + oscillation statistics + engine-vs-scipy agreement).
"""
import json
import os
import subprocess
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
SP = Path(os.environ.get('HR_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/'
                         'scratchpad/wt2/higgs_reheat'))
BIN = SP / 'build' / 'rpreheat'
SRC = HERE / 'rpreheat.c'
B_I, GAM = 3.0, 0.5
MRED, GSTAR = 2.435e18, 106.75
COLS = ['tau', 'a', 'H', 'R', 'NX', 'NY', 'NS', 'Gauss', 'EX', 'eps', 'NYH', 'NYL', 'Nsig', 'absX', 'absY', 'S_over_vev',
        'rho_r', 'rho_I', 'radfrac', 'Hcheck', 'Wcheck', 'w_I', 'm2curv', 'KV', 'Haux', 'rho']


def run_engine(end=80.0, nsamp=16000, H_GeV=1.29e11):
    BIN.parent.mkdir(parents=True, exist_ok=True)
    if not BIN.exists() or BIN.stat().st_mtime < SRC.stat().st_mtime:
        subprocess.run(['gcc', '-O3', '-std=c11', '-Wall', '-o', str(BIN), str(SRC), '-lm'], check=True)
    eta2 = (H_GeV / MRED) ** 2
    g = lambda v: '%.17g' % float(v)
    h = 40.0
    z0 = [0.0, 0.0, 0.0, 0.0, 0.0, -h / np.sqrt(2)]          # X, Y, Sigma (re, im): Sigma = -i h/sqrt2 (vev)
    args = [str(BIN), '0', g(B_I), g(GAM), g(-5 / 6), g(eta2), g(10.0), g(1e-3), g(0.5), g(0.5), g(3.4586847e-3),
            g(30.0), g(h), g(7.159456), g(end), g(1e-11), g(1e-14), str(nsamp)] + [g(v) for v in z0]
    out = subprocess.run(args, capture_output=True, text=True, check=True)
    S = np.array([[float(v) for v in ln.split()[1:]] for ln in out.stdout.split('\n') if ln.startswith('S ')])
    return {c: S[:, i] for i, c in enumerate(COLS)}


def rhs(t, y, b=B_I, gam=GAM):
    lna, I, J, r = y
    H = np.sqrt((0.5 * J * J + 0.5 * b * b * I * I + r) / 3.0)
    return [H, J, -(3 * H + gam) * J - b * b * I, gam * J * J - 4 * H * r]


def background(end=200.0, b=B_I, gam=GAM, rtol=1e-11, atol=1e-14):
    """dense solution; R = 2 b^2 I^2 - J^2 (H_i units), H from Friedmann."""
    y0 = [0.0, np.sqrt(6.0) / b, 0.0, 0.0]
    sol = solve_ivp(rhs, (0, end), y0, method='DOP853', rtol=rtol, atol=atol, dense_output=True, args=(b, gam),
                    max_step=0.02)
    return sol


def fields(sol, t, b=B_I):
    lna, I, J, r = sol.sol(t)
    H = np.sqrt((0.5 * J * J + 0.5 * b * b * I * I + r) / 3.0)
    R = 2 * b * b * I * I - J * J
    Hdot = -(J * J + 4 * r / 3) / 2
    rhoI = 0.5 * J * J + 0.5 * b * b * I * I
    return dict(a=np.exp(lna), H=H, R=R, Hdot=Hdot, rho_r=r, rho_I=rhoI)


def T_GeV(rho_r, H_GeV, gstar=GSTAR):
    rho = rho_r * MRED ** 2 * H_GeV ** 2
    return (30 * rho / (np.pi ** 2 * gstar)) ** 0.25


def main():
    E = run_engine()
    sol = background(end=80.0)
    F = fields(sol, E['tau'])
    dev = {k: float(np.max(np.abs(F[k2] - E[k]) / (np.abs(E[k]).max()))) for k, k2 in
           (('R', 'R'), ('H', 'H'), ('a', 'a'), ('rho_r', 'rho_r'))}
    tau = E['tau']
    # oscillation statistics: local minima/maxima of R
    R = E['R']
    imin = [i for i in range(1, len(R) - 1) if R[i] < R[i - 1] and R[i] <= R[i + 1]]
    imax = [i for i in range(1, len(R) - 1) if R[i] > R[i - 1] and R[i] >= R[i + 1]]
    osc = [dict(tau=float(tau[i]), R_min=float(R[i]), H=float(E['H'][i]), R_min_over_H2=float(R[i] / E['H'][i] ** 2),
                radfrac=float(E['radfrac'][i])) for i in imin[:30]]
    oscM = [dict(tau=float(tau[i]), R_max=float(R[i])) for i in imax[:30]]
    first_neg = float(tau[np.where(R < 0)[0][0]])
    # time-averaged R per oscillation window
    rows = []
    for H_GeV in (1.29e11, 1.51e11):
        Tm = T_GeV(E['rho_r'], H_GeV)
        rows.append(dict(H_GeV=H_GeV, T_max_GeV=float(Tm.max()), tau_Tmax=float(tau[np.argmax(Tm)]),
                         T_at_tau={f'{t}': float(np.interp(t, tau, Tm)) for t in (0.05, 0.1, 0.3, 1, 2, 5, 10, 20, 40, 80)}))
    # radiation-domination time
    i_eq = np.where(E['radfrac'] > 0.5)[0][0]
    step = max(1, len(tau) // 4000)
    table = {c: E[c][::step].tolist() for c in ('tau', 'a', 'H', 'R', 'rho_r', 'rho_I', 'radfrac')}
    out = dict(engine='rpreheat.c (copied unchanged), mode 0, b=3, gamma=0.5, rtol 1e-11',
               engine_vs_scipy_max_rel_dev=dev, R0=float(R[0]), first_R_negative_tau=first_neg,
               R_minima=osc, R_maxima=oscM, radiation_equality_tau=float(tau[i_eq]),
               a_at_equality=float(E['a'][i_eq]), temperatures=rows, table_downsampled=table,
               checks=dict(max_H_consistency='Haux vs constraint: %.2e' % float(np.max(np.abs(E['Hcheck'])))))
    json.dump(out, open(RES / 'background.json', 'w'), indent=None, default=float)
    print(json.dumps({k: v for k, v in out.items() if k not in ('table_downsampled', 'R_maxima')}, indent=1)[:4000])


if __name__ == '__main__':
    main()
