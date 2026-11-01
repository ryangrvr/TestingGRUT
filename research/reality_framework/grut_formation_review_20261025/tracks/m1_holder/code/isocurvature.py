#!/usr/bin/env python3
"""m1_holder track (Claude, own code) -- part 4b: de Sitter suppression of the spectator charge power.

Own mode functions: for a canonically normalised real scalar of constant mass m in exact de Sitter (H = H_i), the
Bunch-Davies mode in the variable S = sqrt(2 k^3) delta r_k / H_i and x = k/(aH) = e^{-N} is exactly
    S(x) = sqrt(pi/2) e^{-pi mu/2} x^{3/2} H^{(1)}_{i mu}(x),   mu = sqrt(m^2/H^2 - 9/4),   dS/dN = -x dS/dx,
normalised so that |S| -> x deep inside the horizon and |S|^2 -> 1 at late times for m = 0 (P = H^2/4 pi^2).
It is evaluated with mpmath (complex order, 40 digits) at every N and cross-checked against an own DOP853 integration of
S_NN + 3 S_N + (e^{-2N} + m^2) S = 0 started at x0 = 1000 and 2000 (Codex's procedure, re-implemented).
Power: P_{delta N/N}(N) = sum_a |L_pos,a S_a(N) + L_vel,a S_N,a(N)|^2 / (4 pi^2) with the late charge response L from
  (i) our own finite differences (results/charge_transfer.json) and (ii) Codex's canonical row [Q perturbations.json]
projected on our eigenvectors.  Codex's table is checked against the e^{-3 Delta N} law.
"""
import json
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
mp.mp.dps = 40
CODEX_TABLE = {0.0: 1.9858166793719456e-05, 2.0: 4.984135128883863e-08, 4.0: 1.2357522432943231e-10, 8.0: 7.592347839932724e-16,
               20.0: 1.7611900904580933e-31, 50.0: 1.4430133116013324e-70, 60.0: 1.3503760669554607e-83}
NEW_PACKET_50 = 3.435e-70            # quoted in BRIEF_M1 for the unpushed packet; not readable from git on 2026-11-01


def S_exact(m2, N):
    """returns (S, S_N) complex at e-fold N after horizon exit (N = ln(aH/k))."""
    mu = mp.sqrt(mp.mpf(m2) - mp.mpf(9) / 4)
    x = mp.e ** (-mp.mpf(N))
    pref = mp.sqrt(mp.pi / 2) * mp.e ** (-mp.pi * mu / 2)
    f = lambda xx: pref * xx ** mp.mpf(1.5) * mp.hankel1(1j * mu, xx)
    S = f(x)
    dSdx = mp.diff(f, x)
    return complex(S), complex(-x * dSdx)


def S_ode(m2s, x0, points):
    """own re-implementation of the integration in N with the exact Hankel start at x0."""
    m2s = np.asarray(m2s, float); n = len(m2s)
    init = np.zeros(2 * n, complex)
    for i, m2 in enumerate(m2s):
        S, SN = S_exact(m2, -np.log(x0)); init[i] = S; init[n + i] = SN
    sol = solve_ivp(lambda N, y: np.r_[y[n:], -3 * y[n:] - (np.exp(-2 * N) + m2s) * y[:n]], (-np.log(x0), max(points)), init,
                    method='DOP853', rtol=1e-11, atol=1e-13, t_eval=points)
    return sol.y[:n], sol.y[n:]


def main():
    ct = json.load(open(RES / 'charge_transfer.json'))
    rr = ct['response_row']
    m2 = np.array(rr['physical_eigenvalues'])
    L_own = np.array(rr['own_richardson']); L_codex = np.array(rr['codex_on_own_basis'])
    out = dict(masses2=m2.tolist(), mu=np.sqrt(m2 - 2.25).tolist())
    # normalisation sanity: massless field -> |S|^2 -> 1
    S0, _ = S_exact(0.0, 10.0)
    out['massless_late_|S|^2'] = abs(S0) ** 2
    # exact vs ODE at N = 2, 4, 8
    pts = np.array([0.0, 2.0, 4.0, 8.0])
    So, Do = S_ode(m2, 1000.0, pts); So2, Do2 = S_ode(m2, 2000.0, pts)
    cmp = []
    for j, N in enumerate(pts):
        for i, mm in enumerate(m2):
            Se, De = S_exact(mm, N)
            cmp.append(dict(N=N, m2=mm, abs_S_exact=abs(Se), abs_S_ode=abs(So[i, j]), rel_S=abs(So[i, j] / Se - 1), rel_SN=abs(Do[i, j] / De - 1),
                            rel_S_x0_2000=abs(So2[i, j] / Se - 1)))
    out['exact_vs_ode'] = cmp
    out['max_rel_diff_exact_vs_ode'] = max(max(c['rel_S'], c['rel_SN']) for c in cmp)
    # power table
    Ns = [0.0, 2.0, 4.0, 8.0, 20.0, 50.0, 60.0]
    table = []
    for N in Ns:
        resp_own = np.zeros(5, complex); resp_codex = np.zeros(5, complex); amp = []
        for i, mm in enumerate(m2):
            S, SN = S_exact(mm, N)
            resp_own[i] = L_own[i] * S + L_own[5 + i] * SN
            resp_codex[i] = L_codex[i] * S + L_codex[5 + i] * SN
            amp.append(abs(S) ** 2)
        P_own = float(np.sum(abs(resp_own) ** 2) / (4 * np.pi ** 2)); P_codex_row = float(np.sum(abs(resp_codex) ** 2) / (4 * np.pi ** 2))
        table.append(dict(N=N, power_own_L=P_own, power_codex_L=P_codex_row, power_codex_table=CODEX_TABLE[N],
                          ratio_own_to_codex_table=P_own / CODEX_TABLE[N], ratio_codexL_to_codex_table=P_codex_row / CODEX_TABLE[N],
                          sqrt_power_own=np.sqrt(P_own), mode_abs_S_sq=amp, mode_contributions_own=(abs(resp_own) ** 2 / (4 * np.pi ** 2)).tolist()))
    out['power_table'] = table
    # e^{-3 Delta N} law on Codex's table and on ours
    law = []
    keys = sorted(CODEX_TABLE)
    for a, b in zip(keys[:-1], keys[1:]):
        r_codex = CODEX_TABLE[a] / CODEX_TABLE[b]
        r_own = [t for t in table if t['N'] == a][0]['power_own_L'] / [t for t in table if t['N'] == b][0]['power_own_L']
        law.append(dict(from_N=a, to_N=b, codex_ratio=r_codex, exp_3dN=float(np.exp(3 * (b - a))), codex_ratio_over_exp=r_codex / np.exp(3 * (b - a)),
                        own_ratio_over_exp=r_own / np.exp(3 * (b - a))))
    out['e_minus_3N_law'] = law
    # late-time analytic: |S|^2 -> e^{-3N} |c|^2 with oscillation amplitude ~ e^{-pi mu} for the lightest modes
    out['oscillation_amplitude_estimate_e^-pi_mu'] = [float(np.exp(-np.pi * np.sqrt(mm - 2.25))) for mm in m2]
    out['power_at_50_own'] = [t for t in table if t['N'] == 50.0][0]['power_own_L']
    out['power_at_50_codex_f8186fa'] = CODEX_TABLE[50.0]
    out['power_at_50_new_packet_quoted'] = NEW_PACKET_50
    out['ratio_new_packet_to_f8186fa'] = NEW_PACKET_50 / CODEX_TABLE[50.0]
    # if the new packet's number came from a different late normalisation, e.g. (S, S_N) -> conjugate pairing or adding the
    # position/velocity powers incoherently, test two candidate conventions
    inco = 0.0; posonly = 0.0
    for i, mm in enumerate(m2):
        S, SN = S_exact(mm, 50.0)
        inco += (L_own[i] * abs(S)) ** 2 + (L_own[5 + i] * abs(SN)) ** 2
        posonly += (L_own[i] * abs(S)) ** 2
    out['candidate_conventions_at_50'] = dict(incoherent_pos_plus_vel=inco / (4 * np.pi ** 2), position_only=posonly / (4 * np.pi ** 2))
    json.dump(out, open(RES / 'isocurvature.json', 'w'), indent=1, default=float)
    print(json.dumps({k: out[k] for k in out if k != 'exact_vs_ode'}, indent=1, default=float))
    print('max rel diff exact vs ode', out['max_rel_diff_exact_vs_ode'])


if __name__ == '__main__':
    main()
