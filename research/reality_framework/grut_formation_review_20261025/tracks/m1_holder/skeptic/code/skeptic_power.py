#!/usr/bin/env python3
"""m1_holder_skeptic (Claude, adversarial reviewer; own code). Independent check of the spectator charge-power table
(Codex f8186fa perturbations.json; track m1_holder isocurvature.py) by two methods that differ from both:
  (1) closed form: small-x expansion of the de Sitter Hankel mode, |S|^2 -> e^{-3N} e^{pi mu}/(2 mu sinh pi mu) *
      (1 - 2 e^{-pi mu} cos psi + e^{-2 pi mu}), S_N -> -(3/2 + i mu) S for the dominant branch, so that
      P(N) = e^{-3N}/(4 pi^2) sum_a [(L_pos,a - 1.5 L_vel,a)^2 + mu_a^2 L_vel,a^2] / (mu_a (1 - e^{-2 pi mu_a}))  (leading)
      and the two-branch form with the Gamma-function phases (next order);
  (2) ODE in the conformal variable x = k/(aH) for chi = S/x: chi'' + (1 + (m^2 - 2)/x^2) chi = 0, WKB-initialised at
      x0 = 2000 (|chi| -> 1), integrated inward with DOP853; neither the N-variable (Codex, track) nor mpmath Hankel
      functions (track) are used.
Rows: Codex's physical position/velocity response [Q perturbations.json] and the track's own FD row [results/charge_transfer.json].
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import loggamma

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
TRACK = HERE.parent.parent / 'results'
CODEX_TABLE = {0.0: 1.9858166793719456e-05, 2.0: 4.984135128883863e-08, 4.0: 1.2357522432943231e-10, 8.0: 7.592347839932724e-16,
               20.0: 1.7611900904580933e-31, 50.0: 1.4430133116013324e-70, 60.0: 1.3503760669554607e-83}
CODEX_M2 = np.array([13.091408829452904, 14.589249978666942, 538.146643427554, 1546.659721973107, 3239.8852511586983])
CODEX_ROW = np.array([-0.015482021001551589, -0.030182507194285172, 0.0018709020686744039, -0.00037602803986460016, -7.819116991468326e-06,
                      -0.006355186321161243, 0.00707974696172499, -6.150069853994815e-05, 6.805619706250102e-07, -3.8916041663512055e-08])
ct = json.load(open(TRACK / 'charge_transfer.json'))['response_row']
TRACK_ROW = np.array(ct['own_richardson'])
TRACK_TABLE = {0.0: 1.9858019408836794e-05, 2.0: 4.984097943550111e-08, 4.0: 1.2357430227682094e-10, 8.0: 7.592291189755589e-16,
               20.0: 1.7611769494668086e-31, 50.0: 1.4430025458437605e-70, 60.0: 1.3503659971971494e-83}
mu = np.sqrt(CODEX_M2 - 2.25)


def P_leading(N, L):
    Lp, Lv = L[:5], L[5:]
    return np.exp(-3 * N) / (4 * np.pi ** 2) * np.sum(((Lp - 1.5 * Lv) ** 2 + mu ** 2 * Lv ** 2) / (mu * (1 - np.exp(-2 * np.pi * mu))))


def S_two_branch(N, m2):
    """small-x two-branch expansion (no O(x^2) corrections): returns S, S_N."""
    mu_ = np.sqrt(m2 - 2.25); x = np.exp(-N)
    g_plus = np.exp(loggamma(1 + 1j * mu_)); g_minus = np.exp(loggamma(1 - 1j * mu_))
    pref = np.sqrt(np.pi / 2) * np.exp(-np.pi * mu_ / 2) * (-1 / np.sinh(np.pi * mu_))
    t1 = (x / 2) ** (-1j * mu_) / g_minus; t2 = np.exp(np.pi * mu_) * (x / 2) ** (1j * mu_) / g_plus
    S = pref * x ** 1.5 * (t1 - t2)
    S_N = -pref * x ** 1.5 * ((1.5 - 1j * mu_) * t1 - (1.5 + 1j * mu_) * t2)      # -x d/dx of x^{3/2 -/+ i mu}
    return S, S_N


def P_two_branch(N, L):
    tot = 0.0
    for a in range(5):
        S, SN = S_two_branch(N, CODEX_M2[a])
        tot += abs(L[a] * S + L[5 + a] * SN) ** 2
    return tot / (4 * np.pi ** 2)


def modes_ode(m2, Ns, x0=2000.0):
    """chi = S/x, chi'' + (1 + (m2-2)/x^2) chi = 0 integrated inward from x0 with WKB start; returns S, S_N at e^{-N}."""
    b = m2 - 2.0
    w = np.sqrt(1 + b / x0 ** 2); wp = -b / (w * x0 ** 3)
    chi0 = w ** -0.5 * (1 + 0j); dchi0 = chi0 * (1j * w - 0.5 * wp / w)
    xs = np.exp(-np.asarray(Ns))
    sol = solve_ivp(lambda x, y: [y[1], -(1 + b / x ** 2) * y[0]], (x0, xs.min()), [chi0, dchi0], method='DOP853', rtol=1e-12, atol=1e-16,
                    t_eval=np.sort(xs)[::-1])
    out = {}
    for x, chi, dchi in zip(sol.t, sol.y[0], sol.y[1]):
        S = x * chi; S_x = chi + x * dchi; out[float(-np.log(x))] = (S, -x * S_x)
    return out


def main():
    out = dict(mu=mu.tolist(), relative_oscillation_amplitude_of_S2=(2 * np.exp(-np.pi * mu)).tolist(),
               late_S2_prefactor_e3N=(np.exp(np.pi * mu) / (2 * mu * np.sinh(np.pi * mu))).tolist())
    rows = []
    for N in sorted(CODEX_TABLE):
        rows.append(dict(N=N, codex_table=CODEX_TABLE[N], track_own_row=TRACK_TABLE[N],
                         leading_closed_form_codex_row=P_leading(N, CODEX_ROW), two_branch_codex_row=P_two_branch(N, CODEX_ROW),
                         leading_closed_form_track_row=P_leading(N, TRACK_ROW), two_branch_track_row=P_two_branch(N, TRACK_ROW)))
        rows[-1]['leading_over_codex_table'] = rows[-1]['leading_closed_form_codex_row'] / CODEX_TABLE[N]
        rows[-1]['two_branch_over_codex_table'] = rows[-1]['two_branch_codex_row'] / CODEX_TABLE[N]
        rows[-1]['two_branch_track_row_over_track_table'] = rows[-1]['two_branch_track_row'] / TRACK_TABLE[N]
    out['closed_form_rows'] = rows
    # ODE in x for N = 0, 2, 4, 8 (full equation, gradient term included)
    Ns = [0.0, 2.0, 4.0, 8.0]
    S = {N: np.zeros(5, complex) for N in Ns}; SN = {N: np.zeros(5, complex) for N in Ns}
    for a in range(5):
        res = modes_ode(CODEX_M2[a], Ns)
        for N in Ns:
            S[N][a], SN[N][a] = res[N]
    ode_rows = []
    for N in Ns:
        Pc = float(np.sum(np.abs(CODEX_ROW[:5] * S[N] + CODEX_ROW[5:] * SN[N]) ** 2) / (4 * np.pi ** 2))
        Pt = float(np.sum(np.abs(TRACK_ROW[:5] * S[N] + TRACK_ROW[5:] * SN[N]) ** 2) / (4 * np.pi ** 2))
        ode_rows.append(dict(N=N, ode_x_codex_row=Pc, codex_table=CODEX_TABLE[N], ratio_codex=Pc / CODEX_TABLE[N],
                             ode_x_track_row=Pt, track_table=TRACK_TABLE[N], ratio_track=Pt / TRACK_TABLE[N],
                             mode_abs_S2=(np.abs(S[N]) ** 2).tolist()))
    out['ode_x_rows'] = ode_rows
    # massless normalisation check of the x-ODE: m2 = 0 -> |S|^2 -> 1 + x^2
    r0 = modes_ode(0.0, [6.0])[6.0]
    out['massless_check_abs_S2_at_N6'] = float(abs(r0[0]) ** 2)
    # how much of the power is velocity response: fraction of the leading closed form from mu^2 L_vel^2 terms
    Lp, Lv = CODEX_ROW[:5], CODEX_ROW[5:]
    terms = ((Lp - 1.5 * Lv) ** 2 + mu ** 2 * Lv ** 2) / (mu * (1 - np.exp(-2 * np.pi * mu)))
    out['velocity_fraction_of_late_power'] = float(np.sum(mu ** 2 * Lv ** 2 / (mu * (1 - np.exp(-2 * np.pi * mu)))) / np.sum(terms))
    out['mode_fraction_of_late_power'] = (terms / terms.sum()).tolist()
    # the 3.435e-70 of the unpushed packet: ratios to candidate quantities at N = 50
    P50 = P_leading(50.0, CODEX_ROW)
    out['unpushed_3p435e70'] = dict(value=3.435e-70, over_f8186fa=3.435e-70 / CODEX_TABLE[50.0], over_leading_closed_form=3.435e-70 / P50,
                                    note='not derivable from f8186fa data; packet absent from origin on 2026-11-01')
    json.dump(out, open(RES / 'skeptic_power.json', 'w'), indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
