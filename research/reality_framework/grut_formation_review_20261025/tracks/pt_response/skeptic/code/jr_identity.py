"""SKEPTIC (pt_response_skeptic), INDEPENDENT code.  Radial action vs lab-frequency statistic of a coherent orbit.

[D] For a quasi-periodic planar orbit psi(s) = sum_n c_n exp(i w_n s), w_n = w_th + n W_r (W_r = 2pi/T_r, w_th = mean
angular velocity), the action-angle identity J_j = <p . dq/dtheta_j> gives
    J_theta = sum_n w_n |c_n|^2 = q (charge),    J_r = sum_n n w_n |c_n|^2 = (<|psi'|^2> - w_th q) / W_r,
and Codex's lab statistic is A_- = sum_{w_n<0} |w_n| |c_n|^2.  When w_0 > 0 > w_{-1} (true here):
    J_r - A_- = sum_{n>=1} n w_n |c_n|^2 + sum_{n<=-2} |n+1| |w_n| |c_n|^2  >= 0,
so the radial action always exceeds the lab statistic, by the weight of the higher harmonics (n >= 1 and n <= -2).
Checks: (1) the reviewed starting orbit (E0, q from COHERENT_MATCHING.json), own DOP853 orbit, three routes to J_r:
loop integral, Fourier identity, <|psi'|^2> formula; compare Codex 1.763725 / 1.737788.  (2) In the medium: J_r per
radial period from FIELD_TRACES (base48) with three weightings (uniform, Hann, Hann^2) and its drift within each window,
against the lab coherent statistic, to size the window/weighting systematic of 'J_r exceeds lab by 0.8-2.7%'.
Usage: python3 -I jr_identity.py <packet dir> <out json>
"""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.signal import argrelmin
P = sys.argv[1].rstrip('/') + '/'; OUT = sys.argv[2]
cm = json.load(open(P + 'COHERENT_MATCHING.json'))
LX = 1e-3; E0 = cm['initial_energy']; q = cm['initial_charge']
# start at r_max on the x axis
f = lambda r: q * q / (2 * r * r) + LX * r ** 4 / 4 - E0
from scipy.optimize import brentq
rmax = brentq(f, 50, 1000); rmin = brentq(f, 1, rmax * .99)
def rhs(t, y): r2 = y[0] ** 2 + y[1] ** 2; return [y[2], y[3], -LX * r2 * y[0], -LX * r2 * y[1]]
y0 = [rmax, 0., 0., q / rmax]
ev = lambda t, y: y[0] * y[2] + y[1] * y[3]; ev.direction = -1
sol = solve_ivp(rhs, (0, 3), y0, rtol=1e-13, atol=1e-10, method='DOP853', events=ev, dense_output=True)
Tr = sol.t_events[0][1] if abs(sol.t_events[0][0]) < 1e-9 else sol.t_events[0][0]
M = 16384; t = np.linspace(0, Tr, M, endpoint=False); Y = sol.sol(t)
psi = Y[0] + 1j * Y[1]; dpsi = Y[2] + 1j * Y[3]
th = np.unwrap(np.angle(psi)); Yend = sol.sol(Tr); dth = np.unwrap([th[-1], np.angle(Yend[0] + 1j * Yend[1])])[1] - th[0]
wth = dth / Tr; Wr = 2 * np.pi / Tr
cn = np.fft.fft(psi * np.exp(-1j * wth * t)) / M; nn = np.fft.fftfreq(M) * M; wn = wth + nn * Wr
P2 = np.abs(cn) ** 2
Jth = float(np.sum(wn * P2)); Jr_f = float(np.sum(nn * wn * P2)); Am = float(np.sum((-wn * P2)[wn < 0])); Ap = float(np.sum((wn * P2)[wn > 0]))
Jr_avg = float((np.mean(np.abs(dpsi) ** 2) - wth * q) / Wr)
r = np.abs(psi); pr = (Y[0] * Y[2] + Y[1] * Y[3]) / r
Jr_loop = float(np.trapezoid(np.r_[pr ** 2, pr[0] ** 2], np.r_[t, Tr]) / (2 * np.pi))
excess_pos = float(np.sum((nn * wn * P2)[nn >= 1])); excess_neg = float(np.sum((np.abs(nn + 1) * np.abs(wn) * P2)[nn <= -2]))
res = dict(starting_orbit=dict(E0=E0, q=q, rmin=rmin, rmax=rmax, radial_period=float(Tr), angular_advance=float(dth),
           J_theta_over_q=Jth / q, Jr_fourier_over_q=Jr_f / q, Jr_average_formula_over_q=Jr_avg / q, Jr_loop_over_q=Jr_loop / q,
           lab_A_minus_over_q=Am / q, lab_A_plus_over_q=Ap / q, Jr_minus_A_over_q=(Jr_f - Am) / q,
           harmonic_excess_over_q=dict(n_ge_1=excess_pos / q, n_le_m2=excess_neg / q),
           codex=dict(radial_action_per_charge=cm['radial_action_per_charge'], lab=cm['lab_negative_frequency_action_per_charge'],
                      radial_period=cm['orbit']['radial_period'], angular_advance=cm['orbit']['angular_advance']),
           top_harmonics=[(int(nn[i]), float(wn[i]), float(P2[i] / P2.max())) for i in np.argsort(P2)[::-1][:6]]))
print(json.dumps(res['starting_orbit'], indent=1))
# (2) medium
F = np.load(P + 'FIELD_TRACES.npz'); names = list(F['names'])
s = F['base48_s']; v = F['base48_values']
g = {k: v[:, names.index(k)] for k in ('q0', 'q1', 'p0', 'p1', 'nQ')}
nQ = float(np.mean(g['nQ']))
win = {}
for lo, hi in ((0., 6.), (6., 12.), (14., 24.), (24., 34.)):
    sel = (s >= lo - 1.5) & (s <= hi + 1.5)
    sp = {k: CubicSpline(s[sel], g[k][sel]) for k in ('q0', 'q1', 'p0', 'p1')}
    sf = np.arange(max(lo - 1, s[0]), min(hi + 1, s[-1]), .001)
    a, b, pa, pb = [sp[k](sf) for k in ('q0', 'q1', 'p0', 'p1')]
    rr = np.hypot(a, b); prr = (a * pa + b * pb) / rr
    mins = argrelmin(rr)[0]
    J = []; mid = []
    for i, j in zip(mins[:-1], mins[1:]):
        J.append(np.trapezoid(prr[i:j + 1] ** 2, sf[i:j + 1]) / (2 * np.pi)); mid.append(.5 * (sf[i] + sf[j]))
    J = np.array(J) / nQ; mid = np.array(mid); inw = (mid > lo) & (mid < hi)
    hann = np.sin(np.pi * (mid[inw] - lo) / (hi - lo)) ** 2
    wd = json.load(open(P + f'WINDOW_base48_{lo:g}_{hi:g}.json'))
    lab = wd['coherent_negative_frequency_action']
    vals = dict(uniform=float(J[inw].mean()), hann=float(np.sum(hann * J[inw]) / hann.sum()), hann2=float(np.sum(hann ** 2 * J[inw]) / np.sum(hann ** 2)))
    win[f'{lo:g}_{hi:g}'] = dict(lab_coherent=lab, Jr=vals, excess_percent={k: 100 * (x / lab - 1) for k, x in vals.items()},
                                 Jr_first_last_period=[float(J[inw][0]), float(J[inw][-1])], periods=int(inw.sum()),
                                 within_window_drift_percent=float(100 * (J[inw][-1] / J[inw][0] - 1)))
    print(lo, hi, 'lab %.4f' % lab, {k: round(x, 4) for k, x in vals.items()}, 'first/last %.3f/%.3f' % tuple(win[f'{lo:g}_{hi:g}']['Jr_first_last_period']))
res['medium_base48'] = win
json.dump(res, open(OUT, 'w'), indent=1)
