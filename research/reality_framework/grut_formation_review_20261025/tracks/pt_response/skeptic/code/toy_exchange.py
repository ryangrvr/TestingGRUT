"""SKEPTIC (pt_response_skeptic), INDEPENDENT toy: is the X <-> Xbar 'conversion' of waves by the static anomalous term
of an eccentric condensate a loss of the Xbar budget, or an exchange with the condensate?

Closed Hamiltonian 3-mode Galerkin truncation of the complex quartic field (lam |psi|^4 / 4, lam = 1e-3):
    psi(x, s) = psi_0(s) + psi_+(s) e^{ikx} + psi_-(s) e^{-ikx},
force projected exactly onto the three modes (8-point FFT, alias-free for cubic terms).  Charge and energy are
exactly conserved by the truncation.  The condensate psi_0 starts on an eccentric orbit with r_max/r_min ~ 3.9
(as measured near s = 20); the waves start as pure Xbar (lab e^{-i Omega s}) with small amplitude eps.
Codex's statistic, A_-^{(j)} = sum_{omega<0} |omega| P_j(omega) (Hann window, length W), is evaluated per mode in
successive windows.  The condensate change is measured against an eps = 0 reference run (same orbit).
If the static anomalous term destroyed Xbar, A_-(waves) + Delta A_-(condensate) would swing with the wave
conversion; if it is the 2<->2 exchange Xbar_w + X_c -> X_w + Xbar_c, the sum stays constant.
Usage: python3 -I toy_exchange.py <out json>
"""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp

OUT = sys.argv[1]
LX = 1e-3; NX = 8
x = 2 * np.pi * np.arange(NX) / NX
E = np.array([np.ones(NX), np.exp(1j * x), np.exp(-1j * x)])  # mode functions

def rhs_factory(k2):
    kk = np.array([0., k2, k2])
    def f(t, y):
        c = y[:3] + 1j * y[3:6]; p = y[6:9] + 1j * y[9:12]
        ps = c @ E
        F = LX * np.abs(ps) ** 2 * ps
        Fm = (E.conj() @ F) / NX
        acc = -kk * c - Fm
        return np.concatenate([p.real, p.imag, acc.real, acc.imag])
    return f

def statistic(sig, ds):
    n = len(sig); h = np.hanning(n); X = np.fft.fft(h * sig); P = np.abs(X) ** 2 / (n * np.sum(h * h)); w = 2 * np.pi * np.fft.fftfreq(n, ds)
    return float(np.sum((-w * P)[w < 0])), float(np.sum((w * P)[w > 0]))

def run(k2, eps, rmax=99.4, rmin=25.5, T=120., ds=.01, phase=(0.3, 1.7)):
    q2 = (LX / 4 * (rmax ** 4 - rmin ** 4)) / (.5 * (1 / rmin ** 2 - 1 / rmax ** 2)); q = np.sqrt(q2)
    Om = np.sqrt(k2 + 2 * LX * .5 * (rmax ** 2 + rmin ** 2))
    c0 = np.array([rmax, eps * np.exp(1j * phase[0]), eps * np.exp(1j * phase[1])], complex)
    p0 = np.array([1j * q / rmax, -1j * Om * c0[1], -1j * Om * c0[2]], complex)  # waves: pure e^{-i Om s} (Xbar)
    y0 = np.concatenate([c0.real, c0.imag, p0.real, p0.imag])
    t = np.arange(0, T + ds / 2, ds)
    sol = solve_ivp(rhs_factory(k2), (0, T), y0, t_eval=t, method='DOP853', rtol=1e-12, atol=1e-9)
    c = sol.y[:3] + 1j * sol.y[3:6]; p = sol.y[6:9] + 1j * sol.y[9:12]
    Q = np.sum((c.conj() * p).imag, axis=0)
    return t, c, p, Q

def windows(t, c, W=10.):
    ds = t[1] - t[0]; out = []
    for lo in np.arange(0, t[-1] - W + 1e-9, W / 2):
        sel = (t >= lo - 1e-9) & (t <= lo + W + 1e-9)
        out.append([lo] + [statistic(c[j, sel], ds) for j in range(3)])
    return out

res = dict(method=__doc__, cases=[])
for k2, eps in ((20., 2.), (20., 5.), (3., 2.), (45., 2.)):
    t, c, p, Q = run(k2, eps); t0, c0_, p0_, Q0 = run(k2, 0.)
    wv = windows(t, c); w0 = windows(t0, c0_)
    rows = []
    for (lo, s0, sp, sm), (_, r0, _, _) in zip(wv, w0):
        A_w = sp[0] + sm[0]; Ap_w = sp[1] + sm[1]
        dA_c = s0[0] - r0[0]; dAp_c = s0[1] - r0[1]
        rows.append(dict(window_start=float(lo), waves_A_minus=A_w, waves_A_plus=Ap_w, cond_dA_minus=dA_c, cond_dA_plus=dAp_c,
                         total_dA_minus_excess=A_w + dA_c))
    Aw = np.array([r['waves_A_minus'] for r in rows]); tot = np.array([r['total_dA_minus_excess'] for r in rows])
    wch = np.array([r['waves_A_plus'] - r['waves_A_minus'] for r in rows])
    case = dict(k2=k2, eps=eps, charge_drift=float(np.ptp(Q) / abs(Q[0])), waves_A_minus_range=[float(Aw.min()), float(Aw.max())],
                total_excess_range=[float(tot.min()), float(tot.max())], waves_charge_range=[float(wch.min()), float(wch.max())], rows=rows)
    res['cases'].append(case)
    print('k2 %.0f eps %.0f : charge drift %.1e ; waves A- %.3f..%.3f ; waves charge %.3f..%.3f ; waves A- + dA-(cond) %.3f..%.3f' % (
        k2, eps, case['charge_drift'], Aw.min(), Aw.max(), wch.min(), wch.max(), tot.min(), tot.max()))
json.dump(res, open(OUT, 'w'), indent=1)
