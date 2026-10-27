"""Independent re-derivation of the homogeneous coherent-matching numbers (own code, Claude pt_rerun).

* J_r by direct quadrature of p_r(r) in r (Gauss-Chebyshev-type endpoint handling via scipy quad with
  algebraic weights), independent of Codex's u-substitution;
* lab-frequency action via a long Hann-windowed time series (different method from Codex's one-period Fourier series);
* slow ramps to a free mass with a different ramp shape (smoothstep) and several final masses m_out;
* sudden switches to several m_out (shows that the sudden range depends on m_out);
* a harmonic check: free orbit psi = A e^{ims} + B e^{-ims} has J_r = m|B|^2 exactly.
Usage: python3 -B coherent_indep.py <early_packet_dir> <published COHERENT_MATCHING.json> <out.json>
"""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

early, pubjson, out = sys.argv[1:4]
LX = 1e-3
d = np.load(f'{early}/MATCH_L12.npz'); y0 = np.r_[d['light_mean_q'][:2], d['light_mean_p'][:2]]
q = y0[0] * y0[3] - y0[1] * y0[2]; E0 = .5 * y0[2:] @ y0[2:] + LX / 4 * (y0[:2] @ y0[:2]) ** 2
pub = json.load(open(pubjson))
res = dict(charge=q, energy=E0, published_radial_action_per_charge=pub['radial_action_per_charge'],
           published_lab_negative_per_charge=pub['lab_negative_frequency_action_per_charge'])


def Jr(E, q, lam, m):
    f = lambda r: 2 * E - q * q / r ** 2 - m * m * r * r - lam * r ** 4 / 2
    # turning points by bracketing around the circular-orbit radius (minimum of the effective potential)
    rc = brentq(lambda r: -q * q / r ** 3 + m * m * r + lam * r ** 3, 1e-6, 1e6)
    assert f(rc) > 0
    r1 = brentq(f, 1e-9 * rc, rc); r2 = brentq(f, rc, 1e6 * rc)
    # p_r = sqrt(f); integrate with the square-root endpoint behaviour removed by r = r1 + (r2-r1) (1-cos t)/2
    g = lambda t: np.sqrt(max(f(r1 + (r2 - r1) * (1 - np.cos(t)) / 2), 0.)) * (r2 - r1) * np.sin(t) / 2
    val = quad(g, 0, np.pi, epsabs=1e-12, epsrel=1e-13, limit=200)[0]
    return val / np.pi, r1, r2


J, r1, r2 = Jr(E0, q, LX, 0.)
res['Jr_per_charge_own'] = J / q; res['Jr_diff_vs_published'] = J / q - pub['radial_action_per_charge']; res['rmin_rmax'] = [r1, r2]

# harmonic sanity: free orbit, J_r = m|B|^2
m = 2.3; A = 1.7 * np.exp(.4j); B = .6 * np.exp(-1.1j)
psi0 = A + B; dpsi0 = 1j * m * (A - B)
qh = (psi0.conjugate() * dpsi0).imag; Eh = .5 * abs(dpsi0) ** 2 + .5 * m * m * abs(psi0) ** 2
res['harmonic_check'] = dict(Jr=Jr(Eh, qh, 0., m)[0], m_absB2=m * abs(B) ** 2, q=qh, m_absA2_minus_m_absB2=m * (abs(A) ** 2 - abs(B) ** 2))

# lab-frequency action from a long Hann window (different method)
def rhs(t, v): return np.r_[v[2:], -LX * (v[:2] @ v[:2]) * v[:2]]
ds = .01; T = 400.; tt = np.arange(0, T + 1e-9, ds)
sol = solve_ivp(rhs, (0, T), y0, t_eval=tt, rtol=1e-12, atol=1e-10, method='DOP853')
psi = sol.y[0] + 1j * sol.y[1]; M = len(psi); h = np.hanning(M); w = 2 * np.pi * np.fft.fftfreq(M, ds)
P = abs(np.fft.fft(h * psi)) ** 2 / (M * np.sum(h * h))
res['lab_negative_per_charge_hann400'] = float(np.sum((-w * P)[w < 0]) / q)
res['lab_positive_minus_negative_per_charge'] = float((np.sum((w * P)[w > 0]) - np.sum((-w * P)[w < 0])) / q)


# ramps: own smoothstep shape f(u)=u^2(3-2u), m(s)=m_out f, lambda(s)=LX(1-f)
def ramp(start, duration, mout, fixed_quartic=False):
    def flow(t, v):
        u = t / duration; f = u * u * (3 - 2 * u); lam = LX if fixed_quartic else LX * (1 - f)
        r2 = v[:2] @ v[:2]
        return np.r_[v[2:], -((mout * f) ** 2 + lam * r2) * v[:2]]
    ss = solve_ivp(flow, (0, duration), start, rtol=1e-11, atol=1e-9, method='DOP853')
    fin = ss.y[:, -1]; r2 = fin[:2] @ fin[:2]
    lam_end = LX if fixed_quartic else 0.
    E = .5 * fin[2:] @ fin[2:] + .5 * mout ** 2 * r2 + lam_end / 4 * r2 * r2
    qf = fin[0] * fin[3] - fin[1] * fin[2]
    nminus_quadratic = ((.5 * fin[2:] @ fin[2:] + .5 * mout ** 2 * r2) / mout - qf) / 2
    J_end = Jr(E, qf, lam_end, mout)[0]
    return dict(nminus_quadratic_per_charge=nminus_quadratic / q, Jr_end_per_charge=J_end / q, charge_rel=abs(qf / q - 1))


# start phases along the orbit
Tr = pub['orbit']['radial_period']
starts = [solve_ivp(rhs, (0, p * Tr), y0, rtol=1e-12, atol=1e-10, method='DOP853').y[:, -1] if p > 0 else y0 for p in np.arange(8) / 8]
slow = []
for mout in [10., 30., 100.]:
    for duration in [80., 200.]:
        vals = [ramp(st, duration, mout) for st in starts[::2]]
        a = np.array([v['nminus_quadratic_per_charge'] for v in vals])
        slow.append(dict(m_out=mout, duration=duration, nminus_min=float(a.min()), nminus_max=float(a.max()),
                         max_dev_from_Jr=float(np.max(abs(a - J / q))), charge_rel_max=max(v['charge_rel'] for v in vals)))
        print('slow', slow[-1], flush=True)
res['slow_smoothstep_ramps'] = slow
fq = [ramp(st, 200., 30., fixed_quartic=True) for st in starts[::2]]
res['fixed_quartic_smoothstep_m30_T200'] = dict(Jr_end_minus_Jr0_max=float(max(abs(v['Jr_end_per_charge'] - J / q) for v in fq)),
                                                 quadratic_projection_range=[float(min(v['nminus_quadratic_per_charge'] for v in fq)), float(max(v['nminus_quadratic_per_charge'] for v in fq))])

# sudden switch to free mass m_out at 64 phases
phases = np.arange(64) / 64
dense = solve_ivp(rhs, (0, Tr), y0, rtol=1e-12, atol=1e-10, method='DOP853', dense_output=True)
sud = []
for mout in [3.68, 5., 10., 30., 100.]:
    vals = []
    for p in phases:
        st = dense.sol(p * Tr); r2 = st[:2] @ st[:2]
        vals.append(((.5 * st[2:] @ st[2:] + .5 * mout ** 2 * r2) / mout - q) / 2 / q)
    vals = np.array(vals)
    sud.append(dict(m_out=mout, min=float(vals.min()), max=float(vals.max()), phase_mean=float(vals.mean())))
res['sudden_switch'] = sud
json.dump(res, open(out, 'w'), indent=1, default=float)
print(json.dumps({k: v for k, v in res.items() if k not in ('slow_smoothstep_ramps',)}, indent=1, default=float))
