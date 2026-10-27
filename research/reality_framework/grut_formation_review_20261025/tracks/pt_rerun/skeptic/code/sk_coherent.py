"""Skeptic of pt_rerun: homogeneous coherent orbit by a third, independent method (own code, Claude).

Different from Codex (sin^2 substitution in r + DOP853/FFT) and pt_rerun (cos substitution in r + long Hann series):
* u = r^2 variable: (du/dt)^2 = P(u) = -2 lam u^3 + 8 E u - 4 q^2 (cubic). Roots u1 < 0 < u2 < u3.
  T_r = 4 K(k) / sqrt(2 lam (u3-u1)), k^2 = (u3-u2)/(u3-u1)   (closed form, mpmath ellipk)
  J_r = (1/pi) int_{u2}^{u3} sqrt(P(u)) / (4u) du             (mpmath tanh-sinh, 30 digits)
  Delta_theta = 2 int_{u2}^{u3} q/u du/sqrt(P)                  (mpmath)
* lab-frequency action from our own fixed-step 8th-order symplectic composition (Kahan-Li s15odr8 of Stormer-Verlet)
  over one radial period, trapezoid-rule Fourier coefficients (spectral for periodic integrands).
* sudden switch to a free mass m: n_-(u) = ((E - lam u^2/4 + m^2 u/2)/m - q)/(2q) depends on r^2 only (energy
  conservation), so its exact range over the orbit is attained at u2, u3 or at the interior stationary point u*=m^2/lam.
* slow matching with a C-infinity ramp f(x)=psi(x)/(psi(x)+psi(1-x)), psi(x)=exp(-1/x) (all derivatives vanish at both
  ends), several durations: the deviation from J_r should fall faster than any power of T.
Usage: python3 -B sk_coherent.py <MATCH_L12.npz> <published COHERENT_MATCHING.json> <out.json>
"""
import sys, json
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp

mp.mp.dps = 30
match, pubjson, out = sys.argv[1:4]
LX = 1e-3
d = np.load(match); y0 = np.r_[d['light_mean_q'][:2], d['light_mean_p'][:2]]
pub = json.load(open(pubjson))
q = y0[0] * y0[3] - y0[1] * y0[2]
E = .5 * y0[2:] @ y0[2:] + LX / 4 * (y0[:2] @ y0[:2]) ** 2
res = dict(q=q, E=E)

# ------------------------------------------------------------------ closed-form / mpmath radial quantities in u = r^2
qm, Em, lm = mp.mpf(q), mp.mpf(E), mp.mpf(LX)
P = lambda u: -2 * lm * u ** 3 + 8 * Em * u - 4 * qm ** 2
roots = sorted([mp.re(r) for r in mp.polyroots([-2 * lm, 0, 8 * Em, -4 * qm ** 2], maxsteps=200, extraprec=60)])
u1, u2, u3 = roots
kk2 = (u3 - u2) / (u3 - u1)
Tr = 4 * mp.ellipk(kk2) / mp.sqrt(2 * lm * (u3 - u1))
Jr = mp.quad(lambda u: mp.sqrt(max(P(u), 0)) / (4 * u), [u2, (u2 + u3) / 2, u3]) / mp.pi
Dth = 2 * mp.quad(lambda u: qm / u / mp.sqrt(max(P(u), mp.mpf(0))), [u2, (u2 + u3) / 2, u3])
res['u_roots'] = [float(u1), float(u2), float(u3)]
res['rmin_rmax'] = [float(mp.sqrt(u2)), float(mp.sqrt(u3))]
res['radial_period'] = float(Tr); res['angular_advance'] = float(Dth)
res['Jr_per_charge'] = float(Jr / qm)
res['Jr_per_charge_minus_published'] = float(Jr / qm) - pub['radial_action_per_charge']
res['Tr_minus_published'] = float(Tr) - pub['orbit']['radial_period']
res['Dtheta_minus_published'] = float(Dth) - pub['orbit']['angular_advance']
res['mean_angular_frequency'] = float(Dth / Tr)

# ------------------------------------------------------------------ own symplectic integrator (Kahan-Li order 8, 15 stages)
KL8 = [0.74167036435061295344822780, -0.40910082580003159399730010, 0.19075471029623837995387626,
       -0.57386247111608226665638773, 0.29906418130365592384446354, 0.33462491824529818378495798,
       0.31529309239676659663205666, -0.79688793935291635401978884]
gam = np.array(KL8 + KL8[-2::-1])  # symmetric 15-stage composition


def verlet(x, p, h, force):
    p = p + .5 * h * force(x); x = x + h * p; p = p + .5 * h * force(x); return x, p


def kl8_step(x, p, h, force):
    for g in gam: x, p = verlet(x, p, g * h, force)
    return x, p


force0 = lambda x: -LX * (x @ x) * x
Trf = float(Tr); Dthf = float(Dth)
lab = {}
for M in (2048, 4096, 8192):
    h = Trf / M; x = y0[:2].copy(); p = y0[2:].copy(); psi = np.empty(M, complex)
    for j in range(M):
        psi[j] = x[0] + 1j * x[1]; x, p = kl8_step(x, p, h, force0)
    end_err = float(np.max(abs(np.r_[x, p] - np.r_[np.cos(Dthf) * y0[0] - np.sin(Dthf) * y0[1], np.sin(Dthf) * y0[0] + np.cos(Dthf) * y0[1],
                                                       np.cos(Dthf) * y0[2] - np.sin(Dthf) * y0[3], np.sin(Dthf) * y0[2] + np.cos(Dthf) * y0[3]])))
    t = h * np.arange(M)
    c = np.fft.fft(psi * np.exp(-1j * Dthf / Trf * t)) / M
    om = Dthf / Trf + 2 * np.pi / Trf * np.fft.fftfreq(M) * M
    neg = float(-np.sum((om * abs(c) ** 2)[om < 0]) / q); pos = float(np.sum((om * abs(c) ** 2)[om > 0]) / q)
    lab[M] = dict(lab_negative_per_charge=neg, charge_check=pos - neg, end_rotation_error=end_err,
                  minus_published=neg - pub['lab_negative_frequency_action_per_charge'])
res['lab_frequency_KL8'] = lab

# ------------------------------------------------------------------ sudden switch: exact range over the orbit
def n_sudden(u, m):
    return ((E - LX * u * u / 4 + m * m * u / 2) / m - q) / (2 * q)


U2, U3 = float(u2), float(u3); sud = []
for m in [3.678, 5., 10., 25.92, 30., 100.]:
    cand = [U2, U3]; us = m * m / LX
    if U2 < us < U3: cand.append(us)
    vals = [n_sudden(u, m) for u in cand]
    # radial-phase mean (time average over the orbit) = <n(u)> with dt = du / sqrt(P)
    num = mp.quad(lambda u: ((Em - lm * u * u / 4 + m * m * u / 2) / m - qm) / (2 * qm) / mp.sqrt(max(P(u), mp.mpf(0))), [u2, (u2 + u3) / 2, u3])
    den = mp.quad(lambda u: 1 / mp.sqrt(max(P(u), mp.mpf(0))), [u2, (u2 + u3) / 2, u3])
    sud.append(dict(m_out=m, exact_min=float(min(vals)), exact_max=float(max(vals)), interior_stationary=bool(U2 < us < U3),
                    time_average_over_radial_period=float(num / den)))
res['sudden_exact'] = sud
res['sudden_published_8phase_m30'] = [pub['summary'][0]['minimum'], pub['summary'][0]['maximum']]
res['m_for_zero_minimum'] = float(q / U2)  # (q/r^2)^2 r^2/... : n_sudden(u2,m) minimal near m = q/u2 when lam u^2 << E

# ------------------------------------------------------------------ slow matching with a C-infinity ramp
def finf(x):
    x = min(max(x, 0.), 1.)
    if x <= 0: return 0.
    if x >= 1: return 1.
    a = np.exp(-1 / x); b = np.exp(-1 / (1 - x)); return a / (a + b)


def ramp(start, T, mout, fixed=False):
    def rhs(t, v):
        f = finf(t / T); lam = LX if fixed else LX * (1 - f); r2 = v[:2] @ v[:2]
        return np.r_[v[2:], -((mout * f) ** 2 + lam * r2) * v[:2]]
    s = solve_ivp(rhs, (0, T), start, rtol=1e-12, atol=1e-10, method='DOP853')
    v = s.y[:, -1]; r2 = v[:2] @ v[:2]; qf = v[0] * v[3] - v[1] * v[2]
    Ef = .5 * v[2:] @ v[2:] + .5 * mout ** 2 * r2
    return ((Ef / mout - qf) / 2) / q, abs(qf / q - 1)


Jq = float(Jr / qm)
dense = solve_ivp(lambda t, v: np.r_[v[2:], -LX * (v[:2] @ v[:2]) * v[:2]], (0, Trf), y0, rtol=1e-13, atol=1e-11, method='DOP853', dense_output=True)
starts = [dense.sol(ph * Trf) for ph in np.arange(8) / 8]
slow = []
for T in [5., 10., 20., 40., 80.]:
    vals = [ramp(st, T, 30.) for st in starts]
    a = np.array([v[0] for v in vals])
    slow.append(dict(T=T, m_out=30., min=float(a.min()), max=float(a.max()), max_abs_dev_from_Jr=float(np.max(abs(a - Jq))),
                     charge_rel=float(max(v[1] for v in vals))))
    print('slow Cinf', slow[-1], flush=True)
res['slow_Cinfinity_ramp'] = slow
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != 'slow_Cinfinity_ramp'}, indent=1))
