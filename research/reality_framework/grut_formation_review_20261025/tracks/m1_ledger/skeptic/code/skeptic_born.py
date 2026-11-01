#!/usr/bin/env python3
"""m1_ledger_skeptic (Claude, own code): the vacuum-Born X^6 pair integral on the M1 and old profiles with a
different numerical method from both Codex (CubicSpline to linspace(0,70,12001), Simpson; q = qmax sin(theta)
Gauss-Legendre) and m1_ledger (CubicSpline to dr = 0.01 up to 80, Simpson; q = qmax(1-u^2) Gauss-Legendre).

Here: (A) the profile's NATIVE grid (no interpolation; M1 dr = 0.05 to 160, old dr = 0.075 to 120), Simpson for
the spherical transform Ft4(q) = 4 pi int r^2 F^4 j0(qr) dr, and scipy.integrate.quad with the algebraic endpoint
weight (qmax - q)^{1/2} for the q integral (beta = sqrt((qmax-q)(qmax+q)/(W^2-q^2))); (B) a trapezoid variant on a
2x refined spline grid; (C) the Parseval identity int_0^inf q^2 |Ft4|^2 dq = 2 pi^2 int r^2 F^8 dr as an absolute
normalisation check of the transform; (D) the rate prefactor re-derived (see NOTES) and the c6 = 1 rate.

Own derivation of the prefactor (same as m1_ledger's and Codex's, done again from scratch):
  L = A X^6 + h.c., X = X0 + chi, X0 = (v/sqrt2) F e^{-i Omega t}.  The chi^2 term is C(6,2) A X0^4 chi^2 = 15 A X0^4 chi^2;
  chi contains b^dagger, so chi^2 creates two anti-X with <k1 k2|chi^2|0> = 2 e^{i(k1+k2).x}/sqrt(2E1 2E2 (2pi)^6).
  Amplitude: 30 A (v/sqrt2)^4 Ft4(q) 2pi delta(E1+E2-W), W = 4 Omega m_X, q = k1+k2.
  Rate = (1/2)_identical int d^3k1 d^3k2/((2pi)^6 2E1 2E2) |30 A (v/sqrt2)^4 Ft4|^2 2pi delta(E1+E2-W)
       = 450 A^2 (v/sqrt2)^8 (2pi/(2pi)^6) int d^3q |Ft4(q)|^2 (pi/2) beta(s),  s = W^2 - q^2, beta = sqrt(1-4m^2/s)
       = 450 A^2 (v/sqrt2)^8 (4 pi^3 / (64 pi^6)) I = (900/(32 pi^3)) A^2 (v/sqrt2)^8 I / m_X^3,  I = int q^2 |Ft4|^2 beta dq.
  The two-body factor: int d^3k1 d^3k2/(2E1 2E2) delta^4(P-k1-k2) = (pi/2) beta  (CM frame: 4pi k^2/(4E^2) x E/(2k)).
  Each event absorbs 4 condensate quanta (core charge -4) and emits two anti-X (charge -2): X-number violation 6 = Z6.
Output: ../results/skeptic_born.json
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad, simpson
from scipy.interpolate import CubicSpline

SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/m1_ledger_skeptic/data')
OUT = Path(__file__).resolve().parents[1] / 'results'
HBAR = 6.582119569e-25; MRED = 2.435e18; V = 0.924089637375362; MX = 0.1 * V; Q = 200000; AGE = 4.3549488e17

PROF = {'M1': dict(csv=SP / 'PROFILE_M1.csv', Omega=0.8817204919524214, codex_I=2868466.5109355855, codex_born=9.471343360338484e-47,
                   codex_events=4.735671680169242e-42, m1_I=2868466.5109466882, m1_born=9.471343360375144e-47),
        'old': dict(csv=SP / 'profile_old.csv', Omega=0.8611951308217225, codex_I=None, codex_born=1.0399701293581432e-46,
                    codex_events=None, m1_I=3149626.589328209, m1_born=1.0399701293584032e-46)}


def simpson_weights(n, h):
    w = np.ones(n); w[1:-1:2] = 4; w[2:-1:2] = 2
    return w * h / 3


def make_transform(r, F4):
    """Ft4(q) on a uniform grid with Simpson weights; returns a scalar-callable and a vector version."""
    h = r[1] - r[0]
    assert len(r) % 2 == 1, 'Simpson needs an odd number of points'
    base = simpson_weights(len(r), h) * r * r * F4
    def ft_scalar(q):
        return 4 * np.pi * float(np.sum(base * np.sinc(q * r / np.pi)))
    def ft_vec(q):
        return 4 * np.pi * np.array([np.sum(base * np.sinc(qq * r / np.pi)) for qq in q])
    return ft_scalar, ft_vec


def born_native(prof):
    a = np.loadtxt(prof['csv'], delimiter=',', skiprows=1)
    r, F = a[:, 0], a[:, 1]
    dr = np.diff(r)
    assert np.allclose(dr, dr[0], rtol=1e-9), 'native grid must be uniform'
    W = 4 * prof['Omega']; qmax = np.sqrt(W * W - 4)
    ft, ftv = make_transform(r, F ** 4)
    # (A) adaptive quad with the algebraic endpoint weight (qmax - q)^{1/2}; integrand = q^2 |Ft4|^2 sqrt((qmax+q)/(W^2-q^2))
    g = lambda q: q * q * ft(q) ** 2 * np.sqrt((qmax + q) / (W * W - q * q))
    I_quad, err = quad(g, 0.0, qmax, weight='alg', wvar=(0.0, 0.5), limit=400, epsabs=0, epsrel=1e-12)
    # (C) Parseval: int d^3q/(2pi)^3 |Ft4|^2 = int d^3x F^8  =>  int_0^inf q^2 |Ft4|^2 dq = 8 pi^3 int r^2 F^8 dr.
    # The native grid aliases above q ~ 1/dr, so the transform for this check is built on a 10x refined cubic-spline
    # grid (dr_native/10, Simpson) and integrated to q = 60 (phase step 0.3 rad per point there); Gaussian control below.
    spl = CubicSpline(r, F); rf = np.linspace(0, r[-1], 10 * (len(r) - 1) + 1); Ff = spl(rf)
    ftf, _ = make_transform(rf, Ff ** 4)
    lhs, lerr = quad(lambda q: q * q * ftf(q) ** 2, 0.0, 60.0, limit=800, epsabs=0, epsrel=1e-10)
    tail = quad(lambda q: q * q * ftf(q) ** 2, 60.0, 120.0, limit=200)[0]
    rhs = 8 * np.pi ** 3 * simpson(rf * rf * Ff ** 8, x=rf)
    pref = 900 / (32 * np.pi ** 3) * (V / np.sqrt(2)) ** 8 / (MRED ** 4 * MX ** 3)      # A = 1/MRED^2 (c6 = 1)
    events = pref * I_quad / HBAR
    return dict(grid=dict(points=len(r), dr=float(dr[0]), rmax=float(r[-1])), W=W, qmax=qmax, I=I_quad, quad_abs_err=err,
                parseval_lhs=lhs, parseval_tail_60_120=tail, parseval_rhs=rhs, parseval_rel=(lhs + tail) / rhs - 1,
                events_per_second_c1=events, born_fractional_core_rate_c1=4 * events / Q, F0=float(F[0]),
                codex_I=prof['codex_I'], rel_I_codex=(I_quad / prof['codex_I'] - 1) if prof['codex_I'] else None,
                m1_I=prof['m1_I'], rel_I_m1=I_quad / prof['m1_I'] - 1,
                rel_born_codex=4 * events / Q / prof['codex_born'] - 1, rel_born_m1=4 * events / Q / prof['m1_born'] - 1)


def born_spline_trapz(prof, dr=0.025, rmax=100.0, nq=4001):
    """(B) a cruder, deliberately different variant: cubic spline, trapezoid in r and in theta (q = qmax sin theta, trapezoid)."""
    a = np.loadtxt(prof['csv'], delimiter=',', skiprows=1)
    spl = CubicSpline(a[:, 0], a[:, 1])
    r = np.arange(0, rmax + dr / 2, dr); F = spl(r); F[r > a[-1, 0]] = 0
    W = 4 * prof['Omega']; qmax = np.sqrt(W * W - 4)
    th = np.linspace(0, np.pi / 2, nq); q = qmax * np.sin(th)
    base = r * r * F ** 4
    ft = np.array([4 * np.pi * np.trapezoid(base * np.sinc(qq * r / np.pi), r) for qq in q])
    beta = np.sqrt(np.maximum(1 - 4 / (W * W - q * q), 0))
    I = np.trapezoid(q * q * ft * ft * beta * qmax * np.cos(th), th)
    return dict(I=I, rel_I_codex=(I / prof['codex_I'] - 1) if prof['codex_I'] else None, rel_I_m1=I / prof['m1_I'] - 1, dr=dr, rmax=rmax, nq=nq)


res = {}
rg = np.linspace(0, 60, 24001); L = 3.0
ftg, _ = make_transform(rg, np.exp(-rg * rg / (2 * L * L)))
qg = np.array([0.0, 0.5, 1.0, 2.0, 3.0])
exact = (2 * np.pi) ** 1.5 * L ** 3 * np.exp(-L * L * qg * qg / 2)
res['gaussian_control'] = dict(L=L, q=qg.tolist(), max_abs_err_over_peak=float(np.max(np.abs(np.array([ftg(q) for q in qg]) - exact)) / exact[0]))
print('gaussian control', res['gaussian_control'])
for k, p in PROF.items():
    nat = born_native(p); var = born_spline_trapz(p)
    res[k] = dict(native_quad=nat, spline_trapezoid=var)
    print(k, 'I_native', nat['I'], 'rel codex', nat['rel_I_codex'], 'rel m1', nat['rel_I_m1'], 'parseval', nat['parseval_rel'],
          'born', nat['born_fractional_core_rate_c1'], 'rel codex', nat['rel_born_codex'], '| trapz I', var['I'], var['rel_I_m1'])
res['profile_factor_native'] = res['M1']['native_quad']['born_fractional_core_rate_c1'] / res['old']['native_quad']['born_fractional_core_rate_c1']
res['m1_ledger_profile_factor'] = 0.910732250186683
# Gamma t0 chain with the native Born rates + Codex's coherent channel [Q] + Codex's A
coh = dict(old=3.91958726382564e-50, M1=6.099550997404576e-51)
A = dict(old=2.108859793121751e-31, M1=2.079605512546138e-31)
G = {k: (coh[k] + res[k]['native_quad']['born_fractional_core_rate_c1']) * (A[k] * MRED ** 2) ** 2 * AGE for k in ('old', 'M1')}
res['Gamma_t0'] = dict(old=G['old'], M1=G['M1'], codex_old=7.08367857848652e-17, codex_M1=6.271629359413795e-17,
                       rel_old=G['old'] / 7.08367857848652e-17 - 1, rel_M1=G['M1'] / 6.271629359413795e-17 - 1, ratio=G['M1'] / G['old'])
print(json.dumps(res['Gamma_t0'], indent=1))
OUT.mkdir(parents=True, exist_ok=True)
(OUT / 'skeptic_born.json').write_text(json.dumps(dict(method=__doc__, results=res), indent=1, default=float))
