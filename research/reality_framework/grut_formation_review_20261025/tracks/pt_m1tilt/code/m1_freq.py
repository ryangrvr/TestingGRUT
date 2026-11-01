"""Track pt_m1tilt diagnostic: frequency locking between the phi oscillation and the X radial oscillation along M.
For one run (orbit, eps, phase) sample, at ~40 log-spaced M between 1.5 and 60 M_c,peri, a window of 24 radial
periods densely (4096 points) and take the dominant Fourier frequency of r^2(s) (f_r) and of chi(s) (f_phi), plus the
power fraction of chi at the harmonics of f_r (slaved part) and the X free-basis n_-/q and N_phi/q over the window.
Usage: python3 m1_freq.py <orbit> <eps> <phase>.  Writes $SP/review/pt_m1tilt/freq_<...>.json.  [C] classical model."""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m1_common import *

name, eps, ph = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
y4, info = orbit(name); st = phase_state(y4, info['T_r'], ph)
q = charge(st)
Mperi = np.sqrt(G2 / LP) * info['rmin']; Mapo = np.sqrt(G2 / LP) * info['rmax']
M0 = 0.5 * Mperi; kappa = eps * info['omega_r']
par = SPar(3, 0, LX, G2, LP, M0, kappa, 0.0, EPST_M1)
chi0, _ = phi_quasistatic(st[0] ** 2 + st[1] ** 2, M0, EPST_M1)
R = Runner(par, [st[0], st[1], chi0, st[2], st[3], 0.0], rtol=1e-11)
s_of = lambda M: np.log(M / M0) / kappa
rows = []
NP, NS = 24, 4096
for Mg in np.geomspace(1.5 * Mperi, 60 * Mperi, 40):
    R.to(s_of(Mg))
    r2 = R.y[0] ** 2 + R.y[1] ** 2
    EX = 0.5 * (R.y[3] ** 2 + R.y[4] ** 2) + 0.25 * LX * r2 * r2 + 0.5 * G2 * Mg * Mg * r2
    Tr = radial_quad(EX, q, LX, G2 * Mg * Mg)['T_r']     # frozen-orbit estimate, only to size the window
    s0 = R.s.value; ts = s0 + np.linspace(0, NP * Tr, NS + 1)
    Y = np.empty((NS + 1, 7)); Y[0] = R.y
    for i, t in enumerate(ts[1:], 1):
        Y[i] = R.to(t)
    dt = ts[1] - ts[0]
    rr = Y[:-1, 0] ** 2 + Y[:-1, 1] ** 2; cc = Y[:-1, 2]
    w = np.hanning(NS)
    Fr = np.abs(np.fft.rfft((rr - rr.mean()) * w)); Fc = np.abs(np.fft.rfft((cc - cc.mean()) * w))
    fr = np.fft.rfftfreq(NS, dt)
    ir = np.argmax(Fr[1:]) + 1; ic = np.argmax(Fc[1:]) + 1
    # refine peaks by parabolic interpolation in log power
    def refine(F, i):
        if 1 <= i < len(F) - 1:
            a, b, c = np.log(F[i - 1] + 1e-300), np.log(F[i] + 1e-300), np.log(F[i + 1] + 1e-300)
            d = 0.5 * (a - c) / (a - 2 * b + c) if (a - 2 * b + c) != 0 else 0.0
            return fr[i] + d * (fr[1] - fr[0])
        return fr[i]
    f_r, f_c = refine(Fr, ir), refine(Fc, ic)
    # power of chi near harmonics n f_r (within +-1.5 bins) vs total
    P = Fc ** 2; tot = P[1:].sum(); harm = 0.0
    for n in range(1, 12):
        k = int(round(n * f_r / (fr[1] - fr[0])))
        lo, hi = max(1, k - 2), min(len(P), k + 3); harm += P[lo:hi].sum()
    nM, nF, Nphi, M = measures(par, ts, Y, q)
    av = lambda v: float(np.trapezoid(v, ts) / (ts[-1] - ts[0]))
    rows.append(dict(M_over_Mperi=float(Mg / Mperi), f_r=float(f_r), f_phi=float(f_c), ratio_phi_over_r=float(f_c / f_r),
                     chi_harmonic_power_fraction=float(harm / tot), nminus=av(nM), Nphi=av(Nphi),
                     dchi_over_M_rms=float(np.sqrt(np.mean((cc / Mg - 1) ** 2))), m_phi_lin_over_2pi=float(np.sqrt(2 * LP) * Mg / (2 * np.pi)),
                     f_r_frozen=float(1 / Tr)))
    r = rows[-1]
    print(f"M/Mperi {r['M_over_Mperi']:7.2f}  f_r {f_r:8.3f} (frozen {1/Tr:8.3f})  f_phi {f_c:8.3f} (lin m_phi/2pi {r['m_phi_lin_over_2pi']:8.3f})  "
          f"f_phi/f_r {r['ratio_phi_over_r']:.4f}  harm.frac {r['chi_harmonic_power_fraction']:.3f}  rms dchi/M {r['dchi_over_M_rms']:.3f}  n-/q {r['nminus']:.4f}  Nphi/q {r['Nphi']:.4f}")
out = os.path.join(SCR, f'freq_{name}_eps{eps:.0e}_ph{ph:g}.json')
json.dump(dict(orbit=name, eps=eps, phase=ph, J0_per_q=info['J_per_q'], rows=rows), open(out, 'w'))
print('wrote', out)
