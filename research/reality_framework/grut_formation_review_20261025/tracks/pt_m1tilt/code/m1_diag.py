"""Track pt_m1tilt diagnostic: where along M does the loss happen under the M1 potential?
For one run (orbit, eps, phase) record on a log grid in M (1 .. 40 M_c,apo): the free-basis n_-/q and N_phi/q averaged
over one X period at that M, the frozen-orbit radial frequency omega_r(M) of the X field (quartic + mass g M, at the
current q and E_X) and the phi frequency m_phi(M) = sqrt(2 lamP M^2 + g^2 <r^2>), so that the 1:1 crossing
omega_r = m_phi can be located.  Usage: python3 m1_diag.py <orbit> <eps> <phase> [epst] [mode]
Writes $SP/review/pt_m1tilt/diag_<orbit>_eps<eps>_ph<phase>_et<epst>.json.  [C] homogeneous classical model."""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m1_common import *

name, eps, ph = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
epst = float(sys.argv[4]) if len(sys.argv) > 4 else EPST_M1
mode = int(sys.argv[5]) if len(sys.argv) > 5 else 3
y4, info = orbit(name); st = phase_state(y4, info['T_r'], ph)
q = charge(st)
Mperi = np.sqrt(G2 / LP) * info['rmin']; Mapo = np.sqrt(G2 / LP) * info['rmax']
M0 = 0.5 * Mperi; kappa = eps * info['omega_r']
par = SPar(mode, 0, LX, G2, LP, M0, kappa, 0.0, epst)
chi0, _ = phi_quasistatic(st[0] ** 2 + st[1] ** 2, M0, epst)
R = Runner(par, [st[0], st[1], chi0, st[2], st[3], 0.0], rtol=1e-11)
s_of = lambda M: np.log(M / M0) / kappa
rows = []
for Mg in np.geomspace(1.0 * Mperi, 40 * Mapo, 120):
    R.to(s_of(Mg))
    # one-X-period window (period of the massive X: 2 pi/(g M) for the harmonic part; use 2 radial periods of the frozen orbit)
    r2 = R.y[0] ** 2 + R.y[1] ** 2
    EX = 0.5 * (R.y[3] ** 2 + R.y[4] ** 2) + 0.25 * LX * r2 * r2 + 0.5 * G2 * Mg * Mg * r2
    try:
        o = radial_quad(EX, q, LX, G2 * Mg * Mg); Tr = o['T_r']; r2m = o['r2mean']; Jr = o['J_r'] / q
    except Exception:
        Tr = 2 * np.pi / (GG * Mg); r2m = r2; Jr = float('nan')
    s0 = R.s.value; ts = s0 + np.linspace(0, 2 * Tr, 201)
    Y = np.empty((201, 7)); Y[0] = R.y
    for i, t in enumerate(ts[1:], 1):
        Y[i] = R.to(t)
    nM, nF, Nphi, M = measures(par, ts, Y, q)
    av = lambda v: float(np.trapezoid(v, ts) / (ts[-1] - ts[0]))
    mphi = np.sqrt(2 * LP * Mg * Mg + G2 * r2m)
    rows.append(dict(M_over_Mperi=float(Mg / Mperi), nminus=av(nM), Nphi=av(Nphi), Jr_frozen=float(Jr), omega_r=2 * np.pi / Tr,
                     m_phi=float(mphi), ratio=float(2 * np.pi / Tr / mphi), chi_over_M_mean=float(np.mean(Y[:, 2]) / Mg),
                     dchi_over_M_max=float(np.max(np.abs(Y[:, 2] - Mg)) / Mg),
                     Ephi_over_q=float(av(0.5 * (Y[:, 5] - kappa * M) ** 2 + 0.25 * LP * (Y[:, 2] ** 2 - M ** 2) ** 2 - epst * M * (Y[:, 2] - M) ** 3) / q)))
    print(f"M/Mperi {Mg/Mperi:8.3f}  n-/q {rows[-1]['nminus']:.4f}  Jr_frozen {Jr:.4f}  Nphi/q {rows[-1]['Nphi']:.4f}  omega_r {2*np.pi/Tr:.3f}  m_phi {mphi:.3f}  ratio {rows[-1]['ratio']:.3f}  "
          f"<chi>/M {rows[-1]['chi_over_M_mean']:.3f}  max|dchi|/M {rows[-1]['dchi_over_M_max']:.3f}  Ephi/q {rows[-1]['Ephi_over_q']:.3f}")
out = os.path.join(SCR, f'diag_{name}_eps{eps:.0e}_ph{ph:g}_et{epst:g}_m{mode}.json')
json.dump(dict(orbit=name, eps=eps, phase=ph, epst=epst, mode=mode, J0_per_q=info['J_per_q'], rows=rows), open(out, 'w'))
print('wrote', out)
