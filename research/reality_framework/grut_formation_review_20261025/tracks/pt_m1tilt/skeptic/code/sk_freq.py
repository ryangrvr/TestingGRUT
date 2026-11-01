"""Mechanism diagnostic (skeptic): along M, in windows of NP radial periods (radial period measured from the
pericentre passages inside the window, not the frozen-orbit estimate), the spectra of chi(s) and r^2(s):
dominant frequencies, the chi power at integer harmonics of f_r, at HALF-integer multiples (k/2 f_r, k odd: the
signature of a 2 omega_phi = 3 omega_r parametric lock), and elsewhere; plus the free-basis n_-, N_phi, the
combination 2 n_- + 3 N_phi, the rms swing and the linear m_phi/omega_r (measured).  Also a continuous record of
the exchange: n_- and N_phi over the window.
Usage: python3 -B sk_freq.py <orbit> <eps> <phase> [epst] [Mlo Mhi n]"""
import sys, os, json, numpy as np
from sk_common import *

name, eps, ph = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
epst = float(sys.argv[4]) if len(sys.argv) > 4 else EPST
Mlo, Mhi, nM = (float(sys.argv[5]), float(sys.argv[6]), int(sys.argv[7])) if len(sys.argv) > 7 else (1.5, 40.0, 48)
y4, info = orbit(name); st = phase_state(y4, info['T_r'], ph); q = charge(st)
par = make_par(info, eps, epst=epst)
chi0 = chi_start(st[0] ** 2 + st[1] ** 2, par.M0, epst)
S = Sim(par, np.r_[st[:2], chi0], np.r_[st[2:], 0.0], nper=32)
NP, NS = 32, 8192
rows = []
for Mg in np.geomspace(Mlo * info['Mperi'], Mhi * info['Mperi'], nM):
    S.to(s_of_M(par, Mg))
    Tr_est = 2 * np.pi / np.sqrt(info['omega_r'] ** 2 + (2 * GG * Mg) ** 2)
    # first pass: 2 estimated periods to measure the actual radial period from pericentre passages
    y0 = S.y.copy(); Y = np.vstack([y0, S.sample(4 * Tr_est / 800, 800)])
    r = np.hypot(Y[:, 0], Y[:, 1]); rd = Y[:, 0] * Y[:, 3] + Y[:, 1] * Y[:, 4]
    idx = np.where((rd[:-1] < 0) & (rd[1:] >= 0))[0]
    Tr = float(np.mean(np.diff(Y[idx, 6]))) if len(idx) >= 3 else 2 * Tr_est   # crude; refined by the FFT below
    ds = NP * Tr / NS
    y0 = S.y.copy(); Y = np.vstack([y0, S.sample(ds, NS)])[:-1]
    ts = Y[:, 6]; rr = Y[:, 0] ** 2 + Y[:, 1] ** 2; cc = Y[:, 2]
    w = np.hanning(NS)
    Fr = np.abs(np.fft.rfft((rr - rr.mean()) * w)) ** 2; Fc = np.abs(np.fft.rfft((cc - cc.mean()) * w)) ** 2
    fr = np.fft.rfftfreq(NS, ds); df = fr[1] - fr[0]
    def refine(F, i):
        if 1 <= i < len(F) - 1:
            a, b, c = np.log(F[i - 1] + 1e-300), np.log(F[i] + 1e-300), np.log(F[i + 1] + 1e-300)
            den = a - 2 * b + c
            return fr[i] + (0.5 * (a - c) / den if den != 0 else 0.0) * df
        return fr[i]
    ir = np.argmax(Fr[1:]) + 1; f_r = refine(Fr, ir)
    ic = np.argmax(Fc[1:]) + 1; f_c = refine(Fc, ic)
    # second and third chi peaks (separated by > 3 bins from the first)
    P = Fc.copy(); P[max(0, ic - 3):ic + 4] = 0; ic2 = np.argmax(P[1:]) + 1; f_c2 = refine(Fc, ic2); P2 = Fc[ic2]
    P[max(0, ic2 - 3):ic2 + 4] = 0; ic3 = np.argmax(P[1:]) + 1; f_c3 = refine(Fc, ic3)
    tot = Fc[1:].sum()
    def band_power(mult):
        k = int(round(mult * f_r / df)); lo, hi = max(1, k - 2), min(len(Fc), k + 3)
        return Fc[lo:hi].sum()
    harm = sum(band_power(n) for n in range(1, 13))
    half = sum(band_power(n + 0.5) for n in range(0, 12))
    m = free_measures(Y, ts, par, q, epst)
    nm, Nf = tav(m['nminus'], ts), tav(m['Nphi'], ts)
    mphi_lin = np.sqrt(2 * LP) * Mg
    rows.append(dict(M_over_Mperi=float(Mg / info['Mperi']), f_r=float(f_r), f_phi=float(f_c), f_phi2=float(f_c2), f_phi3=float(f_c3),
                     P2_over_P1=float(P2 / Fc[ic]), ratio=float(f_c / f_r), ratio2=float(f_c2 / f_r), harm_frac=float(harm / tot),
                     half_frac=float(half / tot), mphi_lin_over_wr=float(mphi_lin / (2 * np.pi * f_r)),
                     nminus=nm, Nphi=Nf, combo=2 * nm + 3 * Nf, dchi_rms=float(np.sqrt(np.mean((cc / Mg - 1) ** 2))),
                     chi_min=float(cc.min() / Mg), chi_max=float(cc.max() / Mg), Tr_meas=float(1 / f_r)))
    r_ = rows[-1]
    print(f"M/Mp {r_['M_over_Mperi']:7.2f} f_r {f_r:7.4f} f_phi {f_c:7.4f} ({r_['ratio']:.4f} f_r) 2nd {f_c2:7.4f} ({r_['ratio2']:.4f} f_r, P2/P1 {r_['P2_over_P1']:.3f}) 3rd {f_c3:7.4f} ({f_c3/f_r:.3f}) | harm {r_['harm_frac']:.3f} half {r_['half_frac']:.3f} | mphi_lin/wr {r_['mphi_lin_over_wr']:.3f} | n- {nm:.4f} Nphi {Nf:.4f} 2n+3N {r_['combo']:.4f} | rms {r_['dchi_rms']:.3f} [{r_['chi_min']:.2f},{r_['chi_max']:.2f}]", flush=True)
out = os.path.join(SCR, f'freq_{name}_eps{eps:.0e}_ph{ph:g}_et{epst:g}.json')
json.dump(dict(orbit=name, eps=eps, phase=ph, epst=epst, J0_per_q=info['J_per_q'], NP=NP, NS=NS, rows=rows), open(out, 'w'))
print('wrote', out)
