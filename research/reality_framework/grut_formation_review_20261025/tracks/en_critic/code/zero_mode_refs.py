"""Own code: lab-frame X-bar content of the coherent zero mode (windowed spectrum of the spatial mean) and the two
Codex reference projections at the tau-800 snapshot, for the X-seed-amplitude controls.

Zero mode: psi(s) = mean chi1 + i mean chi2, resampled uniformly in conformal time s, Hann window. With
psi = sum_j c_j exp(i nu_j s), the line charge is nu_j |c_j|^2, so n0_minus = int_{nu<0} |nu| S(nu) dnu/2pi when the
net charge is positive. The window-averaged zero-mode charge is reconstructed as a check.

usage: python3 -I zero_mode_refs.py <out_json> <label=run_dir> ...
"""
import sys, json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline

LX, G2 = .001, .01


def zero_mode(trace, lo, hi, ds=0.002):
    t = trace; m = (t['tau'] >= lo - 1) & (t['tau'] <= hi + 1)
    s = t['eta'][m]; keep = np.r_[True, np.diff(s) > 0]; s = s[keep]
    q0 = t['q0'][m][keep]; q1 = t['q1'][m][keep]; Qtot = float(np.mean(t['nQ'][m]))
    s_lo = np.interp(lo, t['tau'], t['eta']); s_hi = np.interp(hi, t['tau'], t['eta'])
    grid = np.arange(s_lo, s_hi, ds)
    psi = CubicSpline(s, q0)(grid) + 1j * CubicSpline(s, q1)(grid)
    dpsi = CubicSpline(s, q0)(grid, 1) + 1j * CubicSpline(s, q1)(grid, 1)
    w = np.hanning(len(grid)); norm = np.sum(w * w) * ds
    F = np.fft.fft(w * psi) * ds; nu = 2 * np.pi * np.fft.fftfreq(len(grid), ds)
    S = abs(F) ** 2 / norm; dnu = 2 * np.pi / (len(grid) * ds)
    n_plus = float(np.sum((nu > 0) * nu * S) * dnu / (2 * np.pi))
    n_minus = float(np.sum((nu < 0) * -nu * S) * dnu / (2 * np.pi))
    Q0_direct = float(np.sum(w * w * (psi.conj() * dpsi).imag) / np.sum(w * w))
    sign = np.sign(Qtot)
    nm = n_minus if sign > 0 else n_plus
    return dict(window_tau=[lo, hi], window_eta=float(s_hi - s_lo), Q_total=Qtot,
                zero_mode_nminus_over_Q=nm / abs(Qtot), zero_mode_charge_from_lines_over_Q=(n_plus - n_minus) / Qtot,
                zero_mode_charge_direct_over_Q=Q0_direct / Qtot)


def symbol(N, L):
    a = 2 * np.pi * np.fft.fftfreq(N); dx = L / N
    f = lambda x: 4 * np.sin(x / 2) ** 2 * (1 + np.sin(x / 2) ** 2 / 3) / dx ** 2
    return f(a)[:, None, None] + f(a)[None, :, None] + f(a)[None, None, :]


def references(path):
    with open(path, 'rb') as fh:
        N = int(np.frombuffer(fh.read(4), '<i4')[0]); L = float(np.frombuffer(fh.read(8), '<f8')[0])
        s = float(np.frombuffer(fh.read(8), '<f8')[0]); fh.read(4)
        a = np.fromfile(fh, '<f8').reshape(6, N, N, N)
    q, p = a[:3], a[3:]
    qh = np.fft.fftn(q[:2], axes=(-3, -2, -1), norm='forward'); ph = np.fft.fftn(p[:2], axes=(-3, -2, -1), norm='forward')
    k2 = symbol(N, L); nz = k2 > 0
    Qm = (np.conj(qh[0]) * ph[1] - np.conj(qh[1]) * ph[0]).real
    Qtot = float(np.mean(q[0] * p[1] - q[1] * p[0]))
    r2 = q[0] ** 2 + q[1] ** 2; hart = 2 * LX * np.mean(r2) + G2 * np.mean(q[2] ** 2)
    out = dict(eta=s, Q_total=Qtot, Q0_over_Q=float((q[0].mean() * p[1].mean() - q[1].mean() * p[0].mean()) / Qtot),
               Hartree_mass2=float(hart))
    for name, m2 in [('massless', 0.), ('Hartree', hart), ('m2_12', 12.)]:
        w = np.sqrt(np.where(nz, k2 + m2, 1.))
        nr = (np.sum(abs(ph) ** 2, 0) + w * w * np.sum(abs(qh) ** 2, 0)) / (2 * w)
        out[name + '_nminus_over_Q'] = float(np.sum(((nr - Qm) / 2)[nz]) / Qtot)
    return out


def main():
    out = Path(sys.argv[1]); res = {}
    for arg in sys.argv[2:]:
        label, d = arg.split('=', 1); d = Path(d)
        t = np.genfromtxt(str(d / 'state.tsv'), names=True, delimiter='\t')
        row = dict(zero_mode=[zero_mode(t, lo, hi) for lo, hi in [(340, 798), (547, 800), (700, 800)]])
        snaps = sorted(d.glob('state_s*.bin'))
        if snaps:
            last = snaps[-1]; row['snapshot'] = last.name; row['refs_tau800'] = references(last)
        res[label] = row
        z = row['zero_mode']
        print(label, ['%.4f' % x['zero_mode_nminus_over_Q'] for x in z],
              'charge check', ['%.4f/%.4f' % (x['zero_mode_charge_from_lines_over_Q'], x['zero_mode_charge_direct_over_Q']) for x in z],
              row.get('refs_tau800'))
    out.write_text(json.dumps(res, indent=2) + '\n')


if __name__ == '__main__':
    main()
