"""Own diagnostics on an early3d trace (same definitions as Codex analyze.trace_quantities, re-implemented)."""
import numpy as np

LX, G2, LP = .001, .01, .03


def load(path):
    return np.genfromtxt(str(path), names=True, delimiter='\t')


def diag(t):
    rho = t['rho']; r2bar = t['q0'] ** 2 + t['q1'] ** 2
    EX0 = .5 * (t['p0'] ** 2 + t['p1'] ** 2) + LX / 4 * r2bar ** 2 + G2 / 2 * t['qphi'] ** 2 * r2bar
    EP0 = .5 * t['pphi'] ** 2 + LP / 4 * t['qphi'] ** 4
    Ephi = t['Kphi'] + t['Gphi'] + t['Uphi']
    contrast = np.sqrt(np.maximum(t['charge_RMS'] ** 2 - t['nQ'] ** 2, 0)) / abs(t['nQ'])
    fluct = np.maximum(t['X2'] - r2bar, 0) / t['X2']
    return dict(tau=t['tau'], eta=t['eta'], rho=rho, Ephi=Ephi / rho, phimix=(Ephi + t['Umix']) / rho,
                outside=1 - (EX0 + EP0) / rho, Q0Q=(t['q0'] * t['p1'] - t['q1'] * t['p0']) / t['nQ'],
                contrast=contrast, fluct=fluct, Emin=.75 * LX ** (1 / 3) * abs(t['nQ']) ** (4 / 3) / rho,
                r2bar=r2bar, X2=t['X2'], phi2=t['phi2'])


def first(tau, v, level):
    ix = np.flatnonzero(v >= level)
    return float(tau[ix[0]]) if len(ix) else None


def window(d, key, lo, hi):
    m = (d['tau'] >= lo) & (d['tau'] <= hi); v = d[key][m]; e = d['eta'][m]
    # conformal-time (eta) weighted mean
    return dict(min=float(v.min()), max=float(v.max()), eta_mean=float(np.trapezoid(v, e) / (e[-1] - e[0])))
