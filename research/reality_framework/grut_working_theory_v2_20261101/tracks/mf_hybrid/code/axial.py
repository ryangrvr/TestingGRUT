"""Non-perturbative fermion production from an axial (mass-phase) coupling  -(M/2) Nbar exp(i gamma5 phi/f) N  [C].

Rotating the phase away gives an axial chemical potential mu5 = phidot/(2f). Per helicity h and momentum k the
two-level Hamiltonian is H = [[p, M], [M, -p]] with p = h k - mu5(t). For phi = Sigma cos(m t), one half-period
(turning point to turning point) is integrated exactly for a grid of k; n_k is the occupation in the free basis at the
next turning point (mu5 = 0 there, so the particle number is unambiguous). Units: m = 1.
  A = Sigma/(2 f) = mu5_max/m,  Mh = M/m,  S = Sigma/m.
Energy produced per half-period (one Majorana species; only one helicity crosses per half-period):
  eps = (1/2 pi^2) int k^2 n_k sqrt(k^2 + M^2) dk,  fraction F = eps / (S^2/2).
Landau-Zener check: two crossings with P = exp(-pi M^2/|mu5dot|), |mu5dot| = sqrt(A^2 - k^2) (m = 1) at the crossing,
averaged Stueckelberg phase: n_LZ = 2 P (1 - P).
The produced number does not depend on S at fixed A, so F = eps(A, Mh)/(S^2/2) with S = 2 A f/m.
Also: perturbative pseudoscalar width Gamma = (M/f)^2 m beta/(16 pi), beta = sqrt(1 - 4M^2/m^2) (Majorana).
Usage: python3 -B axial.py  -> ../results/axial_check.json
"""
import math, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def half_period(A, Mh, nk=1500, steps_per_unit=None):
    k = np.linspace(1e-3, 1.25 * A + 3 * Mh + 2, nk)
    h = -1.0                       # the helicity that crosses when mu5 < 0 (phi = S cos t, mu5 = -A sin t)
    wmax = k[-1] + A + Mh
    nsteps = int(max(2000, 40 * wmax * math.pi))
    dt = math.pi / nsteps

    def Hm(t):
        p = h * k + A * math.sin(t)
        return p

    p0 = h * k
    E0 = np.sqrt(p0 ** 2 + Mh ** 2)
    # negative-energy eigenvector of [[p, M],[M,-p]]: (M, -(E+p)) normalised
    a = Mh * np.ones_like(k)
    b = -(E0 + p0)
    nrm = np.sqrt(a * a + b * b)
    psi = np.stack([a / nrm, b / nrm]).astype(complex)

    def f(t, y):
        p = Hm(t)
        return -1j * np.stack([p * y[0] + Mh * y[1], Mh * y[0] - p * y[1]])

    t = 0.0
    for _ in range(nsteps):        # RK4
        k1 = f(t, psi)
        k2 = f(t + dt / 2, psi + dt / 2 * k1)
        k3 = f(t + dt / 2, psi + dt / 2 * k2)
        k4 = f(t + dt, psi + dt * k3)
        psi = psi + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
    pe = h * k
    E = np.sqrt(pe ** 2 + Mh ** 2)
    # positive-energy eigenvector: (M, E - p)
    a, b = Mh * np.ones_like(k), E - pe
    nrm = np.sqrt(a * a + b * b)
    n = np.abs((a * psi[0] + b * psi[1]) / nrm) ** 2
    eps = np.trapezoid(k ** 2 * n * np.sqrt(k ** 2 + Mh ** 2), k) / (2 * math.pi ** 2)
    num = np.trapezoid(k ** 2 * n, k) / (2 * math.pi ** 2)
    # LZ estimate
    kk = k[k < A]
    rate = np.sqrt(np.maximum(A * A - kk * kk, 1e-30))
    P = np.exp(-math.pi * Mh ** 2 / rate)
    nlz = 2 * P * (1 - P)
    eps_lz = np.trapezoid(kk ** 2 * nlz * np.sqrt(kk ** 2 + Mh ** 2), kk) / (2 * math.pi ** 2)
    return dict(A=A, Mh=Mh, eps=float(eps), num=float(num), eps_LZ=float(eps_lz), nmax=float(n.max()),
                norm_err=float(np.max(np.abs(np.sum(np.abs(psi) ** 2, axis=0) - 1))))


def frac(A, Mh, f_over_m, **kw):
    r = half_period(A, Mh, **kw)
    S = 2 * A * f_over_m
    r['S'] = S
    r['F'] = r['eps'] / (0.5 * S * S)
    r['F_LZ'] = r['eps_LZ'] / (0.5 * S * S)
    return r


def table(Mh, f_over_m, As):
    return [frac(A, Mh, f_over_m) for A in As]


def gamma_pert(M, f, m):
    if 2 * M >= m:
        return 0.0
    return (M / f) ** 2 * m * math.sqrt(1 - 4 * M * M / (m * m)) / (16 * math.pi)


if __name__ == '__main__':
    out = dict(label='C', note=__doc__)
    rows = []
    for (A, Mh) in [(5, 0.3), (5, 1.0), (10, 0.42), (10, 1.0), (20, 0.42), (20, 2.0), (40, 1.0), (40, 2.0)]:
        r = frac(A, Mh, 1.0)
        rows.append(r)
        print('A %5.1f M %.2f  eps %.4g  eps_LZ %.4g  ratio %.3f  nmax %.3f  norm %.1e' % (
            A, Mh, r['eps'], r['eps_LZ'], r['eps'] / max(r['eps_LZ'], 1e-300), r['nmax'], r['norm_err']))
    out['checks'] = rows
    json.dump(out, open(os.path.join(HERE, '..', 'results', 'axial_check.json'), 'w'), indent=1)
