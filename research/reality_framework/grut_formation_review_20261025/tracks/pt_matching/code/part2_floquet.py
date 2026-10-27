"""Part 2a: orbit family and Floquet stability of the homogeneous phi = 0 state on a frozen quartic X orbit.
phi'' = -(g^2 r(s)^2 - lambda_phi M^2) phi with r(s) the pure-quartic orbit (phi = 0 exactly does not
disturb X).  M = v a is frozen.  Independent code."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from pm_common import *
from scipy.optimize import brentq


def peri_state(q, r0):
    """Pure-quartic orbit with charge q and pericentre (or apocentre) r0: tangential velocity only."""
    return np.array([r0, 0.0, 0.0, q / r0])


def orbit_info(y4, lam=LX, nsteps=4096):
    q = charge(y4)
    E = 0.5 * (y4[2] ** 2 + y4[3] ** 2) + 0.25 * lam * (y4[0] ** 2 + y4[1] ** 2) ** 2
    orb = orbit_quadrature(E, q, lam)
    Tr, Dth = orb['T_r'], orb['dtheta']
    par = make_par(mode=0, lamramp=0, T=1e-9, mout=0.0, lamX0=lam, g2=0.0, lamP=0.0)
    y = y8(y4[0], y4[1], 0.0, y4[2], y4[3], 0.0, s=1.0)
    out = c_fixed(par, y, Tr / nsteps, nsteps, 8)
    psi = np.r_[y4[0] + 1j * y4[1], out[:-1, 0] + 1j * out[:-1, 1]]
    t = np.arange(nsteps) * Tr / nsteps
    Om, wr = Dth / Tr, 2 * np.pi / Tr
    c = np.fft.fft(psi * np.exp(-1j * Om * t)) / nsteps
    n = np.fft.fftfreq(nsteps, 1.0 / nsteps).round().astype(int)
    w = Om + n * wr; P = np.abs(c) ** 2
    r2 = np.abs(psi) ** 2
    return dict(q=q, E=E, J_per_q=orb['J_r'] / q, A_minus_per_q=-np.sum((w * P)[w < 0]) / q, T_r=Tr, omega_r=wr, Omega=Om,
                rmin=orb['rmin'], rmax=orb['rmax'], ratio=orb['rmax'] / orb['rmin'], r2_mean=r2.mean(), r2_min=r2.min(), r2_max=r2.max())


def floquet(y4, M, Tr, nsteps=4096, eps0=1e-100):
    tr = []
    mono = np.zeros((2, 2))
    for j, (f0, p0) in enumerate([(eps0, 0.0), (0.0, eps0)]):
        par = make_par(mode=1, M0=M, kappa=0.0)
        y = y8(y4[0], y4[1], f0, y4[2], y4[3], p0, s=0.0)
        c_fixed(par, y, Tr / nsteps, nsteps, 8, store=False)
        mono[:, j] = [y[2] / eps0, y[5] / eps0]
    t = np.trace(mono)
    nu = np.arccosh(abs(t) / 2) / Tr if abs(t) > 2 else 0.0
    return t, nu, np.linalg.det(mono)


def scan(y4, info, Mmax_factor=1.3, n=3000):
    Mperi = np.sqrt(G2 / LP) * info['rmin']; Mapo = np.sqrt(G2 / LP) * info['rmax']; Mmean = np.sqrt(G2 / LP * info['r2_mean'])
    Ms = np.linspace(0.0, Mmax_factor * Mapo, n)
    rows = []
    for M in Ms:
        t, nu, det = floquet(y4, M, info['T_r'])
        rows.append((M, t, nu, det))
    rows = np.array(rows)
    unstable = rows[:, 2] > 0
    # band edges
    edges = []
    for i in range(1, n):
        if unstable[i] != unstable[i - 1]:
            edges.append(dict(M=0.5 * (rows[i, 0] + rows[i - 1, 0]), into_unstable=bool(unstable[i])))
    first = rows[unstable][0, 0] if unstable.any() else None
    return dict(M_c_peri=Mperi, M_c_apo=Mapo, M_c_meanr2=Mmean, M_first_unstable=first, band_edges=edges[:60],
                max_det_err=float(np.max(np.abs(rows[:, 3] - 1))),
                table=[dict(M=float(r[0]), M_over_Mperi=float(r[0] / Mperi), trace=float(r[1]), nu=float(r[2]), nu_over_omega_r=float(r[2] / info['omega_r'])) for r in rows[::10]],
                ), rows


def main():
    y4 = tau20_state(); q = charge(y4)
    tau20 = orbit_info(y4)
    print('tau20', tau20)
    # late orbit: pure-quartic orbit with the same q and (a) A_-/q = coherent share per coherent charge in s 24-34,
    # (b) the sampled r_max/r_min of the coherent trace in s 24-34 (base48)
    A_target = 0.4760441973431541 / 1.0090661119866706   # Codex WINDOW_base48_24_34: coherent A_- per nQ, divided by mean Q0/Q
    rc = (q * q / LX) ** (1 / 6)
    fA = lambda r0: orbit_info(peri_state(q, r0), nsteps=2048)['A_minus_per_q'] - A_target
    r0A = brentq(fA, 10.0, 0.999 * rc, xtol=1e-10)
    lateA = orbit_info(peri_state(q, r0A))
    fR = lambda r0: orbit_info(peri_state(q, r0), nsteps=1024)['ratio'] - 3.9646286202357226
    r0R = brentq(fR, 10.0, 0.999 * rc, xtol=1e-10)
    lateR = orbit_info(peri_state(q, r0R))
    print('late (A match)', lateA); print('late (ratio match)', lateR)
    res = dict(tau20=tau20, late_A_match=dict(A_target=A_target, r0=r0A, **lateA), late_ratio_match=dict(ratio_target=3.9646286202357226, r0=r0R, **lateR))
    for name, yy, info in [('tau20', y4, tau20), ('late_A_match', peri_state(q, r0A), lateA)]:
        sc, rows = scan(yy, info)
        res['floquet_' + name] = sc
        np.save(os.path.join(SCRATCH, f'floquet_{name}.npy'), rows)
        print(name, 'M_peri', sc['M_c_peri'], 'M_apo', sc['M_c_apo'], 'M_meanr2', sc['M_c_meanr2'], 'first unstable', sc['M_first_unstable'],
              'det err', sc['max_det_err'])
        for e in sc['band_edges'][:20]: print('   edge', e['M'], e['M'] / sc['M_c_peri'], e['into_unstable'])
    save_json('floquet.json', res)


if __name__ == '__main__':
    main()
