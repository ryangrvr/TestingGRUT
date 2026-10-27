"""Part 1a: starting homogeneous quartic orbit -- charge, energy, radial action (quadrature), exact
quasi-periodic Fourier decomposition, lab negative-frequency statistic, and the spectral (Percival) action
identities  q = sum_n w_n |c_n|^2,  J_r = sum_n n w_n |c_n|^2,  w_n = Omega + n w_r.
Independent code (own C symplectic integrator + mpmath quadrature)."""
import sys, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
from pm_common import *


def fourier_orbit(y4, lam, m2, nsteps, order=8):
    """Integrate exactly one radial period with fixed steps; return Fourier data of psi e^{-i Omega t}."""
    q = charge(y4)
    E = 0.5 * (y4[2] ** 2 + y4[3] ** 2) + 0.5 * m2 * (y4[0] ** 2 + y4[1] ** 2) + 0.25 * lam * (y4[0] ** 2 + y4[1] ** 2) ** 2
    orb = orbit_quadrature(E, q, lam, m2)
    Tr, Dth = orb['T_r'], orb['dtheta']
    par = make_par(mode=0, lamramp=0, T=1e-9, mout=np.sqrt(m2), lamX0=lam, g2=0.0, lamP=0.0)
    y = y8(y4[0], y4[1], 0.0, y4[2], y4[3], 0.0, s=1.0)
    h = Tr / nsteps
    out = c_fixed(par, y, h, nsteps, order)
    psi = np.empty(nsteps, complex)
    psi[0] = y4[0] + 1j * y4[1]
    psi[1:] = out[:-1, 0] + 1j * out[:-1, 1]
    psiT = out[-1, 0] + 1j * out[-1, 1]
    closure = abs(psiT - psi[0] * np.exp(1j * Dth)) / abs(psi[0])
    t = np.arange(nsteps) * h
    Om, wr = Dth / Tr, 2 * np.pi / Tr
    c = np.fft.fft(psi * np.exp(-1j * Om * t)) / nsteps
    n = np.fft.fftfreq(nsteps, 1.0 / nsteps).round().astype(int)
    w = Om + n * wr
    P = np.abs(c) ** 2
    res = dict(q=q, E=E, quad=orb, Omega=Om, omega_r=wr, closure=closure,
               A_minus=float(-np.sum((w * P)[w < 0])), A_plus=float(np.sum((w * P)[w > 0])),
               q_spectral=float(np.sum(w * P)), J_spectral=float(np.sum(n * w * P)),
               r2_mean=float(np.sum(P)))
    # decomposition of J_r - A_- into its harmonic terms
    pos = n >= 1; neg2 = n <= -2
    res['diff_pos_harmonics'] = float(np.sum((n * w * P)[pos]))
    res['diff_neg_harmonics'] = float(np.sum(((np.abs(n) - 1) * np.abs(w) * P)[neg2]))
    order_idx = np.argsort(-P)[:8]
    res['top_harmonics'] = [dict(n=int(n[i]), omega=float(w[i]), abs_c2_over_q=float(P[i] / q), omega_c2_over_q=float(w[i] * P[i] / q)) for i in order_idx]
    return res


def main():
    y4 = tau20_state()
    out = {}
    for ns in [1024, 2048, 4096, 8192]:
        r = fourier_orbit(y4, LX, 0.0, ns)
        out[f'quartic_N{ns}'] = r
        q = r['q']
        print(ns, 'q', q, 'E', r['E'], 'J/q quad', r['quad']['J_r'] / q, 'J/q spec', r['J_spectral'] / q,
              'A-/q', r['A_minus'] / q, 'A+/q', r['A_plus'] / q, 'qspec/q-1', r['q_spectral'] / q - 1, 'closure', r['closure'])
    r = out['quartic_N8192']; q = r['q']
    summary = dict(
        q=q, E=r['E'], rmin=r['quad']['rmin'], rmax=r['quad']['rmax'], T_r=r['quad']['T_r'], dtheta=r['quad']['dtheta'],
        Omega_over_omega_r=r['Omega'] / r['omega_r'],
        J_r_quadrature=r['quad']['J_r'], J_r_per_q=r['quad']['J_r'] / q, J_r_per_orig_net=r['quad']['J_r'] / Q_ORIG_NET,
        J_r_spectral_per_q=r['J_spectral'] / q, spectral_identity_rel=r['J_spectral'] / r['quad']['J_r'] - 1,
        A_minus_per_q=r['A_minus'] / q, A_minus_per_orig_net=r['A_minus'] / Q_ORIG_NET, A_plus_per_q=r['A_plus'] / q,
        charge_identity_rel=r['q_spectral'] / q - 1, closure=r['closure'],
        J_minus_A_per_q=(r['quad']['J_r'] - r['A_minus']) / q,
        J_minus_A_from_pos_harmonics_per_q=r['diff_pos_harmonics'] / q,
        J_minus_A_from_neg_harmonics_per_q=r['diff_neg_harmonics'] / q,
        top_harmonics=r['top_harmonics'],
        codex=dict(J_r_per_q=1.763725279396997, A_minus_per_q=1.7377883531120057, J_r_per_orig_net=1.7646827365883024,
                   A_minus_per_orig_net=1.7387317301645944, T_r=0.814465681774558, dtheta=2.995589266991105,
                   rmin=17.69707411368106, rmax=142.6824536721767, q=8117.881574752555, E=105233.42336942146))
    # harmonic and mixed-potential family at the same q and same r_min/r_max shape scale: show J_r - A_- vs anharmonicity
    fam = []
    rr = summary
    for frac in [0.0, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0]:
        # potential m2 r^2/2 + lam r^4/4 with lam = frac*LX and m2 chosen so the circular-orbit frequency scale matches
        lam = frac * LX
        m2 = (1 - frac) * LX * 0.5 * (rr['rmin'] ** 2 + rr['rmax'] ** 2)
        # same initial position and momentum as the tau-20 state (different energy); keeps q fixed
        rf = fourier_orbit(y4, lam, m2, 8192)
        fam.append(dict(quartic_fraction=frac, m2=m2, lam=lam, J_per_q=rf['quad']['J_r'] / rf['q'], A_minus_per_q=rf['A_minus'] / rf['q'],
                        J_spectral_per_q=rf['J_spectral'] / rf['q'], Omega_over_omega_r=rf['Omega'] / rf['omega_r'],
                        rmax_over_rmin=rf['quad']['rmax'] / rf['quad']['rmin']))
        print('family', fam[-1])
    summary['family'] = fam
    save_json('orbit_spectral.json', dict(summary=summary, convergence={k: dict(J_spectral_per_q=v['J_spectral'] / v['q'], A_minus_per_q=v['A_minus'] / v['q'], closure=v['closure'], charge_identity_rel=v['q_spectral'] / v['q'] - 1) for k, v in out.items()}))
    print({k: v for k, v in summary.items() if k not in ('family', 'top_harmonics')})


if __name__ == '__main__':
    main()
