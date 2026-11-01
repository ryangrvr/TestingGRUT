"""Skeptic check of m1_transport T13: is the core interior 'classically forbidden' to a slow X, and what
is the elastic X-core cross section at k = 0.0684 m_X in the linear three-channel problem?

Own code.  Uses Codex's stationary Q = 2e5 M1 profile (d4747cb PROFILE_M1.csv: rho = m_X r, F, S) and the
canonical radial operator of bf8ab3b DERIVATION.md,
    -y'' + [ M(r) + l(l+1)/r^2 - diag((nu+Omega)^2, (nu-Omega)^2, nu^2) ] y = 0,
    M = [[S^2+2aF^2, aF^2, sqrt2 F S], [aF^2, S^2+2aF^2, sqrt2 F S], [sqrt2 F S, sqrt2 F S, D]],
    D = b(3S^2-1) - 6c(S-1) + F^2,  a = 0.1, b = 3, c = 0.2,  y = r (u, v, sqrt2 w).
For an incoming X (particle) of momentum k: E = sqrt(1+k^2), nu = E - Omega, so the u channel is open
with momentum k, the v (Xbar) channel has energy E - 2 Omega < 1 (closed below the reverse threshold)
and the neutral channel has energy E - Omega < sqrt6 (closed).  Only elastic scattering is open, so
the S matrix is a phase shift per partial wave and sigma_el = 4 pi/k^2 sum (2l+1) sin^2 delta_l.

Method: integrate three regular solutions from the origin with RK4 and periodic QR re-orthonormalisation
(stable against the closed-channel exponential growth), impose decay of the closed channels at r_match,
read the open-channel phase shift from Riccati-Bessel matching.  Validated on (i) the vacuum (delta = 0),
(ii) an analytic s-wave square barrier, (iii) the hard-sphere limit.
"""
import json, math, os
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.special import spherical_jn, spherical_yn

RES = Path(__file__).resolve().parents[1] / 'results'
SP = Path(os.environ.get('SP', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad')) / 'review' / 'm1_transport_skeptic'
OMEGA = 0.88172049195
a_q, b_q, c_q = 0.1, 3.0, 0.2

prof = np.genfromtxt(SP / 'codex' / 'PROFILE_M1.csv', delimiter=',', names=True)
r_p, F_p, S_p = prof['rho'], prof['F'], prof['S']
F_i = PchipInterpolator(r_p, F_p); S_i = PchipInterpolator(r_p, S_p)


def Mmat(r, couple=True, Fscale=1.0):
    F = float(F_i(r)) * Fscale; S = float(S_i(r))
    d = S * S + 2 * a_q * F * F
    D = b_q * (3 * S * S - 1) - 6 * c_q * (S - 1) + F * F
    x = a_q * F * F if couple else 0.0
    y = math.sqrt(2) * F * S if couple else 0.0
    return np.array([[d, x, y], [x, d, y], [y, y, D]])


def channel_energies(k):
    E = math.sqrt(1 + k * k); nu = E - OMEGA
    return E, nu, np.array([(nu + OMEGA)**2, (nu - OMEGA)**2, nu**2])


def ricc(l, x):
    j = x * spherical_jn(l, x); y = x * spherical_yn(l, x)
    jp = spherical_jn(l, x) + x * spherical_jn(l, x, derivative=True)
    yp = spherical_yn(l, x) + x * spherical_yn(l, x, derivative=True)
    return j, y, jp, yp


def solve(k, l, Wfun, r0=0.02, rmax=40.0, hstep=0.004, renorm_every=25):
    """Generic 3-channel solver.  Wfun(r) -> 3x3 matrix W with y'' = W y.  Returns (psi, psi') at rmax."""
    n = 3
    psi = np.eye(n) * r0**(l + 1); dpsi = np.eye(n) * (l + 1) * r0**l
    r = r0; steps = int(round((rmax - r0) / hstep)); hstep = (rmax - r0) / steps
    def f(r, P, dP):
        return dP, Wfun(r) @ P
    for i in range(steps):
        k1p, k1d = f(r, psi, dpsi)
        k2p, k2d = f(r + hstep / 2, psi + hstep / 2 * k1p, dpsi + hstep / 2 * k1d)
        k3p, k3d = f(r + hstep / 2, psi + hstep / 2 * k2p, dpsi + hstep / 2 * k2d)
        k4p, k4d = f(r + hstep, psi + hstep * k3p, dpsi + hstep * k3d)
        psi = psi + hstep / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
        dpsi = dpsi + hstep / 6 * (k1d + 2 * k2d + 2 * k3d + k4d)
        r += hstep
        if (i + 1) % renorm_every == 0:
            Q, _ = np.linalg.qr(np.vstack([psi, dpsi]))
            psi, dpsi = Q[:n], Q[n:]
    return psi, dpsi, r


def phase_shift(k, l, Wfun, closed_kappa, open_idx=0, **kw):
    psi, dpsi, r = solve(k, l, Wfun, **kw)
    closed = [i for i in range(3) if i != open_idx]
    A = np.array([dpsi[i] + closed_kappa[i] * psi[i] for i in closed])      # 2 x 3
    _, _, vh = np.linalg.svd(A)
    cvec = vh[-1]
    yu = float(psi[open_idx] @ cvec); dyu = float(dpsi[open_idx] @ cvec)
    j, y, jp, yp = ricc(l, k * r)
    delta = math.atan2(k * yu * jp - dyu * j, k * yu * yp - dyu * y)
    # closed-channel leakage check: size of the closed components relative to the open one at rmax
    leak = max(abs(float(psi[i] @ cvec)) for i in closed) / max(abs(yu), 1e-300)
    return delta, leak


def xcore_W(k, l, couple=True, Fscale=1.0):
    _, _, w2 = channel_energies(k)
    def W(r):
        return Mmat(r, couple, Fscale) + (l * (l + 1) / (r * r)) * np.eye(3) - np.diag(w2)
    return W


def kappas(k):
    _, _, w2 = channel_energies(k)
    return np.array([0.0, math.sqrt(1 - w2[1]), math.sqrt(6 - w2[2])])   # closed channels decay rates in vacuum


def main():
    out = {}
    k0 = 0.0684
    E, nu, w2 = channel_energies(k0)
    out['channels_at_k0'] = {'k': k0, 'E': E, 'nu': nu, 'omega2': w2.tolist(), 'Xbar_channel_energy': E - 2 * OMEGA,
                             'neutral_channel_energy': E - OMEGA, 'closed_kappas': kappas(k0).tolist()}
    # ---- local analysis: eigenvalues of W(r) (l = 0); negative => locally propagating
    rr = np.arange(0, 12.01, 0.25)
    loc = []
    for r in rr:
        W = Mmat(max(r, 1e-9)) - np.diag(w2)
        ev = np.linalg.eigvalsh(W)
        Wd = Mmat(max(r, 1e-9), couple=False) - np.diag(w2)
        W2 = (Mmat(max(r, 1e-9)) - np.diag(w2))[:2, :2]
        loc.append({'r': float(r), 'diag_u_barrier': float(Wd[0, 0]), 'eig_min_full3x3': float(ev[0]), 'eig_charged2x2_min': float(np.linalg.eigvalsh(W2)[0]),
                    'F': float(F_i(r)), 'S': float(S_i(r))})
    out['local_eigenvalues_l0'] = loc
    # homogeneous Bogoliubov dispersion at the central density: (k^2 + aF^2 - nu^2)^2 = 4 nu^2 Omega^2 + a^2 F^4
    F0 = float(F_i(0.0)); S0 = float(S_i(0.0))
    kin2 = math.sqrt(4 * nu**2 * OMEGA**2 + a_q**2 * F0**4) - a_q * F0**2 + nu**2
    out['bogoliubov_interior'] = {'F0': F0, 'S0': S0, 'diag_S2_2aF2': S0**2 + 2 * a_q * F0**2, 'anomalous_aF2': a_q * F0**2,
                                  'interior_k2_from_uniform_dispersion': kin2, 'interior_k': math.sqrt(max(kin2, 0)),
                                  'phonon_sound_speed2': a_q * F0**2 / (2 * OMEGA**2 + a_q * F0**2),
                                  'note': 'gapless Bogoliubov branch: a slow X propagates inside the condensate as a phonon; the diagonal barrier S^2+2aF^2-1 > 0 does not by itself forbid the interior.'}
    # ---- validation 1: vacuum
    kap = kappas(k0)
    Wvac = lambda r: np.diag([1.0, 1.0, 6.0]) - np.diag(w2)
    d_vac, _ = phase_shift(k0, 0, lambda r: Wvac(r), kap)
    # ---- validation 2: analytic s-wave square barrier V0 for r < R (channel 0 only)
    V0, R = 0.4, 4.0
    def Wsq(r):
        W = np.diag([1.0, 1.0, 6.0]) - np.diag(w2)
        if r < R: W[0, 0] += V0
        return W
    d_sq, _ = phase_shift(k0, 0, Wsq, kap)
    kap_in = math.sqrt(V0 - k0 * k0)
    d_sq_exact = math.atan(k0 / (kap_in / math.tanh(kap_in * R))) - k0 * R
    d_sq_exact = (d_sq_exact + math.pi / 2) % math.pi - math.pi / 2
    # ---- validation 3: hard sphere limit (V0 = 1e4)
    def Whs(r):
        W = np.diag([1.0, 1.0, 6.0]) - np.diag(w2)
        if r < R: W[0, 0] += 1e4
        return W
    d_hs, _ = phase_shift(k0, 0, Whs, kap, hstep=0.001)
    out['validation'] = {'vacuum_delta0': d_vac, 'square_barrier': {'V0': V0, 'R': R, 'delta_numeric': d_sq, 'delta_exact': d_sq_exact},
                         'hard_sphere_R4': {'delta_numeric': d_hs, 'delta_exact_minus_kR': -k0 * R}}
    # ---- the X-core problem at k0: phase shifts, coupled vs uncoupled diagonal
    lmax = 8
    res = {}
    for label, couple in (('coupled_3channel', True), ('uncoupled_diagonal', False)):
        deltas = []; leaks = []
        for l in range(lmax + 1):
            d, lk = phase_shift(k0, l, xcore_W(k0, l, couple), kap)
            deltas.append(d); leaks.append(lk)
        sig = 4 * math.pi / k0**2 * sum((2 * l + 1) * math.sin(d)**2 for l, d in enumerate(deltas))
        res[label] = {'deltas': deltas, 'closed_channel_leak_at_rmax': leaks, 'sigma_el_mX2': sig,
                      'sigma_l0_only': 4 * math.pi / k0**2 * math.sin(deltas[0])**2,
                      'scattering_length_minus_tan_delta0_over_k': -math.tan(deltas[0]) / k0}
    out['xcore_k0'] = res
    # convergence checks on the coupled s-wave
    conv = {}
    for hs, rm in ((0.002, 40.0), (0.004, 60.0), (0.008, 40.0)):
        d, lk = phase_shift(k0, 0, xcore_W(k0, 0, True), kap, hstep=hs, rmax=rm)
        conv[f'h_{hs}_rmax_{rm}'] = {'delta0': d, 'leak': lk}
    out['xcore_k0']['convergence_l0'] = conv
    # hard-sphere reference used by the reviewer
    hsref = {}
    for Rr in (3.6, 4.65):
        hsref[f'R_{Rr}'] = 4 * math.pi / k0**2 * sum((2 * l + 1) * math.sin(math.atan2(spherical_jn(l, k0 * Rr), spherical_yn(l, k0 * Rr)))**2 for l in range(40))
    out['reviewer_hard_sphere_sigma'] = hsref
    # ---- k dependence of the coupled elastic cross section (s-wave + up to l=6)
    scan = {}
    for k in (0.02, 0.04, 0.0684, 0.1, 0.15, 0.2, 0.3, 0.5):
        kp = kappas(k)
        ds = [phase_shift(k, l, xcore_W(k, l, True), kp)[0] for l in range(7)]
        scan[f'{k:g}'] = {'sigma_el_mX2': 4 * math.pi / k**2 * sum((2 * l + 1) * math.sin(d)**2 for l, d in enumerate(ds)), 'delta0': ds[0], 'delta1': ds[1]}
    out['k_scan_coupled'] = scan
    # ---- rates implied (same n_core and H as the reviewer, labelled)
    n_core = 1.2141071422639381e-18; mX = 0.09240896373753621; H2 = 1.7859349018710206e-24
    for label in ('coupled_3channel', 'uncoupled_diagonal'):
        s = res[label]['sigma_el_mX2']
        GH = n_core * s * k0 * mX / H2
        res[label]['Gamma_over_H_2MeV'] = GH
        res[label]['x_dec_if_GH_gt_1'] = math.sqrt(GH) if GH > 1 else None
    out['reviewer_values'] = {'sigma_mX2': [160, 265], 'Gamma_over_H': [6.9e5, 1.1e6], 'x_dec': [830, 1060]}
    RES.mkdir(parents=True, exist_ok=True)
    (RES / 'sk_xcore.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(json.dumps({k: v for k, v in out.items() if k != 'local_eigenvalues_l0'}, indent=1, default=float))
    print('local eigenvalues (r, diag barrier, min eig 3x3, min eig charged 2x2):')
    for row in loc[::2]:
        print(f"{row['r']:5.2f} {row['diag_u_barrier']:+.4f} {row['eig_min_full3x3']:+.4f} {row['eig_charged2x2_min']:+.4f}  F={row['F']:.3f} S={row['S']:.3f}")


if __name__ == '__main__':
    main()
