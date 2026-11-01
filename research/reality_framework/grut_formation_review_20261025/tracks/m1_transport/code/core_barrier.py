"""Track m1_transport, [E]: elastic X-core scattering at k = m_X v_rms, hard-sphere / barrier estimate.

Independent code. Reads Codex's stationary Q=2e5 M1 profile (d4747cb PROFILE_M1.csv, columns rho=m_X r, F, S)
only to locate the radius where the X 'barrier' S^2 + 2 a F^2 - 1 (bf8ab3b DERIVATION: charged normal diagonal
S^2 + 2aF^2, a = 0.1) changes sign, then evaluates the hard-sphere partial-wave cross section
sigma = (4 pi/k^2) sum (2l+1) sin^2 delta_l, tan delta_l = j_l(kR)/y_l(kR), at k = 0.0684 m_X.
This is an order-of-magnitude estimate of a channel Codex's transport omits, not a solved three-channel problem.
"""
import json, math, os
from pathlib import Path
import numpy as np
from scipy.special import spherical_jn, spherical_yn
ROOT = Path(__file__).resolve().parents[1]
SP = Path(os.environ.get('SP', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad')) / 'review' / 'm1_transport'
a_q = 0.1
k = 0.0684
out = {'k_over_mX': k, 'unitarity_s_wave_4pi_over_k2': 4 * math.pi / k**2}
prof = SP / 'codex' / 'PROFILE_M1.csv'
if prof.exists():
    data = np.genfromtxt(prof, delimiter=',', names=True)
    names = data.dtype.names
    r = data['rho']; F = data['F']; S = data['S']
    barrier = S**2 + 2 * a_q * F**2 - 1.0
    # radius where the X barrier drops below k^2 (classically allowed outside)
    idx = np.where(barrier > k**2)[0]
    R_barrier = float(r[idx[-1]]) if len(idx) else None
    half = np.where(F >= 0.5 * F[0])[0]
    out['profile'] = {'columns': names, 'F0': float(F[0]), 'S0': float(S[0]), 'barrier_center': float(barrier[0]),
                      'barrier_max': float(barrier.max()), 'r_barrier_edge_mX^-1': R_barrier, 'r_half_F_mX^-1': float(r[half[-1]]),
                      'r_where_F_below_0.01F0': float(r[np.where(F > 0.01 * F[0])[0][-1]])}
out['hard_sphere'] = {}
for R in (3.6, 4.65, 6., 10., 15., 20.):
    kR = k * R
    sig = 4 * math.pi / k**2 * sum((2 * l + 1) * math.sin(math.atan2(spherical_jn(l, kR), spherical_yn(l, kR)))**2 for l in range(60))
    out['hard_sphere'][f'R_{R:g}'] = {'kR': kR, 'sigma_mX2': sig, 'sigma_over_4piR2': sig / (4 * math.pi * R * R)}
n_core_2MeV = 1.2141071422639381e-18
mX = 0.09240896373753621; H2 = 1.7859349018710206e-24
out['rates'] = {R: {'Gamma_over_H_2MeV': n_core_2MeV * v['sigma_mX2'] * k * mX / H2, 'x_dec_RD': math.sqrt(n_core_2MeV * v['sigma_mX2'] * k * mX / H2)} for R, v in out['hard_sphere'].items()}
out['codex_XX_rate_for_comparison'] = 15.9
out['note'] = ('For X at k = 0.068 m_X the core interior (S^2 + 2 a F^2 > 1 + k^2) is classically forbidden, so the core acts as a '
               'nearly hard sphere of radius ~ r_barrier_edge (3.6/m_X) to r_half_F (4.65/m_X); the field tail extends to ~12/m_X; the anomalous (X Xbar) and neutral couplings add inelastic channels that '
               'Codex found to be ~1 % at k = 0.4-1. Codex evolve.c/elastic.py include only XX self-scattering.')
(ROOT / 'results' / 'core_barrier.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(out, indent=1))
