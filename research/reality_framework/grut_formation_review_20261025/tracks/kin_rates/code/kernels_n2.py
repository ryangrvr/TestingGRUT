#!/usr/bin/env python3
"""kin_rates (Claude, own code): energy-dependent kernels of the n -> 2 (3->2, 4->2) and 2 -> n (2->3, 2->4)
channels by Monte Carlo over the full tree amplitude, with the real-intermediate-state (RIS) subtraction.

Kernel for n -> 2 (MB thermal average factorises on it):
    Gt(s) = int dPhi_n(in; s) int dPhi_2(out; s) |M|^2_sub          [no symmetry factors; GeV^(2n-6)... see units]
    gamma(T) = (1/(S_in S_out)) (T/(8 pi^3)) int ds sqrt s K_1(sqrt s/T) Gt(s)        reactions / volume / time
For 2 -> n:  sigma(s) = (1/S_out) (1/(2 lambda^{1/2}(s,m1^2,m2^2))) int dPhi_n(out) |M|^2_sub.

RIS subtraction (automatic): every internal line S whose propagating field (s for neutral lines, X for charge +-1)
is within |p_S^2 - m^2| < win * m Gamma_reg is decomposed, M = N_S D_S + M_S^off (M is linear in each propagator;
M_S^off is obtained by switching that propagator off), and
    |M|^2_sub = |M|^2 - sum_S 2 (Im D_S)^2 |N_S|^2 + sum_{S1,S2 disjoint} 4 (Im D_1)^2 (Im D_2)^2 |N_12|^2 ,
which replaces |D|^2 -> Re(D^2) in the squared-resonant terms (the on-shell delta pieces are the sequential
processes already in the network) and keeps all interference and double-pole-interference terms.  The result is
regulator independent as Gamma_reg -> 0 (Hadamard finite part); Gamma_reg is taken as a fraction f_reg of the mass
for Monte Carlo variance and the dependence on f_reg is reported.  Importance sampling: incoming X Xbar pairs in
positions (0,1) [and (2,3)] get a mixed flat/Breit-Wigner map of their invariant mass at m_s^2.
"""
import itertools
import json
import sys
import time
import numpy as np
from scipy.special import kve
import model as M
from amp import Amplitude, kallen
import phasespace as PS

MASS = {'X': M.MX, 'Xb': M.MX, 's': M.MS}

# incoming legs first (ordered so that neutral X Xb pairs sit at (0,1) and (2,3)), then outgoing
CHANNELS_N2 = {
    'A_XXXbXb_XXb': dict(legs=['X', 'Xb', 'X', 'Xb', ('X', 'out'), ('Xb', 'out')], Sin=4, Sout=1),
    'B_XXXXb_XX': dict(legs=['X', 'Xb', 'X', 'X', ('X', 'out'), ('X', 'out')], Sin=6, Sout=2),
    'C_XXXbXb_ss': dict(legs=['X', 'Xb', 'X', 'Xb', ('s', 'out'), ('s', 'out')], Sin=4, Sout=2),
    'H_XXbss_XXb': dict(legs=['X', 'Xb', 's', 's', ('X', 'out'), ('Xb', 'out')], Sin=2, Sout=1),
    'J_XXXbs_Xs': dict(legs=['X', 'Xb', 'X', 's', ('X', 'out'), ('s', 'out')], Sin=2, Sout=1),
    'D_XXXb_Xs': dict(legs=['X', 'Xb', 'X', ('X', 'out'), ('s', 'out')], Sin=2, Sout=1),
    'E_XXbs_XXb': dict(legs=['X', 'Xb', 's', ('X', 'out'), ('Xb', 'out')], Sin=1, Sout=1),
    'E2_XXs_XX': dict(legs=['X', 'X', 's', ('X', 'out'), ('X', 'out')], Sin=2, Sout=2),
    'F_sss_ss': dict(legs=['s', 's', 's', ('s', 'out'), ('s', 'out')], Sin=6, Sout=2),
    'G_XXbs_ss': dict(legs=['X', 'Xb', 's', ('s', 'out'), ('s', 'out')], Sin=1, Sout=2),
    'G2_Xss_Xs': dict(legs=['X', 's', 's', ('X', 'out'), ('s', 'out')], Sin=2, Sout=1),
    # 2->2 controls for the engine
    'ctrl_XX_XX': dict(legs=['X', 'X', ('X', 'out'), ('X', 'out')], Sin=2, Sout=2),
    'ctrl_XXb_XXb': dict(legs=['X', 'Xb', ('X', 'out'), ('Xb', 'out')], Sin=1, Sout=1),
    'ctrl_Xs_Xs': dict(legs=['X', 's', ('X', 'out'), ('s', 'out')], Sin=1, Sout=1),
    'ctrl_ss_XXb': dict(legs=['s', 's', ('X', 'out'), ('Xb', 'out')], Sin=2, Sout=1),
}
CHANNELS_2N = {
    'XXb_XXbs': dict(legs=['X', 'Xb', ('X', 'out'), ('Xb', 'out'), ('s', 'out')], Sin=1, Sout=1),
    'XX_XXs': dict(legs=['X', 'X', ('X', 'out'), ('X', 'out'), ('s', 'out')], Sin=2, Sout=2),
    'XXb_sss': dict(legs=['X', 'Xb', ('s', 'out'), ('s', 'out'), ('s', 'out')], Sin=1, Sout=6),
    'Xs_Xss': dict(legs=['X', 's', ('X', 'out'), ('s', 'out'), ('s', 'out')], Sin=1, Sout=2),
    'XXb_XXbXXb': dict(legs=['X', 'Xb', ('X', 'out'), ('Xb', 'out'), ('X', 'out'), ('Xb', 'out')], Sin=1, Sout=4),
    'XX_XXXXb': dict(legs=['X', 'X', ('X', 'out'), ('Xb', 'out'), ('X', 'out'), ('X', 'out')], Sin=2, Sout=6),
    'XXb_ssss': dict(legs=['X', 'Xb', ('s', 'out'), ('s', 'out'), ('s', 'out'), ('s', 'out')], Sin=1, Sout=24),
}

_AMP = {}


def species(leg):
    return leg[0] if isinstance(leg, tuple) else leg


def is_out(leg):
    return isinstance(leg, tuple) and leg[1] == 'out'


def get_amp(legs, freg):
    key = (tuple(legs), freg)
    if key not in _AMP:
        _AMP[key] = Amplitude(legs, widths=(freg * M.MX, max(freg * M.MS, M.GAMMA_S_VAC)))
    return _AMP[key]


def subtracted_M2(A, mom, win=40.0):
    """|M|^2, |M|^2_sub and the subtraction pieces for momenta mom (N, n, 4)."""
    Mfull, inv = A(mom, return_invariants=True)
    M2 = np.abs(Mfull)**2
    N = M2.shape[0]
    sub = np.zeros(N)
    flagged = {}
    for S in A.internal:
        fld = A.propagating(S)
        if fld is None:
            continue
        idx = 0 if fld == 'X' else 2
        m2c = A.m2c[idx]
        m2 = m2c.real; eps = -m2c.imag
        x = inv[S] - m2
        near = np.abs(x) < win * eps
        if not near.any():
            continue
        Moff = A(mom[near], switch={S: 0.0})
        D = 1.0 / (x[near] + 1j * eps)
        Nres2 = np.abs(Mfull[near] - Moff)**2 / np.abs(D)**2
        sub[near] += 2.0 * D.imag**2 * Nres2
        full_off = np.zeros(N, dtype=complex); full_off[near] = Moff
        flagged[S] = (near, full_off, D, eps, m2)
    add = np.zeros(N)
    keys = list(flagged)
    for S1, S2 in itertools.combinations(keys, 2):
        if S1 & S2:
            continue
        n1, off1, _, e1, m21 = flagged[S1]; n2, off2, _, e2, m22 = flagged[S2]
        both = n1 & n2
        if not both.any():
            continue
        Mboth = A(mom[both], switch={S1: 0.0, S2: 0.0})
        D1 = 1.0 / (inv[S1][both] - m21 + 1j * e1); D2 = 1.0 / (inv[S2][both] - m22 + 1j * e2)
        N12 = Mfull[both] - off1[both] - off2[both] + Mboth
        add[both] += 4.0 * D1.imag**2 * D2.imag**2 * np.abs(N12)**2 / (np.abs(D1)**2 * np.abs(D2)**2)
    return M2, M2 - sub + add, sub, add, len(flagged)


def kernel_n2(ch, sqrt_s, N, rng, freg=0.02, win=40.0, mix=0.5, ret_mom=False):
    """Gt(s) for an n->2 channel: dict with estimate, MC error, unsubtracted value, subtraction size."""
    spec = CHANNELS_N2[ch]
    legs = spec['legs']
    m_in = [MASS[species(l)] for l in legs if not is_out(l)]
    m_out = [MASS[species(l)] for l in legs if is_out(l)]
    if sqrt_s <= max(sum(m_in), sum(m_out)):
        return dict(G=0.0, err=0.0, G_full=0.0, sub=0.0)
    A = get_amp(legs, freg)
    sp_in = [species(l) for l in legs if not is_out(l)]
    kw = {}
    if len(sp_in) >= 2 and set(sp_in[:2]) == {'X', 'Xb'}:
        kw.update(peak12=M.MS**2, width12=M.MS * max(freg * M.MS, M.GAMMA_S_VAC) * 3)
    if len(sp_in) == 4 and set(sp_in[2:]) == {'X', 'Xb'}:
        kw.update(peak34=M.MS**2, width34=M.MS * max(freg * M.MS, M.GAMMA_S_VAC) * 3)
    mom_in, w_in = PS.generate(sqrt_s, m_in, N, rng, **kw)
    mom_out, w_out = PS.generate(sqrt_s, m_out, N, rng)
    mom = np.concatenate([mom_in, mom_out], 1)
    M2, M2s, sub, add, nfl = subtracted_M2(A, mom, win)
    w = w_in * w_out
    f = w * M2s
    out = dict(G=float(f.mean()), err=float(f.std() / np.sqrt(N)), G_full=float((w * M2).mean()),
               sub=float((w * sub).mean()), add=float((w * add).mean()), lines_flagged=nfl, N=N, freg=freg)
    if ret_mom:
        out['mom'] = mom
    return out


def sigma_2n(ch, sqrt_s, N, rng, freg=0.02, win=40.0):
    """2->n cross section (GeV^-2) with final-state RIS subtraction."""
    spec = CHANNELS_2N[ch]
    legs = spec['legs']
    m_in = [MASS[species(l)] for l in legs if not is_out(l)]
    m_out = [MASS[species(l)] for l in legs if is_out(l)]
    s = sqrt_s**2
    if sqrt_s <= sum(m_out):
        return dict(sigma=0.0, err=0.0)
    A = get_amp(legs, freg)
    # incoming along z
    p = np.sqrt(kallen(s, m_in[0]**2, m_in[1]**2)) / (2 * sqrt_s)
    E1 = (s + m_in[0]**2 - m_in[1]**2) / (2 * sqrt_s); E2 = sqrt_s - E1
    pin = np.zeros((N, 2, 4)); pin[:, 0] = [E1, 0, 0, p]; pin[:, 1] = [E2, 0, 0, -p]
    kw = {}
    sp_out = [species(l) for l in legs if is_out(l)]
    if len(sp_out) >= 2 and set(sp_out[:2]) == {'X', 'Xb'}:
        kw.update(peak12=M.MS**2, width12=M.MS * max(freg * M.MS, M.GAMMA_S_VAC) * 3)
    if len(sp_out) == 4 and set(sp_out[2:]) == {'X', 'Xb'}:
        kw.update(peak34=M.MS**2, width34=M.MS * max(freg * M.MS, M.GAMMA_S_VAC) * 3)
    mom_out, w_out = PS.generate(sqrt_s, m_out, N, rng, **kw)
    mom = np.concatenate([pin, mom_out], 1)
    M2, M2s, sub, add, nfl = subtracted_M2(A, mom, win)
    flux = 2 * np.sqrt(kallen(s, m_in[0]**2, m_in[1]**2))
    f = w_out * M2s / flux / spec['Sout']
    return dict(sigma=float(f.mean()), err=float(f.std() / np.sqrt(N)), sigma_full=float((w_out * M2).mean() / flux / spec['Sout']),
                sub=float((w_out * sub).mean() / flux / spec['Sout']), lines_flagged=nfl, N=N, freg=freg)


def thermal_from_table(sq_grid, G_grid, T, Sin, Sout, sq_th):
    """gamma(T) = (1/(Sin Sout)) (T/8 pi^3) int ds sqrt s K1(sqrt s/T) Gt(s), Gt interpolated log-log in (sqrt s - sq_th)."""
    x = np.log(np.asarray(sq_grid) - sq_th)
    y = np.log(np.maximum(np.asarray(G_grid), 1e-300))
    sq = np.geomspace(sq_th * (1 + 1e-9), sq_th + 60 * T + 5 * sq_th, 4000)
    # extend beyond the table with the last slope
    lx = np.log(sq - sq_th)
    ly = np.interp(lx, x, y)
    hi = lx > x[-1]
    if hi.any():
        slope = (y[-1] - y[-2]) / (x[-1] - x[-2])
        ly[hi] = y[-1] + slope * (lx[hi] - x[-1])
    G = np.exp(ly)
    integrand = 2 * sq * sq * kve(1, sq / T) * np.exp(-sq / T) * G      # ds = 2 sq d(sq)
    val = np.trapezoid(integrand, sq)
    return val * T / (8 * np.pi**3) / (Sin * Sout)


if __name__ == '__main__':
    rng = np.random.default_rng(5)
    import xsec22 as X
    t0 = time.time()
    # control: XX->XX (no poles)
    for sq in (2.2, 3.0, 6.0):
        s = (sq * M.MX)**2
        r = kernel_n2('ctrl_XX_XX', sq * M.MX, 100000, rng)
        ana = kallen(s, M.MX**2, M.MX**2) * 2 * X.sigma('XX_XX', s, X.Params()) / (4 * np.pi * s)
        print(f'ctrl XX->XX sqrt s={sq} m: MC {r["G"]:.5e} +- {r["err"]:.1e}  analytic {ana:.5e}  ratio {r["G"]/ana:.4f}')
    # control: XXb->XXb at the regulator width, subtracted (same mathematical object by quadrature)
    for freg in (0.05, 0.02):
        Pw = X.Params(Gamma_s=freg * M.MS)
        for sq in (2.3, 2.449, 2.6, 4.0):
            s = (sq * M.MX)**2
            r = kernel_n2('ctrl_XXb_XXb', sq * M.MX, 200000, rng, freg=freg)
            ana = kallen(s, M.MX**2, M.MX**2) * X.sigma('XXb_XXb', s, Pw, mode='sub') / (4 * np.pi * s)
            anaf = kallen(s, M.MX**2, M.MX**2) * X.sigma('XXb_XXb', s, Pw, mode='full') / (4 * np.pi * s)
            print(f'ctrl XXb->XXb f={freg} sqrt s={sq}: MC sub {r["G"]:.5e} +- {r["err"]:.1e} quad sub {ana:.5e} | MC full {r["G_full"]:.5e} quad full {anaf:.5e} flagged {r["lines_flagged"]}')
    # control: Xs->Xs (u-pole inside the angular range for 3.6 < sqrt s/m < 5)
    for freg in (0.05, 0.02):
        Pw = X.Params(Gamma_X_reg=freg)
        for sq in (3.8, 4.5, 6.0):
            s = (sq * M.MX)**2
            r = kernel_n2('ctrl_Xs_Xs', sq * M.MX, 200000, rng, freg=freg)
            ana = kallen(s, M.MX**2, M.MS**2) * X.sigma('Xs_Xs', s, Pw, mode='sub') / (4 * np.pi * s)
            print(f'ctrl Xs->Xs f={freg} sqrt s={sq}: MC sub {r["G"]:.5e} +- {r["err"]:.1e} quad sub {ana:.5e} ratio {r["G"]/ana:.4f}')
    # threshold limit of A and B vs the threshold constants
    for ch, Aold, Sout in (('A_XXXbXb_XXb', 0.0842146549366979, 1), ('B_XXXXb_XX', 0.002361976125302127, 2)):
        sq = 4 * M.MX * (1 + 1e-4)
        r = kernel_n2(ch, sq, 50000, rng)
        mom, w = PS.generate(sq, [M.MX] * 4, 50000, rng)
        phi4 = w.mean()
        print(f'{ch} near threshold: G/Phi4 = {r["G"]/phi4:.5e}; A_old*(2m)^4*Sout = {Aold*(2*M.MX)**4*Sout:.5e}')
    print('time', time.time() - t0)
