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
    'K_XXbss_ss': dict(legs=['X', 'Xb', 's', 's', ('s', 'out'), ('s', 'out')], Sin=2, Sout=2),
    'L_sss_XXb': dict(legs=['s', 's', 's', ('X', 'out'), ('Xb', 'out')], Sin=6, Sout=1),
    'N_ssss_ss': dict(legs=['s', 's', 's', 's', ('s', 'out'), ('s', 'out')], Sin=24, Sout=2),
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
    # reverses of the n->2 channels not covered above (for the table-level detailed-balance check)
    'ss_XXXbXb': dict(legs=['s', 's', ('X', 'out'), ('Xb', 'out'), ('X', 'out'), ('Xb', 'out')], Sin=2, Sout=4),
    'Xs_XXXb': dict(legs=['X', 's', ('X', 'out'), ('X', 'out'), ('Xb', 'out')], Sin=1, Sout=2),
    'ss_sss': dict(legs=['s', 's', ('s', 'out'), ('s', 'out'), ('s', 'out')], Sin=2, Sout=6),
    'ss_XXbs': dict(legs=['s', 's', ('X', 'out'), ('Xb', 'out'), ('s', 'out')], Sin=2, Sout=1),
    'XXb_XXbss': dict(legs=['X', 'Xb', ('X', 'out'), ('Xb', 'out'), ('s', 'out'), ('s', 'out')], Sin=1, Sout=2),
    'Xs_XXXbs': dict(legs=['X', 's', ('X', 'out'), ('X', 'out'), ('Xb', 'out'), ('s', 'out')], Sin=1, Sout=2),
    'ss_ssss': dict(legs=['s', 's', ('s', 'out'), ('s', 'out'), ('s', 'out'), ('s', 'out')], Sin=2, Sout=24),
    'ss_XXbss': dict(legs=['s', 's', ('X', 'out'), ('Xb', 'out'), ('s', 'out'), ('s', 'out')], Sin=2, Sout=2),
}

_AMP = {}


def species(leg):
    return leg[0] if isinstance(leg, tuple) else leg


def is_out(leg):
    return isinstance(leg, tuple) and leg[1] == 'out'


def get_amp(legs, freg, fregX=None):
    """Amplitude with regulated propagators: Gamma_reg,s = max(freg m_s, Gamma_s), Gamma_reg,X = fregX m_X."""
    fregX = freg if fregX is None else fregX
    key = (tuple(legs), freg, fregX)
    if key not in _AMP:
        _AMP[key] = Amplitude(legs, widths=(fregX * M.MX, max(freg * M.MS, M.GAMMA_S_VAC)))
    return _AMP[key]


def pairings(sp, freg):
    """Importance-sampling channels for a list of species sp (3 or 4 bodies): one channel per neutral X Xb pair
    (Breit-Wigner map of its invariant mass at m_s^2); for 4 bodies with two disjoint neutral pairs both are mapped.
    Returns a list of channel dicts for phasespace.multichannel (flat single channel if no neutral pair)."""
    n = len(sp)
    wid = M.MS * max(freg * M.MS, M.GAMMA_S_VAC) * 3
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n) if {sp[i], sp[j]} == {'X', 'Xb'}]
    if not pairs:
        return [dict(perm=tuple(range(n)))]
    chans = []
    if n == 3:
        for (i, j) in pairs:
            k = [x for x in range(3) if x not in (i, j)][0]
            chans.append(dict(perm=(i, j, k), peak12=M.MS**2, width12=wid))
    else:
        done = set()
        for (i, j) in pairs:
            rest = tuple(x for x in range(4) if x not in (i, j))
            key = frozenset([(i, j), rest])
            if key in done:
                continue
            done.add(key)
            ch = dict(perm=(i, j) + rest, peak12=M.MS**2, width12=wid)
            if {sp[rest[0]], sp[rest[1]]} == {'X', 'Xb'}:
                ch.update(peak34=M.MS**2, width34=wid)
            chans.append(ch)
    return chans


def _compatible(S1, S2):
    """Two internal lines can belong to the same tree iff their leg sets are disjoint or nested."""
    inter = S1 & S2
    return inter == 0 or inter == S1 or inter == S2


def subtracted_M2(A, mom, win=25.0, max_order=3):
    """|M|^2, |M|^2_sub and the subtraction pieces for momenta mom (N, n, 4).

    M is multilinear in the propagators D_S of its internal lines.  For every set T of flagged lines that can
    coexist in one tree (pairwise disjoint or nested), the coefficient N_T of prod_{S in T} D_S is obtained by
    inclusion-exclusion over amplitudes with subsets of T switched off, and
        |M|^2_sub = sum_T (-2)^|T| prod_{S in T} (Im D_S)^2 |N_T|^2        (T = empty gives |M|^2),
    which is the replacement |D_S|^2 -> Re(D_S^2) for every flagged line (on-shell delta pieces removed, all
    interference and off-shell pieces kept).  Orders |T| <= max_order are included (3 is enough for the 6-point
    channels here: at most three lines can be simultaneously on shell).  Nested lines are the sequential chains
    (e.g. X Xb -> s -> X Xb*, Xb* X -> s) and must be included: skipping them leaves a 1/eps^2 regulator
    artefact (found and fixed on resumption)."""
    Mfull, inv = A(mom, return_invariants=True)
    M2 = np.abs(Mfull)**2
    N = M2.shape[0]
    corr = np.zeros(N)
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
        if near.any():
            flagged[S] = (near, 1.0 / (x + 1j * eps))
    keys = list(flagged)
    pieces = {}
    cache = {}

    def amp_off(off):
        """Amplitude with the lines `off` switched off, evaluated once on the rows where all of them are near shell
        (a superset of every row set that needs it); returns (rows mask, values)."""
        if off not in cache:
            mask = np.logical_and.reduce([flagged[S][0] for S in off])
            cache[off] = (mask, A(mom[mask], switch={S: 0.0 for S in off}))
        return cache[off]
    for order in range(1, max_order + 1):
        for T in itertools.combinations(keys, order):
            if not all(_compatible(a, b) for a, b in itertools.combinations(T, 2)):
                continue
            both = np.logical_and.reduce([flagged[S][0] for S in T])
            if not both.any():
                continue
            # inclusion-exclusion: N_T prod D = sum_{U subset T} (-1)^{|T|-|U|} M(T \ U off), on the rows `both`
            NTD = np.zeros(N, dtype=complex)
            for k in range(order + 1):
                for U in itertools.combinations(T, k):
                    off = tuple(S for S in T if S not in U)
                    if not off:
                        NTD[both] += (-1)**(order - k) * Mfull[both]
                    else:
                        mask, vals = amp_off(off)
                        full_vals = np.zeros(N, dtype=complex); full_vals[mask] = vals
                        NTD[both] += (-1)**(order - k) * full_vals[both]
            prodD = np.ones(N, dtype=complex); prodIm2 = np.ones(N)
            for S in T:
                D = flagged[S][1]
                prodD *= D; prodIm2 *= D.imag**2
            term = np.where(both, (-2.0)**order * prodIm2 * np.abs(NTD)**2 / np.abs(prodD)**2, 0.0)
            corr += term
            pieces[order] = pieces.get(order, 0.0) + term
    sub = -pieces.get(1, np.zeros(N)); add = pieces.get(2, np.zeros(N)); third = pieces.get(3, np.zeros(N))
    return M2, M2 + corr, sub, add + third, len(flagged)


def kernel_n2(ch, sqrt_s, N, rng, freg=0.02, fregX=0.05, win=25.0, ret_mom=False):
    """Gt(s) for an n->2 channel: dict with estimate, MC error, unsubtracted value, subtraction size."""
    spec = CHANNELS_N2[ch]
    legs = spec['legs']
    m_in = [MASS[species(l)] for l in legs if not is_out(l)]
    m_out = [MASS[species(l)] for l in legs if is_out(l)]
    if sqrt_s <= max(sum(m_in), sum(m_out)):
        return dict(G=0.0, err=0.0, G_full=0.0, sub=0.0)
    A = get_amp(legs, freg, fregX)
    sp_in = [species(l) for l in legs if not is_out(l)]
    if len(m_in) == 2:
        mom_in, w_in = PS.generate(sqrt_s, m_in, N, rng)
    else:
        mom_in, w_in = PS.multichannel(sqrt_s, m_in, N, rng, pairings(sp_in, freg))
    N = mom_in.shape[0]
    mom_out, w_out = PS.generate(sqrt_s, m_out, N, rng)
    mom = np.concatenate([mom_in, mom_out], 1)
    M2, M2s, sub, add, nfl = subtracted_M2(A, mom, win)
    w = w_in * w_out
    f = w * M2s
    out = dict(G=float(f.mean()), err=float(f.std() / np.sqrt(N)), G_full=float((w * M2).mean()),
               err_full=float((w * M2).std() / np.sqrt(N)),
               sub=float((w * sub).mean()), add=float((w * add).mean()), lines_flagged=nfl, N=N, freg=freg, fregX=fregX)
    if ret_mom:
        out['mom'] = mom
    return out


def sigma_2n(ch, sqrt_s, N, rng, freg=0.02, fregX=0.05, win=25.0):
    """2->n cross section (GeV^-2) with final-state RIS subtraction."""
    spec = CHANNELS_2N[ch]
    legs = spec['legs']
    m_in = [MASS[species(l)] for l in legs if not is_out(l)]
    m_out = [MASS[species(l)] for l in legs if is_out(l)]
    s = sqrt_s**2
    if sqrt_s <= sum(m_out):
        return dict(sigma=0.0, err=0.0)
    A = get_amp(legs, freg, fregX)
    sp_out = [species(l) for l in legs if is_out(l)]
    mom_out, w_out = PS.multichannel(sqrt_s, m_out, N, rng, pairings(sp_out, freg))
    N = mom_out.shape[0]
    # incoming along z
    p = np.sqrt(kallen(s, m_in[0]**2, m_in[1]**2)) / (2 * sqrt_s)
    E1 = (s + m_in[0]**2 - m_in[1]**2) / (2 * sqrt_s); E2 = sqrt_s - E1
    pin = np.zeros((N, 2, 4)); pin[:, 0] = [E1, 0, 0, p]; pin[:, 1] = [E2, 0, 0, -p]
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
