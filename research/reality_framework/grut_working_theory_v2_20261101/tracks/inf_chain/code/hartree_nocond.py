#!/usr/bin/env python3
"""inf_chain (Claude, own code): relaxed Hartree X, X-bar, phi gas WITHOUT a condensate, for W above the point where the
circular condensate of the copied solver (dq_entropy/code/hartree_dq.py, imported read-only for its integrals bint, nb, sb
and constants) evaporates (condensed fraction -> 0 near W ~ 5.5e3).  Same Hartree conventions: X mass^2 M2 = 2 lX D + g^2 P,
phi mass^2 g^2 D + 3 lP P, D and P self-consistent, chemical potential mu < M free, fixed comoving N_X and energy W E_c.
Energy: E = lX/4 (8 d^2) + g^2/2 P D + 3 lP/4 P^2 + eX - M2 D/2 + ephi - Mp2 P/2, d = D/2 (the copied solver with rc2 = 0).
Labels: [C]."""
import math, sys
import numpy as np
from scipy.optimize import root
sys.path.insert(0, '/home/user/TestingGRUT/research/reality_framework/grut_working_theory_v2_20261101/tracks/dq_entropy/code')
import hartree_dq as hq


def state(T, mu, D, P):
    LX, G2, LP = hq.LX, hq.G2, hq.LP
    M2 = 2 * LX * D + G2 * P
    Mp2 = G2 * D + 3 * LP * P
    w = lambda k: math.sqrt(k * k + M2)
    W = lambda k: math.sqrt(k * k + max(Mp2, 0.0))
    npl = hq.bint(lambda k: hq.nb((w(k) - mu) / T), T); nmi = hq.bint(lambda k: hq.nb((w(k) + mu) / T), T)
    Dn = hq.bint(lambda k: (hq.nb((w(k) - mu) / T) + hq.nb((w(k) + mu) / T)) / w(k), T)
    Pn = hq.bint(lambda k: hq.nb(W(k) / T) / W(k) if W(k) > 0 else 0.0, T)
    eX = hq.bint(lambda k: w(k) * (hq.nb((w(k) - mu) / T) + hq.nb((w(k) + mu) / T)), T)
    sX = hq.bint(lambda k: hq.sb(hq.nb((w(k) - mu) / T)) + hq.sb(hq.nb((w(k) + mu) / T)), T)
    eP = hq.bint(lambda k: W(k) * hq.nb(W(k) / T), T); sP = hq.bint(lambda k: hq.sb(hq.nb(W(k) / T)), T)
    d = D / 2
    E = LX / 4 * (8 * d * d) + G2 / 2 * P * D + 3 * LP / 4 * P * P + eX - 0.5 * M2 * D + eP - 0.5 * Mp2 * P
    return dict(T=T, mu=mu, D=D, P=P, M2=M2, Q=npl - nmi, E=E, S=sX + sP, n_plus=npl, n_minus=nmi, Dn=Dn, Pn=Pn)


def solve(NX, EX_tot, guess):
    def F(z):
        T, D, P = math.exp(z[0]), math.exp(z[2]), math.exp(z[3])
        M = math.sqrt(2 * hq.LX * D + hq.G2 * P)
        mu = M * math.tanh(z[1])                 # |mu| < M
        s = state(T, mu, D, P)
        return [s['Q'] / NX - 1, s['E'] / EX_tot - 1, s['Dn'] / D - 1, s['Pn'] / P - 1]
    r = root(F, guess, method='hybr', options=dict(xtol=1e-13))
    T, D, P = math.exp(r.x[0]), math.exp(r.x[2]), math.exp(r.x[3])
    M = math.sqrt(2 * hq.LX * D + hq.G2 * P)
    s = state(T, M * math.tanh(r.x[1]), D, P)
    s['resid'] = max(abs(x) for x in F(r.x)); s['x'] = list(r.x); s['mu_over_M'] = math.tanh(r.x[1])
    return s


def scan(NX, EX, Ws_cond, Ws_nocond):
    """condensed branch by continuation (copied solver) up to the last Ws_cond, then the no-condensate branch."""
    hq.set_release(NX, EX)
    Ec = (abs(NX) / 6.977568) ** (4 / 3)
    out = []; g = None; Wp = None
    for W in Ws_cond:
        if Wp is not None:
            for Wm in np.geomspace(Wp, W, int(math.ceil(math.log(W / Wp) / math.log(1.2))) + 1)[1:]:
                gg = list(g); gg[0] += 0.25 * math.log(Wm / Wp)
                st, g = hq.solve(Wm * Ec / EX - 1, True, 0.0, gg); Wp = Wm
        else:
            st, g = hq.solve(W * Ec / EX - 1, True, 0.0, None); Wp = W
        out.append(dict(W=W, branch='condensate', S_per_Q=st['S'] / abs(NX), cond_frac=st['cond_frac'], Xbar_per_Q=st['n_minus'] / abs(NX),
                        resid=max(abs(x) for x in st['resid'])))
    M = math.sqrt(2 * hq.LX * st['D'] + hq.G2 * st['P'])
    guess = [math.log(st['T']), math.atanh(min(0.999999, st['mu'] / M)), math.log(st['D']), math.log(st['P'])]
    Wq = Ws_cond[-1]
    for W in Ws_nocond:
        for Wm in np.geomspace(Wq, W, max(2, int(math.ceil(math.log(max(W, Wq * 1.0001) / Wq) / math.log(1.08))) + 1))[1:]:
            gg = list(guess); f = math.log(Wm / Wq)
            gg[0] += 0.25 * f; gg[2] += 0.5 * f; gg[3] += 0.5 * f; gg[1] = min(gg[1], 2.5)
            s = solve(abs(NX), Wm * Ec, gg); guess = s['x']; Wq = Wm
        out.append(dict(W=W, branch='no condensate', S_per_Q=s['S'] / abs(NX), mu_over_M=s['mu_over_M'], Xbar_per_Q=s['n_minus'] / abs(NX),
                        resid=s['resid']))
    return out
