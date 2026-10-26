#!/usr/bin/env python3
"""kin_rates (Claude, own code): lab-frame Monte Carlo of thermal collision integrals with Bose statistics.

Sampling: for species of mass m at temperature T the momentum magnitude is drawn from p^2 exp(-(E-m)/T) (inverse
CDF on a grid); the measure int d^3p/((2pi)^3 2E) f(E) is then  mean over samples of  w = Z F(E)/(4 pi^2 E),
Z = m^2 T K2(m/T) e^{m/T} (= m^2 T kve(2,m/T)), F(E) = f(E) e^{(E-m)/T}:
    MB: F = e^{(mu-m)/T};   BE: F = e^{(mu-m)/T} / (1 - e^{-(E-mu)/T}).
2->2 rate density: gamma = (1/(S_in S_out)) < w1 w2 Phi2(s) (1+f3)(1+f4) |M|^2_sub >, with the final pair generated
isotropically in the pair rest frame (Phi2 = |k*|/(4 pi sqrt s)) and boosted to the lab.
Outputs: Bose/MB ratios of the thermal rates; energy-exchange coefficients between species at different
temperatures (elastic X s); hard-particle energy-loss and collision rates on a thermal bath; detailed-balance
checks of s <-> X Xb (decay by MC, inverse decay by exact 2-D quadrature) and s s <-> X Xb with Bose factors and
chemical potentials; and the point-kinematics gain/loss identity samples.
"""
import json
import sys
import numpy as np
from scipy.special import kve
from scipy import integrate
import model as M
import xsec22 as X
from amp import kallen
import phasespace as PS

_CDF = {}


def thermal_momenta(m, T, N, rng):
    key = (round(m, 12), round(T, 12))
    if key not in _CDF:
        pmax = 45 * T + 12 * np.sqrt(2 * m * T) + m
        grid = np.linspace(0, pmax, 20001)
        E = np.sqrt(grid * grid + m * m)
        dens = grid * grid * np.exp(-(E - m) / T)
        cdf = np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(grid))])
        _CDF[key] = (grid, cdf / cdf[-1])
    grid, cdf = _CDF[key]
    p = np.interp(rng.uniform(0, 1, N), cdf, grid)
    c = rng.uniform(-1, 1, N); ph = rng.uniform(0, 2 * np.pi, N); sn = np.sqrt(1 - c * c)
    E = np.sqrt(p * p + m * m)
    return np.column_stack([E, p * sn * np.cos(ph), p * sn * np.sin(ph), p * c])


def Z(m, T):
    return m * m * T * kve(2, m / T)


def Fw(E, m, T, mu, stat):
    if stat == 'MB':
        return np.exp((mu - m) / T) * np.ones_like(E)
    return np.exp((mu - m) / T) / (1 - np.exp(-(E - mu) / T))


def fdist(E, T, mu, stat):
    if stat == 'MB':
        return np.exp(-(E - mu) / T)
    return 1.0 / np.expm1((E - mu) / T)


def weights(mom, m, T, mu, stat):
    E = mom[:, 0]
    return Z(m, T) * Fw(E, m, T, mu, stat) / (4 * np.pi**2 * E)


def m2_sub(ch, s, t, u, P, mode='sub'):
    """|M|^2 with Re(D^2) in place of |D|^2 on the physical pole (RIS subtracted), vectorised over s,t,u."""
    if ch == 'XXb_XXb':
        A = -4 * P.lX - P.lam3**2 / (t - P.ms2); B = -P.lam3**2; x = s - P.ms2; eps = P.eps_s
    elif ch == 'Xs_Xs':
        A = -2 * P.g2 - P.lam3**2 / (s - P.mX2) - P.lam3 * P.lam_s3 / (t - P.ms2); B = -P.lam3**2; x = u - P.mX2; eps = P.eps_X
    elif ch == 'XX_XX':
        return (-4 * P.lX - P.lam3**2 * (1 / (t - P.ms2) + 1 / (u - P.ms2)))**2
    elif ch == 'ss_ss':
        return (-P.lam_s4 - P.lam_s3**2 * (1 / (s - P.ms2) + 1 / (t - P.ms2) + 1 / (u - P.ms2)))**2
    elif ch in ('ss_XXb', 'XXb_ss'):
        return (-2 * P.g2 - P.lam_s3 * P.lam3 / (s - P.ms2) - P.lam3**2 * (1 / (t - P.mX2) + 1 / (u - P.mX2)))**2
    else:
        raise KeyError(ch)
    D = 1 / (x + 1j * eps)
    if mode == 'nores':
        return A * A
    if mode == 'full':
        return A * A + 2 * A * B * D.real + B * B * np.abs(D)**2
    return A * A + 2 * A * B * D.real + B * B * (D * D).real


SPEC = {  # (m1,m2,m3,m4) keys, S_in, S_out, species names for the Bose factors (X,Xb,s)
    'XXb_XXb': (('mX', 'mX', 'mX', 'mX'), 1, 1, ('X', 'Xb', 'X', 'Xb')),
    'XX_XX': (('mX', 'mX', 'mX', 'mX'), 2, 2, ('X', 'X', 'X', 'X')),
    'Xs_Xs': (('mX', 'ms', 'mX', 'ms'), 1, 1, ('X', 's', 'X', 's')),
    'ss_ss': (('ms', 'ms', 'ms', 'ms'), 2, 2, ('s', 's', 's', 's')),
    'ss_XXb': (('ms', 'ms', 'mX', 'mX'), 2, 1, ('s', 's', 'X', 'Xb')),
    'XXb_ss': (('mX', 'mX', 'ms', 'ms'), 1, 2, ('X', 'Xb', 's', 's')),
}


def rate_22(ch, T, P, N, rng, stat='MB', mu=None, Ts=None, mode='sub', energy=False):
    """Thermal rate density of a 2->2 channel.  mu: dict of chemical potentials per species; Ts: dict of
    temperatures per species (default T for all).  energy=True also returns <(E4-E2)> flux (energy gained by
    species 2 per volume per time) and <(E3-E1)>."""
    keys, Sin, Sout, sp = SPEC[ch]
    m = [getattr(P, k) for k in keys]
    mu = mu or {}
    Ts = Ts or {}
    T1, T2 = Ts.get(sp[0], T), Ts.get(sp[1], T)
    p1 = thermal_momenta(m[0], T1, N, rng); p2 = thermal_momenta(m[1], T2, N, rng)
    w = weights(p1, m[0], T1, mu.get(sp[0], 0.0), stat) * weights(p2, m[1], T2, mu.get(sp[1], 0.0), stat)
    Ptot = p1 + p2
    s = Ptot[:, 0]**2 - np.sum(Ptot[:, 1:]**2, 1)
    ok = s > (m[2] + m[3])**2 * (1 + 1e-12)
    k1, k2, phi2 = PS.two_body(Ptot, m[2], m[3], rng, N)
    t = (p1 - k1); u = (p1 - k2)
    t = t[:, 0]**2 - np.sum(t[:, 1:]**2, 1); u = u[:, 0]**2 - np.sum(u[:, 1:]**2, 1)
    M2 = m2_sub(ch, s, t, u, P, mode)
    bose = np.ones(N)
    if stat == 'BE':
        bose = (1 + fdist(k1[:, 0], Ts.get(sp[2], T), mu.get(sp[2], 0.0), 'BE')) * (1 + fdist(k2[:, 0], Ts.get(sp[3], T), mu.get(sp[3], 0.0), 'BE'))
    f = np.where(ok, w * phi2 * bose * M2, 0.0) / (Sin * Sout)
    out = dict(gamma=float(f.mean()), err=float(f.std() / np.sqrt(N)))
    if energy:
        g2 = f * (k2[:, 0] - p2[:, 0]); g1 = f * (k1[:, 0] - p1[:, 0])
        out.update(Q2=float(g2.mean()), Q2_err=float(g2.std() / np.sqrt(N)), Q1=float(g1.mean()), Q1_err=float(g1.std() / np.sqrt(N)))
    return out


def decay_bose(T, P, N, rng, mu_s=0.0, mu_p=0.0, mu_m=0.0):
    """gamma_{s->X Xb} with Bose statistics: int dPi_s f_s int dPhi_2 lam3^2 (1+f_X)(1+f_Xb)  (MC)."""
    ps = thermal_momenta(P.ms, T, N, rng)
    w = weights(ps, P.ms, T, mu_s, 'BE')
    k1, k2, phi2 = PS.two_body(ps, P.mX, P.mX, rng, N)
    bose = (1 + fdist(k1[:, 0], T, mu_p, 'BE')) * (1 + fdist(k2[:, 0], T, mu_m, 'BE'))
    f = w * phi2 * P.lam3**2 * bose
    return float(f.mean()), float(f.std() / np.sqrt(N))


def inverse_decay_bose(T, P, mu_s=0.0, mu_p=0.0, mu_m=0.0):
    """gamma_{X Xb -> s} with Bose statistics by exact 2-D quadrature:
    gamma = int dPi1 dPi2 f1 f2 (1+f_s) lam3^2 2pi delta(s - ms^2),  angular delta -> 1/(2 p1 p2) on |c0|<=1."""
    m, ms2 = P.mX, P.ms2

    def inner(p1):
        E1 = np.sqrt(p1 * p1 + m * m)

        def g(p2):
            E2 = np.sqrt(p2 * p2 + m * m)
            c0 = (2 * m * m + 2 * E1 * E2 - ms2) / (2 * p1 * p2)
            if abs(c0) > 1:
                return 0.0
            f1 = 1 / np.expm1((E1 - mu_p) / T); f2 = 1 / np.expm1((E2 - mu_m) / T)
            fs = 1 / np.expm1((E1 + E2 - mu_s) / T)
            return p2 * p2 / (2 * E2) * f1 * f2 * (1 + fs) * 2 * np.pi / (2 * p1 * p2)
        # p2 range where |c0| <= 1
        pmax = 60 * T + 5 * ms2**0.5
        return p1 * p1 / (2 * E1) * integrate.quad(g, 1e-9, pmax, limit=400, points=[m, 3 * T])[0]
    pmax = 60 * T + 5 * ms2**0.5
    val = integrate.quad(inner, 1e-9, pmax, limit=400, points=[m, 3 * T])[0]
    # measure factors: d^3p/((2pi)^3 2E) -> 4 pi p^2 dp /((2pi)^3 2E) for p1; 2 pi p^2 dp dcos/((2pi)^3 2E) for p2
    return val * P.lam3**2 * (4 * np.pi) * (2 * np.pi) / (2 * np.pi)**6


def hard_particle(E, T, P, N, rng, species_bath=('X', 'Xb', 's'), stat='BE'):
    """A test X of lab energy E on a thermal bath: collision rate (1/time) and energy-loss rate per channel."""
    out = {}
    p1 = np.zeros((N, 4)); p1[:, 0] = E; p1[:, 3] = np.sqrt(E * E - P.mX2)
    for b in species_bath:
        ch = {'X': 'XX_XX', 'Xb': 'XXb_XXb', 's': 'Xs_Xs'}[b]
        keys, Sin, Sout, sp = SPEC[ch]
        mb = P.ms if b == 's' else P.mX
        p2 = thermal_momenta(mb, T, N, rng)
        w2 = weights(p2, mb, T, 0.0, stat)
        Ptot = p1 + p2
        s = Ptot[:, 0]**2 - np.sum(Ptot[:, 1:]**2, 1)
        k1, k2, phi2 = PS.two_body(Ptot, P.mX, mb, rng, N)
        t = p1 - k1; u = p1 - k2
        t = t[:, 0]**2 - np.sum(t[:, 1:]**2, 1); u = u[:, 0]**2 - np.sum(u[:, 1:]**2, 1)
        mode = 'nores' if ch == 'XXb_XXb' else 'sub'
        M2 = m2_sub(ch, s, t, u, P, mode)
        bose = np.ones(N)
        if stat == 'BE':
            bose = (1 + fdist(k1[:, 0], T, 0.0, 'BE')) * (1 + fdist(k2[:, 0], T, 0.0, 'BE'))
        # rate per test particle: (1/(2E)) int dPi2 f2 int dPhi2 |M|^2 (1+f)(1+f) / S_out
        f = w2 * phi2 * bose * M2 / (2 * E) / Sout
        dE = E - k1[:, 0]                       # energy lost by the test particle (continuation = particle 3)
        dE_max = E - np.maximum(k1[:, 0], k2[:, 0])
        out[b] = dict(rate=float(f.mean()), rate_err=float(f.std() / np.sqrt(N)),
                      dEdt=float((f * dE).mean()), dEdt_err=float((f * dE).std() / np.sqrt(N)),
                      dEdt_fastest_survives=float((f * dE_max).mean()),
                      mean_s_over_m2=float(np.mean(s) / P.mX2))
    return out


def point_kinematics_checks(P, rng, nsamp=80):
    """Gain/loss identity prod_in f prod_out (1+f) = prod_out f prod_in (1+f) for on-shell, energy-conserving
    configurations with sum mu_in = sum mu_out (Bose).  Channels: s<->XXb (boosted), ss<->XXb, 4->2 A, 3->2 D."""
    res = []
    for _ in range(nsamp):
        T = P.mX * 10**rng.uniform(-0.5, 1.2)
        mu = P.mX * rng.uniform(-0.7, 0.7)
        cases = []
        # s -> X Xb with a boosted parent: energies from a two-body decay in a moving frame
        b = rng.uniform(0, 0.8); g = 1 / np.sqrt(1 - b * b); c = rng.uniform(-1, 1)
        pst = np.sqrt(P.ms2 / 4 - P.mX2)
        Es = g * P.ms; E1 = g * (P.ms / 2 + b * pst * c); E2 = Es - E1
        cases.append(([Es], [2 * mu], [E1, E2], [mu, -mu]))
        # ss -> X Xb: random CM energy and angles
        sq = 2 * P.ms * (1 + rng.uniform(0, 2))
        k = np.sqrt(sq * sq / 4 - P.mX2); ps = np.sqrt(sq * sq / 4 - P.ms2)
        cases.append(([sq / 2 + 0.1 * ps, sq / 2 - 0.1 * ps], [2 * mu, 2 * mu], [sq / 2 + 0.3 * k, sq / 2 - 0.3 * k], [mu, -mu]))
        # X X Xb Xb -> X Xb
        E = P.mX * (1 + rng.uniform(0, 2, 4)); Etot = E.sum()
        cases.append((list(E), [mu, mu, -mu, -mu], [Etot * 0.6, Etot * 0.4], [mu, -mu]))
        # X X Xb -> X s
        E = P.mX * (1 + rng.uniform(0, 3, 3)); Etot = E.sum()
        cases.append((list(E), [mu, mu, -mu], [Etot - P.ms * 1.2, P.ms * 1.2], [mu, 2 * mu]))
        for ein, uin, eout, uout in cases:
            fi = [1 / np.expm1((e - u) / T) for e, u in zip(ein, uin)]
            fo = [1 / np.expm1((e - u) / T) for e, u in zip(eout, uout)]
            fwd = np.prod(fi) * np.prod([1 + f for f in fo]); bwd = np.prod(fo) * np.prod([1 + f for f in fi])
            res.append(abs(fwd - bwd) / (fwd + bwd))
    return dict(samples=len(res), max_residual=float(max(res)), mean_residual=float(np.mean(res)))


def main(outpath):
    rng = np.random.default_rng(2026)
    P = X.Params(Gamma_X_reg=0.05)
    out = dict(scope='Lab-frame MC thermal integrals (own code); Bose statistics; MB references from xsec22 quadrature.')
    N = 400000
    # 1. Bose/MB ratios and detailed balance of the thermal rates
    rows = []
    for t in (0.3, 1.0, 3.0, 12.0):
        T = t * M.MX
        row = dict(T_over_mX=t)
        for ch, mode in (('XXb_XXb', 'nores'), ('XX_XX', 'sub'), ('Xs_Xs', 'sub'), ('ss_ss', 'sub'), ('ss_XXb', 'sub'), ('XXb_ss', 'sub')):
            mb = rate_22(ch, T, P, N, rng, 'MB', mode=mode)
            be = rate_22(ch, T, P, N, rng, 'BE', mode=mode)
            ref = X.thermal_rate(ch, T, P, mode=mode)
            row[ch] = dict(gamma_MB_mc=mb['gamma'], gamma_MB_mc_err=mb['err'], gamma_MB_quad=ref,
                           gamma_BE_mc=be['gamma'], gamma_BE_mc_err=be['err'],
                           bose_over_MB=be['gamma'] / mb['gamma'] if mb['gamma'] else None, mode=mode)
        # Bose detailed balance ss <-> XXb with chemical potentials mu_X = mu_Xb = 0.3 m, mu_s = 0.6 m
        mu = dict(X=0.3 * M.MX, Xb=0.3 * M.MX, s=0.6 * M.MX)
        f = rate_22('ss_XXb', T, P, N, rng, 'BE', mu=mu); r = rate_22('XXb_ss', T, P, N, rng, 'BE', mu=mu)
        row['detailed_balance_ss_XXb_BE_mu'] = dict(forward=f['gamma'], forward_err=f['err'], reverse=r['gamma'], reverse_err=r['err'],
                                                   ratio=f['gamma'] / r['gamma'])
        # s <-> XXb with Bose factors: decay (MC) vs inverse decay (quadrature), mu_s = mu_+ + mu_-
        for (ms_, mp_, mm_) in ((0.0, 0.0, 0.0), (0.5 * M.MX, 0.4 * M.MX, 0.1 * M.MX)):
            d, de = decay_bose(T, P, N, rng, ms_, mp_, mm_)
            inv = inverse_decay_bose(T, P, ms_, mp_, mm_)
            row[f'detailed_balance_s_XXb_BE_mus{ms_/M.MX:.1f}'] = dict(decay_mc=d, decay_err=de, inverse_quad=inv, ratio=d / inv,
                                                                       decay_MB_analytic=X.decay_rate_density(T, P) if ms_ == 0 else None)
        rows.append(row)
        print(json.dumps(row, indent=None, default=float)[:600], flush=True)
    out['thermal_rates'] = rows
    # 2. energy exchange X s elastic at two temperatures (MB): Q_s/(n_X n_s) in GeV^-1 ... (GeV per unit n n time)
    ex = []
    for tX in (0.5, 1.0, 3.0):
        for ts in (0.5, 1.0, 3.0):
            TX, Tsx = tX * M.MX, ts * M.MX
            r = rate_22('Xs_Xs', TX, P, N, rng, 'MB', Ts={'X': TX, 's': Tsx}, energy=True)
            nX = X.n_eq(P.mX, TX); ns = X.n_eq(P.ms, Tsx)
            ex.append(dict(T_X_over_m=tX, T_s_over_m=ts, gamma=r['gamma'], Q_s=r['Q2'], Q_s_err=r['Q2_err'], Q_X=r['Q1'],
                           sum_Q_over_abs=(r['Q1'] + r['Q2']) / (abs(r['Q1']) + abs(r['Q2']) + 1e-300),
                           Q_s_per_nXns_GeV_minus1=r['Q2'] / (nX * ns), gamma_per_nXns=r['gamma'] / (nX * ns)))
    out['energy_exchange_Xs'] = ex
    print('energy exchange done', flush=True)
    # 3. hard particle on a thermal bath
    hp = []
    for t in (0.3, 1.0, 3.0):
        T = t * M.MX
        for Eom in (3.0, 10.0, 100.0, 1e3, 1e4, 7.7e4):
            r = hard_particle(Eom * M.MX, T, P, 200000, rng)
            hp.append(dict(T_over_m=t, E_over_m=Eom, channels=r))
    out['hard_particle'] = hp
    print('hard particle done', flush=True)
    # 4. point kinematics identity samples
    out['point_kinematics'] = point_kinematics_checks(P, rng)
    json.dump(out, open(outpath, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main(sys.argv[1])
