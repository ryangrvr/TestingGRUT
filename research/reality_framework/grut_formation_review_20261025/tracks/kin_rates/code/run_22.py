#!/usr/bin/env python3
"""kin_rates: tables of the 2->2 cross sections, the 1<->2 rates and their MB thermal averages (own quadrature code,
xsec22.py).  usage: python3 -B run_22.py <out.json>
sigma(s) on a geometric grid from threshold to sqrt s = 1e6 m_X; modes: 'sub' (RIS-subtracted finite part),
'nores' (non-resonant A^2 only), for Xs_Xs also 'res' and 'int'; for XXb_XXb the components (sigma0, sigma1, sigma2)
with sigma = sigma0 + sigma1 Re D + sigma2 {Re D^2 | |D|^2}, D = 1/(s - ms^2 + i ms Gamma_s), so that the solver
can rebuild any prescription.  Transport cross sections (weight 1 - cos theta) for the elastic channels.
Thermal averages gamma(T) (reactions per volume per time, MB, mu = 0) and <sigma v> = gamma S_in/(n1 n2) for T/m_X
from 0.03 to 30; the XXb_XXb s-pole integrated with the finite-part quadrature (fp_quad/pv_quad)."""
import json
import sys
import time
import numpy as np
import model as M
import xsec22 as X

P = X.Params()          # physical s width; X regulator 1e-7 m_X (finite part limit)
CH = ['XXb_XXb', 'XX_XX', 'Xs_Xs', 'ss_ss', 'ss_XXb', 'XXb_ss']
TGRID = sorted(set(np.round(np.concatenate([np.geomspace(0.03, 30, 46), [0.3, 1.0, 3.0, 7.344, 8.118, 12.0, 13.06]]), 6)))


def sigma_tables():
    out = {}
    for ch in CH:
        m, Sin, Sout = X.masses_of(ch, P)
        sq_th = max(m[0] + m[1], m[2] + m[3])
        deltas = np.geomspace(1e-4, 1e6 * M.MX / sq_th - 1, 110)
        if ch == 'Xs_Xs':   # resolve the window where the u-pole is inside the angular range (3.606 .. 5.0 m_X)
            extra = np.concatenate([np.linspace(3.55, 3.66, 12), np.linspace(4.9, 5.1, 12)]) * M.MX / sq_th - 1
            deltas = np.sort(np.concatenate([deltas, extra[extra > 0]]))
        rows = []
        for d in deltas:
            sq = sq_th * (1 + d); s = sq * sq
            r = dict(sqrt_s_over_mX=sq / M.MX, s_GeV2=s)
            if ch == 'XXb_XXb':
                s0, s1, s2 = X.sigma_parts(s, P)
                t0, t1, t2 = X.sigma_parts(s, P, transport=True)
                Dr = (s - P.ms2) / ((s - P.ms2)**2 + P.eps_s**2)
                r.update(sigma0=s0, sigma1=s1, sigma2=s2, sigma_nores=s0, sigma_sub_pointwise=s0 + s1 * Dr + s2 * ((s - P.ms2)**2 - P.eps_s**2) / ((s - P.ms2)**2 + P.eps_s**2)**2,
                         sigma_full_pointwise=s0 + s1 * Dr + s2 / ((s - P.ms2)**2 + P.eps_s**2),
                         transport0=t0, transport1=t1, transport2=t2)
            else:
                r['sigma_sub'] = X.sigma(ch, s, P, mode='sub')
                r['sigma_nores'] = X.sigma(ch, s, P, mode='nores')
                if ch == 'Xs_Xs':
                    r['sigma_res'] = X.sigma(ch, s, P, mode='res'); r['sigma_int'] = X.sigma(ch, s, P, mode='int')
                if ch in ('XX_XX', 'Xs_Xs', 'ss_ss'):
                    r['transport_sub'] = X.sigma(ch, s, P, mode='sub', transport=True)
                    r['transport_nores'] = X.sigma(ch, s, P, mode='nores', transport=True)
            rows.append(r)
        # high-energy behaviour: log-log slope over the last decade and the contact constant
        key = 'sigma_nores'
        hi = [(x['sqrt_s_over_mX'], x[key]) for x in rows if x['sqrt_s_over_mX'] > 1e5]
        slope = np.polyfit(np.log([h[0] for h in hi]), np.log([abs(h[1]) for h in hi]), 1)[0]
        Mc = {'XXb_XXb': 4 * P.lX, 'XX_XX': 4 * P.lX, 'Xs_Xs': 2 * P.g2, 'ss_ss': P.lam_s4, 'ss_XXb': 2 * P.g2, 'XXb_ss': 2 * P.g2}[ch]
        out[ch] = dict(masses=m, Sin=Sin, Sout=Sout, sqrt_s_th_over_mX=sq_th / M.MX, rows=rows,
                       high_energy=dict(loglog_slope_in_sqrt_s=float(slope), contact_limit='sigma -> |M_c|^2/(16 pi s S_out)',
                                        M_contact=Mc, sigma_times_s_contact=Mc * Mc / (16 * np.pi * Sout),
                                        sigma_times_s_table_at_1e6=hi[-1][1] * (hi[-1][0] * M.MX)**2))
    return out


def thermal_tables():
    out = {}
    for ch in CH:
        m, Sin, Sout = X.masses_of(ch, P)
        rows = []
        t0 = time.time()
        for t in TGRID:
            T = t * M.MX
            n1, n2 = X.n_eq(m[0], T), X.n_eq(m[1], T)
            r = dict(T_over_mX=t, T_GeV=T, n1_eq=n1, n2_eq=n2)
            for mode in ('sub', 'nores'):
                g = X.thermal_rate(ch, T, P, mode=mode)
                r[f'gamma_{mode}'] = g; r[f'sigmav_{mode}'] = g * Sin / (n1 * n2)
            if ch in ('XX_XX', 'Xs_Xs', 'ss_ss', 'XXb_XXb'):
                for mode in ('sub', 'nores'):
                    g = X.thermal_rate(ch, T, P, mode=mode, transport=True)
                    r[f'gamma_transport_{mode}'] = g; r[f'sigmav_transport_{mode}'] = g * Sin / (n1 * n2)
            if ch == 'XXb_XXb':
                r['gamma_ris'] = X.thermal_rate(ch, T, P, mode='ris')
                r['gamma_full'] = X.thermal_rate(ch, T, P, mode='full')
                r['gamma_decay_MB'] = X.decay_rate_density(T, P)
            rows.append(r)
        out[ch] = dict(Sin=Sin, Sout=Sout, rows=rows, seconds=time.time() - t0)
        print(ch, 'thermal done', time.time() - t0, flush=True)
    return out


def decay_tables():
    rows = []
    for t in TGRID:
        T = t * M.MX
        ns = X.n_eq(P.ms, T)
        g = X.decay_rate_density(T, P)
        rows.append(dict(T_over_mX=t, T_GeV=T, n_s_eq=ns, gamma_decay_MB=g, Gamma_eff_over_Gamma_vac=g / (ns * P.Gamma_s),
                         gamma_inverse_decay_MB_mu0=g))
    return dict(Gamma_s_vac_GeV=P.Gamma_s, formula='gamma_dec = n_s Gamma_s <m_s/E> = m_s^2 T Gamma_s K1(m_s/T)/(2 pi^2); '
                'inverse decay = decay at mu = 0 (MB); with Bose factors and mu see thermal.width_s_medium and mc_checks', rows=rows)


if __name__ == '__main__':
    t0 = time.time()
    res = dict(params=dict(mX=P.mX, ms=P.ms, lam3=P.lam3, lam_s3=P.lam_s3, lam_s4=P.lam_s4, g2=P.g2, lX=P.lX, Gamma_s=P.Gamma_s,
                           eps_s=P.eps_s, eps_X=P.eps_X),
               sigma=sigma_tables(), thermal=thermal_tables(), decay=decay_tables(), seconds=time.time() - t0)
    json.dump(res, open(sys.argv[1], 'w'), indent=1, default=float)
    for ch in CH:
        print(ch, res['sigma'][ch]['high_energy'])
    print('done', time.time() - t0)
