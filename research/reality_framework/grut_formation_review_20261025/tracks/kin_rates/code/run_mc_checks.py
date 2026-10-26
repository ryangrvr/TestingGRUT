#!/usr/bin/env python3
"""kin_rates: Monte Carlo checks with Bose statistics (own code; engine thermal_mc.py).  usage: python3 -B run_mc_checks.py <out.json>
1. Bose/MB ratios of the 2->2 thermal rates at mu = 0 and the MB MC against the quadrature of xsec22.
2. Detailed balance with Bose factors and chemical potentials: s <-> X Xb (decay by MC, inverse decay by exact 2-D
   quadrature; mu_s = mu_+ + mu_-) and s s <-> X Xb (both directions by MC; mu_s = mu_+ + mu_- = 0 is NOT required
   for this pair-number-changing... note: ss -> XXb changes the pair inventory, balance needs 2 mu_s = mu_+ + mu_-).
3. Vanishing internal energy moment of every elastic kernel at a common temperature, and the X-s energy exchange
   coefficient at different species temperatures (the closure input).
4. Hard X on a thermal bath: collision rate and energy-loss rate per channel for E/m_X up to 7.7e4.
5. Point-kinematics gain/loss identity for every forward/reverse pair of the network (>= 240 cases).
"""
import json
import sys
import time
import numpy as np
import model as M
import xsec22 as X
import thermal_mc as TM
import kernels_n2 as K
import phasespace as PS


def TM_mass(b, P):
    return P.ms if b == 's' else P.mX


def point_identity(rng, per_pair=20):
    """prod_in f prod_out (1+f) = prod_out f prod_in (1+f) on random on-shell, energy-conserving configurations
    (boosted n-body + 2-body phase-space points at the same total momentum) with Bose f's at a common T and
    chemical potentials obeying the pair's constraint (mu_s = mu_+ + mu_- and, for reactions that change the
    pair inventory, mu_+ + mu_- = 0 as well)."""
    pairs = {ch: spec for ch, spec in K.CHANNELS_N2.items() if not ch.startswith('ctrl')}
    pairs.update({'decay_s_XXb': dict(legs=['s', ('X', 'out'), ('Xb', 'out')]),
                  'ss_XXb': dict(legs=['s', 's', ('X', 'out'), ('Xb', 'out')]),
                  'XXb_XXb': dict(legs=['X', 'Xb', ('X', 'out'), ('Xb', 'out')]),
                  'Xs_Xs': dict(legs=['X', 's', ('X', 'out'), ('s', 'out')])})
    rows = []
    worst = 0.0
    for ch, spec in pairs.items():
        legs = spec['legs']
        sp_in = [K.species(l) for l in legs if not K.is_out(l)]; sp_out = [K.species(l) for l in legs if K.is_out(l)]
        m_in = [K.MASS[s] for s in sp_in]; m_out = [K.MASS[s] for s in sp_out]
        # pair-inventory change: P = n_Xb + n_s
        dP = (sp_out.count('Xb') + sp_out.count('s')) - (sp_in.count('Xb') + sp_in.count('s'))
        res = []
        for _ in range(per_pair):
            T = M.MX * 10**rng.uniform(-0.5, 1.2)
            mup = M.MX * rng.uniform(-0.6, 0.6)
            mum = -mup if dP != 0 else M.MX * rng.uniform(-0.6, 0.6)
            mu = {'X': mup, 'Xb': mum, 's': mup + mum}
            sq = max(sum(m_in), sum(m_out)) * (1 + rng.uniform(0.05, 3.0))
            if len(m_in) == 1:
                pin = np.array([[[sq, 0, 0, 0.0]]])
            else:
                pin, _ = PS.generate(sq, m_in, 1, rng)
            if len(m_out) == 1:
                pout = np.array([[[sq, 0, 0, 0.0]]])
            else:
                pout, _ = PS.generate(sq, m_out, 1, rng)
            # common random boost
            b = rng.uniform(0, 0.9); c = rng.uniform(-1, 1); ph = rng.uniform(0, 2 * np.pi); sn = np.sqrt(1 - c * c)
            g = 1 / np.sqrt(1 - b * b)
            Pm = np.array([[g * sq, g * sq * b * sn * np.cos(ph), g * sq * b * sn * np.sin(ph), g * sq * b * c]])
            Ein = [PS.boost(pin[:, i, :], Pm)[0, 0] for i in range(len(m_in))]
            Eout = [PS.boost(pout[:, i, :], Pm)[0, 0] for i in range(len(m_out))]
            fi = [1 / np.expm1((E - mu[s]) / T) for E, s in zip(Ein, sp_in)]
            fo = [1 / np.expm1((E - mu[s]) / T) for E, s in zip(Eout, sp_out)]
            fwd = np.prod(fi) * np.prod([1 + f for f in fo]); bwd = np.prod(fo) * np.prod([1 + f for f in fi])
            r = abs(fwd - bwd) / (fwd + bwd)
            res.append(r)
            # energy residual of the configuration
            eres = abs(sum(Ein) - sum(Eout)) / sum(Ein)
            worst = max(worst, eres)
        rows.append(dict(pair=ch, legs=[list(l) if isinstance(l, tuple) else l for l in legs], cases=per_pair,
                         pair_inventory_change=dP, mu_constraint='mu_s = mu_+ + mu_-' + ('; mu_+ + mu_- = 0' if dP else ''),
                         max_residual=float(max(res)), mean_residual=float(np.mean(res))))
    return dict(total_cases=sum(r['cases'] for r in rows), max_residual=float(max(r['max_residual'] for r in rows)),
                max_energy_residual=float(worst), pairs=rows)


def main(outpath):
    rng = np.random.default_rng(2026)
    P = X.Params(Gamma_X_reg=0.05)
    out = dict(scope='Lab-frame MC thermal integrals with Bose statistics (own code); MB references from the xsec22 quadrature. '
                     'The X-line regulator of the Xs_Xs u-pole finite part is Gamma_X_reg = 0.05 m_X here (MC variance); '
                     'the quadrature tables use 1e-7 m_X.')
    N = 200000
    t0 = time.time()
    rows = []
    for t in (0.3, 1.0, 3.0, 7.344, 12.0):
        T = t * M.MX
        row = dict(T_over_mX=t)
        for ch, mode in (('XXb_XXb', 'nores'), ('XX_XX', 'sub'), ('Xs_Xs', 'sub'), ('Xs_Xs', 'nores'), ('ss_ss', 'sub'), ('ss_XXb', 'sub'), ('XXb_ss', 'sub')):
            mb = TM.rate_22(ch, T, P, N, rng, 'MB', mode=mode)
            be = TM.rate_22(ch, T, P, N, rng, 'BE', mode=mode)
            ref = X.thermal_rate(ch, T, P, mode=mode)
            row[ch if mode != 'nores' or ch == 'XXb_XXb' else ch + '_nores'] = dict(gamma_MB_mc=mb['gamma'], gamma_MB_mc_err=mb['err'], gamma_MB_quad=ref,
                           mc_over_quad=mb['gamma'] / ref if ref else None,
                           gamma_BE_mc=be['gamma'], gamma_BE_mc_err=be['err'],
                           bose_over_MB=be['gamma'] / mb['gamma'] if mb['gamma'] else None, mode=mode)
        # ss <-> XXb with Bose factors and chemical potentials obeying 2 mu_s = mu_+ + mu_- (inventory-changing)
        for (mp_, mm_) in ((0.0, 0.0), (0.3 * M.MX, 0.3 * M.MX), (0.5 * M.MX, -0.2 * M.MX)):
            mu = dict(X=mp_, Xb=mm_, s=0.5 * (mp_ + mm_))
            f = TM.rate_22('ss_XXb', T, P, N, rng, 'BE', mu=mu); r = TM.rate_22('XXb_ss', T, P, N, rng, 'BE', mu=mu)
            row[f'db_ss_XXb_BE_mup{mp_/M.MX:.1f}_mum{mm_/M.MX:.1f}'] = dict(mu_over_mX=dict(X=mp_ / M.MX, Xb=mm_ / M.MX, s=mu['s'] / M.MX),
                                                                           forward=f['gamma'], forward_err=f['err'], reverse=r['gamma'], reverse_err=r['err'],
                                                                           ratio=f['gamma'] / r['gamma'], sigma_dev=(f['gamma'] - r['gamma']) / np.hypot(f['err'], r['err']))
        # s <-> XXb: decay (MC) vs inverse decay (quadrature), mu_s = mu_+ + mu_-
        for (mp_, mm_) in ((0.0, 0.0), (0.4 * M.MX, 0.1 * M.MX), (0.6 * M.MX, -0.3 * M.MX)):
            ms_ = mp_ + mm_
            d, de = TM.decay_bose(T, P, N, rng, ms_, mp_, mm_)
            inv = TM.inverse_decay_bose(T, P, ms_, mp_, mm_)
            row[f'db_s_XXb_BE_mup{mp_/M.MX:.1f}_mum{mm_/M.MX:.1f}'] = dict(mu_over_mX=dict(X=mp_ / M.MX, Xb=mm_ / M.MX, s=ms_ / M.MX),
                                                                          decay_mc=d, decay_err=de, inverse_quad=inv, ratio=d / inv, sigma_dev=(d - inv) / de,
                                                                          decay_MB_analytic=X.decay_rate_density(T, P) if ms_ == 0 else None)
        rows.append(row)
        print('T', t, 'done', time.time() - t0, flush=True)
    out['thermal_rates_and_balance'] = rows
    # 3. elastic energy moments
    el = []
    for t in (0.3, 1.0, 3.0, 12.0):
        T = t * M.MX
        for ch in ('XXb_XXb', 'XX_XX', 'Xs_Xs', 'ss_ss'):
            mode = 'nores' if ch in ('XXb_XXb', 'Xs_Xs') else 'sub'
            r = TM.rate_22(ch, T, P, N, rng, 'BE', mode=mode, energy=True)
            el.append(dict(T_over_mX=t, channel=ch, gamma=r['gamma'], Q1=r['Q1'], Q2=r['Q2'],
                           sum_over_abs=(r['Q1'] + r['Q2']) / (abs(r['Q1']) + abs(r['Q2']) + 1e-300), Q1_err=r['Q1_err']))
    out['elastic_energy_moment'] = el
    ex = []
    for tX in (0.5, 1.0, 3.0, 7.344):
        for ts in (0.5, 1.0, 3.0, 7.344):
            TX, Tsx = tX * M.MX, ts * M.MX
            r = TM.rate_22('Xs_Xs', TX, P, N, rng, 'MB', Ts={'X': TX, 's': Tsx}, energy=True, mode='nores')
            nX = X.n_eq(P.mX, TX); ns = X.n_eq(P.ms, Tsx)
            ex.append(dict(T_X_over_m=tX, T_s_over_m=ts, mode='nores', gamma=r['gamma'], Q_s=r['Q2'], Q_s_err=r['Q2_err'], Q_X=r['Q1'],
                           sum_Q_over_abs=(r['Q1'] + r['Q2']) / (abs(r['Q1']) + abs(r['Q2']) + 1e-300),
                           Q_s_per_nXns_GeV_minus1=r['Q2'] / (nX * ns), gamma_per_nXns=r['gamma'] / (nX * ns)))
    out['energy_exchange_Xs_MB'] = ex
    print('energy exchange done', time.time() - t0, flush=True)
    # 4. hard particle on a thermal bath
    hp = []
    for t in (0.3, 1.0, 3.0, 7.344):
        T = t * M.MX
        nbath = {b: TM.n_bose(TM_mass(b, P), T) for b in ('X', 'Xb', 's')}
        for Eom in (3.0, 10.0, 100.0, 1e3, 1e4, 7.72e4):
            r = TM.hard_particle(Eom * M.MX, T, P, 100000, rng)
            for b in r:
                r[b]['n_bath_BE_GeV3'] = nbath[b]
                r[b]['rate_per_density_GeV_minus2'] = r[b]['rate'] / nbath[b]
                r[b]['dEdt_per_density_GeV_minus1'] = r[b]['dEdt'] / nbath[b]
            hp.append(dict(T_over_m=t, E_over_m=Eom, channels=r))
    out['hard_particle'] = hp
    print('hard particle done', time.time() - t0, flush=True)
    out['point_kinematics'] = point_identity(rng, per_pair=20)
    out['seconds'] = time.time() - t0
    json.dump(out, open(outpath, 'w'), indent=1, default=float)
    print(json.dumps(out['point_kinematics'], indent=None, default=float)[:500])


if __name__ == '__main__':
    main(sys.argv[1])
