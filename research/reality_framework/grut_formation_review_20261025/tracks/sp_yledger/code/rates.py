#!/usr/bin/env python3
"""sp_yledger (Claude, own code): order-of-magnitude rate table for the hidden sector after the release
(item b).  All rates are tree-level estimates [E] with the stated couplings; H = H_i/a^2 (radiation,
instant visible reheating at H_i); T_vis = T_RI/a (g* constant, adequate above 200 GeV).

Rates per particle, GeV:
  Y_H -> Y_L + A_L (Goldstone equivalence), sigma -> A_L A_L, Y_L -> 3X, Y -> XX on the condensate,
  A -> X Xbar (time-dilated while relativistic), hard-daughter energy sharing on the soft X population
  (2->2 through lambda_X, sigma = lambda^2/(pi s), s ~ 2 E w), thermal 2->2 and 2->4 in a thermalised hidden
  bath (lambda_X only, and with the light-target X-phi coupling g_X^2 = 0.01 once phi is populated).
Usage: python3 rates.py <point_ledger.json> <out.json>
"""
import json, sys
import numpy as np

H_I = 1.1211167304108865e13
MRED, GSTAR = 2.435e18, 106.75
LX, GX2 = 1e-3, 0.01
ZETA3 = 1.2020569


def main():
    L = json.load(open(sys.argv[1]))
    W = L['widths_GeV']; M = L['masses_GeV']; Nn = L['numbers_late']; En = L['energies_over_EX']
    EX = L['windows']['2000-3000']['EX']; NX = abs(Nn['NX'])
    Ncirc = 2 * np.sqrt(2 * LX) * (EX / (3 * LX))**0.75
    TRI = (90 * MRED**2 * H_I**2 / (np.pi**2 * GSTAR))**0.25
    Xconf = (EX / LX)**0.25
    a_dec = L['a_mean_decay_radiation']['YL']; a_Adec = L['a_mean_decay_radiation']['A']
    RY = En['Y_deposited']; RA = En['A_deposited_range'][0]
    rows = []
    for la in np.arange(1.5, 13.01, 0.5):
        a = 10**la
        Hh = H_I / a**2
        Tv = TRI / a
        # hidden energy (conformal) deposited so far, thermal temperature if thermalised (g_h = 2: X, Xbar)
        Eh = EX * (1 + (RY if a > a_dec else 0) + (RA if a > a_Adec else 0))
        Th = (30 * Eh * H_I**4 / a**4 / (np.pi**2 * 2))**0.25
        # soft X population (condensate or its fragments): number ~ Ncirc, quantum energy w
        n_soft = Ncirc * H_I**3 / a**3
        w = EX / Ncirc * H_I / a
        # daughters (exist after a_dec): energy m_L/3 redshifted
        Ed = M['m_YL'] / 3 * min(1.0, a_dec / a)
        G_share = LX**2 * n_soft / (2 * np.pi * Ed * w)
        # A energy: emitted near a ~ 1e2 with E ~ Delta_M; time dilation factor m_A/E_A
        EA = max(M['m_A'], M['Delta_M'] * 1e2 / a)
        G_A = W['A_to_XXbar'] * M['m_A'] / EA
        # thermal bath rates (if thermalised): n = 2 zeta3/pi^2 T^3, <s> ~ 18 T^2
        n_th = 2 * ZETA3 / np.pi**2 * Th**3
        G22_l = n_th * LX**2 / (np.pi * 18 * Th**2)
        G22_g = n_th * GX2**2 / (np.pi * 18 * Th**2)
        G24_l = G22_l * LX / (16 * np.pi**2)
        G24_g = G22_g * GX2 / (16 * np.pi**2)
        rows.append(dict(a=a, T_vis_GeV=Tv, H_GeV=Hh, T_h_if_thermal_GeV=Th, T_h_over_T_vis=Th / Tv,
                         YH_to_YL_A=W['YH_to_YL_A'] / Hh, sigma_to_AA=W['sigma_to_AA'] / Hh,
                         YL_to_3X=W['YL_3X'] / Hh, Y_to_XX_condensate=W['Y_to_XX_condensate_over_H'],
                         A_to_XXbar=G_A / Hh, daughter_energy_sharing=G_share / Hh if a > a_dec else None,
                         thermal_22_lambda=G22_l / Hh, thermal_22_gX=G22_g / Hh,
                         thermal_24_lambda=G24_l / Hh, thermal_24_gX=G24_g / Hh))

    def first(key, thr=1.0, after=0):
        for r in rows:
            v = r[key]
            if v is not None and r['a'] > after and v >= thr:
                return dict(a=r['a'], T_vis_GeV=r['T_vis_GeV'])
        return None
    summary = dict(
        a_dec_YL=a_dec, a_dec_A=a_Adec, T_vis_A_decay=TRI / a_Adec,
        daughter_sharing_at_a_dec=next(r['daughter_energy_sharing'] for r in rows if r['a'] > a_dec),
        daughter_sharing_reaches_H=first('daughter_energy_sharing', after=a_dec),
        thermal_24_lambda_reaches_H=first('thermal_24_lambda', after=a_dec),
        thermal_24_gX_reaches_H=first('thermal_24_gX', after=a_dec),
        thermal_22_lambda_reaches_H=first('thermal_22_lambda'),
        Y_to_XX_condensate_fraction_before_3body=W['Y_to_XX_condensate_over_H'] * np.log(a_dec / 100),
        quartic_self_resonance=('Floquet exponent ~0.036 per unit of the Lame variable x (lambda phi^4, '
                                'Greene-Kofman-Linde-Starobinsky 1997 [R]); x = sqrt(2 lambda)|X|a: %.2f per unit a, so '
                                'growth ~exp(0.036 * %.2f a): e^30 by a ~ %.0f [E]. Lattice not run.'
                                % (np.sqrt(2 * LX) * Xconf, np.sqrt(2 * LX) * Xconf, 30 / (0.036 * np.sqrt(2 * LX) * Xconf))))
    out = dict(rows=rows, summary=summary, labels='[E] tree-level order-of-magnitude estimates; [D] for widths')
    json.dump(out, open(sys.argv[2], 'w'), indent=1, default=float)
    print(json.dumps(summary, indent=1, default=float))
    for r in rows:
        print(' '.join(f'{k}={v:.3g}' for k, v in r.items() if v is not None))


if __name__ == '__main__':
    main()
