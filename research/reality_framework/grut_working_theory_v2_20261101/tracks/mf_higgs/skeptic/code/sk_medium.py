#!/usr/bin/env python3
"""mf_higgs SKEPTIC: how much of the medium mass is needed?  D at M_t = 173.45 / 174.05 / 174.6 with
m_p^2 = f x (N gas + SM hard gas) of sk_plasma's model, for two energy laws: E fixed at injection (conservative,
m^2 ~ a^-4 x fractions) and E redshifting as 1/a (relativistic quanta; N floored at M).  Own modes (sk_modes runs).
[C on E]  Output: results/sk_medium.json
"""
import sys
sys.dont_write_bytecode = True
import math
import numpy as np
import sk_common as C
import sk_modes as M

MNU, VEW = 0.05e-9, 174.1


def curves(o, hp, redshift):
    d, h = C.load_hist(hp)
    Hi = d['H_i_GeV']; mP = C.MPL / Hi
    t = h['t']; rho = 3 * mP ** 2 * h['H'] ** 2
    ie = int(np.argmin(np.abs(t - d['waterfall']['t_eps1'])))
    E0 = 0.75 * math.sqrt(2 * h['f_chi'][ie] * rho[ie]) / 4000.0
    E_N = np.maximum(E0 * (h['a'][ie] / h['a']) if redshift else E0 + 0 * t, 120.0)
    E_SM = E_N / 2 if not redshift else np.maximum(E0 / 2 * (h['a'][ie] / h['a']), 1.0)
    fSM = np.nan_to_num(h['f_SM']); fN = np.maximum(np.nan_to_num(h['f_rad']) - fSM, 0)
    m2 = MNU * 120 * Hi / VEW ** 2 * fN * rho / E_N ** 2 + 0.0228 * fSM * rho / E_SM ** 2
    return np.interp(o['t'], t, m2)


out = {}
for lab, rn, hp in (('EXP run, MC_TR fractions', 'EXP_base', M.HPATH['MC_TR']),
                    ('MC_P0 @ 2.4e11', 'MC_P0_Hi2.4e+11', M.HPATH['MC_P0']),
                    ('MC_P0 @ 1.3e11', 'MC_P0_Hi1.3e+11', M.HPATH['MC_P0'])):
    o = M.load(rn)
    for rs in (False, True):
        m2 = curves(o, hp, rs)
        row = {}
        for f in (1.0, 1e-2, 1e-3, 3e-4, 1e-4):
            row[f'{f:g}'] = {f'{Mt}': M.danger(o, Mt, mp2=f * m2)['D'] for Mt in (173.45, 174.05, 174.6)}
        # smallest factor (log grid) with D < 0 at all three M_t
        fs = np.logspace(0, -6, 25)
        fmin = None
        for f in fs:
            if all(not np.isfinite(M.danger(o, Mt, mp2=f * m2)['D']) or M.danger(o, Mt, mp2=f * m2)['D'] < 0 for Mt in (174.6,)):
                fmin = float(f)
            else:
                break
        key = f'{lab} | E ' + ('redshifting' if rs else 'fixed (conservative)')
        out[key] = dict(D=row, smallest_factor_safe_at_174p6=fmin,
                        m2_at={f'{x}': float(np.interp(x, o['t'], m2)) for x in (0.2, 0.3, 1.0, 3.0, 10.0)})
        print(key, 'fmin', fmin, {k: {a: (round(b, 1) if np.isfinite(b) else 'no roll') for a, b in v.items()} for k, v in row.items()},
              flush=True)
C.jdump(out, C.RES / 'sk_medium.json')
