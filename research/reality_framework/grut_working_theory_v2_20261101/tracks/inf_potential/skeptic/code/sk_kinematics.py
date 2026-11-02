#!/usr/bin/env python3
"""inf_potential SKEPTIC (Claude): (1) kinematic test of the engine width Gamma_I = 0.5 H_i for the brief's singlet-fermion
channel y phi_I N N (Gamma = N_f y^2 m_I/8 pi): field-dependent mass y|phi_I| vs m_I/2 at the first turning point, and the
matter-era length until the condition holds by redshift alone (A ∝ a^-3/2) [D/E];  (2) the HT6t 'clean window' against
the H_i systematics of WORKING_THEORY section 5 (portals vanishing at M_P: H_i = 1.08e11; history -16%/+21%) from the
track's scan_H.json [C interpolation];  (3) the yield-calibration shift implied by the benchmarks' own smoke release
(Y ∝ N_X H_i^{3/2}/K^{3/2} at fixed dimensionless history) [E].  Output ../results/sk_kinematics.json."""
import json, math
from pathlib import Path
import numpy as np
TRK = Path(__file__).resolve().parents[2]; OUT = Path(__file__).resolve().parent.parent / 'results'
MP = 2.435e18
out = {}
eng = json.load(open(TRK / 'results' / 'engine_test.json'))['smoke']
for n in ('HT6t', 'T1', 'E1'):
    d = json.load(open(TRK / 'results' / ('bench_%s.json' % n)))['by_Hi'][0]
    mI = d['m_I_over_Hi']; Hi = d['H_i_GeV']; A1 = abs(d['osc_gamma0']['turning'][0]['phiI_turn']) * MP / Hi   # in H_i
    r = {}
    for Nf in (1, 3):
        y = math.sqrt(8 * math.pi * 0.5 / (Nf * mI))
        ratio = y * A1 / (mI / 2)
        akin = ratio ** (2 / 3)
        K = d['Hend_over_Hi'] * math.sqrt(akin)          # instant decay at a_kin after a pure matter era from the end
        r['Nf%d' % Nf] = dict(y=y, y_phi_over_half_mI=ratio, a_kin=akin, Nmatter=math.log(akin), dNstar=-math.log(akin) / 4,
                              K_if_decay_at_akin=K, Kratio_vs_track=K / d['K_end_gamma05'],
                              t_kin_Hi=(2 / 3) / d['Hend_over_Hi'] * akin ** 1.5)
    NXr = abs(eng[n]['NX_endnorm']) / 8045.485140302082
    Yr = NXr * (1.98581 / d['K_end_gamma05']) ** 1.5
    out[n] = dict(m_I_over_Hi=mI, A1_over_Hi=A1, kinematics=r, NX_ratio_vs_M1r=NXr, yield_ratio_vs_M1r=Yr,
                  H_i_shift_factor=Yr ** (-2 / 3), H_i_recal_GeV=Hi * Yr ** (-2 / 3))
    print(n, json.dumps(out[n], default=lambda x: float("%.4g" % x)))
sc = json.load(open(TRK / 'results' / 'scan_H.json'))['rows']
print([ (k) for k in sc[0].keys()][:20])
json.dump(out, open(OUT / 'sk_kinematics.json', 'w'), indent=1)
