#!/usr/bin/env python3
"""rp_critic: arithmetic cross-checks of numbers shared between the three tracks and Codex e160356 (own evaluation).
usage: rc_arith.py CODEX_REHEATING_RESULTS.json OUT.json"""
import sys, json, math
M_R = 2.435e18; H_I = 1.3334961939833566e11; G = 106.75; Y_Q = 4.8294009077e-9
out = {}
# 1. Dilution bookkeeping D = (r a^4/3)^{3/4} from Codex's recorded final states, vs yield/(|NX|/2788.18)
d = json.load(open(sys.argv[1]))
rows = {}
for r in d['rows']:
    a = r['final_a']; rr = r['final_state'][14]
    D = (rr * a ** 4 / 3) ** 0.75
    ratio = abs(r['late_NX']) / 2788.18
    # yield directly in GeV units from Codex's late NX and final radiation density
    rho = rr * M_R ** 2 * H_I ** 2; T = (30 * rho / (math.pi ** 2 * G)) ** 0.25; s = 2 * math.pi ** 2 * G * T ** 3 / 45
    Y = abs(r['late_NX']) * H_I ** 3 / (a ** 3 * s)
    rows[r['tag']] = dict(D=D, charge_ratio=ratio, ratio_over_D=ratio / D, codex_yield=abs(r['yield_over_YQ']), own_Y_over_YQ=Y / Y_Q,
                          T_RH_est_GeV=0.54 * math.sqrt(r['Gamma_over_H'] * H_I * M_R) if not r['rd_control'] else None,
                          Tvis_final_GeV=r['Tvis_final_GeV'])
out['dilution_rows'] = rows
# 2. T_RH prefactor and values on the yield-restoring line
pref = (90 / (math.pi ** 2 * G)) ** 0.25
out['T_RH'] = dict(prefactor=pref, **{f'Gamma_{g}H_i': pref * math.sqrt(g * H_I * M_R) for g in [0.03, 0.1, 0.3, 0.5, 1.0, 3.0, 10.0]})
# 3. Higgs-portal source arithmetic
b = 5e-11; vH = 246.0; mh = 125.0; ms = 0.2264; mX = 0.1 * 0.924089637
beta = math.sqrt(1 - 4 * ms ** 2 / mh ** 2)
Ghss = b ** 2 * vH ** 2 * beta / (8 * math.pi * mh)
Hmh = math.sqrt(math.pi ** 2 * G / 90) * mh ** 2 / M_R
Ypair = 2 * 135 * Ghss / (8 * math.pi ** 3 * G * Hmh)
frest = 2 * mX * Ypair / 4.4e-10
# conformal energy of the prompt daughters: J = <E_h>/T * Y_h... each h gives 2 pairs; prompt bound J_h = 5 Y_hss (Codex/rp_portals convention)
Yhss = Ypair / 2; J_h = 5 * Yhss
EXN = 38386.250633; F_R = 4 * (M_R / H_I) ** 2 / EXN
out['higgs'] = dict(Gamma_hss=Ghss, H_mh=Hmh, Y_pair=Ypair, f_rest=frest, J_h_prompt=J_h, R_h=F_R * J_h, F_R=F_R,
                    H_at_mX_equals_H_T_GeV=pref * math.sqrt(mX * M_R), H_at_T60GeV=math.sqrt(math.pi ** 2 * G / 90) * 60 ** 2 / M_R,
                    H_at_mA=math.sqrt(math.pi ** 2 * G / 90) * (3.2003908656e6) ** 2 / M_R)
# 4. widths
gY = 0.37; mA = 3.2003908656e6; eX = 1e-7; eps = 1e-7
GSM = 41 * eps ** 2 * gY ** 2 * mA / (96 * math.pi); GX = eX ** 2 * mA / (48 * math.pi) * (1 - 4 * mX ** 2 / mA ** 2) ** 1.5
out['widths'] = dict(Gamma_SM=GSM, Gamma_X=GX, B=GX / (GX + GSM), sumY2_weyl=3 * (6 / 36 + 3 * 4 / 9 + 3 / 9 + 2 / 4 + 1))
# 5. toy-model exact constants
out['boundary_exact'] = dict(T_amp=12 * 0.04 * 2 * math.sqrt(2), E_amp=2 * 0.04 * 2 * math.sqrt(2))
json.dump(out, open(sys.argv[2], 'w'), indent=1)
for k, v in rows.items():
    print(f"{k:14s} D={v['D']:.3f} ratio={v['charge_ratio']:.4f} ratio/D={v['ratio_over_D']:.5f} codex={v['codex_yield']:.5f} ownY={v['own_Y_over_YQ']:.5f} T_RH={v['T_RH_est_GeV']}")
print('T_RH', out['T_RH']); print('higgs', out['higgs']); print('widths', out['widths'])
