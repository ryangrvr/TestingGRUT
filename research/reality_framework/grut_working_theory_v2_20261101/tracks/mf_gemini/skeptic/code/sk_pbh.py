"""mf_gemini SKEPTIC: PBH bookkeeping on this skeptic's waterfall runs [E].

  * f_PBH from entropy conservation: rho_PBH,0 = beta (3/4) T_f s_0 (radiation at formation, g_* = g_s = 106.75),
    f = rho_PBH,0/rho_DM,0 (Omega_c h^2 = 0.120).  The track's (T_f/T_eq)(g_f/g_eq)^{1/3} form omits the factor
    (g_f/g_eq)(g_sf/g_seq)^{-4/3} = 0.389 from the non-a^{-4} radiation density (it overstates f by 2.6).
  * Formation epoch: instant radiation (track) vs the slow reheating the track itself finds (perturbative cap
    Gamma = 1.8e-8 H_end, matter era N_MD = (2/3) ln(H_end/Gamma)); the peak modes (N_b >= 16) re-enter after
    reheating, at H_entry = H_end exp(-2 N_b - N_MD/2).
  * Patches with no completed waterfall (concave V_SR): beta >= P_no_waterfall (they keep inflating).
Run: python3 -B sk_pbh.py -> ../results/sk_pbh.json
"""
import json
import os
import numpy as np

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
MP = 2.435e18
GEV_G = 1.78266192e-24
S0 = 2891.2 * (1.9733e-14)**3            # GeV^3
RHO_DM = 0.120 * 1.05371e-5 * (1.9733e-14)**3   # GeV^4
G = 106.75
GAMMA = 0.2


def f_pbh(beta, H_f):
    T_f = (90 * H_f**2 * MP**2 / (np.pi**2 * G))**0.25
    return beta * 0.75 * T_f * S0 / RHO_DM, T_f


def main():
    with open(os.path.join(RES, 'sk_waterfall.json')) as fh:
        wf = json.load(fh)
    H0 = wf['H0_GeV']
    rows = []
    for c in wf['cases']:
        r = c['results']['q1']
        Hend = r['H_end'] * H0
        Nb = r['N_before_end_peak']
        # instant radiation
        Hf = Hend * np.exp(-2 * Nb)
        M = GAMMA * 4 * np.pi * MP**2 / Hf * GEV_G
        fr = {}
        for z in ['0.5', '0.67', '1.0']:
            fr[z] = f_pbh(r['P_tail'][z], Hf)[0]
        # slow reheating (perturbative cap)
        Gam = 1.8e-8 * Hend
        NMD = (2.0 / 3.0) * np.log(Hend / Gam)
        HfM = Hend * np.exp(-2 * Nb - NMD / 2)
        MM = GAMMA * 4 * np.pi * MP**2 / HfM * GEV_G
        fM = f_pbh(r['P_tail']['0.67'], HfM)[0]
        fnw = f_pbh(r['P_no_waterfall'], Hf)[0] if r['P_no_waterfall'] > 0 else 0.0
        rows.append(dict(case=c['name'], N_b=Nb, M_PBH_g_instantRD=M, f_PBH_instantRD=fr,
                         track_formula_factor=1 / 0.389, N_MD_slow=NMD, M_PBH_g_slowRH=MM, f_PBH_slowRH_z067=fM,
                         P_no_waterfall=r['P_no_waterfall'], f_PBH_from_no_waterfall_patches_lower=fnw,
                         Hawking_survives=bool(min(M, MM) > 5e14)))
    out = dict(rows=rows, ranges=dict(
        M_PBH_g=[min(x['M_PBH_g_instantRD'] for x in rows), max(x['M_PBH_g_instantRD'] for x in rows)],
        f_instant_z067=[min(x['f_PBH_instantRD']['0.67'] for x in rows), max(x['f_PBH_instantRD']['0.67'] for x in rows)],
        f_instant_z1=[min(x['f_PBH_instantRD']['1.0'] for x in rows), max(x['f_PBH_instantRD']['1.0'] for x in rows)],
        f_slowRH=[min(x['f_PBH_slowRH_z067'] for x in rows), max(x['f_PBH_slowRH_z067'] for x in rows)],
        M_slowRH_g=[min(x['M_PBH_g_slowRH'] for x in rows), max(x['M_PBH_g_slowRH'] for x in rows)]))
    with open(os.path.join(RES, 'sk_pbh.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    for x in rows:
        print(x['case'], '%.2e %.2e %.2e | slow %.2e %.2e | nowf %.2e' % (x['M_PBH_g_instantRD'], x['f_PBH_instantRD']['0.67'],
              x['f_PBH_instantRD']['1.0'], x['M_PBH_g_slowRH'], x['f_PBH_slowRH_z067'], x['f_PBH_from_no_waterfall_patches_lower']))
    print(out['ranges'])


if __name__ == '__main__':
    main()
