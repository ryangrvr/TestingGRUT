#!/usr/bin/env python3
"""Assemble the structure-screen table and classification (Claude review, track `conversion`).

The classification is approximate and retrospective.  For each model it records:
  * dA (Murgia area criterion, my PKMAX = 120 h/Mpc pipeline) against my own calibrators;
  * the Delta chi^2 against the DESI DR1 compressed P1D constraint (Delta^2_*, n_*), relative
    to LCDM;
  * the mixed-WDM bound for the equivalent mass (Diamanti et al. 2017 Fig. 5; satellites:
    Tan et al. 2025).
Labels:
  allowed   dA < dA(5.7 keV pure) and Delta chi^2_DESI < 4 and f below the mixed bound
  marginal  dA between dA(5.7 keV) and dA(3.5 keV), or 4 <= Delta chi^2_DESI < 9, or f within
            a factor 1.5 of the mixed bound
  excluded (under these approximations)  dA > dA(3.5 keV) or Delta chi^2_DESI >= 9 or f above the bound
"""
import json
from pathlib import Path

RES = Path(__file__).resolve().parents[1] / 'results'
cal = {r['label']: r for r in json.loads((RES / 'class_screen_calibration.json').read_text())['rows']}
dA57 = cal['thermal_5.7keV_full']['dA']; dA35 = cal['thermal_3.5keV_full']['dA']
TNU_MEAN = 3.151 * 0.71611 * 2.7255 * 8.617333262e-5     # <p> of FD at T_nu, eV


def mixed_bound(m_eq):
    # Diamanti et al. 2017, 2-sigma fermion, per mass decade (read from Fig. 5 / text)
    if m_eq < 10: return 0.0085
    if m_eq < 100: return 0.04
    if m_eq < 1000: return 0.16
    if m_eq < 1e4: return 0.29
    return 0.43


def classify(r, f, u):
    m_eq = TNU_MEAN / u
    b = mixed_bound(m_eq)
    dchi = r['chi2_DESI_DR1'] - r['chi2_DESI_DR1_LCDM']
    flags = []
    if r['dA'] > dA35 or dchi >= 9 or f > b:
        lab = 'excluded (approx.)'
    elif r['dA'] > dA57 or dchi >= 4 or f > b / 1.5:
        lab = 'marginal'
    else:
        lab = 'allowed'
    return lab, m_eq, b, dchi


def main():
    rows = []
    for fname, kind in [('class_screen_model.json', 'collisionless'), ('class_screen_elastic.json', 'elastic q_eff')]:
        p = RES / fname
        if not p.exists():
            continue
        models = json.loads(p.read_text())['models']
        conv_u = {m['label'].split('_converted')[0]: m.get('ubar') for m in models if 'converted' in m['label']}
        for r in models:
            allx = 'allfreeX' in r['label']
            if allx and 'q_eff' not in r and 'ubar' not in r:
                r['ubar'] = conv_u.get(r['label'].split('_allfreeX')[0])
            f = 0.0672 if allx else 0.0336
            u = r.get('q_eff', r.get('ubar'))
            if u is None:
                continue
            lab, m_eq, b, dchi = classify(r, f, u)   # all-free-X: total f against the warmer component's bound (conservative)
            rows.append(dict(model=r['label'], kind=kind, a_form=r.get('a_form'),
                             z_form=(1 / r['a_form'] - 1) if r.get('a_form') else None,
                             f=f, q_over_m=u, m_eq_eV=m_eq, dA=r['dA'], sigma8_ratio=r['sigma8_ratio'],
                             P_ratio_z3={k: r[f'P_ratio_z3_k{k}'] for k in ['1.0', '3.0', '5.0', '10.0', '20.0']},
                             D2star_ratio=r['D2star_ratio'], dnstar=r['dnstar'],
                             dchi2_DESI=dchi, dchi2_eBOSS=r['chi2_eBOSS'] - r['chi2_eBOSS_LCDM'],
                             mixed_bound=b, label=lab))
    out = dict(calibrators=dict(dA_5p7keV=dA57, dA_5p3keV=cal['thermal_5.3keV_full']['dA'], dA_3p5keV=dA35,
                                dA_mixed_1keV_f016_Tfull=cal['mixed_1.0keV_f0.16_Tfull']['dA'],
                                dA_mixed_1keV_f016_Town=cal['mixed_1.0keV_f0.16_Town']['dA'],
                                dA_mixed_2keV_f035_Tfull=cal['mixed_2.0keV_f0.35_Tfull']['dA'],
                                dA_mixed_2keV_f035_Town=cal['mixed_2.0keV_f0.35_Town']['dA']),
               rows=rows)
    (RES / 'structure_screen_summary.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out['calibrators'], indent=1))
    for r in sorted(rows, key=lambda x: (x['kind'], 'allfreeX' in x['model'], -(x['z_form'] or 0))):
        print(f"{r['kind']:14s} {r['model'][:38]:38s} z={r['z_form'] or 0:9.3g} f={r['f']:.4f} m_eq={r['m_eq_eV']:9.3g} "
              f"dA={r['dA']:.3f} s8={r['sigma8_ratio']:.4f} P(k=1,5,20)={r['P_ratio_z3']['1.0']:.3f},{r['P_ratio_z3']['5.0']:.3f},{r['P_ratio_z3']['20.0']:.3f} "
              f"D2r={r['D2star_ratio']:.3f} dn={r['dnstar']:+.4f} dchiDESI={r['dchi2_DESI']:+.2f} dchiEB={r['dchi2_eBOSS']:+.1f} -> {r['label']}")


if __name__ == '__main__':
    main()
