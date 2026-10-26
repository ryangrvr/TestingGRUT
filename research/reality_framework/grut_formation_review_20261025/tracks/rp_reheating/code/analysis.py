#!/usr/bin/env python3
"""rp_reheating track (Claude, own code): tables and summary from the run campaign.
Codex's numbers below are quoted [Q] from e160356 REHEATING_RESULTS.json / INDEPENDENT_REHEATING.json (read through git
archive into the scratchpad).  hw_verify numbers are quoted [Q, own track].  The decay ledger is our hw_verify/hw_gates
ledger reused with the reheating correction (Rr/3)^(1/4) on <a> [Q, own tracks; D].
Usage: python3 -B analysis.py ../results      -> results/summary.json and markdown tables on stdout
"""
import json, sys
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import rp

CODEX = {  # tag: (late eps, late NX, late NYH, late NYL, yield/YQ, final a, final H)  [Q e160356]
    'b2_g0.03': (-0.19127890229398642, -12425.475159668482, 2.015855148825169e-05, 0.0035042871115903477, 0.3326972180316388, 96.30425251103159, 0.0006081416606183946),
    'b2_g0.1': (-0.1987690891357234, -12263.839011189242, 2.015543925229144e-05, 0.0035064832732232365, 0.559555155674355, 44.19324828042198, 0.0020242476283894614),
    'b2_g0.3': (-0.21868352842307429, -12041.877784202039, 2.0149540459676422e-05, 0.003511623835289129, 0.7842705752843221, 22.69572716931241, 0.006054112027336305),
    'b10_g0.03': (-0.18147886412418332, -3309.409672295826, 0.0032447285564058876, 1.7397749997223992, 0.18472958785025992, 75.39108631266005, 0.0006080731670753927),
    'b10_g0.1': (-0.18152781879915347, -3314.5903615678008, 0.0031974802403411636, 1.707133532726863, 0.3250488273535754, 34.25137903002546, 0.002023412882428783),
    'b10_g0.3': (-0.17971470443740217, -3339.2056835367794, 0.003068098059885127, 1.6142280972003942, 0.5181727753504077, 17.009124254006217, 0.006042290573268437),
    'long_500': (-0.198329404135299, -12262.117232715218, 2.0155424005444842e-05, 0.003506483642845795, 0.5594765970385078, 62.68783548936551, 0.0010060253811141116),
    'xi_zero_control': (-0.15981961468542605, -2595.7420178295683, 0.017230338964075588, 0.33814045369557455, 0.1184344337543063, 44.19324828092177, 0.0020242476286285098),
}
CODEX_CHECKS = dict(max_Friedmann=1.2685176254738516e-09, Radau_state_error=3.2061945360371023e-12,
                    long_NX=0.00014039473874805086, long_eps=0.002212039117028719, phase_response=4.3523062931427603e-10)
OLD_LATE = dict(eps=-0.14571744282991028, NX=-2788.1767940228406, NYH=0.01725303422930909, NYL=0.3344649635088846)  # [Q parent BACKGROUND]
HW_VERIFY = {0.1: (-0.3804, 5416, 0.545), 0.3: (-0.3456, None, 0.835), 1.0: (-0.2791, None, 1.16), 3.0: (-0.2135, None, 1.30), 10.0: (-0.1720, None, 1.24)}  # [Q own]

# ledger constants at the spot point [D, own hw_verify/hw_gates]
H = rp.H0; M, h, mu, kap, e = 30 * H, 40 * H, rp.P['mu'] * H, rp.P['kap'], 1e-7
mL = np.sqrt(M**2 - np.sqrt(2) * mu * h); mA = 6 * e * h
GL = 3 * kap**2 * mL / (256 * np.pi**3); GA = e**2 * mA / (48 * np.pi)
DEP_L = (mL / H) * np.sqrt(np.pi * H / (2 * GL)); DEP_A = (mA / H) * np.sqrt(np.pi * H / (2 * GA))


def ledger(row):
    L = row['late']; corr = (row['Rr_end'] / 3) ** 0.25
    RY = (L['NYH'] + L['NYL']) * DEP_L * corr / L['EX']; RC = L['NYH'] * DEP_A * corr / L['EX']
    RS = 2 * L['Nsig'] * DEP_A * corr / L['EX']
    return dict(R_Y=RY, R_cascade=RC, R_Sigma=RS, R_central=RY + RC + RS, reheat_factor=corr)


def load(d, name):
    return {r['tag']: r for r in json.loads((d / (name + '.json')).read_text())['rows']}


def main():
    d = Path(sys.argv[1])
    rep, gs, bs, md, ex = load(d, 'reproduce'), load(d, 'gamma_scan'), load(d, 'b_scan'), load(d, 'models'), load(d, 'extra')
    summary = {}
    print('## Reproduction (own C engine, constraint-reconstructed H, DP5(4)) vs Codex [Q]')
    print('| run | eps ours / Codex | N_X ours / Codex | N_YH | N_YL ours / Codex | yield/Y_Q ours / Codex | charge ratio | dilution D | a_f ours / Codex |')
    print('|---|---|---|---|---|---|---|---|---|')
    comp = {}
    for tag, cx in CODEX.items():
        r = rep[tag]; L = r['late']
        rel = dict(eps=abs(L['eps'] / cx[0] - 1), NX=abs(L['NX'] / cx[1] - 1), NYH=abs(L['NYH'] / cx[2] - 1), NYL=abs(L['NYL'] / cx[3] - 1),
                   yld=abs(r['yield_over_YQ'] / cx[4] - 1), a=abs(r['final_a'] / cx[5] - 1), H=abs(r['final_H'] / cx[6] - 1))
        comp[tag] = dict(ours=dict(eps=L['eps'], NX=L['NX'], NYH=L['NYH'], NYL=L['NYL'], yield_over_YQ=r['yield_over_YQ'], final_a=r['final_a'], final_H=r['final_H'],
                                   charge_ratio=r['charge_ratio_to_instant'], dilution=r['dilution'], Rr=r['Rr_end'], checks=r['checks']),
                         codex=dict(eps=cx[0], NX=cx[1], NYH=cx[2], NYL=cx[3], yield_over_YQ=cx[4], final_a=cx[5], final_H=cx[6]), relative_differences=rel)
        print('| %s | %.5f / %.5f | %.1f / %.1f | %.3e | %.4e / %.4e | %.4f / %.4f | %.3f | %.3f | %.3f / %.3f |'
              % (tag, L['eps'], cx[0], L['NX'], cx[1], L['NYH'], L['NYL'], cx[3], r['yield_over_YQ'], cx[4], r['charge_ratio_to_instant'], r['dilution'], r['final_a'], cx[5]))
    worst = max(max(v['relative_differences'].values()) for v in comp.values())
    print('worst relative difference over all compared quantities: %.2e' % worst)
    summary['reproduction'] = comp; summary['worst_relative_difference'] = worst
    # controls
    c, t, lo = rep['b2_g0.1'], rep['central_tight'], rep['central_loose']
    tol = max(abs(c['late'][k] / t['late'][k] - 1) for k in ('NX', 'EX', 'eps')); tol_lo = max(abs(c['late'][k] / lo['late'][k] - 1) for k in ('NX', 'EX', 'eps'))
    ph = max(abs(rep[k]['late']['NX'] / c['late']['NX'] - 1) for k in ('phase_plus', 'phase_minus'))
    slope = (rep['phase_plus']['late']['NX'] - rep['phase_minus']['late']['NX']) / (2e-8 * abs(c['late']['NX']))
    long5, long15 = rep['long_500'], rep['long_1500']
    rd = ex['RD_reference_dS']; rda = rd['a3040']
    rd_err = dict(eps=abs(rda['eps'] / OLD_LATE['eps'] - 1), NX=abs(rda['NX'] / OLD_LATE['NX'] - 1), NYH=abs(rda['NYH'] / OLD_LATE['NYH'] - 1), NYL=abs(rda['NYL'] / OLD_LATE['NYL'] - 1))
    summary['controls'] = dict(
        tight_rtol_1e12_relative_change=tol, loose_rtol_1e8_relative_change=tol_lo, phase_1e8_response=ph, charge_derivative_wrt_phase=slope,
        long_500=dict(NX=abs(long5['late']['NX'] / c['late']['NX'] - 1), eps=abs(long5['late']['eps'] / c['late']['eps'] - 1), eps_value=long5['late']['eps']),
        long_1500=dict(NX=abs(long15['late']['NX'] / c['late']['NX'] - 1), eps=abs(long15['late']['eps'] / c['late']['eps'] - 1), eps_value=long15['late']['eps']),
        RD_reference_a30_40=dict(ours=rda, old_late=OLD_LATE, relative_errors=rd_err),
        max_H_consistency_all_production=max(rep[k]['checks']['max_H_consistency'] for k in CODEX if k in rep),
        max_gauss_all=max(r['checks']['max_gauss'] for r in rep.values()), max_hidden_work_all=max(r['checks']['max_hidden_work'] for r in rep.values() if r['tag'] != 'RD_reference'),
        codex=CODEX_CHECKS)
    print('\ncontrols:', json.dumps(summary['controls'], indent=1, default=float))

    print('\n## Yield deficit decomposition: yield/Y_Q = (|N_X|/2788.18) / D,  D = (rho_r a^4 / 3)^(3/4)')
    print('| m_I/H_i | Gamma/H_i | Gamma/m_I | charge ratio | D | D analytic matter-era | yield/Y_Q | eps | N_YL |')
    print('|---|---|---|---|---|---|---|---|---|')
    dec = []
    allruns = {**gs, **{k: v for k, v in ex.items() if k.startswith('yr_')}}
    for tag in sorted(allruns, key=lambda k: (allruns[k]['b'], allruns[k]['gamma'])):
        r = allruns[tag]
        if r['gamma'] > 0.3 * r['b'] and r['gamma'] > 0.5:   # beyond perturbative decay (Gamma > 0.3 m_I): flagged, kept in JSON
            flag = ' (Gamma > 0.3 m_I)'
        else:
            flag = ''
        Dan = (0.9027452929509336 * (1.5 / r['gamma'])**(2 / 3))**0.75
        dec.append(dict(tag=tag, b=r['b'], gamma=r['gamma'], gamma_over_mI=r['gamma'] / r['b'], charge_ratio=r['charge_ratio_to_instant'], D=r['dilution'],
                        D_analytic=Dan, yield_over_YQ=r['yield_over_YQ'], eps=r['late']['eps'], NYL=r['late']['NYL'], NYH=r['late']['NYH'], Rr=r['Rr_end'], flag=flag.strip()))
        if abs(r['yield_over_YQ']) < 10:
            print('| %g | %g | %.2f | %.3f | %.3f | %.3f | %.3f | %.4f | %.3e |%s' % (r['b'], r['gamma'], r['gamma'] / r['b'], r['charge_ratio_to_instant'], r['dilution'], Dan,
                                                                                    r['yield_over_YQ'], r['late']['eps'], r['late']['NYL'], flag))
    summary['deficit_decomposition'] = dec

    print('\n## m_I scan at Gamma/H_i = 0.1: free Y occupation (resonance with 2 m_I)')
    print('| m_I/H_i | 2 m_I/H_i | N_YL | N_YH | eps | yield/Y_Q |')
    print('|---|---|---|---|---|---|')
    scan = []
    for tag in sorted(bs, key=lambda k: bs[k]['b']):
        r = bs[tag]; scan.append(dict(b=r['b'], NYL=r['late']['NYL'], NYH=r['late']['NYH'], eps=r['late']['eps'], yield_over_YQ=r['yield_over_YQ'], NX=r['late']['NX']))
        print('| %g | %g | %.4e | %.3e | %.4f | %.4f |' % (r['b'], 2 * r['b'], r['late']['NYL'], r['late']['NYH'], r['late']['eps'], r['yield_over_YQ']))
    wL = np.sqrt(900 - 2 * rp.P['mu'] * rp.P['h'] / np.sqrt(2)); wH = np.sqrt(900 + 2 * rp.P['mu'] * rp.P['h'] / np.sqrt(2))
    peakL = max(scan, key=lambda s: s['NYL']); peakH = max(scan, key=lambda s: s['NYH'])
    summary['mI_scan'] = dict(rows=scan, omega_YL_over_H=float(wL), omega_YH_over_H=float(wH), peak_NYL=peakL, peak_NYH=peakH)
    print('omega_YL = %.3f H, omega_YH = %.3f H; N_YL peak %.3f at m_I = %g (2 m_I = %g); N_YH peak %.3f at m_I = %g (2 m_I = %g)'
          % (wL, wH, peakL['NYL'], peakL['b'], 2 * peakL['b'], peakH['NYH'], peakH['b'], 2 * peakH['b']))

    print('\n## Model variants at m_I/H_i = 2 preparation (de Sitter conditional minimum), fluid inflaton backgrounds')
    print('| Gamma/H_i | oscillator (Codex model) eps / yield | fluid, xi<R> = -2.5 H^2 eps / yield | fluid, -c_H H^2 (hw_verify model) eps / yield | hw_verify [Q] eps / diluted yield |')
    print('|---|---|---|---|---|')
    mv = []
    for g in (0.03, 0.1, 0.3, 1.0, 3.0, 10.0):
        o = gs.get('gs_b2_g%g' % g); f1 = ex['fluidxiR_dS_g%g' % g]; f2 = md['fluidcH_g%g' % g]; hw = HW_VERIFY.get(g)
        mv.append(dict(gamma=g, oscillator=dict(eps=o['late']['eps'], yld=o['yield_over_YQ'], NX=o['late']['NX'], NYL=o['late']['NYL']) if o else None,
                       fluid_xiR=dict(eps=f1['late']['eps'], yld=f1['yield_over_YQ'], NX=f1['late']['NX'], NYL=f1['late']['NYL']),
                       fluid_cH=dict(eps=f2['late']['eps'], yld=f2['yield_over_YQ'], NX=f2['late']['NX'], NYL=f2['late']['NYL']), hw_verify=hw))
        print('| %g | %s | %.4f / %.3f | %.4f / %.3f | %s |' % (g, ('%.4f / %.3f' % (o['late']['eps'], o['yield_over_YQ'])) if o else '-', f1['late']['eps'], f1['yield_over_YQ'],
                                                              f2['late']['eps'], f2['yield_over_YQ'], ('%.4f / %.3f' % (hw[0], hw[2])) if hw else '-'))
    summary['model_variants'] = mv

    print('\n## Mediator debris ledger (reused hw_verify ledger with (Rr/3)^(1/4) on <a>) [E]')
    print('| run | N_Y | E_X | reheat factor | R_Y | R_cascade | R_central |')
    print('|---|---|---|---|---|---|---|')
    led = {}
    for tag in ('b2_g0.1', 'b2_g0.3', 'b10_g0.03', 'b10_g0.1', 'b10_g0.3'):
        r = rep[tag]; L = ledger(r); led[tag] = L
        print('| %s | %.4f | %.0f | %.3f | %.3f | %.3f | %.3f |' % (tag, r['late']['NYH'] + r['late']['NYL'], r['late']['EX'], L['reheat_factor'], L['R_Y'], L['R_cascade'], L['R_central']))
    for tag in ('yr_b3_g0.5', 'yr_b5_g1', 'yr_b10_g3', 'yr_b15_g4.5', 'yr_b20_g6', 'yr_b50_g10'):
        src = ex if tag in ex else gs
        if tag not in src:
            tag2 = tag.replace('yr_', 'gs_'); src = gs; r = src[tag2]
        else:
            r = src[tag]
        L = ledger(r); led[tag] = L
        print('| %s | %.4f | %.0f | %.3f | %.3f | %.3f | %.3f |' % (tag, r['late']['NYH'] + r['late']['NYL'], r['late']['EX'], L['reheat_factor'], L['R_Y'], L['R_cascade'], L['R_central']))
    summary['ledger'] = dict(rows=led, deposit_per_Y_quantum_in_EX_units=float(DEP_L), deposit_per_heavy_cascade=float(DEP_A), mL_GeV=float(mL), mA_GeV=float(mA),
                             GammaL_GeV=float(GL), GammaA_GeV=float(GA), note='R relative to E_X; isentropic proxy R_L = 81.5-82 at eps 0.146 and Y_Q; energy-conserving proxy 7-8 [Q own tracks]')
    (d / 'summary.json').write_text(json.dumps(summary, indent=1, default=float))


if __name__ == '__main__':
    main()
