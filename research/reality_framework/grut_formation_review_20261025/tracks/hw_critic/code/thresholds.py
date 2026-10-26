#!/usr/bin/env python3
"""hw_critic (Claude, own code): combine the tracks' debris pieces with every threshold variant on record, and an
independent GMR-formula estimate of the inflationary longitudinal-vector number.

Inputs [Q, own tracks]: mediator ledger at the spot point (hw_verify/hw_gates/hw_envelope agree): R_Y = 0.805,
R_Sigma = 1.2e-7, R_casc = 0.339 at e = 1e-7; vector piece A: N_A = 2.468 H_i^3 (own Proca modes, hw_gates) or
3.26 (GMR formula), R_A = 48.47 / 64.0 at e = 1e-7, deposit per vector 19.64 E_X, E_X = 38398, N_X = 2788.
Scalings [D, hw_gates]: R_A ∝ 1/e (N_A ∝ e^{-1/2}, deposit ∝ e^{-1/2}), R_Sigma and R_casc ∝ e^{-1/2}, R_Y fixed.
Thresholds at eps = 0.1457, normalised yield: sp_critic table 87.85 [Q]; direct proxy 81.5 (hw_verify, sp_yledger
code, endpoint 1e-6 GeV) and 81.9-82.2 (ls_critic, lsrelic code); Codex calibration ~9 % below our table (~80) [Q];
energy-conserving (ke) proxy [ls_critic]: mu_c = E_c = 0.861: f(7) = 1.4e-36, f(8) = 0.0057, f(8.9) = 0.041-0.045
-> crossing 8.2 (log interpolation between 8 and 8.9); mu_c = 0.923: f(5) = 3.7e-4, f(8.9) = 0.249 -> ~6.3;
ke + on-shell s route: f(15) = 4.6e-14, f(20) = 0.093, f(25) = 0.24 -> ~19.
"""
import json, math, sys
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hwc

E0 = 1e-7
PIECES = dict(R_Y=0.805, R_Sigma=1.24e-7, R_casc=0.3387, R_A_modes=48.47, R_A_GMR=64.0)
NX, EX = 2788.2, 38398.0


def crossing(pts):
    """log-linear crossing of f = 0.01 between the bracketing points (R, f)."""
    pts = sorted(pts)
    for (r1, f1), (r2, f2) in zip(pts, pts[1:]):
        if f1 < 0.01 <= f2:
            return r1 + (r2 - r1) * (math.log(0.01) - math.log(f1)) / (math.log(f2) - math.log(f1))
    return None


THR = {
    'table_sp_critic': 87.85,
    'direct_proxy_hw_verify': 81.5,
    'direct_proxy_ls_critic': 82.03,
    'codex_calibration_est': 80.0,
    'ke_muc0.861': crossing([(7.0, 1.4e-36), (8.0, 0.005685), (8.9, 0.0412)]),
    'ke_muc0.923': crossing([(5.0, 3.716e-4), (8.9, 0.2493)]),
    'ke_plus_s': crossing([(15.0, 4.57e-14), (20.0, 0.0931), (25.0, 0.241)]),
}


def pieces(e, vec='modes'):
    s = math.sqrt(E0 / e)
    RA = (PIECES['R_A_modes'] if vec == 'modes' else PIECES['R_A_GMR']) * (E0 / e)
    RS = PIECES['R_Sigma'] * s; RC = PIECES['R_casc'] * s
    return dict(e=e, R_Y=PIECES['R_Y'], R_Sigma=RS, R_casc=RC, R_med=PIECES['R_Y'] + RS + RC, R_A=RA,
                R_total=PIECES['R_Y'] + RS + RC + RA)


def gmr_estimate(H=hwc.H_SPOT, mA_GeV=6 * 1e-7 * 40 * hwc.H_SPOT):
    """Graham-Mardon-Rajendran: Omega_A/Omega_DM = sqrt(m/6e-6 eV) (H_I/1e14 GeV)^2 [R]; rho_DM/s = 4.4e-10 GeV [Q];
    converted to a comoving number per H_i^3 with the prompt-radiation entropy."""
    ratio = math.sqrt(mA_GeV / 6e-15) * (H / 1e14) ** 2
    n_over_s = ratio * 4.4e-10 / mA_GeV
    s, T = hwc.entropy_RH(H)
    N_A = n_over_s * s / H ** 3
    return dict(Omega_A_over_DM=ratio, n_A_over_s=n_over_s, N_A_per_Hi3=N_A, T_R=T, mA_GeV=mA_GeV,
                R_A_from_formula=N_A * 19.64, pairs_per_net_charge=N_A / NX)


if __name__ == '__main__':
    out = {'thresholds': THR, 'gmr': gmr_estimate(), 'rows': []}
    for e in (1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 1e-4, 1e-3):
        for vec in ('modes', 'GMR'):
            r = pieces(e, vec); r['vec'] = vec
            r['ratio'] = {k: r['R_total'] / v for k, v in THR.items()}
            r['mediator_only_ratio'] = {k: r['R_med'] / v for k, v in THR.items()}
            out['rows'].append(r)
    # headroom per net charge for the vector under each threshold (R_L - R_med)/(19.64 N_X)
    out['vector_headroom_per_net_charge'] = {k: (v - pieces(1e-7)['R_med']) / (19.64 * NX) for k, v in THR.items()}
    out['pieceA_pairs_per_net_charge'] = dict(modes=2.468 / NX, GMR=3.26 / NX)
    Path(sys.argv[1]).write_text(json.dumps(out, indent=1))
    print('thresholds', {k: round(v, 2) for k, v in THR.items()})
    print('GMR estimate', {k: (round(v, 4) if isinstance(v, float) and abs(v) > 1e-3 else v) for k, v in out['gmr'].items()})
    print('headroom per net charge', {k: '%.2e' % v for k, v in out['vector_headroom_per_net_charge'].items()})
    print('%-8s %-5s %8s %8s | %s' % ('e', 'vec', 'R_med', 'R_tot', '  '.join('%12s' % k[:12] for k in THR)))
    for r in out['rows']:
        print('%-8.0e %-5s %8.3f %8.2f | %s' % (r['e'], r['vec'], r['R_med'], r['R_total'],
                                                 '  '.join('%12.3f' % r['ratio'][k] for k in THR)))
