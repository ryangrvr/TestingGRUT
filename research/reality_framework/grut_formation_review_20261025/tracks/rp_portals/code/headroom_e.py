"""rp_portals / headroom_e.py (Claude, own code).

Allowed eps_BX and b_Hphi when the hidden gauge coupling is raised to e = 3e-6 or 1e-5 (the hw_critic remedy for the
inflationary-vector gate under the energy-conserving thresholds).  Existing debris at e [Q own hw_critic/hw_gates]:
R_A(e) = R_A(1e-7) x 1e-7/e (modes 48.47, GMR 64.0), R_med(e) = R_Y 0.805 + R_casc 0.3387 sqrt(1e-7/e) (R_Sigma ~ 0).
Thermal source at e from the momentum-resolved solver (eps grid), interpolated in log-log.

Usage: python3 -B headroom_e.py <higgs_source.json> <out.json>
"""
import json, math, sys, pathlib
import numpy as np
from freezein import momentum_solver

def main(hg_path, out):
    hg = json.load(open(hg_path)); Rh0 = hg['R_h_gstar106p75']; b0 = 5e-11
    root = pathlib.Path(__file__).resolve().parents[2]
    thr = json.load(open(root/'hw_critic/results/thresholds.json'))['thresholds']
    res = []
    for e in [1e-7, 3e-6, 1e-5]:
        grid = [1e-11, 1e-10, 1e-9, 1e-8, 3e-8, 1e-7, 3e-7]
        rows = [(eps, momentum_solver(e=e, eps=eps, rtol=1e-9)) for eps in grid]
        le = np.log([r[0] for r in rows]); lR = np.log([r[1]['R_thermal'] for r in rows])
        # keep the monotone branch
        imax = int(np.argmax(lR)); le, lR = le[:imax+1], lR[:imax+1]
        Rmed = 0.805 + 0.3387*math.sqrt(1e-7/e)
        for vec, RA0 in [('modes', 48.47), ('GMR', 64.0)]:
            RA = RA0*1e-7/e
            for tname, RL in thr.items():
                head = RL - Rmed - RA
                eps_max = float(math.exp(np.interp(math.log(head), lR, le))) if head > 0 and math.log(head) < lR[-1] else (0.0 if head <= 0 else None)
                res.append(dict(e=e, vector=vec, R_A=RA, R_med=Rmed, threshold=tname, R_L=RL, headroom=head,
                                eps_BX_max=eps_max, eps_over_e_max=(eps_max/e if eps_max else eps_max),
                                b_Hphi_max=(b0*math.sqrt(head/Rh0) if head > 0 else 0.0)))
        print('e', e, 'R_thermal grid', [(r[0], round(r[1]['R_thermal'], 3)) for r in rows], flush=True)
    json.dump(res, open(out, 'w'), indent=1)
    for r in res:
        if r['threshold'] in ('direct_proxy_hw_verify', 'ke_muc0.861', 'ke_plus_s'):
            print(f"e {r['e']:.0e} {r['vector']:6s} {r['threshold']:24s} head {r['headroom']:8.2f} eps_max {r['eps_BX_max']} b_max {r['b_Hphi_max']:.3g}")

if __name__ == '__main__':
    main(*sys.argv[1:3])
