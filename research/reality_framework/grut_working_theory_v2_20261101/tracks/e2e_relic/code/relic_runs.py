#!/usr/bin/env python3
"""e2e_relic item 2 (Claude, own driver over READ-ONLY qm_micro + qm_micro-skeptic code): dependence of the Q-matter relic
on the charge yield Y_Q = n_Q/s (equivalently the hidden-to-visible entropy ratio r_hv = (S/Q) Y_Q) at fixed S/Q.

Pipeline per Y_Q factor f (all files in this track's scratch, SK_SCR = $SCR/yq<f>):
  1. qm_micro/code/qm_vapour.py <hist_3.0.json> 3.0 with QM_MASS = 0.785/f, i.e. Y_Q = 5.07941423e-9 * 0.9374/0.785 * f
     (two-phase equilibrium history on the visible clock);
  2. qm_micro/skeptic/code/coarsen_heat.py, case bose_heat_sphi1 (the skeptic-corrected ripening: Bose factor + latent-heat
     factor F_h, s_phi = 1) -> charge-weighted <Q> at capture freeze-out;
  3. qm_micro/skeptic/code/sk_freeze.py Model.run_energy (energy-conserving two-moment vapour on zero-heat-capacity objects)
     with the objects of step 2 -> f_v; perfect trapping and qm_micro's Beliaev trapping factor at the objects' radius.
Nothing in qm_micro is modified; only the env (QM_MASS, SK_SCR) and the output paths differ.  [C on E models]
Usage: SK_SCR=<dir> python3 -B relic_runs.py <f> <out.json>
"""
import json, math, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
QS = os.path.join(HERE, '..', '..', 'qm_micro', 'skeptic', 'code')
sys.path.insert(0, QS)
import sk_freeze as SF        # noqa: E402  read-only (reads SK_SCR at import)

if __name__ == '__main__':
    f = float(sys.argv[1]); out = sys.argv[2]
    cj = os.path.join(SF.SCR, 'cheat_3.0.json')
    d = json.load(open(cj)); run = [r for r in d['runs'] if r['case'] == 'bose_heat_sphi1'][0]
    res = dict(Y_Q_factor=f, Y_Q=5.07941423e-9 * 0.9374 / 0.785 * f, S_over_Q=3.0,
               coarsen=dict(freeze=run['freeze'], final_f_v_imposedT=run['final']['f_v']), energy={})
    for trap in (False, True):
        t0 = time.time()
        SF.TRAP_R = None
        M = SF.Model('3.0', cj, 'bose_heat_sphi1')
        if trap:
            SF.TRAP_R = float(SF.ch.radius(M.Qeff0))
        re = M.run_energy()
        fo = next((r for r in re['rows'] if r['Gamma_a_over_H'] < 1), None)
        res['energy']['beliaev' if trap else 'perfect'] = dict(final_f_v=re['final']['f_v'], Gamma_a_eq_H=fo, Qeff0=M.Qeff0,
                                                               alpha0=M.alpha0, seconds=time.time() - t0)
        print('f=%g trap=%s: f_v=%.4f  (Gamma_a=H at T_eq=%.4g, T_vis=%.3g MeV)  %.0fs' % (
            f, trap, re['final']['f_v'], fo['T_eq'] if fo else float('nan'), fo['T_vis_GeV'] * 1e3 if fo else float('nan'),
            time.time() - t0), flush=True)
    json.dump(res, open(out, 'w'), indent=1, default=float)
