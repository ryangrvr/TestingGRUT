"""qm_transfer validation 1 (own driver): the kinetic code on the EXACT vapour history of tracks/fx_observables
(history_B, sudden capture freeze-out, T_v = T_eq until a_f), against the skeptic's two-fluid integration
(tracks/fx_observables/skeptic/code/two_fluid_check.py -> skeptic/results/two_fluid_transfer.json).  The track's and
the skeptic's modules are imported read-only (their main() is not run); nothing of theirs is modified.
Same cosmology as theirs (h = 0.674, omega_c = 0.12, omega_b = 0.0224, N_eff = 3.044).  [C]
usage: val_twofluid.py OUT.json [case indices]
"""
import importlib.util, json, math, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from qt_bg import Cosmo
from qt_hidden import Vapour
from qt_boltz import Mode, CDM

FXO = os.path.abspath(os.path.join(HERE, '..', '..', 'fx_observables'))
spec = importlib.util.spec_from_file_location('tfc', os.path.join(FXO, 'skeptic', 'code', 'two_fluid_check.py'))
tfc = importlib.util.module_from_spec(spec); spec.loader.exec_module(tfc)
vt = tfc.vt

CASES = [(1e14, 1.235, 1.0), (1e12, 1.235, 1.0), (1e10, 1.235, 1.0), (1e8, 1.235, 1.0), (2e5, 1.235, 1e-2),
         (1e12, 2.07, 1.0), (1e8, 2.07, 1.0)]
KS = [1.0, 3.0, 10.0, 20.0, 50.0, 100.0, 300.0]


def main():
    out_fn = sys.argv[1]
    sel = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else range(len(CASES))
    c = Cosmo(h=0.674, omb=0.0224, omc=0.12, Neff=3.044)
    ref = {}
    for r_ in json.load(open(os.path.join(FXO, 'skeptic', 'results', 'two_fluid_transfer.json')))['rows']:
        ref[(r_['Q'], r_['S_over_Q'], r_['sticking'])] = r_
    res = json.load(open(out_fn)) if os.path.exists(out_fn) else {}
    cdm_cache = {}
    for i in sel:
        Q, SQ, s = CASES[i]
        t0 = time.time()
        Hf = tfc.history_full(Q, SQ, s, sigma_m=332.0 if Q == 2e5 else None)
        p = vt.obj_props(Q)
        R = p['R']
        sig_m = 332.0 if Q == 2e5 else math.pi * R * R
        fxh = dict(a=Hf['a'], Nv=Hf['Nv'], Tv=Hf['Tv'], af=Hf['af'], am=Hf['track']['a_m'], sig_c=s * math.pi * R * R,
                   sig_m=sig_m, Q=Q, NQ0=vt.NQ0)
        v = Vapour(c, SQ=SQ, Q=Q, s=s, mode='fx', fx_hist=fxh, NQ0=vt.NQ0, Delta=1 - p['Omega'], Nq=28, Lv=10,
                   sig_a=s * math.pi * R * R, sig_el=sig_m - s * math.pi * R * R)
        row = dict(Q=Q, S_over_Q=SQ, sticking=s, summary={k: (float(x) if isinstance(x, (float, np.floating)) else x)
                                                          for k, x in v.summary().items()}, P={})
        rr = ref.get((Q, SQ, s))
        for k in KS:
            m = Mode(c, v, k, K_rsa=200., K_end=400., rtol=1e-7)
            d1 = m.run()[1.0]
            if k not in cdm_cache:
                cdm_cache[k] = Mode(c, CDM(), k, K_rsa=200., K_end=400., rtol=1e-7).run()[1.0]
            P = (d1 / cdm_cache[k])**2
            tf = rr['P'][str(k)] if rr else None
            row['P'][str(k)] = dict(kinetic=P, two_fluid=tf['two_fluid'] if tf else None,
                                    two_fluid_visc=tf['two_fluid_visc'] if tf else None,
                                    track_single_fluid=tf['track_single_fluid'] if tf else None)
            print('Q=%.0e SQ=%.3f s=%g k=%g: kinetic %.4f | two-fluid %.4f | +visc %.4f | track %.4f   (%.0fs)' % (
                Q, SQ, s, k, P, tf['two_fluid'] if tf else -1, tf['two_fluid_visc'] if tf else -1,
                (tf['track_single_fluid'] or -1) if tf else -1, time.time() - t0), flush=True)
        res['%g_%g_%g' % (Q, SQ, s)] = row
        json.dump(res, open(out_fn, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
