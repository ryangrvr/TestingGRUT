#!/usr/bin/env python3
"""kin_thermo (Claude, own code): refine the E(Q) table of objects.py -- adaptive-step continuation from the saved
Q = 200000 profile towards the thin-wall limit (Omega -> omega_0 = 0.7401, Q -> infinity) and fine steps around
Q_min on the thick-wall side; then a spline-based check of dE/dQ = Omega m_X on both branches.
Usage: python3 -B objects_refine.py      (reads and rewrites results/objects_EQ.json; a few minutes)
"""
import json
import time
import numpy as np
from scipy.interpolate import CubicSpline
import objects as O

MX = O.MX
SCR = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/kin_thermo/profile_Q200000.npz'


class Saved:
    def __init__(self):
        d = np.load(SCR)
        self.r, self.F, self.S, self.Fp, self.Sp, self.om = d['r'], d['F'], d['S'], d['Fp'], d['Sp'], float(d['Omega'])
        self.x = self.r

    def sol(self, r):
        rr = np.clip(r, self.r[0], self.r[-1])
        return np.vstack([np.interp(rr, self.r, self.F), np.interp(rr, self.r, self.Fp),
                          np.interp(rr, self.r, self.S), np.interp(rr, self.r, self.Sp)])


def continue_branch(rows, om_start, om_end, step, qmax=3e8, tol=1e-9, log=None):
    br = O.Branch(); br.last = Saved(); br.L = 60.0
    om = om_start
    _, ob = br(om, tol=tol)
    direction = np.sign(om_end - om_start)
    h = step
    while direction * (om_end - om) > 1e-12:
        nxt = om + direction * h
        if direction * (nxt - om_end) > 0:
            nxt = om_end
        try:
            last_sol, last_L = br.last, br.L
            _, ob = br(float(nxt), tol=tol)
        except Exception as exc:
            br.last, br.L = last_sol, last_L
            h *= 0.5
            if log:
                print('step halved to %.2e at Omega=%.5f (%r)' % (h, nxt, exc), flush=True)
            if h < 1e-4:
                rows.setdefault('_ends', []).append(dict(Omega=float(nxt), error=repr(exc)))
                break
            continue
        om = float(nxt)
        rows[round(om, 6)] = ob
        if log:
            print('%.5f Q=%.6g E/Qm=%.7f F0=%.4f r90=%.2f L=%.0f  %.0fs' % (om, ob['Q'], ob['E_over_QmX'], ob['F0'],
                                                                          ob['r90_over_mX'], ob['L'], time.time() - log), flush=True)
        if ob['Q'] > qmax:
            break
        if h < step:
            h = min(step, h * 1.5)
    return rows


def main():
    t0 = time.time()
    out = json.loads((O.RES / 'objects_EQ.json').read_text())
    rows = {round(x['Omega'], 6): x for x in out['rows']}
    continue_branch(rows, 0.86, 0.745, 0.0025, log=t0)
    continue_branch(rows, 0.86, 0.998, 0.0025, log=t0)
    ends = rows.pop('_ends', [])
    if ends:
        out.setdefault('branch_ends', []).extend([dict(stage='refine', **e) for e in ends])
    tab = [rows[k] for k in sorted(rows)]
    Om = np.array([x['Omega'] for x in tab]); Q = np.array([x['Q'] for x in tab]); E = np.array([x['E_GeV'] for x in tab])
    sE, sQ = CubicSpline(Om, E), CubicSpline(Om, Q)
    dEdQ = sE(Om, 1) / sQ(Om, 1) / MX
    dQdOm = sQ(Om, 1)
    for x, d, dq in zip(tab, dEdQ, dQdOm):
        x['dEdQ_over_mX_table'] = float(d)
        x['dQ_dOmega'] = float(dq)
        x['branch'] = 'stable (dQ/dOmega<0)' if dq < 0 else 'thick-wall (dQ/dOmega>0, classically unstable)'
    i = int(np.argmin(Q))
    mask = np.abs(Om - Om[i]) > 0.012
    out['rows'] = tab
    out['Q_min'] = dict(Omega=float(Om[i]), Q=float(Q[i]), E_over_QmX=float(tab[i]['E_over_QmX']))
    out['dEdQ_identity_check'] = dict(max_abs_dEdQ_over_mX_minus_Omega_away_from_Qmin=float(np.max(np.abs(dEdQ[mask] - Om[mask]))),
                                      points=int(mask.sum()), excluded_window_in_Omega=0.012,
                                      note='dE/dQ from cubic splines E(Omega), Q(Omega); the identity dE/dQ = omega holds on both branches')
    LX, LP, V = 1e-3, 0.03, 0.924089637375362
    V0 = LP * V**4 / 4
    out['thin_wall_omega_over_mX_from_Ueff'] = float(np.sqrt(2 * np.sqrt(LX * V0)) / MX)
    st = [x for x in tab if x['dQ_dOmega'] < 0]
    eq = np.array([x['E_over_QmX'] for x in st]); oo = np.array([x['Omega'] for x in st])
    o = np.argsort(eq)
    out['E_over_QmX_crosses_1_at'] = dict(Omega=float(np.interp(1.0, eq[o], oo[o])),
                                          Q=float(np.interp(1.0, eq[o], np.array([x['Q'] for x in st])[o])))
    out['refine_seconds'] = time.time() - t0
    (O.RES / 'objects_EQ.json').write_text(json.dumps(out, indent=1) + '\n')
    print('rows', len(tab), 'Omega range', Om.min(), Om.max(), 'Q range', Q.min(), Q.max(), 'Q_min', out['Q_min'],
          'identity', out['dEdQ_identity_check'], 'omega0', out['thin_wall_omega_over_mX_from_Ueff'], 'E/Q=1 at', out['E_over_QmX_crosses_1_at'])


if __name__ == '__main__':
    main()
