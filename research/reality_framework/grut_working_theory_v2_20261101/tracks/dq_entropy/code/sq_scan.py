#!/usr/bin/env python3
"""dq_entropy item 2 (Claude, own code): entropy per charge S/Q(R) of the relaxed pre-onset Hartree plasma for the M1r
release computed in item 1 (frozen N_X, E_X of rpdq.c), for R in [0, 8] (X-sector conformal energy (1 + R) E_X at fixed
N_X); the release-invariant check (the relaxed state depends only on E/E_c(N), E_c = (N/c)^{4/3}); R_crit at the survival
thresholds of tracks/fx_chain/skeptic (survival_threshold.json); the onset a_on, T_h/m_X, n_Q; the same for the release with
the induced portals vanishing at M_P (fx_holder SKEPTIC N1).
Solver: hartree_dq.py = tracks/fx_chain/skeptic/code/hartree_indep.py plus set_release() (no other change).
Usage: sq_scan.py ../results/freeze.json ../results/sq_scan.json
"""
import json, math, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hartree_dq as HQ

WT = Path(__file__).resolve().parents[2]
SURV = json.load(open(WT / 'fx_chain' / 'skeptic' / 'results' / 'survival_threshold.json'))['thresholds']
SIG = {k: SURV[k]['sigma_min'] for k in ('bulk_Q_inf', 'Q_3e7', 'Q_2e5', 'Q_bind_71263')}
SPOT = dict(NX=HQ.NX, EX=HQ.EX, H=HQ.H, V=HQ.V)
LX = 1e-3
CC = 2 * math.sqrt(2 * LX) * (3 * LX) ** -0.75            # N_c = CC E^{3/4} (circular orbit, conformal quartic)


def e_over_ec(N, E):
    return E / (N / CC) ** (4.0 / 3.0)


class Release:
    def __init__(self, name, N, E, H_onset):
        self.name, self.N, self.E, self.H = name, abs(N), E, H_onset
        self.w0 = e_over_ec(self.N, self.E)
        self._g = None

    def solve(self, R, guess=None):
        HQ.set_release(self.N, self.E, self.H, SPOT['V'])
        st, g = HQ.solve(R, True, 0.0, guess or self._g)
        if (not st['success']) or max(abs(x) for x in st['resid']) > 1e-9:
            st, g = HQ.solve(R, True, 0.0, self.spot_guess(R))
        if (not st['success']) or max(abs(x) for x in st['resid']) > 1e-9:
            raise RuntimeError('%s: Hartree failed at R = %g' % (self.name, R))
        self._g = g
        a = HQ.onset_a(st); mc = HQ.MX * a / HQ.H
        return dict(R=R, S_per_Q=st['S'] / self.N, T_conf=st['T'], Xbar_per_Q=st['n_minus'] / self.N, X_per_Q=st['n_plus'] / self.N,
                    phi_per_Q=st['n_phi'] / self.N, cond_frac=st['cond_frac'], a_on=a, T_h_over_mX_on=st['T'] / mc,
                    nQ_on_mX3=self.N / mc ** 3, ncond_on_mX3=st['cond_frac'] * self.N / mc ** 3,
                    E_per_Q_on_mX=self.E * (1 + R) / self.N / mc, resid=max(abs(x) for x in st['resid']), _st=st, _g=g)

    def spot_guess(self, R):
        """solve at the spot normalisation at the same E/E_c and scale the solution (scale invariance)."""
        HQ.set_release(SPOT['NX'], SPOT['EX'], SPOT['H'], SPOT['V'])
        w_spot = e_over_ec(SPOT['NX'], SPOT['EX'])
        Reff = self.w0 * (1 + R) / w_spot - 1
        st, g = HQ.solve(Reff, True, 0.0, None)
        lam = (self.N / SPOT['NX']) ** (1 / 3)
        HQ.set_release(self.N, self.E, self.H, SPOT['V'])
        return [g[0] + math.log(lam), g[1] + 2 * math.log(lam), g[2] + 2 * math.log(lam), g[3] + 2 * math.log(lam)]


def spot_point(Reff, H_onset):
    HQ.set_release(SPOT['NX'], SPOT['EX'], H_onset, SPOT['V'])
    st, g = HQ.solve(Reff, True, 0.0, None)
    a = HQ.onset_a(st); mc = HQ.MX * a / HQ.H
    return dict(R_eff_spot=Reff, S_per_Q=st['S'] / SPOT['NX'], T_h_over_mX_on=st['T'] / mc, nQ_on_mX3=SPOT['NX'] / mc ** 3,
                Xbar_per_Q=st['n_minus'] / SPOT['NX'], cond_frac=st['cond_frac'], resid=max(abs(x) for x in st['resid']))


def strip(r):
    return {k: v for k, v in r.items() if not k.startswith('_')}


def main():
    t0 = time.time()
    fr = json.load(open(sys.argv[1]))
    H_freeX = 1.293e11            # free-X calibration (wt_critic; recomputed in summary.py); a_on scales as H_i
    rel = dict(benchmark_cMY0=Release('benchmark (portals 0 at M_Y)', fr['base_1e5']['frozen']['NX']['value'], fr['base_1e5']['E_X_late'], H_freeX),
               planck_bc=Release('portals 0 at M_P', fr['planck_1e5']['frozen']['NX']['value'], fr['planck_1e5']['E_X_late'], 1.124e11 * 0.957))
    out = dict(sigma_min=SIG, CC=CC, spot=dict(**SPOT, E_over_Ec=e_over_ec(SPOT['NX'], SPOT['EX'])), releases={})
    Rgrid = [round(x, 4) for x in np.arange(0.0, 8.0001, 0.25)]
    for key, r in rel.items():
        rows = []
        r._g = r.spot_guess(0.0)
        for R in Rgrid:
            rows.append(r.solve(R))
        sq = np.array([x['S_per_Q'] for x in rows])
        print('%s: E/E_c = %.4f; S/Q(R=0) = %.5f, S/Q(4) = %.5f, S/Q(8) = %.5f  (%.0f s)' % (r.name, r.w0, sq[0], sq[16], sq[-1], time.time() - t0), flush=True)
        # R_crit at each threshold
        Rc = {}
        for k, s in SIG.items():
            if sq[0] >= s:
                Rc[k] = dict(sigma=s, R_crit=None, note='S/Q(R = 0) already above the threshold')
                continue
            if sq[-1] <= s:
                Rc[k] = dict(sigma=s, R_crit=None, note='above R = 8')
                continue
            j = int(np.where(sq > s)[0][0])
            Rx = brentq(lambda R: r.solve(R)['S_per_Q'] - s, Rgrid[j - 1], Rgrid[j], xtol=1e-10)
            pt = r.solve(Rx)
            Rc[k] = dict(sigma=s, R_crit=Rx, W_crit_Ec_units=r.w0 * (1 + Rx), onset=strip(pt))
        # release-invariance check: direct M1r solve against the spot-normalised solver at the same E/E_c
        inv = []
        w_spot = e_over_ec(SPOT['NX'], SPOT['EX'])
        for R in (0.0, 1.0, 2.5, 5.0, 8.0):
            d = r.solve(R); s = spot_point(r.w0 * (1 + R) / w_spot - 1, r.H)
            inv.append(dict(R=R, direct=dict(S_per_Q=d['S_per_Q'], T_h_over_mX_on=d['T_h_over_mX_on'], nQ_on_mX3=d['nQ_on_mX3'], Xbar_per_Q=d['Xbar_per_Q']),
                            spot_normalised=s, dS_rel=s['S_per_Q'] / d['S_per_Q'] - 1, dT_rel=s['T_h_over_mX_on'] / d['T_h_over_mX_on'] - 1,
                            dn_rel=s['nQ_on_mX3'] / d['nQ_on_mX3'] - 1))
        # self-consistent onset (tachyonic -kappa a^2 inside the phi Hartree mass) at a few R
        sc = []
        for R in (0.0, 2.0, 4.0):
            d = r.solve(R); g = d['_g']
            HQ.set_release(r.N, r.E, r.H, SPOT['V'])
            def gap(a):
                s2, _ = HQ.solve(R, True, HQ.KAPPA * a * a, g)
                return s2['Mp2'] / (HQ.KAPPA * a * a)
            try:
                asc = brentq(gap, d['a_on'] * 0.999, d['a_on'] * 1.2, xtol=1e-6 * d['a_on'])
                s2, _ = HQ.solve(R, True, HQ.KAPPA * asc * asc, g); mcs = HQ.MX * asc / HQ.H
                sc.append(dict(R=R, a_sc_over_a_on=asc / d['a_on'], T_h_over_mX_sc=s2['T'] / mcs, nQ_sc=r.N / mcs ** 3,
                               S_per_Q_sc=s2['S'] / r.N, S_per_Q=d['S_per_Q']))
            except ValueError as e:
                sc.append(dict(R=R, error=repr(e)))
        out['releases'][key] = dict(name=r.name, N_X=r.N, E_X=r.E, E_over_Ec=r.w0, H_onset_GeV=r.H, scan=[strip(x) for x in rows],
                                    R_crit=Rc, invariance=inv, selfconsistent_onset=sc)
        print('   R_crit:', {k: (round(v['R_crit'], 5) if v['R_crit'] is not None else v['note']) for k, v in Rc.items()})
        print('   invariance max |dS/S| = %.2e, |dT/T| = %.2e, |dn/n| = %.2e' % (max(abs(x['dS_rel']) for x in inv), max(abs(x['dT_rel']) for x in inv),
                                                                            max(abs(x['dn_rel']) for x in inv)))
        print('   self-consistent onset:', sc, flush=True)
    out['seconds'] = time.time() - t0
    json.dump(out, open(sys.argv[2], 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
