"""COPY of fx_formation/code/eq_history.py (qm_micro, 2 Nov 2026).  Changes: imports thermo_copy; RES points to
fx_formation/results (read-only); T_vis_of and hubble take an optional Y_Q (the Q-matter branch has Y_Q = D/(M m_X)).

fx_formation item 1c (Claude, 1 Nov 2026): equilibrium charge fraction in the condensed phase along the history.

At fixed hidden charge-to-entropy ratio Y_h = n_Q/s_h (reversible equilibrium, hidden sector decoupled), for each
T_h find the equilibrium state:
  two-phase:   vapour at mu_c(T) (bulk Q-matter coexistence; thermo.mu_coex), liquid volume fraction phi:
               n_Q = (1-phi) n_v + phi n_B,  s = (1-phi) s_v + phi s_liq,  n_Q = Y_h s  ->  phi;  eta = phi n_B / n_Q
  vapour-only: if phi <= 0: solve n_net(mu, T) = Y_h s_vap(mu, T) for mu < mu_c  (eta = 0)
  liquid-only: phi >= 1 (homogeneous condensate; the T=0 bulk density is used, so this is indicative only).
Also the finite-Q Saha version (objects Q_min..Q_max with g_int = 1e-10, 1, 1e10) for T <= 3, and the visible clock:
s_vis = n_Q/Y_Q (n_Q in GeV^3), T_vis from s_vis = (2 pi^2/45) g_s(T_vis) T_vis^3, H = sqrt(pi^2 g*/90) T_vis^2/M_P.
Y_h for each label is read from tracks/fx_chain/results/heat_onset.json (onset and spinodal crossing, massless 3-dof
entropy: an [E] bracket).
"""
import json, math, os, sys, time
import numpy as np
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import thermo_copy as th

MX = 0.0924089637375362
YQ = 5.07941423e-9
MPL = 2.435e18
FXC = os.path.join(HERE, '..', '..', 'fx_chain', 'results', 'heat_onset.json')

# SM g_*(T) (GeV), approximate [R]: Husdal 2016 / Saikawa-Shirai 2018 shapes, smoothed in log T
_GT = np.array([[1e-4, 3.36], [1e-3, 10.75], [0.02, 10.75], [0.06, 14.0], [0.1, 17.3], [0.15, 25.0], [0.2, 45.0],
                [0.3, 61.75], [1.0, 72.0], [2.0, 77.0], [5.0, 84.0], [10.0, 86.25], [30.0, 86.5], [60.0, 90.0],
                [100.0, 96.0], [200.0, 103.0], [1e3, 106.75], [1e6, 106.75]])


def gstar(T):
    return float(np.interp(math.log(T), np.log(_GT[:, 0]), _GT[:, 1]))


def T_vis_of(nQ_mx3, YQ=YQ):
    s = nQ_mx3 * MX ** 3 / YQ
    f = lambda lt: math.log(2 * math.pi ** 2 / 45 * gstar(math.exp(lt)) * math.exp(3 * lt)) - math.log(s)
    return math.exp(brentq(f, math.log(1e-15), math.log(1e8)))


T_EQ = 0.8e-9   # GeV, matter-radiation equality (approximate, [R]); H^2 = H_rad^2 (1 + T_eq/T) below it


def hubble(Tv):
    return math.sqrt(math.pi ** 2 * gstar(Tv) / 90) * Tv ** 2 / MPL * math.sqrt(1 + T_EQ / Tv)


def yh_labels():
    d = json.load(open(FXC))
    out = {}
    hist = {h['label']: h for h in d['history_M1']}
    for c in d['chain']:
        if not c['phi_in_equilibrium']:
            continue
        o = c['onset_untilted']
        y_on = o['nQ_over_mX3'] / (2 * math.pi ** 2 / 45 * 3 * o['T_h_over_mX'] ** 3)
        y_sp = None
        if c['label'] in hist:
            sp = hist[c['label']]['crossings']['spinodal_233']
            y_sp = sp['n'] / (2 * math.pi ** 2 / 45 * 3 * sp['T_h_over_mX'] ** 3)
        out[c['label']] = dict(R=c['R'], Y_onset=y_on, Y_spinodal=y_sp, T_h_onset=o['T_h_over_mX'],
                               T_h_spinodal=hist[c['label']]['crossings']['spinodal_233']['T_h_over_mX'] if c['label'] in hist else None)
    return out


def state(T, Yh, thermal_liquid=True):
    mc = th.mu_coex(T, thermal_liquid)
    v = th.vapour(mc, T)
    pl, sl = th.liquid_thermal(T)
    nB = th.N_B
    num = Yh * v['s'] - v['n_net']
    den = nB - v['n_net'] - Yh * (sl - v['s'])
    phi = num / den
    if 0 < phi < 1:
        nQ = (1 - phi) * v['n_net'] + phi * nB
        s = (1 - phi) * v['s'] + phi * sl
        return dict(T=T, phase='two-phase', mu=mc, phi=phi, eta=phi * nB / nQ, nQ=nQ, s=s, n_vap=v['n_net'],
                    vap_ratio=v['n_net'] / v['s'], Xbar_per_charge=(1 - phi) * v['n_Xbar'] / nQ)
    if phi >= 1:
        # homogeneous condensate denser than the T = 0 bulk: indicative only (n_Q = Y_h s_liq, s from the bulk spectrum)
        return dict(T=T, phase='liquid-only', mu=None, phi=1.0, eta=1.0, nQ=Yh * sl, s=sl, n_vap=0.0, vap_ratio=v['n_net'] / v['s'],
                    Xbar_per_charge=0.0)
    # vapour only
    g = lambda mu: th.vapour(mu, T)['n_net'] - Yh * th.vapour(mu, T)['s']
    mu = brentq(g, 1e-9, mc, xtol=1e-13)
    w = th.vapour(mu, T)
    return dict(T=T, phase='vapour-only', mu=mu, phi=0.0, eta=0.0, nQ=w['n_net'], s=w['s'], n_vap=w['n_net'],
                vap_ratio=v['n_net'] / v['s'], Xbar_per_charge=w['n_Xbar'] / w['n_net'], mu_c=mc)


def saha_state(fam, T, Yh, Qmax, gint):
    """Finite-Q Saha at fixed Y_h (objects' own entropy neglected; dilute objects): solve
    n_vap(mu) + Sum Q n(Q; mu) = Y_h s_vap(mu) for mu."""
    def resid(mu):
        v = th.vapour(mu, T)
        lc, _ = th.saha_objects(fam, mu, T, Qmax, gint)
        return math.log(v['n_net'] + math.exp(min(lc, 700))) - math.log(Yh * v['s'])
    lo, hi = 1e-6, 0.99999
    if resid(lo) > 0 or resid(hi) < 0:
        return None
    mu = brentq(resid, lo, hi, xtol=1e-14)
    v = th.vapour(mu, T)
    lc, ln = th.saha_objects(fam, mu, T, Qmax, gint)
    nQ = v['n_net'] + math.exp(lc)
    return dict(T=T, mu=mu, eta=math.exp(lc) / nQ, mean_Q=math.exp(lc - ln), nQ=nQ)


def main(out, part):
    t0 = time.time()
    labs = yh_labels()
    for k, v in labs.items():
        print('%-36s R=%7.3f Y_h(onset)=%.4f Y_h(spinodal)=%s' % (k, v['R'], v['Y_onset'], '%.4f' % v['Y_spinodal'] if v['Y_spinodal'] else '-'))
    # the vapour's maximal charge-to-entropy ratio at coexistence (objects survive iff Y_h exceeds it)
    Ts = np.geomspace(0.01, 12.0, 120)
    ratio = []
    for T in Ts:
        mc = th.mu_coex(T)
        v = th.vapour(mc, T)
        ratio.append(dict(T=float(T), mu_c=mc, n_vap=v['n_net'], s_vap=v['s'], ratio=v['n_net'] / v['s'],
                          Xbar_frac=v['n_Xbar'] / v['n_X']))
    rr = np.array([r['ratio'] for r in ratio]); i = int(np.argmax(rr))
    ycrit = dict(T_at_max=float(Ts[i]), Y_crit=float(rr[i]))
    print('vapour n/s at mu_c: max %.4f at T = %.3f' % (rr[i], Ts[i]))
    # same for finite reference objects: coexistence with an object of charge Q (mu = Omega(Q)) -- Gibbs-Thomson
    fam = th.Family()
    ycrit_Q = {}
    for Q, Om in ((2e5, 0.8817204919524271), (1e6, 0.83429668), (1e8, 0.78777238)):
        best = 0; Tb = None
        for T in np.geomspace(0.02, 5.0, 80):
            v = th.vapour(Om, T); r = v['n_net'] / v['s']
            if r > best:
                best, Tb = r, T
        ycrit_Q['%g' % Q] = dict(Omega=Om, Y_crit=best, T_at_max=float(Tb))
        print('  coexistence with Q=%g objects (mu = %.4f): Y_crit = %.4f at T = %.3f' % (Q, Om, best, Tb))
    Yset = {'R0_heavy': labs['heavy_release:R0'], 'med_e1e6_heavy': labs['heavy_release:med_e1e6'],
            'tot_e1e6_band_heavy': labs['heavy_release:tot_e1e6_band_4.97'],
            'tot_e1e6_band_pilot': labs['pilot_state:tot_e1e6_band_4.97'],
            'tot_e1e7_band_heavy': labs['heavy_release:tot_e1e7_band_41.7'], 'R0_pilot': labs['pilot_state:R0']}
    hist = {}
    if part != 'hist':
        Yset_iter = {}
    else:
        Yset_iter = Yset
    Tg = np.geomspace(0.001, 14.0, 70)[::-1]
    for name, L in Yset_iter.items():
        for which in ('Y_onset', 'Y_spinodal'):
            Yh = L[which]
            if Yh is None:
                continue
            rows = []
            for T in Tg:
                st = state(float(T), Yh)
                if not (st['nQ'] > 0):
                    print('bad state', name, which, T, st); continue
                try:
                    st['T_vis_GeV'] = T_vis_of(st['nQ']); st['H_over_mX'] = hubble(st['T_vis_GeV']) / MX
                except ValueError:
                    print('T_vis failed', name, which, T, st, flush=True); continue
                rows.append(st)
            # T_ev (two-phase -> vapour-only going down) and T_s (vapour-only -> two-phase)
            ph = [r['phase'] for r in rows]
            T_ev = T_s = None
            for j in range(1, len(rows)):
                if ph[j - 1] == 'two-phase' and ph[j] == 'vapour-only' and T_ev is None:
                    T_ev = math.sqrt(rows[j - 1]['T'] * rows[j]['T'])
                if ph[j - 1] == 'vapour-only' and ph[j] == 'two-phase':
                    T_s = math.sqrt(rows[j - 1]['T'] * rows[j]['T'])
            etamin = min(r['eta'] for r in rows if r['T'] < 6)
            hist['%s|%s' % (name, which)] = dict(Y_h=Yh, R=L['R'], T_ev=T_ev, T_s=T_s, eta_min_below_T6=etamin, rows=rows)
            print('%-26s %-11s Y_h=%.4f  T_ev=%s T_s=%s  eta_min=%.4f  t=%.0fs' % (name, which, Yh, T_ev, T_s, etamin, time.time() - t0), flush=True)
            for r in rows[::10]:
                print('     T=%7.4f %-11s eta=%.4f mu=%s nQ=%.3e Tvis=%.3g GeV H/m=%.2e' % (r['T'], r['phase'], r['eta'], '%.4f' % r['mu'] if r['mu'] else '-', r['nQ'], r['T_vis_GeV'], r['H_over_mX']))
    # finite-Q Saha histories for the e1e6 band (heavy, onset) and R0, T <= 3
    saha = {}
    for name in (('tot_e1e6_band_heavy', 'R0_heavy', 'med_e1e6_heavy') if part == 'saha' else ()):
        Yh = Yset[name]['Y_onset']
        for Qmax, gint in ((1e8, 1e-10), (1e8, 1.0), (1e8, 1e10), (1e6, 1.0)):
            if True:
                rows = []
                for T in np.geomspace(0.01, 3.0, 25)[::-1]:
                    r = saha_state(fam, float(T), Yh, Qmax, gint)
                    if r is not None:
                        rows.append(r)
                saha['%s|Qmax%g|g%g' % (name, Qmax, gint)] = rows
                # crossing eta = 1/2
                cr = []
                for j in range(1, len(rows)):
                    if (rows[j - 1]['eta'] - 0.5) * (rows[j]['eta'] - 0.5) < 0:
                        cr.append(math.sqrt(rows[j - 1]['T'] * rows[j]['T']))
                print('Saha %-22s Qmax=%g g=%g: eta=1/2 crossings at T = %s; eta(T=0.01)=%.4f eta(T=3)=%.4f t=%.0fs' % (
                    name, Qmax, gint, ['%.4f' % c for c in cr], rows[-1]['eta'], rows[0]['eta'], time.time() - t0), flush=True)
    json.dump(dict(Y_labels=labs, vapour_ratio_at_coexistence=ratio, Y_crit=ycrit, Y_crit_finite_Q=ycrit_Q,
                   histories=hist, saha=saha, runtime_s=time.time() - t0), open(out, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
