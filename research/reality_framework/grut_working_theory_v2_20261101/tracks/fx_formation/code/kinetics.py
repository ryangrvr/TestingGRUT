"""fx_formation item 2-3 (Claude, 1 Nov 2026): is equilibrium maintained, and what sets the object size?  [E]

Reads results/eq_history.json (equilibrium histories at fixed Y_h, with T_vis and H) and results/family.json.
All rates in m_X units.  Ingredients (each an estimate):
  * object radius: bulk-density sphere R = (3Q/4 pi n_B)^{1/3} for Q > 1e8, r90 from the family below; geometric
    absorption cross section pi R^2 times a sticking probability s_stick (1 and 0.01);
  * vapour X mean speed vbar = <p/E> of the Bose distribution (thermo.vapour gives n; vbar from sqrt(8T/pi) NR,
    capped at 1);
  * exchange rate per object charge:  G_ex(Q) = pi R^2 vbar n_X s_stick / Q  (charge turnover of one object);
  * capture rate per vapour X:  G_c = n_obj pi R^2 vbar s_stick, n_obj = eta n_Q / Q;
  * evaporation time of an object at saturation: t_ev = Q / (pi R^2 vbar n_X s_stick);
  * coagulation: Fuchs kernel K = K_c K_fm / (K_c + K_fm), K_c = 16 pi R D, D = T/zeta, zeta = (4 pi/3) R^2 (e+P)_vap vbar,
    K_fm = pi (2R)^2 sqrt(2) v_th, v_th = sqrt(8 T / (pi M)); mean-field dQ/dt = (1/2) s_m K eta n_Q;
  * Ostwald ripening (Wagner 1961, interface-limited, ballistic vapour): d<R^2>/dt = (64/81) k_W,
    k_W = (s_stick vbar/4)(n_X/n_B)(2 sigma/(n_B T)), sigma = 98.8 m_X^3 (from the E(Q) fit);
  * slow passage through the spinodal: Gamma_peak ~ c1 delta, k_peak ~ c2 sqrt(delta) (fitted here from the
    cold-fluid cubic, coldfluid_m1_copy.py), delta = 3 H t, amplification A = c1 delta^2/(6 H) reaches A_* in 10..60;
    Q_frag = n_sp (2 pi / k_peak)^3 at delta_*;
  * cavitation (thin wall) barrier in the stretched liquid between 233 and 449 m_X^3: 16 pi sigma^3 / (3 P^2).
  * nucleation of new objects from a supersaturated vapour: barrier E(Q*) - mu Q*, Omega(Q*) = mu (critical object);
    for mu > Omega(Q_min) the barrier is at least E(Q_min) - mu Q_min.
"""
import json, math, os, sys, time
import numpy as np
from scipy.optimize import brentq, minimize_scalar
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import thermo as th
import coldfluid_m1_copy as cf
import eq_history as eh

RES = os.path.join(HERE, '..', 'results')
NB = th.N_B
SIGMA = 98.8


def radius(fam, Q):
    return float(fam.r90(Q)) if Q <= 1e8 else (3 * Q / (4 * math.pi * NB)) ** (1 / 3)


def vbar(T, m=1.0):
    return min(1.0, math.sqrt(8 * T / (math.pi * m)) if T < 0.5 else math.sqrt(8 * T / (math.pi * m)) / math.sqrt(1 + 8 * T / (math.pi * m)))


def hist_arrays(H):
    rows = [r for r in H['rows'] if 'H_over_mX' in r]
    rows.sort(key=lambda r: -r['T'])
    T = np.array([r['T'] for r in rows]); nQ = np.array([r['nQ'] for r in rows]); eta = np.array([r['eta'] for r in rows])
    Hh = np.array([r['H_over_mX'] for r in rows]); ph = [r['phase'] for r in rows]
    mu = [r['mu'] for r in rows]
    nX, eP = [], []
    for r in rows:
        if r['mu'] is None:
            nX.append(0.0); eP.append(0.0); continue
        v = th.vapour(r['mu'], r['T'])
        nX.append(v['n_X']); eP.append(v['e'] + v['P'])
    lna = -np.log(nQ) / 3
    return dict(T=T, nQ=nQ, eta=eta, H=Hh, phase=ph, mu=mu, nX=np.array(nX), eP=np.array(eP), lna=lna,
                Tvis=np.array([r['T_vis_GeV'] for r in rows]))


def spinodal_passage():
    """c1, c2 near the T=0 spinodal and Q_frag(A_*) with H at n_Q = 233.2 (T_vis from Y_Q)."""
    Fsp = minimize_scalar(lambda F: math.sqrt(cf.S_of_F(F) ** 2 + cf.a * F * F), bracket=(1.7, 1.85, 2.0), tol=1e-13).x
    def st(F):
        S = cf.S_of_F(F); Om = math.sqrt(S * S + cf.a * F * F); D = cf.Uphi_pp(S) + F * F
        return S, Om, D, 100 * Om * F * F
    nsp = st(Fsp)[3]
    rows = []
    for dl in (1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2, 0.4, 0.6):
        F = brentq(lambda F: st(F)[3] - nsp * (1 - dl), 0.05, Fsp, xtol=1e-15)
        S, Om, D, n = st(F)
        ks = np.geomspace(1e-6, 2.0, 500)
        g = np.array([cf.growth(k, F, S, Om, D)[0] for k in ks])
        i = int(np.argmax(g))
        r = minimize_scalar(lambda lk: -cf.growth(math.exp(lk), F, S, Om, D)[0],
                            bracket=(math.log(ks[max(i - 1, 0)]), math.log(ks[i]), math.log(ks[min(i + 1, 499)])), tol=1e-10)
        kp = math.exp(r.x); gp = cf.growth(kp, F, S, Om, D)[0]
        rows.append(dict(delta=dl, n=n, F=F, k_peak=kp, Gamma_peak=gp, Q_cube=n * (2 * math.pi / kp) ** 3))
    small = [r for r in rows if r['delta'] <= 1e-4]
    c1 = float(np.exp(np.mean([math.log(r['Gamma_peak'] / r['delta']) for r in small])))
    c2 = float(np.exp(np.mean([math.log(r['k_peak'] / math.sqrt(r['delta'])) for r in small])))
    Tv = eh.T_vis_of(nsp); Hsp = eh.hubble(Tv) / th.MX if hasattr(th, 'MX') else eh.hubble(Tv) / eh.MX
    cases = []
    for A in (10.0, 20.0, 60.0):
        dstar = math.sqrt(6 * Hsp * A / c1)
        kst = c2 * math.sqrt(dstar)
        cases.append(dict(A=A, delta_star=dstar, k_star=kst, Gamma_star=c1 * dstar, Gamma_over_H=c1 * dstar / Hsp,
                          Q_cube=nsp * (2 * math.pi / kst) ** 3, Q_half_cube=nsp * (math.pi / kst) ** 3))
    return dict(F_sp=Fsp, n_sp=nsp, T_vis_GeV=Tv, H_over_mX=Hsp, c1=c1, c2=c2, table=rows, cases=cases)


def cavitation():
    """Thin-wall barrier for a vacuum (vapour) bubble in the stretched T=0 liquid at n between n_sp and n_B."""
    out = []
    for n in (440., 420., 400., 350., 300., 260., 240., 234.):
        F = brentq(lambda F: 100 * math.sqrt(cf.S_of_F(F) ** 2 + cf.a * F * F) * F * F - n, 1.0, 3.0, xtol=1e-14)
        S = cf.S_of_F(F); Om = math.sqrt(S * S + cf.a * F * F)
        e = 100 * (0.5 * Om * Om * F * F + th.Uphi(S) + 0.5 * S * S * F * F + cf.a / 4 * F ** 4)
        P = n * Om - e
        dF = 16 * math.pi * SIGMA ** 3 / (3 * P * P) if P < 0 else float('inf')
        out.append(dict(n=n, P=P, R_crit=2 * SIGMA / abs(P) if P < 0 else None, barrier_thin_wall=dF))
    return out


def nucleation_barrier(fam, mu):
    """Barrier for forming a critical object from vapour at chemical potential mu."""
    Qg = np.geomspace(fam.Qmin * 1.000001, 1e12, 4000)
    E = fam.e_over_q(Qg) * Qg
    f = E - mu * Qg
    i = int(np.argmax(f[: np.argmin(f) + 1])) if np.any(np.diff(f) < 0) else len(f) - 1
    # for mu > Omega(Q_min) f decreases from Q_min on: barrier = f(Q_min)
    return float(np.max(f[: max(i, 1) + 1])), float(Qg[i])


def main(out):
    t0 = time.time()
    fam = th.Family()
    EH = json.load(open(os.path.join(RES, 'eq_history.json')))
    res = dict(notes='all [E]; m_X units; see module docstring')
    res['spinodal_passage'] = sp = spinodal_passage()
    print('spinodal: n=%.3f c1=%.4f c2=%.4f T_vis=%.1f GeV H/m=%.3e' % (sp['n_sp'], sp['c1'], sp['c2'], sp['T_vis_GeV'], sp['H_over_mX']))
    for c in sp['cases']:
        print('   A=%g delta*=%.2e k*=%.3e Gamma*/H=%.2e Q_cube=%.2e Q_half=%.2e' % (c['A'], c['delta_star'], c['k_star'], c['Gamma_over_H'], c['Q_cube'], c['Q_half_cube']))
    for r in sp['table']:
        print('   delta=%.0e n=%.2f k_peak=%.4e Gamma=%.4e Q_cube=%.3e' % (r['delta'], r['n'], r['k_peak'], r['Gamma_peak'], r['Q_cube']))
    res['cavitation'] = cav = cavitation()
    for c in cav:
        print('   cavitation n=%.0f P=%.2f R_c=%s barrier=%.3g' % (c['n'], c['P'], '%.2f' % c['R_crit'] if c['R_crit'] else '-', c['barrier_thin_wall']))
    Qlist = [3e5, 1e8, 1e14, 1e20, 1e26, 1e32, 1e38]
    res['histories'] = {}
    for key in ('R0_heavy|Y_onset', 'med_e1e6_heavy|Y_onset', 'tot_e1e6_band_heavy|Y_onset', 'tot_e1e6_band_pilot|Y_onset'):
        A = hist_arrays(EH['histories'][key])
        T, nQ, eta, Hh, lna = A['T'], A['nQ'], A['eta'], A['H'], A['lna']
        two = np.array([p == 'two-phase' for p in A['phase']])
        out_h = dict(Y_h=EH['histories'][key]['Y_h'], T_ev=EH['histories'][key]['T_ev'], T_s=EH['histories'][key]['T_s'])
        # (1) exchange / evaporation rates at representative T
        ex = []
        for j in range(len(T)):
            if A['mu'][j] is None:
                continue
            vb = vbar(T[j])
            row = dict(T=T[j], T_vis_GeV=A['Tvis'][j], H=Hh[j], eta=eta[j], phase=A['phase'][j], n_X=A['nX'][j])
            for Q in Qlist:
                R = radius(fam, Q)
                Gex = math.pi * R * R * vb * A['nX'][j] / Q
                Gc = (eta[j] * nQ[j] / Q) * math.pi * R * R * vb
                row['Gex_over_H_Q%.0e' % Q] = Gex / Hh[j]
                row['Gcapt_over_H_Q%.0e' % Q] = Gc / Hh[j]
            ex.append(row)
        out_h['rates'] = ex
        # (2) capture freeze-out of the vapour onto objects of fixed <Q> (two-phase branch below T = 0.5):
        #     d N_v / d ln a = -(G_c/H) (N_v - N_eq);  N in units of the comoving charge; N_eq = 1 - eta_eq
        fr = {}
        jlo = [j for j in range(len(T)) if T[j] < 0.6 and two[j]]
        if jlo:
            for s_st in (1.0, 0.01):
                for Q in Qlist:
                    R = radius(fam, Q)
                    Nv = 1 - eta[jlo[0]]
                    Tf = None
                    for j0, j1 in zip(jlo[:-1], jlo[1:]):
                        Gc_H = (eta[j1] * nQ[j1] / Q) * math.pi * R * R * vbar(T[j1]) * s_st / Hh[j1]
                        dl = lna[j1] - lna[j0]
                        Neq = 1 - eta[j1]
                        # exact relaxation over the step with frozen rate
                        Nv = Neq + (Nv - Neq) * math.exp(-Gc_H * dl) if Gc_H * dl < 700 else Neq
                        if Tf is None and Gc_H < 1:
                            Tf = T[j1]
                    fr['s%g_Q%.0e' % (s_st, Q)] = dict(one_minus_eta_final=Nv, T_freeze=Tf, T_end=float(T[jlo[-1]]))
        out_h['capture_freezeout'] = fr
        # (3) evaporation feasibility (only when a vapour-only window exists): t_ev/t_H at the last two-phase point
        if out_h['T_ev']:
            jev = max(j for j in range(len(T)) if T[j] >= out_h['T_ev'] and two[j])
            Qev = {}
            for Q in Qlist:
                R = radius(fam, Q)
                tev = Q / (math.pi * R * R * vbar(T[jev]) * A['nX'][jev])
                Qev['%.0e' % Q] = dict(t_ev_H=tev * Hh[jev], t_ev_H_stick1e_2=tev * Hh[jev] * 100)
            # largest Q evaporating within one Hubble time at T_ev (s = 1): solve t_ev H = 1 with bulk radius
            # t_ev = (4/3) R n_B / (vbar n_X)
            Rmax = Hh[jev] ** -1 * vbar(T[jev]) * A['nX'][jev] * 3 / (4 * NB)
            Qev['Q_max_evaporating_in_1_over_H'] = NB * 4 * math.pi / 3 * Rmax ** 3
            Qev['charge_per_Hubble_volume'] = nQ[jev] * (1 / Hh[jev]) ** 3
            Qev['T'] = T[jev]
            out_h['evaporation'] = Qev
        # (4) growth after formation: coagulation (Fuchs) and Ostwald (Wagner), two-phase stretch from the spinodal down
        jtw = [j for j in range(len(T)) if two[j]]
        grow = {}
        for Q0 in (3e5, 1e12):
            for s_st, s_m in ((1.0, 1.0), (0.01, 0.1)):
                Qc = Q0; R2 = radius(fam, Q0) ** 2
                hist_g = []
                for j0, j1 in zip(jtw[:-1], jtw[1:]):
                    if A['phase'][j0] != 'two-phase' or A['phase'][j1] != 'two-phase':
                        continue
                    dt = (lna[j1] - lna[j0]) / Hh[j1]
                    Tt = T[j1]; nX = A['nX'][j1]; eP = A['eP'][j1]
                    # coagulation
                    R = radius(fam, Qc); M = 0.8 * Qc
                    vth = math.sqrt(8 * Tt / (math.pi * M))
                    zeta = 4 * math.pi / 3 * R * R * eP * vbar(Tt)
                    Dd = Tt / zeta
                    Kc = 16 * math.pi * R * Dd; Kfm = math.pi * (2 * R) ** 2 * math.sqrt(2) * vth
                    K = Kc * Kfm / (Kc + Kfm)
                    # integrate dQ/dt = 0.5 s_m K eta nQ with K ~ Q^p locally: use small substeps in log
                    rate = 0.5 * s_m * K * eta[j1] * nQ[j1]
                    Qc = Qc + rate * dt if rate * dt < 10 * Qc else Qc * (1 + rate * dt / Qc) ** 0.9
                    # Ostwald
                    kW = s_st * vbar(Tt) / 4 * (nX / NB) * (2 * SIGMA / (NB * Tt))
                    R2 = R2 + 64 / 81 * kW * dt
                    hist_g.append(dict(T=Tt, Q_coag=Qc, Q_ostwald=NB * 4 * math.pi / 3 * R2 ** 1.5, K=K, Kc=Kc, Kfm=Kfm,
                                       coag_rate_over_H=rate / Qc / Hh[j1], kW=kW))
                grow['Q0_%.0e_s%g_sm%g' % (Q0, s_st, s_m)] = hist_g
        out_h['growth'] = grow
        res['histories'][key] = out_h
        print('%s: Y_h=%.4f T_ev=%s T_s=%s' % (key, out_h['Y_h'], out_h['T_ev'], out_h['T_s']))
        for r in ex[::8]:
            print('   T=%.4f Tvis=%.3g GeV %s eta=%.3f nX=%.2e  Gex/H(3e5,1e20,1e32)=%.1e %.1e %.1e  Gc/H(3e5,1e20,1e32)=%.1e %.1e %.1e' % (
                r['T'], r['T_vis_GeV'], r['phase'], r['eta'], r['n_X'], r['Gex_over_H_Q3e+05'], r['Gex_over_H_Q1e+20'], r['Gex_over_H_Q1e+32'],
                r['Gcapt_over_H_Q3e+05'], r['Gcapt_over_H_Q1e+20'], r['Gcapt_over_H_Q1e+32']))
        for k, v in fr.items():
            print('   freeze %-14s 1-eta=%.4f T_f=%s' % (k, v['one_minus_eta_final'], v['T_freeze']))
        if 'evaporation' in out_h:
            e = out_h['evaporation']
            print('   evaporation at T=%.3f: Q evaporating within 1/H up to %.2e; charge per Hubble volume %.2e' % (e['T'], e['Q_max_evaporating_in_1_over_H'], e['charge_per_Hubble_volume']))
        for k, g in grow.items():
            if g:
                print('   growth %-22s end T=%.3f Q_coag=%.2e Q_ostwald=%.2e (coag rate/H at end %.1e)' % (k, g[-1]['T'], g[-1]['Q_coag'], g[-1]['Q_ostwald'], g[-1]['coag_rate_over_H']))
    # (5) nucleation from a supersaturated vapour (the e1e6 band history continued below T_s with no objects)
    nuc = []
    for key in ('tot_e1e6_band_heavy|Y_onset', 'tot_e1e6_band_pilot|Y_onset'):
        Yh = EH['histories'][key]['Y_h']
        for T in (0.08, 0.05, 0.03, 0.02, 0.01, 0.005, 0.002):
            g = lambda mu: th.vapour(mu, T)['n_net'] - Yh * th.vapour(mu, T)['s']
            try:
                mu = brentq(g, 1e-6, 0.999999, xtol=1e-14)
            except ValueError:
                nuc.append(dict(history=key, T=T, note='no vapour solution (BEC)')); continue
            v = th.vapour(mu, T)
            bar, Qs = nucleation_barrier(fam, mu)
            fug = math.exp((mu - 1) / T)
            nuc.append(dict(history=key, T=T, mu=mu, fugacity=fug, barrier=bar, Q_crit=Qs, barrier_over_T=bar / T,
                            supersaturation_ln=(mu - th.mu_coex(T)) / T))
            print('   supersaturated vapour %s T=%.3f mu=%.5f fugacity e^{(mu-m)/T}=%.3e  barrier=%.1f m_X (Q*=%.3g)  barrier/T=%.2e' % (
                key.split('|')[0], T, mu, fug, bar, Qs, bar / T))
    res['nucleation_from_vapour'] = nuc
    res['runtime_s'] = time.time() - t0
    json.dump(res, open(out, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main(sys.argv[1])
