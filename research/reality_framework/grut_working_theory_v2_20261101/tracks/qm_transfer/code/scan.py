"""qm_transfer production scan (own driver, Claude 2 Nov 2026).  [C]

Families (Planck 2018 background, K_rsa = 200, K_end = 400, Nq = 24, L = 10):
  P : physical M1r two-phase model, g = 1, 'thermo' kinetics, s = 1:  S/Q in {0.5, 1.235, 3.4} x zeta grid
  S : the same with 'sudden' capture freeze-out (fx_formation picture), S/Q = 1.235
  E : s = 0.01 (elastic (1-s) pi R^2), 'thermo', S/Q = 1.235, zeta <= 2.5e-4 (Q >= Q_min)
  X : parametric: target vapour charge fraction f_c at freeze-out, target freeze-out temperature T_f, zeta (s = 1):
      the equilibrium history is deformed by a saturation factor g and an adiabat S/Q chosen so that the
      two-phase equilibrium gives N_v = f_c, T = T_f where Gamma_a(vbar)/H = 1 (then the full kinetics runs).
usage: scan.py OUTDIR worker_index n_workers [time_budget_s] [families]
Writes OUTDIR/model_<name>.json per model (restartable) and OUTDIR/cdm_ref.json.
"""
import json, math, os, sys, time
import numpy as np
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from qt_bg import Cosmo
from qt_hidden import Vapour, Frozen, KGEO, MX, RHOC_GEV4_H2, Omega_fit
from qt_boltz import Mode, CDM

KS = [1.5, 3.0, 5.0, 8.0, 12.0, 17.0, 24.0, 34.0, 48.0, 68.0, 96.0, 136.0]
A_OUT = (1.0, 0.25)
ZETAS = [1e-7, 1e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2]
KW = dict(K_rsa=200.0, K_end=400.0, a_out=A_OUT, rtol=1e-7)


def specs():
    out = []
    for SQ in (1.235, 0.5, 3.4):
        for z in ZETAS:
            out.append(dict(name='P_SQ%g_z%.0e' % (SQ, z), fam='P', mode='thermo', SQ=SQ, zeta=z, s=1.0, g=1.0))
    for z in ZETAS:
        out.append(dict(name='S_SQ1.235_z%.0e' % z, fam='S', mode='sudden', SQ=1.235, zeta=z, s=1.0, g=1.0))
    for z in (1e-7, 1e-6, 1e-5, 3e-5, 1e-4, 2.5e-4):
        out.append(dict(name='E_s0.01_z%.0e' % z, fam='E', mode='thermo', SQ=1.235, zeta=z, s=0.01, g=1.0))
    for Tf in (0.004, 0.010):
        for fc in (0.03, 0.1, 0.35):
            for z in ZETAS:
                out.append(dict(name='X_fc%g_Tf%g_z%.0e' % (fc, Tf, z), fam='X', mode='thermo', fc=fc, Tf=Tf, zeta=z, s=1.0))
    for m in (5.3, 6.5, 3.5, 5.7, 10.0, 20.0):
        out.append(dict(name='W_m%g' % m, fam='W', m_keV=m))
    # tracks/qm_micro evaluation points (NOTES.md section 3 table): base = computed sticking s(E_k), perfect trapping;
    # bel = Beliaev-trapping variant, approximated by a constant trapping factor on s(E_k) (s_eff(T_f) ~ 0.04-0.07)
    for SQ, Q in ((1.5, 5.5e19), (2.5, 1.2e20), (3.2, 2.7e20)):
        out.append(dict(name='M_base_SQ%g' % SQ, fam='M', mode='thermo', SQ=SQ, Q=Q, s_scale=1.0))
    out.append(dict(name='M_base_SQ3.2_sudden', fam='M', mode='sudden', SQ=3.2, Q=2.7e20, s_scale=1.0))
    for SQ, Q, sc in ((1.5, 2.3e20, 0.08), (3.2, 1.6e21, 0.15)):
        out.append(dict(name='M_bel_SQ%g' % SQ, fam='M', mode='thermo', SQ=SQ, Q=Q, s_scale=sc))
    # reorder: M before the long X family
    xs = [o for o in out if o['fam'] == 'X']; rest = [o for o in out if o['fam'] != 'X']
    return rest + xs


def s_of_E_qm_micro():
    """angle-averaged sticking s(E_k) of tracks/qm_micro/results/stick_thermal.json (planar M1 interface, tree level);
    below the table: threshold law s = c sqrt(E_k) with the table's c_angle_avg."""
    d = json.load(open(os.path.join(HERE, '..', '..', 'qm_micro', 'results', 'stick_thermal.json')))
    E = np.array(d['T_kin']); S = np.array(d['s_X_angle_avg'])
    cth = d['threshold_law']['c_angle_avg']
    lE, lS = np.log(E), np.log(np.maximum(S, 1e-300))
    def f(Ek):
        Ek = np.maximum(np.asarray(Ek, float), 1e-30)
        out = np.exp(np.interp(np.log(Ek), lE, lS))
        lo = Ek < E[0]
        out = np.where(lo, cth * np.sqrt(Ek), out)
        return np.where(Ek > E[-1], S[-1], out)
    return f


def invert(c, fc, Tf, zeta, s=1.0):
    """(f_c, T_f, zeta) -> (S/Q, g): the deformed two-phase equilibrium freezes out at N_v = f_c, T = T_f."""
    Q = (s / zeta)**3
    D = 1 - Omega_fit(Q)
    NQ0 = c.omc * RHOC_GEV4_H2 / (MX**4 * ((1 - D) + fc * D))
    vb = math.sqrt(8 * Tf / math.pi)
    from qt_bg import GEV_TO_INVMPC
    def F(la):
        a = math.exp(la)
        H = c.calH(a) / a / (MX * GEV_TO_INVMPC)
        return math.log(zeta * KGEO * NQ0 / a**3 * (1 - fc)**(2 / 3) * vb / H)
    la = brentq(F, math.log(1e-16), 0.0)
    af = math.exp(la)
    nQ = NQ0 / af**3
    g = fc * nQ / ((Tf / (2 * math.pi))**1.5 * math.exp(-D / Tf))
    SQ = fc * (2.5 + D / Tf)
    return dict(Q=Q, Delta=D, a_f_target=af, g=g, SQ=SQ, NQ0=NQ0)


def build(c, sp):
    if sp['fam'] == 'W':     # pure thermal relic (FD, g = 2), own-code reference curve
        Tnu = (4 / 11)**(1 / 3) * c.T0 * 8.617333262e-5
        Tq = Tnu * (c.omc * 94.0 / (sp['m_keV'] * 1e3))**(1 / 3) / (sp['m_keV'] * 1e3)
        fr = Frozen(c, fraction=0.999, dist='FD', Tq=Tq, Nq=28, Lv=14, a_lo=1e-11)
        cW = Cosmo(h=c.h, omb=c.omb, omc=c.omc, Neff=c.Neff, dm_efac=fr.efac)
        fr.c = cW
        return fr, dict(Tq=Tq, cos=cW)
    if sp['fam'] == 'M':
        v = Vapour(c, SQ=sp['SQ'], Q=sp['Q'], s=1.0, mode=sp['mode'], Nq=24, Lv=10, s_of_E=s_of_E_qm_micro(), s_scale=sp['s_scale'])
        return v, {}
    s = sp['s']
    Q = (s / sp['zeta'])**3
    extra = {}
    if sp['fam'] == 'X':
        inv = invert(c, sp['fc'], sp['Tf'], sp['zeta'], s)
        extra = inv
        v = Vapour(c, SQ=inv['SQ'], Q=Q, s=s, g=inv['g'], mode=sp['mode'], Nq=24, Lv=10)
    else:
        v = Vapour(c, SQ=sp['SQ'], Q=Q, s=s, g=sp['g'], mode=sp['mode'], Nq=24, Lv=10)
    return v, extra


def main():
    outdir = sys.argv[1]; wi = int(sys.argv[2]); nw = int(sys.argv[3])
    budget = float(sys.argv[4]) if len(sys.argv) > 4 else 560.0
    fams = sys.argv[5].split(',') if len(sys.argv) > 5 else ['P', 'S', 'E', 'W', 'M', 'X']
    t_start = time.time()
    os.makedirs(outdir, exist_ok=True)
    c = Cosmo()
    ref_fn = os.path.join(outdir, 'cdm_ref.json')
    if os.path.exists(ref_fn):
        ref = {float(k): {float(a): v for a, v in d.items()} for k, d in json.load(open(ref_fn)).items()}
    else:
        ref = {}
        for k in KS:
            ref[k] = Mode(c, CDM(), k, **KW).run()
        if wi == 0:
            json.dump({str(k): {str(a): v for a, v in d.items()} for k, d in ref.items()}, open(ref_fn, 'w'), indent=1)
    todo = [sp for sp in specs() if sp['fam'] in fams]
    todo = [sp for i, sp in enumerate(todo) if i % nw == wi]
    for sp in todo:
        fn = os.path.join(outdir, 'model_%s.json' % sp['name'])
        if os.path.exists(fn):
            continue
        if time.time() - t_start > budget - 100:
            print('budget reached', flush=True); break
        t0 = time.time()
        try:
            v, extra = build(c, sp)
        except Exception as e:
            json.dump(dict(spec=sp, error='background: %r' % e), open(fn, 'w'), indent=1); print(sp['name'], 'ERR', e); continue
        cm = extra.pop('cos', c)
        summ = ({k: (float(x) if isinstance(x, (float, np.floating, int)) and not isinstance(x, bool) else x) for k, x in v.summary().items()}
                if sp['fam'] != 'W' else dict(sigma_p=v.sigma_p, f_v_mass=v.f_v_mass))
        rows = []
        for k in KS:
            try:
                m = Mode(cm, v, k, **KW)
                r = m.run()
                rows.append(dict(k=k, P={str(a): (r[a] / ref[k][a])**2 for a in A_OUT}, nsteps=sum(l['nsteps'] for l in m.log)))
            except Exception as e:
                rows.append(dict(k=k, error=repr(e)))
        out = dict(spec=sp, inversion=extra, summary=summ, rows=rows, time_s=time.time() - t0,
                   T_f_kin=float(np.interp(math.log(v.a_f), v.k_lna, v.k_T)) if getattr(v, 'a_f', None) else None)
        json.dump(out, open(fn, 'w'), indent=1, default=float)
        print('%s: fv=%.4f sigp=%.3e af=%.3e  P(z=0):' % (sp['name'], v.f_v_mass, v.sigma_p, getattr(v, 'a_f', 0.0)),
              ' '.join('%g:%.3f' % (r['k'], r['P']['1.0']) if 'P' in r else '%g:ERR' % r['k'] for r in rows),
              '(%.0fs)' % (time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
