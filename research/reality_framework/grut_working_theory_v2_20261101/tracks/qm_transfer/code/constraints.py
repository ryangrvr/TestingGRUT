"""qm_transfer: from P_mix/P_LCDM(k) to small-scale-structure constraints (own code, Claude, 2 Nov 2026).  [C]/[E]

Criteria (all applied identically to the model and to thermal-WDM reference curves at the published bound masses):
 (1) Lyman-alpha, area criterion (Murgia, Merle, Viel, Totzauer & Schneider 2017, arXiv:1704.07838) [Q, read]:
     P_1D(k) = (1/2pi) int_k^inf k' P(k') dk',  r = P_1D/P_1D^LCDM,  A = int_{0.5 h}^{20 h /Mpc} r dk,
     dA = (A_LCDM - A)/A_LCDM;  excluded if dA > dA(thermal WDM at the bound).  Bounds: 5.3 keV (Irsic et al. 2017,
     arXiv:1702.01764), 5.7 keV (Irsic et al. 2024, arXiv:2309.04533), 3.5 keV conservative (1702.01764) [R].
 (2) MW satellites proxy [E]: cumulative halo abundance N(>M) at z = 0, sharp-k filter (c = 2.5) + Sheth-Tormen
     (Schneider, Smith & Reed 2013; Benson et al. 2013 [R]), at M = 1e8 M_sun; bound: thermal WDM 6.5 keV
     (Nadler et al. 2021, arXiv:2008.00022) [R].
 (3) strong-lensing proxy [E]: differential abundance dn/dlnM at M = 10^7.5 M_sun; bounds 5.2 keV (Gilman et al. 2020,
     arXiv:1908.06983: m_hm < 10^7.8) and 6.1 keV (Keeley et al. 2024, arXiv:2405.01620) [R].
 (4) future reach [R forecasts, used as thresholds only]: the same proxies at 10 and 20 keV (Rubin-era satellites
     10-20 keV, Nadler et al. 2024 forecast; lensing with larger samples ~10 keV), and dA at 10 keV for a
     next-generation Lyman-alpha (DESI-like high-resolution follow-up) [E].
Thermal-WDM reference curves: (a) Viel et al. 2005 fit (astro-ph/0501562), the form behind the published bounds;
(b) the same relic through our own Boltzmann code (scan family W), which is colder than the Viel fit, as is CLASS
(Vogel & Abazajian 2022, arXiv:2210.10753).  Both are reported.  P_LCDM shape: Eisenstein & Hu 1998 no-wiggle
(astro-ph/9709112) with n_s = 0.9649, sigma_8 = 0.811 (only ratios and the Sheth-Tormen nu depend on it).
"""
import json, math, os, sys, glob
import numpy as np
from scipy.integrate import quad, cumulative_trapezoid

H = 0.6736
OM = 0.3153
NS = 0.9649
SIG8 = 0.811
RHO_M = OM * 2.775e11 * H**2                   # M_sun / Mpc^3
KMIN_LY, KMAX_LY = 0.5 * H, 20.0 * H


def eh_nw(k, h=H, om=OM, omb=0.02237, T0=2.7255):
    omh2 = om * h * h
    fb = omb / omh2
    th = T0 / 2.7
    s = 44.5 * np.log(9.83 / omh2) / np.sqrt(1 + 10 * omb**0.75)
    ag = 1 - 0.328 * np.log(431 * omh2) * fb + 0.38 * np.log(22.3 * omh2) * fb**2
    G = om * h * (ag + (1 - ag) / (1 + (0.43 * k * s)**4))
    q = k / h * th**2 / G
    L0 = np.log(2 * np.e + 1.8 * q); C0 = 14.2 + 731 / (1 + 62.5 * q)
    return L0 / (L0 + C0 * q * q)


KG = np.exp(np.linspace(np.log(1e-4), np.log(3e3), 6000))
P_LCDM = KG**NS * eh_nw(KG)**2
# sigma_8 normalisation (top hat R = 8/h)
_R8 = 8.0 / H
_W = lambda x: 3 * (np.sin(x) - x * np.cos(x)) / x**3
_s8 = np.sqrt(np.trapezoid(KG**2 * P_LCDM * _W(KG * _R8)**2, KG) / (2 * np.pi**2))
P_LCDM *= (SIG8 / _s8)**2


def viel_T2(m_keV, k, h=H, Om_dm=0.1200 / H**2):
    al = 0.049 * m_keV**(-1.11) * (Om_dm / 0.25)**0.11 * (h / 0.7)**1.22 / h
    return (1 + (al * k)**2.24)**(-10 / 1.12)


def ratio_on_grid(ks, Ps):
    """interpolate P/P_LCDM in ln k on KG; below the first point a power-law deficit, above the last a constant."""
    ks = np.asarray(ks, float); Ps = np.asarray(Ps, float)
    lk = np.log(KG)
    r = np.interp(lk, np.log(ks), Ps)
    d0, d1 = 1 - Ps[0], 1 - Ps[1]
    if d0 > 1e-5 and d1 > d0:
        n = math.log(d1 / d0) / math.log(ks[1] / ks[0])
        n = min(max(n, 0.5), 4.0)
    else:
        n = 2.0
    lo = KG < ks[0]
    r[lo] = 1 - max(d0, 0.0) * (KG[lo] / ks[0])**n
    return np.clip(r, 0.0, None)


KUP_LY = 75.0 * H     # upper limit of the P_1D integral; 75 h/Mpc reproduces Murgia et al.'s dA_REF = 0.38 (3.5 keV)
                      # and 0.21 (5.3 keV) for the Viel fit (their P(k) output limit) [C, see NOTES]; None = infinity


def dA_of(r, kup=KUP_LY):
    P = P_LCDM * r
    sel0 = KG <= (kup if kup else KG[-1])
    K = KG[sel0]
    # P_1D(k) = (1/2pi) int_k^kup k' P dk' via reverse cumulative integral
    I = cumulative_trapezoid((K * P[sel0])[::-1], -K[::-1], initial=0.0)[::-1] / (2 * np.pi)
    I0 = cumulative_trapezoid((K * P_LCDM[sel0])[::-1], -K[::-1], initial=0.0)[::-1] / (2 * np.pi)
    sel = (K >= KMIN_LY) & (K <= KMAX_LY)
    rr = I[sel] / I0[sel]
    A = np.trapezoid(rr, K[sel])
    return 1 - A / (KMAX_LY - KMIN_LY)


def hmf(r, Ms=(10**6.5, 1e7, 10**7.5, 1e8, 10**8.5, 1e9)):
    """sharp-k + Sheth-Tormen: returns (dn/dlnM, N(>M)) at Ms (z = 0)."""
    P = P_LCDM * r
    c = 2.5
    Mg = np.exp(np.linspace(np.log(1e5), np.log(1e13), 800))
    Rg = (3 * Mg / (4 * np.pi * RHO_M))**(1 / 3) / c
    cum = cumulative_trapezoid(KG**2 * P, KG, initial=0.0) / (2 * np.pi**2)
    s2 = np.interp(1 / Rg, KG, cum)
    dlns2_dlnR = -np.interp(1 / Rg, KG, P) / (2 * np.pi**2 * Rg**3 * s2)
    nu = 1.686**2 / s2
    A, a, p = 0.3222, 0.707, 0.3
    f = A * np.sqrt(2 * a * nu / np.pi) * (1 + (a * nu)**(-p)) * np.exp(-a * nu / 2)
    dndlnM = RHO_M / Mg * f * np.abs(dlns2_dlnR) / 6.0      # dlnR/dlnM = 1/3, d ln nu = -dlns2 -> |.|/2 * 1/3
    Ncum = cumulative_trapezoid(dndlnM[::-1], -np.log(Mg)[::-1], initial=0.0)[::-1]
    return (np.interp(np.log(Ms), np.log(Mg), dndlnM), np.interp(np.log(Ms), np.log(Mg), Ncum))


_R1 = np.ones_like(KG)
_DN0, _NC0 = hmf(_R1)


def metrics(r):
    dn, Nc = hmf(r)
    return dict(dA=float(dA_of(r)), dA_inf=float(dA_of(r, None)), N_gt_1e8=float(Nc[3] / _NC0[3]), N_gt_1e7=float(Nc[1] / _NC0[1]),
                dn_1e7p5=float(dn[2] / _DN0[2]), dn_1e8=float(dn[3] / _DN0[3]), dn_1e7=float(dn[1] / _DN0[1]))


def half_mode(ks, Ps):
    """k where P/P_LCDM = 1/4 (T = 1/2), or None"""
    Ps = np.asarray(Ps); ks = np.asarray(ks)
    j = np.where(Ps < 0.25)[0]
    if len(j) == 0:
        return None
    j = j[0]
    if j == 0:
        return float(ks[0])
    x = np.interp(0.25, [Ps[j], Ps[j - 1]], [np.log(ks[j]), np.log(ks[j - 1])])
    return float(math.exp(x))


def m_equiv_keV(sigma_p):
    """thermal-relic mass with the same 1D comoving velocity dispersion today (FD, g = 2): sigma = 2.08 Tq,
    Tq(1 keV) = 3.762e-8, Tq ~ m^-4/3."""
    return (2.0806 * 3.762e-8 / sigma_p)**0.75


def references(Wdir=None):
    refs = {}
    for m in (3.5, 5.2, 5.3, 5.7, 6.1, 6.5, 10.0, 20.0):
        refs['viel_%g' % m] = metrics(viel_T2(m, KG))
    if Wdir:
        for fn in glob.glob(os.path.join(Wdir, 'model_W_m*.json')):
            d = json.load(open(fn))
            ks = [r['k'] for r in d['rows'] if 'P' in r]; Ps = [r['P']['1.0'] for r in d['rows'] if 'P' in r]
            refs['own_%g' % d['spec']['m_keV']] = metrics(ratio_on_grid(ks, Ps))
            refs['own_%g' % d['spec']['m_keV']]['P'] = dict(zip(map(str, ks), Ps))
    return refs


def judge(mt, refs, which='viel'):
    g = lambda m, key: refs['%s_%g' % (which, m)][key] if ('%s_%g' % (which, m)) in refs else refs['viel_%g' % m][key]
    out = dict(
        lya_5p3=mt['dA'] <= g(5.3, 'dA'), lya_5p7=mt['dA'] <= g(5.7, 'dA'), lya_3p5=mt['dA'] <= g(3.5, 'dA'),
        sat_6p5=mt['N_gt_1e8'] >= g(6.5, 'N_gt_1e8'),
        lens_5p2=mt['dn_1e7p5'] >= refs['viel_5.2']['dn_1e7p5'], lens_6p1=mt['dn_1e7p5'] >= refs['viel_6.1']['dn_1e7p5'],
        fut_sat_10=mt['N_gt_1e8'] < g(10.0, 'N_gt_1e8'), fut_sat_20=mt['N_gt_1e8'] < g(20.0, 'N_gt_1e8'),
        fut_lens_10=mt['dn_1e7p5'] < g(10.0, 'dn_1e7p5'), fut_lens_20=mt['dn_1e7p5'] < g(20.0, 'dn_1e7p5'),
        fut_lya_10=mt['dA'] > g(10.0, 'dA'))
    out['allowed_now'] = bool(out['lya_5p7'] and out['sat_6p5'] and out['lens_6p1'])
    out['allowed_conservative'] = bool(out['lya_3p5'] and out['lens_5p2'])
    out['band_10keV'] = bool(out['allowed_now'] and (out['fut_sat_10'] or out['fut_lens_10'] or out['fut_lya_10']))
    out['band_20keV'] = bool(out['allowed_now'] and (out['fut_sat_20'] or out['fut_lens_20']))
    return out


def mixed_dm_crosscheck(f_v, sigma_p):
    """Baur et al. 2017 (arXiv:1706.03118) as quoted in later work: f_WDM <~ 0.14 (m/keV)^1.37 [R];
    Boyarsky et al. 2009 (arXiv:0812.0010): any mass allowed for f_WDM < 0.35 (95%, Bayesian) [R]."""
    m = m_equiv_keV(sigma_p)
    return dict(m_equiv_keV=m, baur_ok=bool(f_v <= 0.14 * m**1.37 or m >= 4.2), boyarsky_any_mass=bool(f_v < 0.35))


def main():
    scan_dir = sys.argv[1]; out_fn = sys.argv[2]
    refs = references(scan_dir)
    models = []
    for fn in sorted(glob.glob(os.path.join(scan_dir, 'model_*.json'))):
        d = json.load(open(fn))
        if d['spec']['fam'] == 'W' or 'rows' not in d:
            continue
        ks = [r['k'] for r in d['rows'] if 'P' in r]; Ps = [r['P']['1.0'] for r in d['rows'] if 'P' in r]
        Ps3 = [r['P']['0.25'] for r in d['rows'] if 'P' in r]
        if len(ks) < len(d['rows']):
            print('incomplete', fn)
        r = ratio_on_grid(ks, Ps)
        mt = metrics(r)
        S = d['summary']
        row = dict(name=d['spec']['name'], spec=d['spec'], f_v_mass=S['f_v_mass'], f_v_charge=S.get('f_v_charge'),
                   sigma_p=S['sigma_p'], a_f=S.get('a_f_kin') or S.get('a_f_eq'), a_coll_end=S.get('a_coll_end'),
                   T_f_eq=S.get('T_f_eq'), T_f_kin=d.get('T_f_kin'), zeta=d['spec'].get('zeta'),
                   P_z0=dict(zip(map(str, ks), Ps)), P_z3=dict(zip(map(str, ks), Ps3)), metrics=mt,
                   k_half=half_mode(ks, Ps), judge_viel=judge(mt, refs, 'viel'), judge_own=judge(mt, refs, 'own'),
                   mixed=mixed_dm_crosscheck(S['f_v_mass'], S['sigma_p']))
        models.append(row)
    json.dump(dict(scope=__doc__.strip().splitlines()[0], references=refs, models=models), open(out_fn, 'w'), indent=1,
              default=float)
    for m in models:
        j = m['judge_viel']
        print('%-28s fv=%.3f sigp=%.2e dA=%.3f N8=%.3f dn7.5=%.3f  now=%d band10=%d band20=%d | own: now=%d' % (
            m['name'], m['f_v_mass'], m['sigma_p'], m['metrics']['dA'], m['metrics']['N_gt_1e8'], m['metrics']['dn_1e7p5'],
            j['allowed_now'], j['band_10keV'], j['band_20keV'], m['judge_own']['allowed_now']))
    for k in sorted(refs):
        print(k, {x: round(y, 4) for x, y in refs[k].items() if x != 'P'})


if __name__ == '__main__':
    main()
