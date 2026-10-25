#!/usr/bin/env python3
"""First structure screen (retrospective, approximate) -- Claude review, track `conversion`.

Linear CLASS (Lesgourgues et al.; classy 3.3.2.0) power spectra for CDM plus
a non-cold fraction with the computed momentum distribution of the converted X
(from my Monte Carlo escape times, q/m = u_esc a_esc), optionally plus the
original positive daughters (q/m = a_f sqrt(0.5)).  Compared with LCDM at equal
total omega_dm.  Small-scale screen: the Murgia et al. (2017) 1D area criterion,
dA = 1 - A/A_LCDM, A = int_{0.5}^{20 h/Mpc} P1D/P1D_LCDM dk, P1D(k) =
(1/2pi) int_k^inf k' P(k') dk', linear z = 0, with dA_ref = 0.38 (3.5 keV
thermal, conservative) and 0.21 (5.3 keV).  My pipeline is calibrated by
recomputing dA for pure thermal WDM at 3.5 and 5.3 keV, and by computing dA for
1 keV thermal mixed at f = 0.16 (the 2-sigma bound of Garcia-Gallego et al. 2025).
This is a screen, not a likelihood.  Run with the scratch venv python.
"""
import json
import sys
from pathlib import Path

import numpy as np
from classy import Class

SCR = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/conversion')
RES = Path(__file__).resolve().parents[1] / 'results'
H = 0.674
OMB = 0.0224
OMDM = 0.315 * H * H - OMB            # 0.12070
TCMB = 2.7255
TNC = 0.71611                        # T_ncdm / T_cmb
KT = TNC * TCMB * 8.617333262e-5     # eV, today
UESC = 2.532075973720193
X0 = 3.0                             # psd peak position in units of T_ncdm
EV3_TO_CM3 = 1 / (1.973269804e-5) ** 3
RHOC_H2 = 10537.4986                 # eV cm^-3 (rho_crit / h^2)

BASE = {'h': H, 'omega_b': OMB, 'N_ur': 3.046, 'A_s': 2.1e-9, 'n_s': 0.965, 'tau_reio': 0.054,
        'output': 'mPk', 'z_pk': '0,3', 'YHe': 0.245, 'l_max_ncdm': 12}
CLS = {'output': 'tCl,pCl,lCl', 'lensing': 'yes', 'l_max_scalars': 2500}


NBINS = 60
PKMAX = 120.
DO_CLS = True
CLS_TAGS = {'T1keV_codex_1', 'epoch_codex_z1e06', 'epoch_codex_z3e05'}


def psd_file(u, tag, min_rel=0.10, width_scale=1.0):
    """psd in x = q/T_ncdm with the sample mean at X0.  Narrow samples (relative
    spread < min_rel): Gaussian of relative width min_rel (P(k) was checked to be
    insensitive to this width: 3.6% vs 10.7% changes P(k) by < 3e-4).  Broad
    samples: Gaussian-kernel density of the sample (bandwidth 5% of X0).  Always
    strictly positive so CLASS's dln f/dln q stays finite."""
    ubar = float(np.mean(u))
    x = X0 * u / ubar
    rel = float(np.std(u) / ubar)
    if rel < min_rel:
        sig = min_rel * X0 * width_scale
        xs = np.linspace(X0 - 7 * sig, X0 + 7 * sig, 500)
        f = np.exp(-0.5 * ((xs - X0) / sig) ** 2) / xs ** 2
    else:
        bw = 0.05 * X0 * width_scale
        xs = np.linspace(max(0.02, x.min() - 6 * bw), x.max() + 6 * bw, 800)
        sub = x if x.size <= 20000 else np.random.default_rng(0).choice(x, 20000, replace=False)
        dens = np.zeros_like(xs)
        for chunk in np.array_split(sub, max(1, sub.size // 4000)):
            dens += np.exp(-0.5 * ((xs[:, None] - chunk[None, :]) / bw) ** 2).sum(axis=1)
        dens = np.maximum(dens / dens.max(), 1e-200)
        f = dens / xs ** 2
    path = SCR / 'classrun' / f'psd_{tag}.dat'
    np.savetxt(path, np.c_[xs, f])
    norm = np.trapezoid(f * xs ** 2, xs)
    return path, ubar, norm


def ncdm_entry(u, frac, tag, width_scale=1.0):
    path, ubar, norm = psd_file(u, tag, width_scale=width_scale)
    m = X0 * KT / ubar                                   # eV
    n_cm3 = (KT ** 3) / (2 * np.pi ** 2) * norm * EV3_TO_CM3   # per dof, factor 1
    om_unit = m * n_cm3 / RHOC_H2                        # omega_ncdm for psd factor 1
    factor = frac * OMDM / om_unit
    xmax = float(np.loadtxt(path)[-1, 0])
    return dict(file=str(path), m=m, factor=factor, ubar=ubar, omega=frac * OMDM, qmax=xmax)


def run(ncdm_list, label, cls=None):
    cls = DO_CLS if cls is None else cls
    pars = dict(BASE)
    pars['P_k_max_h/Mpc'] = PKMAX
    om_nc = sum(e['omega'] for e in ncdm_list)
    pars['omega_cdm'] = OMDM - om_nc
    if ncdm_list:
        pars.update({'N_ncdm': len(ncdm_list),
                     'use_ncdm_psd_files': ','.join(['1' if 'file' in e else '0' for e in ncdm_list]),
                     'm_ncdm': ','.join(f"{e['m']:.10e}" for e in ncdm_list),
                     'T_ncdm': ','.join(f"{e.get('T', TNC):.8f}" for e in ncdm_list),
                     'omega_ncdm': ','.join(f"{e['omega']:.10e}" for e in ncdm_list),
                     'deg_ncdm': ','.join('1' for _ in ncdm_list)})
        files = [e['file'] for e in ncdm_list if 'file' in e]
        if files:
            pars['ncdm_psd_filenames'] = ','.join(files)
            # fixed trapezoidal momentum grid for tabulated (narrow) distributions
            pars['ncdm_quadrature_strategy'] = ','.join('2' if 'file' in e else '0' for e in ncdm_list)
            pars['ncdm_maximum_q'] = ','.join(f"{e.get('qmax', 15.0):.4f}" for e in ncdm_list)
            pars['ncdm_N_momentum_bins'] = ','.join(str(NBINS) if 'file' in e else '150' for e in ncdm_list)
    c = Class(); c.set(pars); c.compute()
    k = np.geomspace(1e-4, PKMAX * 0.995, 1500)          # h/Mpc
    pk0 = np.array([c.pk_lin(kk * H, 0.0) for kk in k]) * H ** 3
    pk3 = np.array([c.pk_lin(kk * H, 3.0) for kk in k]) * H ** 3
    s8 = c.sigma(8 / H, 0.0)
    om_ncdm = c.Omega_nu * H * H if ncdm_list else 0.0
    c.struct_cleanup(); c.empty()
    # separate run for CMB spectra (keeps the Cl k-sampling away from large k)
    clpp = np.ones(4); cltt = np.ones(3)
    if cls:
        pc = {kk: vv for kk, vv in pars.items() if kk not in ('output', 'P_k_max_h/Mpc', 'z_pk')}
        pc.update(CLS)
        c = Class(); c.set(pc); c.compute()
        cl = c.lensed_cl(2500)
        clpp = cl['pp'][[100, 400, 1000, 2000]]
        cltt = cl['tt'][[200, 1000, 2000]]
        c.struct_cleanup(); c.empty()
    return dict(label=label, k=k, pk0=pk0, pk3=pk3, sigma8=s8, omega_ncdm_class=om_ncdm, clpp=clpp, cltt=cltt)


def p1d(k, pk):
    # (1/2pi) int_k^inf k' P dk' with a power-law tail beyond the last point
    integrand = k * pk
    cum = np.zeros_like(k)
    seg = 0.5 * (integrand[1:] + integrand[:-1]) * np.diff(k)
    cum[:-1] = np.cumsum(seg[::-1])[::-1]
    slope = np.log(integrand[-1] / integrand[-20]) / np.log(k[-1] / k[-20])
    tail = -integrand[-1] * k[-1] / (slope + 1) if slope < -1 else 0.0
    return (cum + tail) / (2 * np.pi)


def delta_area(ref, mod, kmin=0.5, kmax=20.):
    r = p1d(mod['k'], mod['pk0']) / p1d(ref['k'], ref['pk0'])
    sel = (ref['k'] >= kmin) & (ref['k'] <= kmax)
    kk = ref['k'][sel]
    A = np.trapezoid(r[sel], kk)
    return float(1 - A / (kk[-1] - kk[0]))


def thermal(m_kev, frac):
    """thermal relic, g = 1.5 (deg 1 FD), temperature fixed by its own abundance."""
    om = frac * OMDM
    T_over_Tnu = (om * 93.14 / (m_kev * 1000)) ** (1 / 3)
    return dict(m=m_kev * 1000, T=TNC * T_over_Tnu, omega=om)


def thermal_full_convention(m_kev, frac):
    """mixed thermal relic whose temperature is that of an all-DM relic of the same mass."""
    T_over_Tnu = (OMDM * 93.14 / (m_kev * 1000)) ** (1 / 3)
    return dict(m=m_kev * 1000, T=TNC * T_over_Tnu, omega=frac * OMDM, factor=frac)


KSTAR_SKM = 0.009


def kstar_hmpc(z=3.0):
    E = np.sqrt(0.315 * (1 + z) ** 3 + 0.685)
    return KSTAR_SKM * 100 * E / (1 + z)        # h/Mpc (H/(1+z) in h km/s/Mpc)


def compressed(mod, z=3.0):
    ks = kstar_hmpc(z)
    k = mod['k']; pk = mod['pk3']
    lk, lp = np.log(k), np.log(pk)
    i = np.searchsorted(k, ks)
    sl = slice(max(i - 6, 0), i + 6)
    c = np.polyfit(lk[sl] - np.log(ks), lp[sl], 2)
    P = np.exp(c[2]); n = c[1]
    return ks ** 3 * P / (2 * np.pi ** 2), n


# Published compressed Lya constraints at k* = 0.009 s/km, z = 3 [Q]:
DESI_DR1 = dict(D2=0.379, sD2=0.032, n=-2.309, sn=0.019, rho=0.0)    # arXiv:2601.21432 (correlation not used)
EBOSS = dict(D2=0.310, sD2=0.020, n=-2.340, sn=0.006, rho=0.512)     # Rogers & Poulin, arXiv:2311.16377 (eBOSS DR14)


def chi2(D2, n, d):
    x = (D2 - d['D2']) / d['sD2']; y = (n - d['n']) / d['sn']; r = d['rho']
    return float((x * x + y * y - 2 * r * x * y) / (1 - r * r))


def summarize(ref, mod):
    k = ref['k']
    out = dict(label=mod['label'], sigma8=mod['sigma8'], sigma8_ratio=mod['sigma8'] / ref['sigma8'],
               omega_ncdm_class=mod['omega_ncdm_class'], dA=delta_area(ref, mod),
               Clpp_ratio_l100_400_1000_2000=(mod['clpp'] / ref['clpp']).tolist(),
               ClTT_lensed_ratio_l200_1000_2000=(mod['cltt'] / ref['cltt']).tolist())
    D2r, nr = compressed(ref); D2m, nm = compressed(mod)
    # apply the model/LCDM shift to my LCDM baseline
    out.update(D2star_LCDM=D2r, nstar_LCDM=nr, D2star=D2m, nstar=nm,
               D2star_ratio=D2m / D2r, dnstar=nm - nr,
               chi2_DESI_DR1=chi2(D2m, nm, DESI_DR1), chi2_DESI_DR1_LCDM=chi2(D2r, nr, DESI_DR1),
               chi2_eBOSS=chi2(D2m, nm, EBOSS), chi2_eBOSS_LCDM=chi2(D2r, nr, EBOSS))
    for kk in [0.1, 0.3, 1., 3., 5., 10., 20., 50.]:
        out[f'P_ratio_z0_k{kk}'] = float(np.interp(kk, k, mod['pk0'] / ref['pk0']))
        out[f'P_ratio_z3_k{kk}'] = float(np.interp(kk, k, mod['pk3'] / ref['pk3']))
    return out


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    out = {}
    ref = run([], 'LCDM')
    out['LCDM_sigma8'] = ref['sigma8']
    cal = []
    if which == 'elastic':
        global DO_CLS
        DO_CLS = False
    if which in ('all', 'cal'):
        for m in [3.5, 5.3, 5.7]:
            cal.append(summarize(ref, run([dict(thermal(m, 1.0))], f'thermal_{m}keV_full')))
            print(cal[-1], flush=True)
        for m, f in [(1.0, 0.16), (2.0, 0.35)]:
            cal.append(summarize(ref, run([dict(thermal_full_convention(m, f))], f'mixed_{m}keV_f{f}_Tfull')))
            print(cal[-1], flush=True)
            cal.append(summarize(ref, run([dict(thermal(m, f))], f'mixed_{m}keV_f{f}_Town')))
            print(cal[-1], flush=True)
        out['calibration'] = cal
    rows = []
    if which in ('all', 'model'):
        mcfiles = [x for x in [SCR / 'mc_aesc_T1keV_codex_1.npy'] if x.exists()] + sorted(SCR.glob('mc_aesc_epoch_codex_*.npy'))
        for f in mcfiles:
            tag = f.stem.replace('mc_aesc_', '')
            ae = np.load(f)
            u = UESC * ae
            meta = json.loads((f.with_suffix('.json')).read_text()) if f.with_suffix('.json').exists() else {}
            af = meta.get('a_form')
            if af is None:
                af = {'z1e6_codex_1': 1 / 1000001., 'T1keV_codex_1': 2.34865418e-07}.get(tag)
            conv = ncdm_entry(u, 0.033588, tag + '_conv')
            r1 = summarize(ref, run([conv], tag + '_converted_f0.0336', cls=tag in CLS_TAGS))
            r1.update(ubar=conv['ubar'], m_equiv_eV_at_Tnu=conv['m'], a_form=af)
            rows.append(r1); print(r1, flush=True)
            if af is not None:
                ud = np.full(20000, af * np.sqrt(0.5)) * (1 + 0.002 * np.random.default_rng(0).standard_normal(20000))
                dau = ncdm_entry(ud, 0.033588, tag + '_dau')
                r2 = summarize(ref, run([conv, dau], tag + '_allfreeX_f0.0672', cls=False))
                r2.update(ubar_daughters=dau['ubar'], a_form=af)
                rows.append(r2); print(r2, flush=True)
        out['models'] = rows
    if which == 'elastic':
        el = json.loads((RES / 'elastic_epochs.json').read_text())
        for e in el:
            qe = e['q_eff_converted']; qd = e['q_eff_daughters']
            conv = ncdm_entry(np.full(20000, qe), 0.033588, e['tag'] + '_conv_el')
            r1 = summarize(ref, run([conv], e['tag'] + '_converted_elastic_qeff'))
            r1.update(a_form=e['a_form'], q_eff=qe, q_collisionless=e['q_over_m'])
            rows.append(r1); print(r1, flush=True)
            dau = ncdm_entry(np.full(20000, qd), 0.033588, e['tag'] + '_dau_el')
            r2 = summarize(ref, run([conv, dau], e['tag'] + '_allfreeX_elastic_qeff'))
            r2.update(a_form=e['a_form'], q_eff=qe, q_eff_daughters=qd)
            rows.append(r2); print(r2, flush=True)
        out['models'] = rows
    (RES / f'class_screen_{which}.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
