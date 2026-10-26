#!/usr/bin/env python3
"""kin_rates finalisation pass (own code): independent spot recomputations against results/rates_interface.json and
consistency summaries of results/mc_checks.json and results/balance.json in the layout of NOTES.md.
Writes results/finalise_checks.json.  usage: python3 -B finalise_checks.py
Recomputations (all closed-form or direct quadrature, no reuse of the track's modules):
  (a) tree masses, trilinear and vacuum width;  (b) MB decay-rate density m_s^2 T Gamma K1(m_s/T)/(2 pi^2) and <m/E>;
  (c) Bose/MB ratio of the s -> X Xb rate density at mu = 0 by 2-D quadrature (E_s and decay angle);
  (d) <sigma v>(T) by the Gondolo-Gelmini quadrature of the tabulated sigma(s) for X X -> X X and s s -> X Xb;
  (e) contact limit |M_c|^2/(16 pi s S_out) at the tabulated point nearest sqrt s = 1e5 m_X for the six 2->2 channels."""
import json
import numpy as np
from scipy.special import kv, kve
from scipy.integrate import quad

RES = '../results'
I = json.load(open(f'{RES}/rates_interface.json'))
MC = json.load(open(f'{RES}/mc_checks.json'))
B = json.load(open(f'{RES}/balance.json'))
P = I['meta']['parameters']
out = {}

# (a) tree level
g, lX, lP, v = P['gX'], P['lambdaX'], P['lambdaPhi'], P['v_GeV']
mX = g * v; ms = np.sqrt(2 * lP) * v; lam3 = 2 * g * g * v; beta = np.sqrt(1 - 4 * mX**2 / ms**2)
Gam = lam3**2 * beta / (16 * np.pi * ms)
out['a_tree'] = dict(mX=mX, ms=ms, lam3=lam3, Gamma_vac=Gam, ms_over_mX=ms / mX,
                     rel_dev=dict(mX=mX / P['mX_GeV'] - 1, ms=ms / P['ms_GeV'] - 1, lam3=lam3 / P['lam3_sXXb_GeV'] - 1, Gamma=Gam / P['Gamma_s_vac_GeV'] - 1))
print('(a) m_X %.7f  m_s %.7f  lam3 %.7f  Gamma_vac %.5e  (rel dev vs interface: %s)' % (mX, ms, lam3, Gam, {k: '%.1e' % x for k, x in out['a_tree']['rel_dev'].items()}))

# (b) MB decay rate density
dec = I['decay_s_XXb']['MB']
TT = np.array(dec['T_over_mX']); gam_tab = np.array(dec['gamma_decay_MB']); eff_tab = np.array(dec['Gamma_eff_over_Gamma_vac'])
rows = []
for t in (0.3, 1.0, 3.0, 7.344, 12.0):
    T = t * mX
    gam = ms**2 * T * Gam * kv(1, ms / T) / (2 * np.pi**2)
    mE = kv(1, ms / T) / kv(2, ms / T)                       # <m/E> for MB: n Gamma <m/E> / (n Gamma)
    i = int(np.argmin(abs(TT - t)))
    rows.append(dict(T_over_mX=t, gamma_decay_MB=gam, table=gam_tab[i], rel=gam / gam_tab[i] - 1, m_over_E=mE, table_eff=eff_tab[i]))
    print('(b) T=%6.3f gamma_dec %.4e table %.4e (rel %.1e)  <m/E> %.3f table %.3f' % (t, gam, gam_tab[i], gam / gam_tab[i] - 1, mE, eff_tab[i]))
out['b_decay_MB'] = rows


# (c) Bose/MB ratio of the decay rate density, mu = 0
def fB(E, T):
    return 1.0 / np.expm1(E / T)


def bose_ratio(T):
    def ang(E, p):
        f = lambda c: fB(E / 2 + p * beta * c / 2, T)
        return 1 + quad(f, -1, 1, epsabs=0, epsrel=1e-9)[0]      # 1 + <f1 + f2> = 1 + 2 <f1>, <f1> = int dc/2 f
    def integrand(E, bose):
        p = np.sqrt(E * E - ms * ms)
        occ = fB(E, T) if bose else np.exp(-E / T)
        fac = ang(E, p) if bose else 1.0
        return p * E * occ * (ms / E) * fac                        # d^3p/(2pi)^3 -> p E dE/(2 pi^2)
    be = quad(integrand, ms, ms + 60 * T, args=(True,), epsabs=0, epsrel=1e-8, limit=400)[0]
    mb = quad(integrand, ms, ms + 60 * T, args=(False,), epsabs=0, epsrel=1e-8, limit=400)[0]
    return be / mb, mb * Gam / (2 * np.pi**2)


rows = []
for t in (0.3, 1.0, 3.0, 7.344, 12.0):
    r, gmb = bose_ratio(t * mX)
    tab = I['decay_s_XXb']['bose_over_MB_mu0'][str(t)]
    rows.append(dict(T_over_mX=t, bose_over_MB_quadrature=r, mc_decay=tab['decay_BE_over_MB'], inverse_quadrature=tab['inverse_BE_quadrature_over_MB'],
                     rel_mc=tab['decay_BE_over_MB'] / r - 1, rel_inv=tab['inverse_BE_quadrature_over_MB'] / r - 1, gamma_MB_check=gmb))
    print('(c) T=%6.3f Bose/MB decay: own quadrature %.4f  track MC %.4f (%.2e)  track inverse-decay quadrature %.4f (%.2e)' % (t, r, tab['decay_BE_over_MB'], rows[-1]['rel_mc'], tab['inverse_BE_quadrature_over_MB'], rows[-1]['rel_inv']))
out['c_decay_bose'] = rows


# (d) Gondolo-Gelmini thermal averages from the tabulated sigma(s)
def n_eq(m, T):
    return m * m * T * kve(2, m / T) * np.exp(-m / T) / (2 * np.pi**2)


def sigmav_GG(ch, col, T, m1, m2):
    c = I['channels_22'][ch]
    sq = np.array(c['sigma_table']['sqrt_s_over_mX']) * mX; sig = np.array(c['sigma_table'][col])
    ok = sig > 0
    lsq, lsig = np.log(sq[ok]), np.log(sig[ok])
    sth = (m1 + m2)**2
    def sigma_of(s):
        x = np.log(np.sqrt(s))
        return np.exp(np.interp(x, lsq, lsig, left=lsig[0], right=lsig[-1] - 2 * (x - lsq[-1])))
    def integrand(u):   # u = sqrt s
        s = u * u
        lam = (s - (m1 + m2)**2) * (s - (m1 - m2)**2)
        return 2 * u * sigma_of(s) * lam * kve(1, u / T) * np.exp(-u / T) / u
    top = np.sqrt(sth) + 80 * T
    val = quad(integrand, np.sqrt(sth), top, epsabs=0, epsrel=1e-7, limit=800, points=[np.sqrt(sth) * (1 + 1e-3), np.sqrt(sth) + T, np.sqrt(sth) + 10 * T])[0]
    return T / (32 * np.pi**4) * val / (n_eq(m1, T) * n_eq(m2, T))


rows = []
for ch, col, t, m1, m2 in (('XX_XX', 'sigma_sub', 1.0, mX, mX), ('ss_XXb', 'sigma_sub', 3.0, ms, ms), ('XX_XX', 'sigma_sub', 7.344, mX, mX), ('Xs_Xs', 'sigma_nores', 1.0, mX, ms)):
    sv = sigmav_GG(ch, col, t * mX, m1, m2)
    th = I['channels_22'][ch]['thermal_MB']
    i = int(np.argmin(abs(np.array(th['T_over_mX']) - t)))
    key = 'sigmav_sub' if col == 'sigma_sub' else 'sigmav_nores'
    rows.append(dict(channel=ch, column=col, T_over_mX=t, sigmav_GG=sv, table=th[key][i], table_T=th['T_over_mX'][i], rel=sv / th[key][i] - 1))
    print('(d) %s %s T=%.3f: GG quadrature of the table %.4e  track %.4e  (rel %.1e)' % (ch, col, t, sv, th[key][i], sv / th[key][i] - 1))
out['d_sigmav_GG'] = rows

# (e) contact limits at the tabulated point nearest 1e5 m_X
rows = []
for ch, c in I['channels_22'].items():
    sq = np.array(c['sigma_table']['sqrt_s_over_mX']); col = 'sigma_nores' if 'sigma_nores' in c['sigma_table'] and ch in ('XXb_XXb', 'Xs_Xs') else 'sigma_sub'
    sig = np.array(c['sigma_table'][col]); i = int(np.argmin(abs(np.log(sq) - np.log(1e5))))
    s = (sq[i] * mX)**2; Mc = c['high_energy']['M_contact']; contact = Mc**2 / (16 * np.pi * s * c['Sout'])
    rows.append(dict(channel=ch, column=col, sqrt_s_over_mX=sq[i], sigma_table=sig[i], contact=contact, rel=sig[i] / contact - 1))
    print('(e) %-8s %-12s sqrt s = %.4e m_X: table %.4e  contact %.4e  (rel %.1e)' % (ch, col, sq[i], sig[i], contact, sig[i] / contact - 1))
out['e_contact'] = rows

# high-T estimate of T_c (closed form) and the m_s/2m_X ratio at 12 m_X from the table
Tc_hiT = v * np.sqrt(lP / (lP / 4 + g * g / 6))
ol = I['masses_T']['oneloop']; i12 = int(np.argmin(abs(np.array(ol['T_over_mX']) - 12)))
out['f_thermal'] = dict(Tc_hiT_over_mX=Tc_hiT / mX, Tc_table_oneloop=ol['T_c_over_mX'], ms_over_2mX_at_12=ol['ms_over_2mX'][i12], T_close=ol['T_close_s_decay_over_mX'])
print('(f) high-T T_c = %.3f m_X (table one-loop %.3f); m_s/2m_X at T=%.2f: %.4f; T_close %.3f' % (Tc_hiT / mX, ol['T_c_over_mX'], ol['T_over_mX'][i12], ol['ms_over_2mX'][i12], ol['T_close_s_decay_over_mX']))

# ---- consistency summaries from mc_checks / balance in the NOTES layout ----
print('\n(g) Bose/MB factors (mc_checks): T | XXb nores | XX sub | Xs nores | ss sub | ss->XXb sub | XXb->ss sub')
bose = {}
for row in MC['thermal_rates_and_balance']:
    t = row['T_over_mX']
    vals = [row['XXb_XXb']['bose_over_MB'], row['XX_XX']['bose_over_MB'], row['Xs_Xs_nores']['bose_over_MB'], row['ss_ss']['bose_over_MB'], row['ss_XXb']['bose_over_MB'], row['XXb_ss']['bose_over_MB']]
    bose[str(t)] = vals
    print('   %6.3f | ' % t + ' | '.join('%.2f' % x for x in vals), '  mc/quad MB:', ', '.join('%.4f' % row[k]['mc_over_quad'] for k in ('XXb_XXb', 'XX_XX', 'Xs_Xs_nores', 'ss_ss', 'ss_XXb', 'XXb_ss')))
out['g_bose_factors'] = bose
print('\n(h) Bose detailed balance, integrated (mc_checks): s<->XXb ratio decay/inverse and sigma; ss<->XXb forward/reverse ratio and sigma')
db = {}
for row in MC['thermal_rates_and_balance']:
    t = row['T_over_mX']; d = {}
    for k, vv in row.items():
        if k.startswith('db_'):
            d[k] = (vv['ratio'], vv['sigma_dev'])
    db[str(t)] = d
    print('   T=%6.3f ' % t + '  '.join('%s: %.4f (%.2f sig)' % (k[3:].replace('_BE_mup', ' mu+=').replace('_mum', ' mu-='), r, s) for k, (r, s) in d.items()))
out['h_bose_db'] = db
rs = [abs(x[0] - 1) for d in db.values() for k, x in d.items() if k.startswith('db_s_')]; ss = [abs(x[1]) for d in db.values() for k, x in d.items() if k.startswith('db_s_')]
rs2 = [x[0] for d in db.values() for k, x in d.items() if k.startswith('db_ss_')]; ss2 = [abs(x[1]) for d in db.values() for k, x in d.items() if k.startswith('db_ss_')]
print('   s<->XXb: max |ratio-1| %.4f, max |sigma| %.2f over %d cases; ss<->XXb: ratio range %.3f-%.3f, max |sigma| %.2f over %d cases' % (max(rs), max(ss), len(rs), min(rs2), max(rs2), max(ss2), len(rs2)))
out['h_summary'] = dict(s_XXb_max_abs_ratio_dev=max(rs), s_XXb_max_sigma=max(ss), ss_XXb_ratio_min=min(rs2), ss_XXb_ratio_max=max(rs2), ss_XXb_max_sigma=max(ss2))
pk = MC['point_kinematics']; em = MC['elastic_energy_moment']
print('(i) point kinematics: %d cases, max residual %.2e, energy %.2e; elastic energy moment max |(Q1+Q2)/(|Q1|+|Q2|)| = %.2e over %d cases' % (pk['total_cases'], pk['max_residual'], pk['max_energy_residual'], max(abs(r['sum_over_abs']) for r in em), len(em)))
out['i_point_kin'] = dict(total_cases=pk['total_cases'], max_residual=pk['max_residual'], max_energy_residual=pk['max_energy_residual'], elastic_moment_max=max(abs(r['sum_over_abs']) for r in em), n_elastic=len(em))
print('(j) hard particle on bath, rows with E/m = 7.72e4:')
for r in MC['hard_particle']:
    if r['E_over_m'] > 1e4:
        print('   T=%.3f E/m=%.3g: ' % (r['T_over_m'], r['E_over_m']) + '; '.join('%s rate/n %.2e dEdt/n %.2e' % (k, c['rate_per_density_GeV_minus2'], c['dEdt_per_density_GeV_minus1']) for k, c in r['channels'].items()))
print('(k) old s route: T, 1/ratio (old/new):', ', '.join('%.3g: %.3g' % (r['T_over_mX'], 1 / r['ratio_sub_to_old']) for r in B['old_s_route'] if abs(r['T_over_mX'] - round(r['T_over_mX'], 3)) < 1e-9 and r['T_over_mX'] in (0.03, 0.3, 1.0, 3.0, 7.344, 12.0)))
print('(l) nr_limits:')
for ch, d in B['nr_limits'].items():
    print('   %s A_old %.5g: ' % (ch, d['A_old']) + ', '.join('T=%.3g sub %.3g full %.3g' % (r['T_over_mX'], r['ratio_sub_to_old'], r['ratio_fullreg_to_old']) for r in d['rows']))
# table-level detailed balance summary
dbt = B['detailed_balance']
polefree = ('F_sss_ss', 'N_ssss_ss', 'L_sss_XXb', 'E2_XXs_XX')
pf = [r for r in dbt if r['channel'] in polefree]; res = [r for r in dbt if r['channel'] not in polefree]
print('(m) table-level MB detailed balance: %d rows; pole-free channels max |ratio_full-1| %.3f; resonant channels ratio_full range %.2f-%.2f (T>=1: %.2f-%.2f)' % (
    len(dbt), max(abs(r['ratio_full'] - 1) for r in pf), min(r['ratio_full'] for r in res), max(r['ratio_full'] for r in res),
    min(r['ratio_full'] for r in res if r['T_over_mX'] >= 1), max(r['ratio_full'] for r in res if r['T_over_mX'] >= 1)))
cons = [r for r in res if r['T_over_mX'] >= 1]
ok = sum(1 for r in cons if abs(r['gamma_n2_sub'] - r['gamma_2n_sub']) <= 2 * np.hypot(r['gamma_n2_sub_err'], r['gamma_2n_sub_err']) + 0.1 * (abs(r['gamma_n2_sub']) + abs(r['gamma_2n_sub'])))
print('    subtracted remainders at T >= 1 m_X consistent within 2 sigma (+10%% interpolation): %d of %d' % (ok, len(cons)))
out['m_table_db'] = dict(rows=len(dbt), polefree_max_dev=max(abs(r['ratio_full'] - 1) for r in pf), resonant_full_range=[min(r['ratio_full'] for r in res), max(r['ratio_full'] for r in res)], sub_consistent=[ok, len(cons)])
# grid sizes
out['n_grid'] = dict(T_22=len(I['channels_22']['XX_XX']['thermal_MB']['T_over_mX']), T_n2=len(I['channels_n2']['A_XXXbXb_XXb']['thermal_MB']['T_over_mX']),
                     sigma_22=len(I['channels_22']['XX_XX']['sigma_table']['sqrt_s_over_mX']), Gt_A=len(I['channels_n2']['A_XXXbXb_XXb']['Gt_table_reg0_02']['sqrt_s_over_mX']),
                     n2_channels=len(I['channels_n2']), twon_channels=len(I['channels_2n']))
print('(n) grids:', out['n_grid'])
json.dump(out, open(f'{RES}/finalise_checks.json', 'w'), indent=1, default=float)
