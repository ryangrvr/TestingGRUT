"""Track pt_m1tilt: aggregate every run in $SP/review/pt_m1tilt/runs into results/m1_transition.json, the frequency
diagnostics into results/m1_frequency_lock.json, and print the tables used in NOTES.md.
[C] homogeneous classical model; not a cosmological result."""
import sys, os, json, glob, re, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m1_common import *

RUNS = os.path.join(SCR, 'runs')
RPHI = np.sqrt(2 * LP) / GG     # m_phi / m_X = 2.449 (M1 and untilted)
SKEPTIC = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'pt_matching', 'skeptic', 'results', 'skC_transition.json')


def st(a):
    a = np.asarray(a, float)
    if len(a) == 0: return None
    return dict(n=len(a), mean=float(a.mean()), std=float(a.std(ddof=1)) if len(a) > 1 else 0.0,
                sem=float(a.std(ddof=1) / np.sqrt(len(a))) if len(a) > 1 else 0.0, min=float(a.min()), max=float(a.max()))


def key(r):
    return (r['orbit'], r['mode'], r['epst'], r['beta'], str(r['phi0']), r['eps'], r.get('rtol', 1e-11), r.get('mend', 40.0),
            r.get('pert', 0.0), r.get('M_start_fac', 0.5))


rows = []
for fn in sorted(glob.glob(os.path.join(RUNS, '*.json'))):
    r = json.load(open(fn)); b = os.path.basename(fn)
    r['file'] = b
    r['rtol'] = float(re.search(r'_rt([0-9]e-[0-9]+)', b).group(1))
    r['mend'] = float(re.search(r'_mend([0-9.]+)', b).group(1)) if '_mend' in b else 40.0
    r['design'] = 'pair' if (abs(r['phase'] * 16 - round(r['phase'] * 16)) < 1e-9 or abs((r['phase'] - 1e-3) * 16 - round((r['phase'] - 1e-3) * 16)) < 1e-9
                             or abs(r['phase'] * 8 - round(r['phase'] * 8)) < 1e-9 or abs((r['phase'] - 1e-3) * 8 - round((r['phase'] - 1e-3) * 8)) < 1e-9) else 'rand'
    rows.append(r)
print(len(rows), 'runs')

orbits = {}
for nm in ['tau20', 'late']:
    y4, info = orbit(nm)
    orbits[nm] = dict(q=info['q'], J0_per_q=info['J_per_q'], T_r=info['T_r'], omega_r=info['omega_r'], rmin=info['rmin'], rmax=info['rmax'],
                      r2mean=info['r2mean'], Mperi=np.sqrt(G2 / LP) * info['rmin'], Mapo=np.sqrt(G2 / LP) * info['rmax'])
J0 = {k: v['J0_per_q'] for k, v in orbits.items()}


def block(rs):
    orb = rs[0]['orbit']
    n = np.array([r['nminus'] for r in rs]); N = np.array([r['Nphi'] for r in rs])
    out = dict(nminus=st(n), Nphi=st(N), sum=st(n + N), nminus_over_J0=float(n.mean() / J0[orb]), Nphi_over_J0=float(N.mean() / J0[orb]),
               one_minus_nminus_over_J0=float(1 - n.mean() / J0[orb]), E_final_over_mXq=st(2 * n + RPHI * N),
               max_qrel=max(r['qrel'] for r in rs), max_work_rel=max(r['work_rel'] for r in rs),
               max_work_rel_dyn=max(r.get('work_rel_dyn', r['work_rel']) for r in rs),
               M_escape_over_Mperi=st([r['M_escape_over_Mperi'] for r in rs if r['M_escape_over_Mperi'] is not None]),
               M_half_over_Mperi=st([r['M_half_over_Mperi'] for r in rs if r.get('M_half_over_Mperi') is not None]),
               M_settle_over_Mperi=st([r['M_settle_over_Mperi'] for r in rs if r.get('M_settle_over_Mperi') is not None]),
               mean_runtime_s=float(np.mean([r['runtime'] for r in rs])), mean_steps=float(np.mean([r['steps'] for r in rs])))
    w = [abs(r['windows']['10']['nminus'] - r['windows']['40']['nminus']) for r in rs if '10' in r['windows'] and '40' in r['windows']]
    out['max_abs_n10_minus_n40'] = float(max(w)) if w else None
    out['max_abs_actualmass_minus_gM'] = float(max(abs(r['windows'][max(r['windows'], key=float)]['nminus_actualmass'] - r['nminus']) for r in rs))
    return out


out = dict(meta=dict(model='homogeneous conformal X (two real fields) + phi; V = lamX/4 r^4 + g2/2 chi^2 r^2 + V_phi',
                     V_phi_M1='lamP/4 (chi^2 - M^2)^2 - epst M (chi - M)^3  (Codex GRUT-M1, d4747cb; epst = 0.002, lamP = 0.03)',
                     V_phi_mode1='lamP/4 chi^4 - lamP/2 M^2 chi^2 - beta lamP M^3 chi  (pt_matching skeptic model; beta = 3 epst/lamP = 0.2 is the linear piece of M1)',
                     constants=dict(lamX=LX, g2=G2, g=GG, lamP=LP, epst=EPST_M1, beta_equiv=3 * EPST_M1 / LP, m_phi_over_m_X=RPHI),
                     growth='M = M0 exp(kappa s), kappa = eps omega_r(initial orbit), M0 = M_start_fac * M_c,peri (0.5 unless stated), M_c = g r/sqrt(lamP)',
                     measure='n_-/q = <[(pX^2/2 + g^2 M^2 r^2/2)/(g M) - q]/(2q)>, N_phi/q = <(dchi_dot^2/2 + m_phi^2 dchi^2/2)/m_phi>/q, dchi = chi - M, m_phi^2 = 2 lamP M^2, window 60 X periods at M_end = 40 M_c,apo',
                     integrator='own adaptive DOP853 in C (m1_dop.c, copied from the skeptic sk_dop.c), rtol 1e-11 unless stated; work integral as extra ODE component',
                     label='[C] computed in a homogeneous classical model; not a cosmological result',
                     runs_dir='$SP/review/pt_m1tilt/runs'),
           orbits=orbits, groups={}, main_table={}, checks={}, controls={}, skeptic_reference={})

keys = sorted(set(key(r) for r in rows), key=str)
for k in keys:
    rs = [r for r in rows if key(r) == k]
    orb, mode, epst, beta, phi0, eps, rtol, mend, pert, m0 = k
    name = f'{orb}|mode{mode}|epst{epst:g}|beta{beta:g}|phi0{phi0}|eps{eps:.0e}|rtol{rtol:.0e}|mend{mend:g}|pert{pert:g}|M0fac{m0:g}'
    out['groups'][name] = block(rs)

# main table: M1 potential, phi0 = qs, rtol 1e-11, mend 40, no perturbation, M0 = 0.5 (rand + pair designs pooled; deterministic)
for orb in ['tau20', 'late']:
    tab = []
    for eps in sorted(set(r['eps'] for r in rows if r['orbit'] == orb and r['mode'] == 3 and r['epst'] == EPST_M1), reverse=True):
        rs = [r for r in rows if key(r) == (orb, 3, EPST_M1, 0.0, 'qs', eps, 1e-11, 40.0, 0.0, 0.5)]
        if rs:
            b = block(rs); tab.append(dict(eps=eps, **b))
    out['main_table'][orb] = tab

# phase smoothness (pair design: k/N vs k/N + 1e-3), seed perturbation, rtol, phi0, M0
chk = out['checks']
for orb in ['tau20', 'late']:
    for eps in sorted(set(r['eps'] for r in rows if r['orbit'] == orb), reverse=True):
        base = {round(r['phase'], 6): r for r in rows if key(r) == (orb, 3, EPST_M1, 0.0, 'qs', eps, 1e-11, 40.0, 0.0, 0.5)}
        d = []
        for ph, r in base.items():
            nb = round(ph + 1e-3, 6)
            if nb in base: d.append(abs(r['nminus'] - base[nb]['nminus']))
        if d: chk[f'{orb}_eps{eps:.0e}_phase_pairs_dphase1e-3'] = dict(n_pairs=len(d), max_abs_dnminus=float(max(d)), median_abs_dnminus=float(np.median(d)),
                                                                     nminus_at_grid_phases=[base[p]['nminus'] for p in sorted(base) if abs(p * 16 - round(p * 16)) < 1e-9])
        for pert in [1e-9]:
            pr = {round(r['phase'], 6): r for r in rows if key(r) == (orb, 3, EPST_M1, 0.0, 'qs', eps, 1e-11, 40.0, pert, 0.5)}
            d = [(ph, abs(pr[ph]['nminus'] - base[ph]['nminus']), abs(pr[ph]['Nphi'] - base[ph]['Nphi'])) for ph in pr if ph in base]
            if d: chk[f'{orb}_eps{eps:.0e}_seed_perturbation_rel{pert:g}'] = dict(n=len(d), max_abs_dnminus=float(max(x[1] for x in d)), max_abs_dNphi=float(max(x[2] for x in d)),
                                                                                 pairs=[dict(phase=x[0], dn=x[1]) for x in d])
        for rt in [1e-12, 1e-10]:
            pr = {round(r['phase'], 6): r for r in rows if key(r) == (orb, 3, EPST_M1, 0.0, 'qs', eps, rt, 40.0, 0.0, 0.5)}
            d = [abs(pr[ph]['nminus'] - base[ph]['nminus']) for ph in pr if ph in base]
            if d: chk[f'{orb}_eps{eps:.0e}_rtol{rt:.0e}_vs_1e-11'] = dict(n=len(d), max_abs_dnminus=float(max(d)))
        for p0 in ['0', '1e-06']:
            pr = {round(r['phase'], 6): r for r in rows if key(r) == (orb, 3, EPST_M1, 0.0, p0, eps, 1e-11, 40.0, 0.0, 0.5)}
            d = [abs(pr[ph]['nminus'] - base[ph]['nminus']) for ph in pr if ph in base]
            if d: chk[f'{orb}_eps{eps:.0e}_phi0_{p0}_vs_qs'] = dict(n=len(d), max_abs_dnminus=float(max(d)),
                                                                   chi0_over_Mperi_qs=[base[ph]['chi0_over_Mperi'] for ph in sorted(pr) if ph in base][:4])
        for m0 in [0.25]:
            pr = {round(r['phase'], 6): r for r in rows if key(r) == (orb, 3, EPST_M1, 0.0, 'qs', eps, 1e-11, 40.0, 0.0, m0)}
            if pr:
                d = [(base[ph]['nminus'], pr[ph]['nminus']) for ph in pr if ph in base]
                chk[f'{orb}_eps{eps:.0e}_M0fac{m0:g}'] = dict(group=block(list(pr.values())), runbyrun_M0half_vs_quarter=d)

# controls
ctl = out['controls']
for orb in ['tau20', 'late']:
    for eps in sorted(set(r['eps'] for r in rows), reverse=True):
        rs = [r for r in rows if key(r) == (orb, 3, 0.0, 0.0, '1e-06', eps, 1e-11, 40.0, 0.0, 0.5)]
        if rs: ctl[f'{orb}_eps{eps:.0e}_epst0_untilted'] = block(rs)
        rs = [r for r in rows if key(r) == (orb, 1, 0.0, 0.2, '1e-06', eps, 1e-11, 40.0, 0.0, 0.5)]
        if rs: ctl[f'{orb}_eps{eps:.0e}_mode1_linear_tilt_beta0.2'] = block(rs)
ctl['J_r_reference_per_q'] = dict(tau20=J0['tau20'], late=J0['late'])
# cross-check of the mode-3 code at epst = 0 against the skeptic's own deterministic eps = 1e-2 runs (same phase & seed)
SKR = SP + '/review/pt_matching_skeptic/runs'
cross = []
for nm, ph, sd in [('tau20', 0.0625, 1e-6), ('late', 0.0, 1e-2)]:
    fn = f'{SKR}/{nm}_g0_b0_eps1e-02_ph{ph:.6f}_seed{sd:.3e}_rt1e-11.json'
    if os.path.exists(fn):
        ref = json.load(open(fn)); y4, info = orbit(nm); stt = phase_state(y4, info['T_r'], ph)
        r3 = run_measure(stt, info, 1e-2, sd, mode=3, epst=0.0)
        cross.append(dict(orbit=nm, phase=ph, seed=sd, skeptic_nminus=ref['nminus'], skeptic_Nphi=ref['Nphi'], mode3_epst0_nminus=r3['nminus'], mode3_epst0_Nphi=r3['Nphi'],
                          dn=r3['nminus'] - ref['nminus'], dN=r3['Nphi'] - ref['Nphi']))
ctl['crosscheck_mode3_epst0_vs_skeptic_eps1e-2'] = cross

# skeptic reference numbers
if os.path.exists(SKEPTIC):
    sk = json.load(open(SKEPTIC))
    out['skeptic_reference'] = dict(pooled={k: dict(nminus=v['nminus'], Nphi=v['Nphi']) for k, v in sk['pooled'].items()},
                                    bias={k: dict(nminus=v['nminus'], Nphi=v['Nphi']) for k, v in sk['bias'].items()})

save_json('m1_transition.json', out)

# frequency diagnostics
fq = {}
for fn in sorted(glob.glob(os.path.join(SCR, 'freq_*.json'))):
    d = json.load(open(fn)); fq[os.path.basename(fn).replace('.json', '')] = d
dg = {}
for fn in sorted(glob.glob(os.path.join(SCR, 'diag_*.json'))):
    d = json.load(open(fn)); dg[os.path.basename(fn).replace('.json', '')] = d
save_json('m1_frequency_lock.json', dict(meta=dict(description='m1_freq.py: dominant Fourier frequency of r^2 and chi over 24-radial-period windows along M; '
                                                                'harmonic_power_fraction = share of chi power at harmonics of f_r (slaved response). m1_diag.py: free-basis '
                                                                'n_-/q, N_phi/q, frozen-orbit omega_r and linear m_phi along M.'), freq=fq, diag=dg))

# ---- print tables
print('\n== main table (M1, phi0 = qs, rtol 1e-11, M_end 40 M_apo, M0 = 0.5 M_peri)')
for orb in ['tau20', 'late']:
    print(orb, 'J0/q =', round(J0[orb], 5))
    for t in out['main_table'][orb]:
        n, N = t['nminus'], t['Nphi']
        print(f"  eps {t['eps']:.0e} n={n['n']:3d}  n-/q {n['mean']:.5f} sd {n['std']:.5f} [{n['min']:.5f}, {n['max']:.5f}]  1-n/J0 {t['one_minus_nminus_over_J0']:.2e}  "
              f"Nphi/q {N['mean']:.2e} sd {N['std']:.1e} [{N['min']:.1e}, {N['max']:.1e}]  sum {t['sum']['mean']:.4f}  W {t['max_work_rel']:.1e} q {t['max_qrel']:.1e}  "
              f"esc {t['M_escape_over_Mperi']['mean']:.3f} half {t['M_half_over_Mperi']['mean']:.3f} settle {t['M_settle_over_Mperi']['mean']:.3f}  t {t['mean_runtime_s']:.0f}s")
print('\n== checks')
for k, v in chk.items():
    if 'M0fac' in k:
        g = v['group']; print(f"  {k}: n-/q {g['nminus']['mean']:.4f} sd {g['nminus']['std']:.4f} (n {g['nminus']['n']})  Nphi {g['Nphi']['mean']:.4f}")
    else:
        print(f"  {k}: " + ', '.join(f'{a}={b:.2e}' if isinstance(b, float) else f'{a}={b}' for a, b in v.items() if a not in ('pairs', 'nminus_at_grid_phases', 'chi0_over_Mperi_qs')))
print('\n== controls')
for k, v in ctl.items():
    if isinstance(v, dict) and 'nminus' in v:
        print(f"  {k}: n-/q {v['nminus']['mean']:.4f} sd {v['nminus']['std']:.4f} [{v['nminus']['min']:.4f}, {v['nminus']['max']:.4f}] (n {v['nminus']['n']})  Nphi {v['Nphi']['mean']:.4f} sd {v['Nphi']['std']:.4f}  W {v['max_work_rel']:.1e}")
    else:
        print(' ', k, v)
