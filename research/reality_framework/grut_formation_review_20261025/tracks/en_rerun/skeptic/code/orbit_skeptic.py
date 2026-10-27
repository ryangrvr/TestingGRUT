"""en_rerun_skeptic -- own code, Codex modules NOT imported.
Re-derive the orbit-phase statements (en_rerun N1, 1a/1b qualifications, claim-2 window means) with a different method:
 * traces: en_rerun's reproduced state.tsv (Codex engine) -- first checked against Codex's PUBLISHED FIELD_CURVES.npz;
 * cycles: apocentre-to-apocentre AND pericentre-to-pericentre, extrema refined by parabolic interpolation;
 * averages: conformal-time (eta) trapezoid weighting, not sample (tau) means;
 * late window: an integer number of whole orbits ending before tau 800, compared with the plain tau 700-800 mean.
usage: python3 -I orbit_skeptic.py <published_FIELD_CURVES.npz> <early_workdir> <out_json> tag [tag...]"""
import sys, json, os
import numpy as np

LX, G2, LP = .001, .01, .03
curves_path, work, outp = sys.argv[1:4]
tags = sys.argv[4:]
C = np.load(curves_path)


def derived(t):
    rho = t['rho']
    mr2 = t['q0']**2 + t['q1']**2
    EX0 = .5*(t['p0']**2 + t['p1']**2) + LX/4*mr2**2 + G2/2*t['qphi']**2*mr2
    EP0 = .5*t['pphi']**2 + LP/4*t['qphi']**4
    return dict(
        mr2=mr2,
        outside=1 - (EX0 + EP0)/rho,
        phimix=(t['Kphi'] + t['Gphi'] + t['Uphi'] + t['Umix'])/rho,
        ephi=(t['Kphi'] + t['Gphi'] + t['Uphi'])/rho,
        meanQ=(t['q0']*t['p1'] - t['q1']*t['p0'])/t['nQ'],
        contrast=np.sqrt(np.maximum(t['charge_RMS']**2 - t['nQ']**2, 0))/abs(t['nQ']),
        fluct=np.maximum(t['X2'] - mr2, 0)/t['X2'],
        dX=t['X2'] - mr2, X2=t['X2'])


def extrema(y, kind):
    """indices of local minima/maxima of a smooth sampled curve (strict, window 1) after a 3-sample
    smoothing guard; keep only extrema that are global within +-w samples where w = 1/4 of the local period."""
    s = -y if kind == 'max' else y
    i = np.flatnonzero((s[1:-1] < s[:-2]) & (s[1:-1] <= s[2:])) + 1
    # discard shallow wiggles: require the extremum to be the extreme over +-60 samples (3 tau at dtau .05)
    keep = [j for j in i if s[j] == s[max(0, j-60):j+61].min()]
    return np.array(keep, int)


def refine(x, y, j):
    if j <= 0 or j >= len(y)-1:
        return x[j], y[j]
    y0, y1, y2 = y[j-1], y[j], y[j+1]
    den = y0 - 2*y1 + y2
    d = .5*(y0 - y2)/den if den != 0 else 0.
    return x[j] + d*(x[j+1]-x[j]), y1 - .25*(y0-y2)*d


def eta_mean(eta, v, a, b):
    m = (eta >= eta[a]) & (eta <= eta[b])
    return float(np.trapezoid(v[m], eta[m])/(eta[m][-1] - eta[m][0]))


out = {}
for tag in tags:
    d = os.path.realpath(f'{work}/{tag}')
    t = np.genfromtxt(d + '/state.tsv', names=True, delimiter='\t')
    if tag == 'halfstep48':  # keep only the .05 grid for comparability of extremum windows
        pass
    D = derived(t)
    tau, eta = t['tau'], t['eta']
    # 1. validate the reproduced trace against Codex's published curves
    pt = C[tag + '_tau']
    stride = 1 if tag == 'N48' else (2 if tag == 'halfstep48' else 10)
    idx = np.searchsorted(tau, pt - 1e-7)
    idx = np.clip(idx, 0, len(tau)-1)
    ok = abs(tau[idx] - pt) < 1e-6
    val = {}
    for key, mine in [('energy_outside_spatial_mean_fraction', D['outside']),
                      ('phi_plus_interaction_energy_fraction', D['phimix']),
                      ('spatial_mean_charge_fraction', D['meanQ']),
                      ('charge_contrast_RMS', D['contrast']),
                      ('X_field_fluctuation_variance_fraction', D['fluct'])]:
        pub = C[tag + '_' + key]
        diff = np.max(abs(mine[idx[ok]] - pub[ok]))/max(np.max(abs(pub[ok])), 1e-300)
        val[key] = float(diff)
    # 2. orbit extrema of |chibar|^2
    w = 120 if tag == 'halfstep48' else 60
    mr2 = D['mr2']
    per = extrema(mr2, 'min')
    apo = extrema(mr2, 'max')
    peri_tau = [refine(tau, mr2, j) for j in per]
    apo_tau = [refine(tau, mr2, j) for j in apo]
    # 3. cycles (both conventions), eta-weighted averages
    def cycles(ext):
        rows = []
        for a, b in zip(ext[:-1], ext[1:]):
            rows.append(dict(tau0=float(tau[a]), tau1=float(tau[b]),
                             fluct_ratio_of_means=eta_mean(eta, D['dX'], a, b)/eta_mean(eta, D['X2'], a, b),
                             fluct_mean_of_ratio=eta_mean(eta, D['fluct'], a, b),
                             contrast_mean=eta_mean(eta, D['contrast'], a, b),
                             outside_mean=eta_mean(eta, D['outside'], a, b),
                             meanQ_mean=eta_mean(eta, D['meanQ'], a, b),
                             phimix_mean=eta_mean(eta, D['phimix'], a, b),
                             mr2_min=float(mr2[a:b+1].min()), mr2_max=float(mr2[a:b+1].max())))
        return rows
    cp, ca = cycles(per), cycles(apo)
    first_half = float(tau[np.argmax(D['fluct'] >= .5)])
    first_c1 = float(tau[np.argmax(D['contrast'] >= 1)])
    def first_cycle(rows, key, lev):
        r = [c for c in rows if c[key] >= lev]
        return r[0] if r else None
    # 4. late windows
    w78 = (tau >= 700) & (tau <= 800)
    i700, i800 = np.flatnonzero(w78)[[0, -1]]
    whole_p = [c for c in cp if c['tau1'] <= 800.0001]
    whole_a = [c for c in ca if c['tau1'] <= 800.0001]
    last2_p = whole_p[-2:]
    last2_a = whole_a[-2:]
    def span_mean(rows, key):
        a = np.searchsorted(tau, rows[0]['tau0']-1e-9); b = np.searchsorted(tau, rows[-1]['tau1']-1e-9)
        return eta_mean(eta, D[key], a, b)
    late = dict(
        tau700_800_sample_mean=dict(outside=float(D['outside'][w78].mean()), meanQ=float(D['meanQ'][w78].mean()),
                                    phimix=float(D['phimix'][w78].mean())),
        tau700_800_eta_mean=dict(outside=eta_mean(eta, D['outside'], i700, i800), meanQ=eta_mean(eta, D['meanQ'], i700, i800),
                                 phimix=eta_mean(eta, D['phimix'], i700, i800)),
        last_two_whole_orbits_peri=dict(tau0=last2_p[0]['tau0'], tau1=last2_p[-1]['tau1'],
                                        outside=span_mean(last2_p, 'outside'), meanQ=span_mean(last2_p, 'meanQ'),
                                        phimix=span_mean(last2_p, 'phimix')),
        last_two_whole_orbits_apo=dict(tau0=last2_a[0]['tau0'], tau1=last2_a[-1]['tau1'],
                                       outside=span_mean(last2_a, 'outside'), meanQ=span_mean(last2_a, 'meanQ'),
                                       phimix=span_mean(last2_a, 'phimix')),
        last_whole_orbit_mr2_range=[last2_p[-1]['mr2_min'], last2_p[-1]['mr2_max']],
        final=dict(tau=float(tau[-1]), outside=float(D['outside'][-1]), meanQ=float(D['meanQ'][-1]),
                   phimix=float(D['phimix'][-1]), ephi=float(D['ephi'][-1]), mr2=float(mr2[-1]),
                   fluct=float(D['fluct'][-1])),
        since_last_pericentre=dict(tau0=float(tau[per[-1]]), outside=[float(D['outside'][per[-1]:].min()), float(D['outside'][per[-1]:].max())],
                                   meanQ=[float(D['meanQ'][per[-1]:].min()), float(D['meanQ'][per[-1]:].max())]),
        orbit_tau700_800=dict(outside=[float(D['outside'][w78].min()), float(D['outside'][w78].max())],
                              meanQ=[float(D['meanQ'][w78].min()), float(D['meanQ'][w78].max())],
                              phimix=[float(D['phimix'][w78].min()), float(D['phimix'][w78].max())]))
    out[tag] = dict(
        validation_vs_published_FIELD_CURVES_maxrel=val,
        pericentres_first8=[(round(float(a), 3), round(float(b), 2)) for a, b in peri_tau[:8]],
        apocentres_first8=[(round(float(a), 3), round(float(b), 1)) for a, b in apo_tau[:8]],
        pericentres_last4=[(round(float(a), 3), round(float(b), 1)) for a, b in peri_tau[-4:]],
        apocentres_last4=[(round(float(a), 3), round(float(b), 1)) for a, b in apo_tau[-4:]],
        n_peri=len(per), n_apo=len(apo),
        first_instantaneous_X_variance_half=first_half,
        first_instantaneous_contrast_one=first_c1,
        cycle_peri_containing_first_half=[c for c in cp if c['tau0'] <= first_half < c['tau1']],
        cycle_apo_containing_first_half=[c for c in ca if c['tau0'] <= first_half < c['tau1']],
        first_peri_cycle_contrast_mean_ge1=first_cycle(cp, 'contrast_mean', 1.),
        first_apo_cycle_contrast_mean_ge1=first_cycle(ca, 'contrast_mean', 1.),
        max_cycle_fluct_ratio_of_means_peri=max(c['fluct_ratio_of_means'] for c in cp),
        max_cycle_fluct_mean_of_ratio_peri=max(c['fluct_mean_of_ratio'] for c in cp),
        first_peri_cycle_fluct_mean_of_ratio_ge_half=first_cycle(cp, 'fluct_mean_of_ratio', .5),
        late=late,
        cycles_peri=cp, cycles_apo=ca)
    o = out[tag]
    print(tag, 'valid', {k: f'{v:.1e}' for k, v in val.items()})
    print('  peri', o['pericentres_first8'][:6], 'last', o['pericentres_last4'])
    print('  apo ', o['apocentres_first8'][:5], 'last', o['apocentres_last4'])
    print('  firsthalf', first_half, 'cyc', [(round(c['tau0'], 2), round(c['tau1'], 2), round(c['fluct_ratio_of_means'], 4), round(c['fluct_mean_of_ratio'], 4)) for c in o['cycle_peri_containing_first_half']],
          'apo', [(round(c['tau0'], 2), round(c['tau1'], 2), round(c['fluct_ratio_of_means'], 4), round(c['fluct_mean_of_ratio'], 4)) for c in o['cycle_apo_containing_first_half']])
    fc = o['first_peri_cycle_contrast_mean_ge1']; fa = o['first_apo_cycle_contrast_mean_ge1']
    print('  contrast>=1 peri-cycle', fc and (round(fc['tau0'], 2), round(fc['tau1'], 2), round(fc['contrast_mean'], 3)),
          'apo-cycle', fa and (round(fa['tau0'], 2), round(fa['tau1'], 2), round(fa['contrast_mean'], 3)))
    print('  max cycle fluct: ratio-of-means', round(o['max_cycle_fluct_ratio_of_means_peri'], 4), 'mean-of-ratio', round(o['max_cycle_fluct_mean_of_ratio_peri'], 4),
          'first cyc mean-of-ratio>=.5', o['first_peri_cycle_fluct_mean_of_ratio_ge_half'] and o['first_peri_cycle_fluct_mean_of_ratio_ge_half']['tau0'])
    print('  late', json.dumps({k: v for k, v in late.items()}, default=float))
json.dump(out, open(outp, 'w'), indent=1, default=float)
