"""Compare own independent N48 history (indep3d.py, BM6) with Codex's published N48 curves and RESULTS.json.
Writes results/en_indep_results.json. Own code; Codex files are read as data only."""
import sys, os, json, struct, hashlib
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import indep3d as I

SP = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad'
PKT = f'{SP}/review/en_indep/pkt/grut_early_nonlinear_20261027'      # pristine copy of the published packet
RUNS = f'{SP}/review/en_indep/runs'
OUT = sys.argv[1]

def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()

curves = np.load(f'{PKT}/FIELD_CURVES.npz')
res = json.load(open(f'{PKT}/RESULTS.json'))['final']
raw = {r['path']: r['sha256'] for r in json.load(open(f'{PKT}/RAW_FILES.json'))['files']}
CK = dict(charge_contrast='charge_contrast_RMS', Xvar_fraction='X_field_fluctuation_variance_fraction',
          outside_mean='energy_outside_spatial_mean_fraction', Ephi_fraction='Ephi_fraction',
          phi_plus_mix_fraction='phi_plus_interaction_energy_fraction', mean_charge_fraction='spatial_mean_charge_fraction')
ct = curves['N48_tau'][:-1]          # last entry duplicates tau=800
cv = {k: curves['N48_' + c][:-1] for k, c in CK.items()}
cv['E'] = curves['N48_rho'][:-1]; cv['nQ'] = curves['N48_nQ'][:-1]

runs = {t: I.derived(np.load(f'{RUNS}/{t}/records.npy')) for t in ('bm6_h004', 'bm6_h002', 'twin_eps1e-6')}
A, B, W = runs['bm6_h004'], runs['bm6_h002'], runs['twin_eps1e-6']
out = {}

# ---- provenance ----
init = f'{SP}/review/en_indep/ework/N48_initial.bin'
out['input'] = dict(initial_file=init, initial_sha256=sha(init), published_RAW_FILES_sha256=raw['N48/initial.bin'],
                    hash_match=sha(init) == raw['N48/initial.bin'],
                    MATCH_L12_sha256_used=sha(f'{PKT}/MATCH_L12.npz'),
                    note='initial field built by Codex initial.py run in a pristine packet copy (reproduction of the INPUT only)')

# ---- clock ----
m = np.load(f'{PKT}/MATCH_L12.npz'); bg = m['background']
a_bg = float(np.exp(bg[-1, 15])); H_bg = float(bg[-1, 16])
cs = curves['N48_eta'][:-1]
out['clock'] = dict(a_m=float(m['a_match']), C=float(m['C_match']), a_from_background=a_bg, H_from_background=H_bg,
                    C_check_a2H=a_bg * a_bg * H_bg, C_rel_diff=abs(a_bg * a_bg * H_bg / float(m['C_match']) - 1),
                    radiation_origin_tau=20 - 1 / (2 * H_bg),
                    codex_tau_vs_own_tau_of_s_maxabs=float(np.max(abs(I.tau_of_s(cs, float(m['a_match']), float(m['C_match'])) - ct))),
                    derivation='a(s)=a_m+C s (radiation, a\'\'=0); dtau=a ds => tau=20+a_m s+C s^2/2; C=da/ds=a^2 (dlna/dtau)=a^2 H')

# ---- conservation ----
def ledger(D, upto):
    """trapezoid exchange-work identity E_X(s)-E_X(0) = int g2 phi pi_phi r ds (own records)."""
    r = np.load(f'{RUNS}/{D}/records.npy'); d = {c: r[:, i] for i, c in enumerate(I.COLS)}
    sel = d['tau'] <= upto + 1e-9
    EX = (d['E'] - d['Kphi'] - d['Gphi'] - d['Uphi'])[sel]; w = d['work'][sel]; s = d['s'][sel]
    integ = np.r_[0, np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(s))]
    return float(np.max(abs(EX - EX[0] - integ)) / d['E'][0])
out['conservation'] = {t: dict(max_energy_rel=float(np.max(abs(D['energy_drift']))), max_charge_rel=float(np.max(abs(D['charge_drift']))),
                               work_ledger_rel_fine_part_tau_le_150=ledger(t, 150.0))
                       for t, D in runs.items()}
out['conservation']['codex_N48_curve_max_energy_rel'] = float(np.max(abs(cv['E'] / cv['E'][0] - 1)))

# ---- crossing times ----
def crossing(tau, v, level):
    i = np.flatnonzero(v >= level)
    if not len(i): return None, None
    i = i[0]; t0, t1, v0, v1 = tau[i - 1], tau[i], v[i - 1], v[i]
    return float(tau[i]), float(t0 + (level - v0) * (t1 - t0) / (v1 - v0))
cr = {}
for lab, D in (('own_h004', A), ('own_h002', B)):
    cr[lab] = dict(charge_contrast_one=crossing(D['tau'], D['charge_contrast'], 1.0), X_variance_half=crossing(D['tau'], D['Xvar_fraction'], 0.5))
cr['codex_N48_curve'] = dict(charge_contrast_one=crossing(ct, cv['charge_contrast'], 1.0), X_variance_half=crossing(ct, cv['Xvar_fraction'], 0.5))
cr['codex_reported'] = {t: dict(charge_contrast_one=r['first_charge_contrast_one'], X_variance_half=r['first_X_variance_half']) for t, r in res.items()}
cr['format'] = '[first recorded sample >= level, linear interpolation between samples]'
out['crossings'] = cr

# ---- values at fixed tau ----
def at(D, tau, key):
    i = np.argmin(abs(D['tau'] - tau)); assert abs(D['tau'][i] - tau) < 1e-6; return float(D[key][i])
def cat(tau, key):
    i = np.argmin(abs(ct - tau)); assert abs(ct[i] - tau) < 1e-6; return float(cv[key][i])
vals = {}
for tau in (60, 70, 100, 150, 200, 400, 800):
    row = {}
    for k in CK:
        o, o2, c = at(A, tau, k), at(B, tau, k), cat(tau, k)
        row[k] = dict(own_h004=o, own_h002=o2, codex_N48=c, rel_own_vs_codex=abs(o / c - 1), rel_own_h004_vs_h002=abs(o / o2 - 1))
    vals[str(tau)] = row
out['values_at_tau'] = vals
out['codex_reported_final_N48'] = {k: res['N48'][k] for k in ('energy_outside_spatial_mean', 'Ephi', 'Ephi_plus_mix', 'spatial_mean_charge_over_net')}

# ---- pointwise divergence (diagnostics) ----
def divergence(D, ref_tau, ref):
    i = np.searchsorted(ref_tau, D['tau'] - 1e-9); assert np.allclose(ref_tau[i], D['tau'], atol=1e-7)
    o = {}
    for k in CK:
        r = abs(D[k] - ref[k][i]) / abs(ref[k][i])
        segs = {}
        for lo, hi in ((20, 70), (70, 150), (150, 400), (400, 800)):
            sel = (D['tau'] >= lo) & (D['tau'] <= hi)
            if sel.any(): segs[f'{lo}-{hi}'] = float(r[sel].max())
        first = {}
        for e in (1e-6, 1e-4, 1e-2, 1e-1):
            j = np.flatnonzero(r > e); first[f'{e:g}'] = float(D['tau'][j[0]]) if len(j) else None
        o[k] = dict(max_rel_by_segment=segs, first_tau_rel_exceeds=first)
    return o
out['diagnostic_divergence_own_h004_vs_codex_N48'] = divergence(A, ct, cv)
out['diagnostic_divergence_own_h002_vs_codex_N48'] = divergence(B, ct, cv)
out['diagnostic_divergence_own_h004_vs_own_h002'] = divergence(A, B['tau'], B)
out['diagnostic_divergence_twin_vs_own_h004_from_tau100'] = divergence(W, A['tau'], A)

# ---- late statistics ----
late = {}
for lo, hi in ((400, 800), (700, 800)):
    seg = {}
    for k in CK:
        sa = (A['tau'] >= lo) & (A['tau'] <= hi)
        i = np.searchsorted(ct, A['tau'][sa] - 1e-9)
        seg[k] = dict(own_mean=float(A[k][sa].mean()), own_std=float(A[k][sa].std()), own_min=float(A[k][sa].min()), own_max=float(A[k][sa].max()),
                      codex_mean_same_grid=float(cv[k][i].mean()), codex_std_same_grid=float(cv[k][i].std()))
    late[f'{lo}-{hi}'] = seg
spread = {k: [min(r[k] for r in res.values()), max(r[k] for r in res.values())] for k in ('energy_outside_spatial_mean', 'Ephi_plus_mix', 'spatial_mean_charge_over_net', 'Ephi')}
late['codex_control_spread_at_800'] = spread
out['late_statistics'] = late

# ---- own reference projections at tau=800 (full complex FFT, own code) ----
from refproj import refproj
z = np.load(f'{RUNS}/bm6_h004/snap_tau800.00.npz')
out['reference_projections_tau800_own'] = refproj(z['q'], z['p'])
out['reference_projections_tau800_codex_reported'] = res['N48']['raw_nonzero_reference_nminus_over_net']

json.dump(out, open(OUT, 'w'), indent=1, default=float)
print(json.dumps({k: out[k] for k in ('input', 'clock', 'conservation', 'crossings', 'reference_projections_tau800_own', 'reference_projections_tau800_codex_reported')}, indent=1, default=float))
for tau in ('100', '200', '400', '800'):
    print(tau, {k: (round(v['own_h004'], 7), round(v['codex_N48'], 7), f"{v['rel_own_vs_codex']:.1e}") for k, v in vals[tau].items()})
print(json.dumps(late, indent=1))
for k in ('diagnostic_divergence_own_h004_vs_codex_N48', 'diagnostic_divergence_own_h004_vs_own_h002', 'diagnostic_divergence_twin_vs_own_h004_from_tau100'):
    print(k); [print(' ', kk, json.dumps(vv)) for kk, vv in out[k].items()]
