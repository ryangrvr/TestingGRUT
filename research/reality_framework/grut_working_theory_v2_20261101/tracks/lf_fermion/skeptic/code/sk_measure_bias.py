"""Skeptic S3: energy-measure and grid-refinement bias of the line fermion energy used by lf_fermion (Q = e_proj_cut
minus the t < 0.15 mean). Claude, lf_fermion skeptic (issue #13 context). Read-only on lf_fermion's scratch outputs;
the skeptic's own no-refinement runs (sk_line_drv.py: Nc = Nf, no zero-padding at t_ref) are in the skeptic scratch.
Clean estimator [model choice]: chiral measure e_chir_cut minus a vacuum residual fitted on t < 0.15 as C / a^2
(the residual of the discrete-time map scales as k_Nyq^2 / a^2; checked by the early-time drift and by the grid ratio).
Usage: python3 -B sk_measure_bias.py -> ../results/sk_measure_bias.json"""
import json, math
from pathlib import Path
import numpy as np
A = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2')
LQ, SQ = A / 'lf_fermion/lfq/REC_P0', A / 'lf_fermion_skeptic/lfq/REC_P0'
V = 367046.6168416639
TS = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]


def series(f, nm):
    m = json.load(open(f))['meas']
    t = np.array([r['t'] for r in m]); a = np.array([r['a'] for r in m])
    ep = np.array([r[nm]['e_proj_cut'] for r in m]) / V ** 2
    ec = np.array([r[nm]['e_chir_cut'] for r in m]) / V ** 2
    early = t < 0.15
    q_lf = ep - ep[early].mean()                                       # lf_fermion's Q
    C = np.sum(ec[early] / a[early] ** 2) / np.sum(1 / a[early] ** 4)   # least squares ec = C / a^2
    q_cl = ec - C / a ** 2
    resid = float(np.std(ec[early] - C / a[early] ** 2) / np.mean(ec[early]))
    pick = lambda y: {tt: float(y[np.argmin(abs(t - tt))]) for tt in TS}
    return dict(q_lf=pick(q_lf), q_clean=pick(q_cl), C=float(C), fit_rel_resid=resid)


out = {}
cases = {'N1': ('l3_sig', ['fine8k', 'fine8ks2', 'fine16k']), 'N2': ('l1_chi', ['fine8k', 'fine8ks2']),
         'N3': ('l1_chi', ['fine8k', 'fine8ks2']), 'N4': ('l1_chi', ['fine8k', 'fine8ks2'])}
for nm, (ln, clean_runs) in cases.items():
    base = series(LQ / ('base_%s.json' % ln), nm)
    cl = {r: series(SQ / ('%s_%s.json' % (r, ln)), nm) for r in clean_runs}
    if nm == 'N1':
        cl['nf16k(refined)'] = series(SQ / ('nf16k_%s.json' % ln), nm)
    else:
        cl['nf16k(refined)'] = series(LQ / ('nf16k_%s.json' % ln), nm)
    cmean = {tt: float(np.mean([cl[r]['q_clean'][tt] for r in clean_runs])) for tt in TS}
    out[nm] = dict(line=ln, base_Q_lf=base['q_lf'], clean_runs={r: cl[r] for r in cl}, clean_mean=cmean,
                   bias_base_over_clean={tt: base['q_lf'][tt] / cmean[tt] for tt in TS})
    print(nm, ln, 'C(base,8k,..)=', '%.2f' % base['C'], ' '.join('%s %.2f' % (r, cl[r]['C']) for r in cl))
    print('   lf Q (base proj):   ', ' '.join('%.1f' % base['q_lf'][tt] for tt in TS))
    for r in cl:
        print('   clean %-15s' % r, ' '.join('%.1f' % cl[r]['q_clean'][tt] for tt in TS))
    print('   bias base/clean-mean', ' '.join('t%.1f %.3f' % (tt, out[nm]['bias_base_over_clean'][tt]) for tt in TS))
Path(__file__).resolve().parents[1].joinpath('results', 'sk_measure_bias.json').write_text(
    json.dumps(dict(label='C', note=__doc__, cases=out), indent=1))
