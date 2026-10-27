"""st_rerun_skeptic analysis (own code) on OUR engine outputs (scratch npz) and OUR background (bg_own).

1. tolerance control: k=2.13, rtol 1e-12 vs 1e-10 (envelope-relative difference, crossing shift);
2. pointwise agreement of our responses with Codex's published quiet3/quiet5.npz arrays (read-only) at Codex's samples;
3. gauge sensitivity: Codex e=1e-7 (quiet) vs e=1e-6 (gauge) vs our e=0, at k=0.02, 2.13, 5.8;
4. growth of ln(envelope) per conformal time eta (own eta = int dtau/a), max over k and k=2.13 alone; [E] 1e-9 extrapolation;
5. decay clocks from OUR widths formula (rp.P couplings) and OUR a(tau), C, tau*; closed-form tail;
6. full-wavelength cube charge with OUR N_X normalisation.
Usage: python3 -B analyse.py <scratch_dir> <codex_packet_dir> <out_json>
"""
import sys, json, math
from pathlib import Path
import numpy as np
from scipy.special import gammaincc, gamma as G
L, packet, outp = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3] / 'rp_reheating' / 'code'))
import rp  # noqa: E402
P = rp.P
bgj = json.loads((HERE.parents[1] / 'results' / 'bg_own.json').read_text())
out = {}


def env_rel(x, ref):
    run = np.maximum.accumulate(np.abs(ref)); run = np.maximum(run, 1e-300)
    return float(np.max(np.abs(x - ref) / run))


def first(t, y, lev):
    w = np.flatnonzero(np.abs(y) >= lev)
    if not len(w): return None
    i = w[0]; y1, y2 = abs(y[i - 1]), abs(y[i]); return float(t[i - 1] + (lev - y1) / (y2 - y1) * (t[i] - t[i - 1]))


# widths, own formulas [D]
mL = math.sqrt(P['M'] ** 2 - math.sqrt(2) * P['mu'] * P['h'])
gL = 3 * P['kap'] ** 2 * mL / (256 * math.pi ** 3)
H_GeV = rp.H0
gam = {'YL': gL}
for e, lab in [(1e-6, 'A_e1e-6'), (1e-7, 'A_e1e-7')]:
    mA = 6 * e * P['h']
    gam[lab] = e * e * mA / (48 * math.pi)  # m_X/m_A -> 0; pinned (1-4mX^2/mA^2)^1.5 differs from 1 by <1e-12 relative
out['widths_own'] = dict(mL_over_H=mL, gamma=gam)

for tag in ('quiet3', 'quiet5'):
    row = {}
    m = np.load(L / f'{tag}_main_800_r1e-10_d0.005.npz'); t = m['t']; a = m['a']; amp = m['amp']; chg = m['charge']; ks = m['ks']
    j213 = int(np.argmin(abs(ks - 2.13)))
    # 1. tolerance control
    tt = np.load(L / f'{tag}_k213_800_r1e-12_d0.005.npz')
    assert np.allclose(tt['t'], t)
    row['tolerance_k213_r12_vs_r10'] = dict(amp_env_rel=env_rel(m['amp'][:, j213], tt['amp'][:, 0]),
                                            charge_env_rel=env_rel(chg[:, j213], tt['charge'][:, 0]),
                                            cross_003_r10=first(t, .03 * amp[:, j213], 1.), cross_003_r12=first(t, .03 * tt['amp'][:, 0], 1.),
                                            cross_1e6_r10=first(t, 1e-6 * amp[:, j213], 1.), cross_1e6_r12=first(t, 1e-6 * tt['amp'][:, 0], 1.))
    # 2. agreement with Codex's published arrays
    z = np.load(packet / f'{tag}.npz'); tc = z['times']; kc = z['ks']
    idx = np.round(tc * 200).astype(int); ok = np.abs(t[idx] - tc) < 1e-9
    agree = {}
    for k in (.02, .4, 2., 2.13, 2.4, 5.8, 8.):
        jc = int(np.argmin(abs(kc - k))); jo = int(np.argmin(abs(ks - k)))
        A_c = z['amplitude_per_radian'][ok, jc]; A_o = amp[idx[ok], jo]
        C_c = z['charge_fraction_per_radian'][ok, jc]; C_o = chg[idx[ok], jo]
        late = tc[ok] >= 80
        agree[str(k)] = dict(amp_env_rel=env_rel(A_o, A_c), charge_env_rel=env_rel(C_o, C_c),
                             amp_env_rel_tau_ge80=env_rel(A_o[late], A_c[late]),
                             peak_ratio=float(np.abs(A_o).max() / np.abs(A_c).max()))
    row['agreement_with_codex_published'] = agree
    # 3. gauge sensitivity (Codex quiet e=1e-7 vs gauge e=1e-6) and ours (e=0)
    gz = np.load(packet / ('gauge3.npz' if tag == 'quiet3' else 'gauge5.npz'))
    gs = {}
    for jg, k in enumerate(gz['ks']):
        jc = int(np.argmin(abs(kc - k)))
        gs[str(float(k))] = dict(codex_e1e6_vs_e1e7=env_rel(gz['amplitude_per_radian'][:, jg], z['amplitude_per_radian'][:, jc]))
    row['gauge_sensitivity_codex'] = gs
    # 4. growth per conformal time
    eta = np.r_[0, np.cumsum(np.diff(t) * .5 * (1 / a[1:] + 1 / a[:-1]))]
    fits = {}
    for name, series in [('max_over_k', np.abs(amp).max(axis=1)), ('k2.13', np.abs(amp[:, j213]))]:
        for lo, hi in [(150, 400), (400, 800), (150, 800)]:
            px, py = [], []
            for e0 in np.arange(lo, hi, 20.):
                w = (t >= e0) & (t < e0 + 20)
                i = np.argmax(series[w]); px.append(eta[w][i]); py.append(math.log(series[w][i]))
            fits[f'{name}_{lo}-{hi}'] = float(np.polyfit(px, py, 1)[0])
    s = fits['max_over_k_400-800']; late = t >= 700; env = np.abs(amp).max(axis=1)
    il = int(np.argmax(env[late])); lnA = math.log(env[late][il]); eta_r = float(eta[late][il]); a_r = float(a[late][il])
    C = bgj[tag]['C_own']; off = bgj[tag]['tau_star_own']
    deta = (math.log(1e9) - lnA) / s; a_c = a_r + C * deta
    row['growth'] = dict(eta_800=float(eta[-1]), slopes=fits,
                         E_1e9_unity=dict(tau=off + a_c ** 2 / (2 * C), a=a_c, label='[E] if the 400-800 envelope slope persisted; not computed'))
    # 5. clocks
    bg = np.load(L / f'bg_{tag}.npz'); tb = bg['tau']; ab = bg['a']
    cl = {}
    for name, g in gam.items():
        def early(tau):
            w = tb <= tau + 1e-12; return float(np.trapezoid(g * ab[w] * np.exp(-g * tb[w]), tb[w]))
        tail = lambda tau: math.sqrt(2 * C / g) * math.exp(-g * off) * G(1.5) * gammaincc(1.5, g * (tau - off))
        total = early(800.) + tail(800.)
        med = math.log(2) / g
        at = {}
        for lab, tau in [('codex_unity_sample', 134. if tag == 'quiet3' else 147.5), ('own_first_unity', None), ('own_10pct', None),
                         ('codex_1e6_unity_sample', 772. if tag == 'quiet3' else 782.)]:
            if tau is None:
                key = 'amp1' if lab == 'own_first_unity' else 'amp10'
                tau = json.loads((HERE.parents[1] / 'results' / f'tangent_{tag}_main_800_r1e-10_d0.005.json').read_text())['crossings']['18k_0.03'][key]['tau']
            at[lab] = dict(tau=tau, number_fraction=-math.expm1(-g * tau), energy_fraction=early(tau) / total)
        cl[name] = dict(gamma=g, median_tau=med, median_a=math.sqrt(2 * C * (med - off)), mean_a_per_decay=total,
                        asymptotic_mean_a=math.sqrt(math.pi * C / (2 * g)), at=at)
    row['clocks_own'] = cl
    norm = abs(bgj[tag]['NX_mean_64_80_dense'])
    row['Q_cube'] = {str(k): norm * (2 * math.pi / k) ** 3 for k in (2.13, 5.8)}
    out[tag] = row
    print(tag, json.dumps(row, indent=1, default=float)[:6000])
outp.write_text(json.dumps(out, indent=1, default=float) + '\n')
