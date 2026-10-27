"""st_rerun: DOP853 tolerance audit with Codex's own C integrator (reproduction-type check).

(1) r10 vs r12 dense-time runs (saved by dense_crossings.py): max relative change of the per-radian amplitude and
    charge responses over tau 0..800, all 18 k.
(2) Linear-rescaling test: the tangent equation is linear, so multiplying the seed by s must multiply the response by s.
    Differences measure how much the absolute floor atol = 1e-3 rtol (vs relative control) affects the result.
(3) Tangent growth factor |v(tau)|/|v(0)| of the archived k=2.13 vector.
Usage: python3 -B tolerance_scaling.py <codex_code_dir> <published_dir> <dense_dir> <out_json>
"""
import sys, json
from pathlib import Path
sys.path.insert(0, sys.argv[1])
import timing_common as TC
from timing_common import np
pub, dense, outp = Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
out = {}
for b in ('b3', 'b5'):
    x = np.load(dense / f'time_{b}_r10.npz'); y = np.load(dense / f'time_{b}_r12.npz')
    A, B = x['amp'], y['amp']; t = x['times']
    scale = np.maximum(1., np.abs(B).max(axis=0))  # per-k running scale
    rel = np.abs(A - B) / (1 + np.abs(B))
    relpeak = np.abs(A - B).max(axis=0) / np.abs(B).max(axis=0)
    out[f'dense_{b}_r10_vs_r12'] = dict(max_scaled_amp=float(rel.max()), max_rel_to_peak_amp_over_k=float(relpeak.max()),
        max_rel_to_peak_charge_over_k=float((np.abs(x['charge'] - y['charge']).max(axis=0) / np.abs(y['charge']).max(axis=0)).max()))
    print(b, out[f'dense_{b}_r10_vs_r12'], flush=True)

exe = TC.build(); times = np.unique(np.r_[np.load(TC.M.ROOT / 'MODE_PROFILES.npz')['times'], np.arange(80., 800.00001, .5)])
orig = TC.phase_seed
ref = {}
for s in [1., 1e-6, 1e6]:
    TC.phase_seed = (lambda beta, step=1e-5, s=s: orig(beta, step) * s)
    bg, vec, ctl = TC.integrate(exe, f'scale_{s:g}', 3., .5, np.array([.02, 2.13, 5.8]), times, e=1e-7, rtol=2e-10)
    amp, dq, nx = TC.responses(3., bg, vec)
    ref[s] = (amp / s, dq / s, ctl['check'])
TC.phase_seed = orig
for s in [1e-6, 1e6]:
    a0, q0, _ = ref[1.]; a1, q1, chk = ref[s]
    out[f'seed_scale_{s:g}'] = dict(check=chk, max_rel_to_peak_amp=float((np.abs(a1 - a0).max(0) / np.abs(a0).max(0)).max()),
                                    max_rel_to_peak_charge=float((np.abs(q1 - q0).max(0) / np.abs(q0).max(0)).max()))
    print(s, out[f'seed_scale_{s:g}'], flush=True)
out['seed_scale_1_check'] = ref[1.][2]
z = np.load(pub / 'quiet3.npz'); sv = z['selected_vectors']; tt = z['times']
j = int(np.argmin(abs(z['selected_vector_ks'] - 2.13)))
norm = np.linalg.norm(sv[:, j], axis=1)
out['tangent_growth_k2.13_quiet3'] = dict(log_growth_at=[dict(tau=float(tt[i]), ln_norm_ratio=float(np.log(norm[i] / norm[0])))
                                                         for i in [np.argmin(abs(tt - x)) for x in (80, 134, 400, 772, 800)]])
print(out['tangent_growth_k2.13_quiet3'])
outp.write_text(json.dumps(out, indent=1) + '\n')
