"""st_physics step 4b: what the supplied 'phase' seed is, and its growing-mode weight relative to a genuine X-phase
rotation.  Reads the unit-column transfer output written by s4_vacuum_seed.py (Codex engine, own seeds) from the
GRUT_TIMING_WORKDIR; own post-processing.

usage: python3 -B s4b_seed_projection.py <codex timing code dir> <out.json>   (GRUT_* env as in BRIEF_TIMING.md)
"""
import sys, json
import numpy as np

code_dir, out = sys.argv[1:3]
sys.path.insert(0, code_dir)
import timing_common as T  # noqa: E402

KV = np.array([1.8, 2.0, 2.13, 2.4, 2.8])
NAMES = ['qXr', 'qXi', 'qYr', 'qYi', 'qSr', 'qSi', 'pXr', 'pXi', 'pYr', 'pYi', 'pSr', 'pSi']
res = dict(label='engine output: Codex response.c with own unit seeds; analysis own')
for tag, beta in (('quiet3', 3.), ('quiet5', 5.)):
    pub = np.load(T.ROOT / f'{tag}.npz'); times = pub['times']
    nk = len(KV) * 13
    data = np.fromfile(T.WORK / (tag + '_stphys.bin'), dtype='<f8').reshape(len(times), 20 + 12 * nk)
    bg = data[:, 1:20]; vec = data[:, 20:].reshape(len(times), nk, 12)
    amp, dq, nx = T.responses(beta, bg, vec)
    v, _ = T.M.prepared(beta)
    chi = T.M.realvec(v[2:8:2] + 1j * v[3:8:2])
    seed = T.phase_seed(beta)
    rX = float(np.linalg.norm(chi[:2]))
    comp = dict(abs_chi_X=rX, dlnabsX_per_rad=float(chi[:2] @ seed[:2] / rX ** 2),
                dtheta_X_per_rad=float((chi[0] * seed[1] - chi[1] * seed[0]) / rX ** 2),
                dtheta_Y_per_rad=float((chi[2] * seed[3] - chi[3] * seed[2]) / (chi[2:4] @ chi[2:4])),
                dtheta_Sigma_per_rad=float((chi[4] * seed[5] - chi[5] * seed[4]) / (chi[4:6] @ chi[4:6])),
                note='seed = d(prepared state)/d(preparation phase): Sigma rotated, X kept real, Y and |X| readjusted')
    late = times >= 600
    per_k = []
    for i, k in enumerate(KV):
        A = amp[:, 13 * i:13 * i + 13]
        env = np.max(abs(A[late]), axis=0)
        xphase = rX * (A[:, 1] + A[:, 7])   # local X rotation by 1 rad: dq = dp = |chi| e_perp (chi' = chi at a=1)
        xrad = rX * (A[:, 0] + A[:, 6])     # d|X|/|X| = 1 with the same velocity rule
        per_k.append(dict(k=float(k), unit_column_gain_vs_supplied=dict(zip(NAMES, (env[:12] / env[12]).tolist())),
                          Xphase_1rad_gain_vs_supplied_1rad=float(np.max(abs(xphase[late])) / env[12]),
                          Xradial_unit_gain_vs_supplied_1rad=float(np.max(abs(xrad[late])) / env[12])))
    res[tag] = dict(seed_components=comp, per_k=per_k)
    print(tag, comp, [(p['k'], round(p['Xphase_1rad_gain_vs_supplied_1rad'], 1)) for p in per_k], flush=True)
json.dump(res, open(out, 'w'), indent=1)
