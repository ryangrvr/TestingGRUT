"""SKEPTIC (pt_response_skeptic) probe runner: low-noise conjugate/normal shell traces at further orbital phases.

ENGINE: Codex's fieldstep tangent (exact derivative of the discrete Yoshida map), imported from our own copy of the
5f5513b tree ($PKC/code).  Background states: the regenerated base48 field (bitwise equal to FIELD_TRACES).
PROBE DESIGN, PROJECTIONS AND ESTIMATORS: ours, different from both response.py (one white-noise probe over all
shells) and pt_response/probe_study.py (Hutchinson over shells 1-4 jointly):
  family A: M random-phase unit-modulus probes supported on shell 1 ONLY (no cross-shell leakage noise in shell 1),
  family B: M random-phase unit-modulus probes supported on shells 2-4.
Each probe kicks pi_1 and pi_2 separately (2 tangents per probe).  All 4M tangents ride on one background.
Delays 0..6 step 0.04 (151 samples), hmax 0.004 as response.py.  Chunked with restart files.
Usage: python3 -B probes_phase.py <s0> <M> <outdir> <max_samples_this_call> [seed] [families: AB (default) | A | B]
"""
import sys, os, time, pickle
sys.path.insert(0, '.')
from pt_common import *  # noqa  (Codex module from our copy; engine only)

s0 = float(sys.argv[1]); M = int(sys.argv[2]); outdir = Path(sys.argv[3]); budget = int(sys.argv[4])
seed = int(sys.argv[5]) if len(sys.argv) > 5 else 424242
FAMS = sys.argv[6] if len(sys.argv) > 6 else 'AB'
outdir.mkdir(parents=True, exist_ok=True)
N, L = 48, 12.
tag = f'skp_s{s0:g}_M{M}_seed{seed}' + ('' if FAMS == 'AB' else f'_{FAMS}')
ck = outdir / f'ckpt_{tag}.pkl'
mm = np.fft.fftfreq(N) * N
nx, ny, nz = np.meshgrid(mm, mm, mm, indexing='ij')
kc = 2 * np.pi / L * np.sqrt(nx ** 2 + ny ** 2 + nz ** 2)
shell = np.floor(kc / (2 * np.pi / L) + .5).astype(int).ravel()
allmodes = np.flatnonzero((shell >= 1) & (shell <= 4))

def pos_half(i, j, k):
    v = (mm[i], mm[j], mm[k]); first = next((x for x in v if x != 0), 0); return first > 0

if not ck.exists():
    z = np.load(WORK / 'base48' / f'state_s{s0:g}.npy').copy()
    rng = np.random.default_rng([seed, int(round(s0 * 100)), M])
    pats = []; fams = []
    for fam, shells_ in [x for x in (('A', (1,)), ('B', (2, 3, 4))) if x[0] in FAMS]:
        for j in range(M):
            coef = np.zeros((N, N, N), complex)
            for ix in np.flatnonzero(np.isin(shell, shells_)):
                i, jj, k = np.unravel_index(ix, (N, N, N))
                if not pos_half(i, jj, k): continue
                ph = rng.uniform(0, 2 * np.pi)
                coef[i, jj, k] = np.exp(1j * ph); coef[(-i) % N, (-jj) % N, (-k) % N] = np.exp(-1j * ph)
            pat = ifftn(coef, norm='forward'); assert np.max(abs(pat.imag)) < 1e-9
            pats.append(np.ascontiguousarray(pat.real)); fams.append(fam)
    nt = 2 * len(pats)
    tangent = np.zeros((nt, 6, N, N, N))
    for p, pat in enumerate(pats):
        tangent[2 * p, 3] = pat; tangent[2 * p + 1, 4] = pat
    srcs = np.array([fftn(p, norm='forward').ravel()[allmodes] for p in pats])
    state = dict(z=z, tangent=tangent, j=0, prev=0., fams=fams, srcs=srcs, proj=[])
else:
    state = pickle.load(open(ck, 'rb'))
set_box(L); lib = library()
times = np.arange(151) * .04; t0 = time.monotonic(); done = 0
z = state['z']; tangent = state['tangent']
while state['j'] < len(times) and done < budget:
    t = times[state['j']]
    if state['j']: step(lib, z, t - state['prev'], .004, tangent=tangent, threads=2)
    state['prev'] = t
    f = fftn(tangent[:, :2], axes=(-3, -2, -1), norm='forward', workers=2).reshape(tangent.shape[0], 2, -1)[:, :, allmodes]
    state['proj'].append(f.astype(np.complex128))
    state['j'] += 1; done += 1
pickle.dump(state, open(ck, 'wb'), protocol=4)
print(tag, 'at sample', state['j'], 'of 151; seconds', round(time.monotonic() - t0, 1), flush=True)
if state['j'] == len(times):
    proj = np.array(state['proj'])  # [t, tangent, a, mode]
    np.savez_compressed(outdir / f'PROJ_{tag}.npz', proj=proj, srcs=state['srcs'], allmodes=allmodes,
                        shell=shell[allmodes], fams=np.array(state['fams']), s0=s0)
    os.remove(ck); print('FINISHED', tag)
