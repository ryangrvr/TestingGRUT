"""Low-noise retarded-response probes on the regenerated base48 field.
ENGINE: Codex's fieldstep tangent (exact derivative of the discrete Yoshida map; copied 5f5513b tree, run from $PKC/code).
PROBES, PROJECTIONS, ESTIMATORS: ours (independent of response.py).

Two probe families at a given phase s0 (each kicks one real momentum pi_1 or pi_2 with a real spatial pattern):
  'exact1'  : cos(k.x), sin(k.x) for each of the 9 +-pairs of shell 1 (k = 2pi/L |n|, |n|=1 or sqrt2 rounding to 1)
              -> 36 tangents; gives the full within-shell 2x2-per-mode response G_ab(k',k) for shell 1, hence the
              exact shell trace with no stochastic-estimator noise, and the off-diagonal (scattering) part.
  'hutch'   : M random-phase, unit-modulus probes over all modes of shells 1-4: src_k = e^{i phi_k}, |src_k| = 1;
              estimator R_ab(shell) = mean_{k in shell} q_a,k^(b) conj(src_k)  (Hutchinson trace; 2M tangents).
Delays 0..6 in steps of 0.04 (151), hmax = 0.004 as in response.py.  Chunked with restart files.
Usage: python3 -B probe_study.py <s0> <family> <M or 0> <outdir> <max_samples_this_call>
"""
import sys, os, time, pickle
sys.path.insert(0, '.')
from pt_common import *  # noqa (Codex module from the copy)

s0 = float(sys.argv[1]); fam = sys.argv[2]; M = int(sys.argv[3]); outdir = Path(sys.argv[4]); budget = int(sys.argv[5])
outdir.mkdir(parents=True, exist_ok=True)
N, L = 48, 12.
tag = f'{fam}_s{s0:g}' + (f'_M{M}' if fam == 'hutch' else '')
ck = outdir / f'ckpt_{tag}.pkl'
mm = np.fft.fftfreq(N) * N
nx, ny, nz = np.meshgrid(mm, mm, mm, indexing='ij')
kc = 2 * np.pi / L * np.sqrt(nx ** 2 + ny ** 2 + nz ** 2)
shell = np.floor(kc / (2 * np.pi / L) + .5).astype(int)
x = np.arange(N) * L / N
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')

def halfspace(ix):
    n = np.array([nx.ravel()[ix], ny.ravel()[ix], nz.ravel()[ix]])
    first = next((v for v in n if v != 0), 0)
    return first > 0

modes = {s: np.flatnonzero((shell.ravel() == s)) for s in range(1, 5)}
if not ck.exists():
    z = np.load(WORK / 'base48' / f'state_s{s0:g}.npy').copy()
    pats = []; meta = []
    if fam == 'exact1':
        for ix in modes[1]:
            if not halfspace(ix): continue
            kv = 2 * np.pi / L * np.array([nx.ravel()[ix], ny.ravel()[ix], nz.ravel()[ix]])
            ph = kv[0] * X + kv[1] * Y + kv[2] * Z
            pats += [np.cos(ph), np.sin(ph)]; meta += [(int(ix), 'cos'), (int(ix), 'sin')]
    else:
        rng = np.random.default_rng(1000 + M + int(s0 * 100))
        for j in range(M):
            coef = np.zeros((N, N, N), complex)
            for s in range(1, 5):
                for ix in modes[s]:
                    if not halfspace(ix): continue
                    i, jj, k = np.unravel_index(ix, (N, N, N)); phi = rng.uniform(0, 2 * np.pi)
                    coef[i, jj, k] = np.exp(1j * phi); coef[(-i) % N, (-jj) % N, (-k) % N] = np.exp(-1j * phi)
            pat = ifftn(coef, norm='forward'); assert np.max(abs(pat.imag)) < 1e-10
            pats.append(pat.real); meta.append((j, 'hutch'))
    nt = 2 * len(pats)
    tangent = np.zeros((nt, 6, N, N, N))
    for p, pat in enumerate(pats):
        tangent[2 * p, 3] = pat; tangent[2 * p + 1, 4] = pat
    srcs = np.array([fftn(p, norm='forward').ravel()[np.concatenate([modes[s] for s in range(1, 5)])] for p in pats])
    state = dict(z=z, tangent=tangent, j=0, prev=0., meta=meta, srcs=srcs, proj=[])
else:
    state = pickle.load(open(ck, 'rb'))
set_box(L); lib = library()
allmodes = np.concatenate([modes[s] for s in range(1, 5)])
times = np.arange(151) * .04; t0 = time.monotonic(); done = 0
z = state['z']; tangent = state['tangent']
while state['j'] < len(times) and done < budget:
    t = times[state['j']]
    if state['j']: step(lib, z, t - state['prev'], .004, tangent=tangent, threads=2)
    state['prev'] = t
    f = fftn(tangent[:, :2], axes=(-3, -2, -1), norm='forward').reshape(tangent.shape[0], 2, -1)[:, :, allmodes]
    state['proj'].append(f.astype(np.complex128))
    state['j'] += 1; done += 1
pickle.dump(state, open(ck, 'wb'), protocol=4)
print(tag, 'at sample', state['j'], 'of 151; seconds', round(time.monotonic() - t0, 1), flush=True)
if state['j'] == len(times):
    proj = np.array(state['proj'])  # [t, tangent, a, mode]
    np.savez_compressed(outdir / f'PROJ_{tag}.npz', proj=proj, srcs=state['srcs'], allmodes=allmodes,
                        shell=shell.ravel()[allmodes], meta=np.array([str(m) for m in state['meta']]), s0=s0)
    os.remove(ck); print('FINISHED', tag)
