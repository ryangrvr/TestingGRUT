"""X-seed-amplitude control on Codex's unmodified early3d engine (reproduction-type engine, new initial data).

Builds the published N48/L12/seed1/cutoff8 tau-20 field exactly as Codex's initial.py does (re-implemented here so
that no Codex module with sibling-packet side effects is imported), checks the f=1 file against the published
SHA-256, then multiplies ONLY the nonzero Fourier modes of (chi1,chi2,pi1,pi2) by f. The X zero mode and the whole
phi seed are unchanged. Runs Codex's engine with Codex's run.py arguments except threads=2 and snapshot spacing.

usage: python3 -I xseed_scan.py <packet_dir> <engine> <workdir> <f> <end_tau> [N]
"""
import sys, json, struct, hashlib, subprocess, time
from pathlib import Path
import numpy as np

LX, G2, LP = .001, .01, .03
PUBLISHED_N48_INITIAL_SHA = 'abe46838e082d40215b00f88c5280ebc932b484c5cffbb55c56486cacecacb3d'


def coefficients(pkt, N=48, L=12., seed=1, cutoff=8.):
    data = np.load(Path(pkt) / f'MATCH_L{L:g}.npz')
    keys = data['keys']; lookup = {int(v): i for i, v in enumerate(keys)}
    coeff = np.zeros((6, N, N, N), complex); nmax = int(np.ceil(cutoff * L / (2 * np.pi)))
    for i in range(-nmax, nmax + 1):
        for j in range(-nmax, nmax + 1):
            for k in range(-nmax, nmax + 1):
                n = (i, j, k); first = next((x for x in n if x), 0)
                if first <= 0: continue
                nn = i * i + j * j + k * k; mom = 2 * np.pi / L * np.sqrt(nn)
                if not 1.8 <= mom <= cutoff: continue
                ix = lookup[nn]
                rng = np.random.default_rng(np.random.SeedSequence([seed, int(L * 100), i + 128, j + 128, k + 128]))
                z = (rng.normal(size=14) + 1j * rng.normal(size=14)) / np.sqrt(2)
                full = data['evolved_scalar'][ix] @ z[:12]; ph = data['evolved_phi'][ix] @ z[12:]
                light = np.r_[full[:2], ph[0], full[6:8], ph[1]] / np.sqrt(L ** 3)
                coeff[(slice(None), i % N, j % N, k % N)] = light
                coeff[(slice(None), (-i) % N, (-j) % N, (-k) % N)] = light.conjugate()
    coeff[:, 0, 0, 0] = np.r_[data['light_mean_q'], data['light_mean_p']]
    clock = dict(tau0=20., a0=float(data['a_match']), C=float(data['C_match']))
    return coeff, clock


def fields(coeff, f):
    c = coeff.copy()
    if f != 1.:
        zero = c[:, 0, 0, 0].copy()
        c[[0, 1, 3, 4]] *= f
        c[:, 0, 0, 0] = zero
    x = np.fft.ifftn(c, axes=(-3, -2, -1), norm='forward')
    assert np.max(abs(x.imag)) < 1e-10
    return x.real[:3].copy(), x.real[3:].copy()


def write_initial(path, q, p, L):
    with open(path, 'wb') as fh:
        fh.write(b'GRUT3D1\n'); fh.write(struct.pack('<id4d', q.shape[-1], L, 0., G2, LX, LP))
        np.asarray(q, dtype='<f8').tofile(fh); np.asarray(p, dtype='<f8').tofile(fh)


def main():
    pkt, exe, work, f, end_tau = sys.argv[1], sys.argv[2], Path(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
    N = int(sys.argv[6]) if len(sys.argv) > 6 else 48; L = 12.
    work.mkdir(parents=True, exist_ok=True)
    coeff, clock = coefficients(pkt, N=N, L=L)
    check = None
    if N == 48:
        q1, p1 = fields(coeff, 1.); write_initial(work / 'initial_f1.bin', q1, p1, L)
        sha = hashlib.sha256((work / 'initial_f1.bin').read_bytes()).hexdigest()
        check = dict(f1_initial_sha256=sha, matches_published=sha == PUBLISHED_N48_INITIAL_SHA)
        (work / 'initial_f1.bin').unlink(); del q1, p1
        assert check['matches_published'], check
    q, p = fields(coeff, f)
    Xfl = float(np.mean((q[0] - q[0].mean()) ** 2 + (q[1] - q[1].mean()) ** 2))
    init = work / 'initial.bin'; write_initial(init, q, p, L); del q, p
    tau0, a0, C = clock['tau0'], clock['a0'], clock['C']
    duration = 2 * (end_tau - tau0) / (a0 + np.sqrt(a0 * a0 + 2 * C * (end_tau - tau0)))
    prefix = work / 'state'
    args = [exe, str(init), str(prefix), '0.004', str(duration), '0.05', '100', '4', '2', '64', str(tau0), str(a0), str(C)]
    t0 = time.monotonic()
    with open(work / 'run.log', 'w') as fh: subprocess.run(args, check=True, stdout=fh, stderr=subprocess.STDOUT)
    info = dict(f=f, N=N, L=L, end_tau=end_tau, clock=clock, initial_X_fluctuation_variance=Xfl, check=check,
                seconds=time.monotonic() - t0, args=args[3:])
    (work / 'info.json').write_text(json.dumps(info, indent=2) + '\n')
    print(json.dumps(info))


if __name__ == '__main__':
    main()
