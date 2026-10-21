"""Own re-implementation of the Gaussian preparation recipe stated in Codex's wrapper record:
'NumPy PCG64 standard_normal; rfftn, exp(-k^2/2), zero k=0, unit realization rms' (three fields drawn
sequentially from default_rng(seed), each an (N,N,N) standard-normal array).

Modes:
  native   N dx seed out           -> the recipe exactly at (N,dx)  [same as Codex's wrapper would produce]
  restrict N0 dx0 N1 seed out      -> draw at (N0,dx0), Fourier-restrict to N1 points on the SAME box
                                      L=N0*dx0 (dx1=L/N1), then re-standardize. Because exp(-k^2/2) is
                                      < 1e-5 at the coarse Nyquist for dx1=0.6, this is the same continuum
                                      preparation sampled on a coarser grid [P].
Prints sha256, means, rms, covariance of the written little-endian f64 file.
"""
import hashlib
import json
import sys

import numpy as np


def filt_field(white, dx):
    N = white.shape[0]
    k = 2 * np.pi * np.fft.fftfreq(N, dx)
    kr = 2 * np.pi * np.fft.rfftfreq(N, dx)
    k2 = k[:, None, None] ** 2 + k[None, :, None] ** 2 + kr[None, None, :] ** 2
    f = np.fft.rfftn(white)
    f *= np.exp(-0.5 * k2)
    f.flat[0] = 0
    return f


def standardize(q):
    q = q - q.mean()
    return q / q.std()


def restrict(fhat, N0, N1):
    """Keep |integer wavenumber| < N1/2 (drop the coarse Nyquist plane), rescale for irfftn norm."""
    h = N1 // 2
    idx = np.r_[0:h, N0 - h + 1:N0]          # -h+1..h-1 along full axes
    tidx = np.r_[0:h, N1 - h + 1:N1]
    out = np.zeros((N1, N1, h + 1), dtype=complex)
    sub = fhat[np.ix_(idx, idx, np.arange(h))]
    out[np.ix_(tidx, tidx, np.arange(h))] = sub
    return np.fft.irfftn(out, s=(N1, N1, N1), axes=(0, 1, 2)) * (N1 ** 3 / N0 ** 3)


def main():
    mode = sys.argv[1]
    if mode == 'native':
        N, dx, seed, out = int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
        rng = np.random.default_rng(seed)
        fields = []
        for _ in range(3):
            f = filt_field(rng.standard_normal((N, N, N)), dx)
            fields.append(standardize(np.fft.irfftn(f, s=(N, N, N), axes=(0, 1, 2))))
        info = {'mode': 'native', 'N': N, 'dx': dx, 'seed': seed}
    elif mode == 'restrict':
        N0, dx0, N1, seed, out = int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
        rng = np.random.default_rng(seed)
        fields = []
        maxdiff = []
        for _ in range(3):
            f = filt_field(rng.standard_normal((N0, N0, N0)), dx0)
            fine = np.fft.irfftn(f, s=(N0, N0, N0), axes=(0, 1, 2))
            mu, sd = fine.mean(), fine.std()
            coarse = restrict(f, N0, N1)
            # pointwise check against the fine field on coincident points (only when N0 % N1 == 0
            # would points coincide; otherwise check sd ratio only)
            fields.append((coarse - coarse.mean()) / sd)  # same normalization as the fine draw
            maxdiff.append({'fine_sd': float(sd), 'coarse_sd_over_fine_sd': float(coarse.std() / sd)})
        info = {'mode': 'restrict', 'N0': N0, 'dx0': dx0, 'N1': N1, 'dx1': N0 * dx0 / N1, 'seed': seed,
                'checks': maxdiff}
    else:
        raise SystemExit('mode?')
    with open(out, 'wb') as fo:
        for q in fields:
            q.astype('<f8', copy=False).tofile(fo)
    arr = np.fromfile(out, dtype='<f8').reshape(3, -1)
    info['means'] = [float(x) for x in arr.mean(axis=1)]
    info['rms'] = [float(np.sqrt(np.mean(x * x))) for x in arr]
    info['covariance'] = [[float(np.mean(arr[i] * arr[j])) for j in range(3)] for i in range(3)]
    info['sha256'] = hashlib.sha256(open(out, 'rb').read()).hexdigest()
    info['numpy'] = np.__version__
    print(json.dumps(info))


if __name__ == '__main__':
    main()
