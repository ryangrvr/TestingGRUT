"""Dense Gaussian preparations: independent check of the 'formation imbalance = linear-field statistical imbalance' claim.

Usage: python3 -I dense_check.py <codex_raw_dir> <out_json>
1. Regenerate Codex's noise (gaussian_run.py noise(), transcribed) and compare covariances with Codex's json [P].
2. Evolve S linearly per lattice khat^2 with the solver's KDK leapfrog (dt=0.06), chi'(0)=chi(0)/40 [E].
3. At total rms 0.5 and at the recorded a=1.57/1.60 samples: <sign S>, raw m_05 = (f+ - f-)/(f+ + f-) with |S|>0.5.
4. Ensemble sd of <sign S> by the arcsine law: Var = mean_r (2/pi) arcsin(rho(r)), rho from the ensemble spectrum
   exp(-k^2) g(khat)^2 (independent of the realization) [E].
"""
import json, sys
from pathlib import Path
import numpy as np

raw, outp = Path(sys.argv[1]), Path(sys.argv[2])


def codex_noise(N, dx, seed):
    rng = np.random.default_rng(seed)
    k = 2 * np.pi * np.fft.fftfreq(N, dx); kr = 2 * np.pi * np.fft.rfftfreq(N, dx)
    k2 = k[:, None, None] ** 2 + k[None, :, None] ** 2 + kr[None, None, :] ** 2
    filt = np.exp(-.5 * k2); out = []
    for field in range(3):
        f = np.fft.rfftn(rng.standard_normal((N, N, N))); f *= filt; f.flat[0] = 0
        q = np.fft.irfftn(f, s=(N, N, N), axes=(0, 1, 2)); q -= q.mean(); q /= q.std(); out.append(q)
    return out


def growth_table(kmax2, kap=2.0, dt=0.06, T=30.0, nk=3001):
    k2 = np.linspace(0, kmax2, nk); F2 = kap ** (2 / 3)
    x = np.ones(nk); p = np.full(nk, 1 / 40); steps = int(round(T / dt)); G = [x.copy()]
    for it in range(1, steps + 1):
        a0, a1 = 1 + (it - 1) * dt / 40, 1 + it * dt / 40
        p += 0.5 * dt * (a0 * a0 - F2 - k2) * x; x = x + dt * p; p += 0.5 * dt * (a1 * a1 - F2 - k2) * x; G.append(x.copy())
    return k2, np.array(G)


res = {}
for name, N, seed in (('gaussian_s101', 192, 101), ('gaussian_s202', 192, 202), ('gaussian_large_s101', 288, 101)):
    dx = 0.4; L = N * dx
    u, w, S = codex_noise(N, dx, seed)
    cov = [[float(np.mean(a * b)) for b in (u, w, S)] for a in (u, w, S)]
    ref = json.loads((raw / f'gaussian_N{N}_s{seed}.noise.json').read_text())['covariance']
    covdev = float(np.max(np.abs(np.array(cov) - np.array(ref))))
    del u, w
    Sk = np.fft.rfftn(S); del S
    kk = np.fft.fftfreq(N) ; kr = np.fft.rfftfreq(N)
    kh2 = (4 * np.sin(np.pi * kk[:, None, None]) ** 2 + 4 * np.sin(np.pi * kk[None, :, None]) ** 2 + 4 * np.sin(np.pi * kr[None, None, :]) ** 2) / dx ** 2
    kc2 = ((2 * np.pi / dx) ** 2) * (kk[:, None, None] ** 2 + kk[None, :, None] ** 2 + kr[None, None, :] ** 2)
    ktab, G = growth_table(float(kh2.max()) + 1e-9)
    rec = np.genfromtxt(raw / f'{name}.csv', delimiter=',', names=True)
    out = {'N': N, 'L': L, 'cov_dev_vs_codex': covdev, 'C_uS': cov[0][2], 'C_wS': cov[1][2]}
    # rms of linear field vs time (Parseval) to locate rms 0.5
    w8 = np.full(Sk.shape, 2.0); w8[..., 0] = 1.0
    if N % 2 == 0: w8[..., -1] = 1.0
    P = (np.abs(Sk) ** 2) * w8 / N ** 6
    sel = np.argsort(kh2.ravel()); kh2s = kh2.ravel()[sel]; Ps = P.ravel()[sel]
    times = np.arange(G.shape[0]) * 0.06
    marks = {}
    for i in range(0, G.shape[0], 5):
        a = 1 + times[i] / 40
        gi = np.interp(kh2s, ktab, G[i])
        rms = 0.01 * np.sqrt(np.sum(Ps * gi ** 2)) / a
        if rms >= 0.5 and 'rms0.5' not in marks: marks['rms0.5'] = i
    for a_ in (1.57, 1.60):
        marks[f'a{a_}'] = int(round((a_ - 1) * 40 / 0.06))
    for lab, i in marks.items():
        a = 1 + times[i] / 40
        gk = np.interp(kh2, ktab, G[i])
        f = 0.01 * np.fft.irfftn(Sk * gk, s=(N, N, N), axes=(0, 1, 2)) / a
        fp, fm = float(np.mean(f > 0.5)), float(np.mean(f < -0.5))
        # ensemble sd via arcsine law with the ensemble spectrum exp(-kc^2) g^2 on the full rfft grid
        Pe = np.exp(-kc2) * gk ** 2; Pe.flat[0] = 0
        corr = np.fft.irfftn(Pe, s=(N, N, N), axes=(0, 1, 2)); corr /= corr.flat[0]
        sd_sign = float(np.sqrt(np.mean((2 / np.pi) * np.arcsin(np.clip(corr, -1, 1)))))
        j = int(np.argmin(np.abs(rec['a'] - a)))
        ob = (rec['fplus_05'][j] - rec['fminus_05'][j]) / max(rec['fplus_05'][j] + rec['fminus_05'][j], 1e-300)
        out[lab] = {'a': float(a), 'rms_lin': float(f.std()), 'mean_sign_lin': float(np.mean(np.sign(f))),
                    'm05_lin': (fp - fm) / max(fp + fm, 1e-300), 'f05_lin_total': fp + fm, 'ensemble_sd_mean_sign': sd_sign,
                    'recorded_a': float(rec['a'][j]), 'recorded_m05': float(ob), 'recorded_f05_total': float(rec['fplus_05'][j] + rec['fminus_05'][j]),
                    'recorded_mean_S': float(rec['mean_S'][j])}
        del f, corr
    res[name] = out
    print(name, json.dumps(out), flush=True)
outp.write_text(json.dumps(res, indent=1))
