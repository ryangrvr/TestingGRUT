"""en_neutral_skeptic item 3: own reading of the reproduction-run snapshots (32-bit, every 2 tau) to check
en_neutral's N3 numbers (matched-band phi energy, share above k_lat = 8) and to measure the classical occupation of
the unseeded UV modes against the half-quantum (vacuum) level.
Quadratic phi energy per mode e_k = 1/2(|pi_k|^2 + w_k^2 |phi_k|^2), w_k^2 = k_lat^2 + g2 <|chi_X|^2> (+3 lamP <phi^2>),
FFT norm 1/N^3 so sum_k e_k is an energy density. Occupation-type ratio e_k V / w_k: 1/2 for a vacuum-like mode.
OWN code (sk_lib FD4 symbol, own reader).
usage: python3 -I snap_check.py <rundir> <out_json>
"""
import sys, json, struct
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sk_lib as K

rd, out = Path(sys.argv[1]), Path(sys.argv[2])
tr = K.read_trace(rd / 'state.tsv'); rho0 = float(tr['rho'][0])
a0 = float(tr['a'][0]); C = float((tr['a'][-1] - tr['a'][0]) / tr['eta'][-1])
snaps = sorted(rd.glob('state_s*.bin'))
N, L = 48, 12.0; V = L ** 3
K2 = K.lattice_k2_grid(N, L); kl = np.sqrt(K2)
band = (kl >= 6.74) & (kl <= 7.27)
rows = []
for path in snaps:
    with open(path, 'rb') as f:
        n_, L_, s_, bits = struct.unpack('<iddi', f.read(24))
        a = np.fromfile(f, dtype='<f4' if bits == 32 else '<f8').reshape(6, n_, n_, n_).astype(np.float64)
    tau = 20 + a0 * s_ + 0.5 * C * s_ * s_
    j = int(np.argmin(abs(tr['eta'] - s_)))
    X2 = float(np.mean(a[0] ** 2 + a[1] ** 2)); P2 = float(np.mean(a[2] ** 2))
    ph = np.fft.fftn(a[2]) / N ** 3; pp = np.fft.fftn(a[5]) / N ** 3
    w2 = K2 + K.G2 * X2 + 3 * K.LP * P2
    e = 0.5 * (abs(pp) ** 2 + w2 * abs(ph) ** 2); e[0, 0, 0] = 0.0
    w = np.sqrt(w2); occ = e * V / w
    r = dict(tau=tau, X2_snapshot=X2, X2_trace=float(tr['X2'][j]),
             band_energy=float(e[band].sum() / rho0), total_quadratic=float(e.sum() / rho0),
             frac_k_gt_8=float(e[kl > 8].sum() / e.sum()), frac_k_lt_1p8=float(e[(kl < 1.8) & (kl > 0)].sum() / e.sum()),
             mean_occ_k_gt_8=float(occ[kl > 8].mean()), mean_occ_seeded_1p8_8=float(occ[(kl >= 1.8) & (kl <= 8)].mean()),
             mean_occ_band=float(occ[band].mean()), median_occ_k_gt_8=float(np.median(occ[kl > 8])),
             frac_modes_k_gt_8_below_half=float((occ[kl > 8] < 0.5).mean()),
             vacuum_energy_k_gt_8_over_rho0=float(0.5 * w[kl > 8].sum() / V / rho0),
             classical_energy_k_gt_8_over_rho0=float(e[kl > 8].sum() / rho0))
    rows.append(r)
tau = np.array([r['tau'] for r in rows]); be = np.array([r['band_energy'] for r in rows])
rm = np.convolve(be, np.ones(5) / 5, mode='same')
pk = int(np.argmax(rm[2:-2])) + 2
def at(x, key): return rows[int(np.argmin(abs(tau - x)))][key]
res = dict(label='OWN analysis of en_neutral reproduction snapshots (Codex engine N48, 32-bit, every 2 tau)', rho0=rho0,
           band_energy_runmean10tau={str(x): float(rm[int(np.argmin(abs(tau - x)))]) for x in (24, 40, 60, 66, 70, 76, 80, 90, 100, 200, 400, 600, 796)},
           band_energy_peak=dict(tau=float(tau[pk]), value=float(rm[pk])),
           late={str(x): {k: at(x, k) for k in ('frac_k_gt_8', 'mean_occ_k_gt_8', 'median_occ_k_gt_8', 'frac_modes_k_gt_8_below_half', 'mean_occ_seeded_1p8_8',
                                                 'mean_occ_band', 'vacuum_energy_k_gt_8_over_rho0', 'classical_energy_k_gt_8_over_rho0', 'total_quadratic')}
                 for x in (20, 140, 400, 600, 790, 800)},
           frac_k_gt_8_last5_mean=float(np.mean([r['frac_k_gt_8'] for r in rows[-5:]])),
           mean_occ_k_gt_8_last5_mean=float(np.mean([r['mean_occ_k_gt_8'] for r in rows[-5:]])))
out.write_text(json.dumps(res, indent=1) + '\n')
print(json.dumps(res, indent=1))
