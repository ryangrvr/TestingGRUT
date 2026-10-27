"""Checks of MATHEMATICS.md items no other track covered (own code).

1. Fixed-charge energy bound E >= (3/4) lambda^(1/3) |n_Q|^(4/3) V: pointwise |pi_X|^2 >= j^2/|chi_X|^2, FD4 symbol >= 0,
   phi terms >= 0, Jensen on |j|^(4/3)  ->  a bound over ALL lattice states at fixed charge, not only homogeneous ones.
   Numerical: published curves (E0 vs the bound at Q0) and a direct test on the tau-800 N48 snapshot.
2. Frozen-Sigma static quadratic-Y estimate 3 kappa^2 |chi|^2 / (4 lambda_X (M_Y^2 - sqrt2 mu h) a^2): solve the static
   quadratic Y equation exactly for random X phases / Sigma phases and compare with the bound.
3. Gauge omission (e=0 in matching and evolution): Coulomb energy of the measured charge fluctuations for e_X=1e-7,1e-6.

usage: GRUT_REHEATING_WORKDIR=... python3 -B -I math_checks.py <reality_framework_dir> <packet_dir> <N48_trace> <N48_tau800_snapshot|-> <out>
"""
import sys, json
from pathlib import Path
import numpy as np

RF, PKT, TRACE, SNAP, OUT = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4], Path(sys.argv[5])
sys.path.insert(0, str(RF / 'grut_reheating_preparation_20261026' / 'code'))
import reheating as R  # noqa: E402

LX, G2, LP = .001, .01, .03
res = {}

# 1. energy bound on published curves
c = np.load(PKT / 'FIELD_CURVES.npz')
tags = sorted(set(k[:-4] for k in c.files if k.endswith('_tau')))
rows = {}
for t in tags:
    E0 = 1 - c[t + '_energy_outside_spatial_mean_fraction']; Em = c[t + '_minimum_energy_at_total_charge_fraction']
    QQ = c[t + '_spatial_mean_charge_fraction']
    rows[t] = dict(Emin_over_E=float(Em[0]), Emin_drift=float(np.ptp(Em)),
                   min_margin_E0_minus_bound_at_Q0=float(np.min(E0 - Em * abs(QQ) ** (4 / 3))))
res['energy_bound_published_curves'] = rows
if SNAP != '-':
    with open(SNAP, 'rb') as fh:
        N = int(np.frombuffer(fh.read(4), '<i4')[0]); L = float(np.frombuffer(fh.read(8), '<f8')[0]); fh.read(12)
        a = np.fromfile(fh, '<f8').reshape(6, N, N, N)
    q, p = a[:3], a[3:]; r2 = q[0] ** 2 + q[1] ** 2; j = q[0] * p[1] - q[1] * p[0]
    kin = .5 * (p[0] ** 2 + p[1] ** 2)
    res['energy_bound_snapshot_tau800'] = dict(
        min_pointwise_kinetic_minus_j2_over_2r2=float(np.min(kin - j * j / (2 * r2))),
        local_bound_over_global_bound=float(np.mean(.75 * LX ** (1 / 3) * abs(j) ** (4 / 3)) / (.75 * LX ** (1 / 3) * abs(j.mean()) ** (4 / 3))),
        note='>=0 and >=1 respectively, as the pointwise/Jensen steps require')

# 2. frozen-Sigma static quadratic-Y estimate
rng = np.random.default_rng(7); M2 = R.M ** 2; h = R.HH; mu = R.MU; kap = R.KAP
ratios = []
for _ in range(2000):
    x = np.exp(1j * rng.uniform(0, 2 * np.pi)) * rng.uniform(1, 100)
    s = h / np.sqrt(2) * np.exp(1j * rng.uniform(0, 2 * np.pi))
    # static quadratic Y: M2 y + 2 mu s conj(y) = -kap x^3  (real 2x2 system in (Re y, Im y))
    b = -kap * x ** 3
    A = np.array([[M2 + 2 * mu * s.real, 2 * mu * s.imag], [2 * mu * s.imag, M2 - 2 * mu * s.real]])
    yr, yi = np.linalg.solve(A, [b.real, b.imag]); y = yr + 1j * yi
    dF = abs(3 * kap * y * np.conj(x) ** 2); F = abs(2 * LX * abs(x) ** 2 * x)
    bound = 3 * kap ** 2 * abs(x) ** 2 / (2 * LX * (M2 - np.sqrt(2) * mu * h))
    ratios.append(dF / F / bound)
res['frozen_Y_estimate'] = dict(kappa=kap, mu=mu, M_Y=R.M, h=h, smallest_Y_mass2=M2 - np.sqrt(2) * mu * h,
                                max_exact_over_bound=float(max(ratios)), min_exact_over_bound=float(min(ratios)),
                                conformal_form='|x|^2=|chi|^2/(2a^2) turns 3k^2|x|^2/(2 lX m^2) into 3k^2|chi|^2/(4 lX m^2 a^2)')

# 3. gauge omission
tr = np.genfromtxt(TRACE, names=True, delimiter='\t')
var = np.maximum(tr['charge_RMS'] ** 2 - tr['nQ'] ** 2, 0); kmin = 2 * np.pi / 12
res['gauge_estimate'] = {f'e={e:g}': dict(coulomb_upper_over_rho_max=float(np.max(.5 * e * e * var / kmin ** 2 / tr['rho'])),
                                          IR_scale_k_where_e2chi2_over_k2_is_1=float(e * np.sqrt(tr['X2'].max())))
                         for e in [1e-7, 1e-6]}
OUT.write_text(json.dumps(res, indent=2) + '\n'); print(json.dumps(res, indent=1))
