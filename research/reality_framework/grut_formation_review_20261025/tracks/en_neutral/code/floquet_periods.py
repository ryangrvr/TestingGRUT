"""en_neutral item 3b: period-by-period Floquet exponents of delta-phi along the nonlinear N48 history. OWN code.

The first-pass finite-time Lyapunov estimate (mechanism.py, ln sigma_max / W) is dominated by the non-normal
intra-period amplitude swing (mu ~ 0.4-0.7 even for stable modes) and is NOT used for conclusions.  Here each
radial period of the coherent zero mode (pericentre to pericentre, minima of |chibar|^2 in the dense trace) defines
a one-period monodromy M_i(k) for phi'' = -(k^2 + m^2(s)) phi, and mu_i(k) = ln max|eig M_i| / T_i
(the Floquet exponent of a pump that repeated period i exactly).  Pumps:
  V1 free quartic orbit with the zero mode's own (E0, Q0) at the period start (exact periodic orbit; depletion test)
  V2 measured zero-mode pump  g2 |chibar(s)|^2
  V3 measured homogeneous Hartree mass  g2 <|chi_X|^2>(s) + 3 lamP <phi^2>(s)
  V4 V2 + constant g2 <S> + 3 lamP <phi^2> (period mean)  -> constant-shift detuning only
  V5 V1 orbit but in the Hartree potential (mean feels + 2 lamX S + g2<phi^2>) and phi mass + g2 S + 3 lamP<phi^2>
usage: python3 -I floquet_periods.py <rundir> <out_json>
"""
import sys, json, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import en_lib as E

rundir, out = Path(sys.argv[1]), Path(sys.argv[2])
t0 = time.time()
tr = np.genfromtxt(rundir / 'state.tsv', names=True, delimiter='\t')
s, tau = tr['eta'], tr['tau']
rb2 = tr['q0'] ** 2 + tr['q1'] ** 2
E0 = 0.5 * (tr['p0'] ** 2 + tr['p1'] ** 2) + E.LX / 4 * rb2 ** 2
Q0 = tr['q0'] * tr['p1'] - tr['q1'] * tr['p0']
S = tr['X2'] - rb2
kgrid = np.arange(0.02, 14.0001, 0.02); k2 = kgrid ** 2


def monodromy_series(sg, m2):
    nk = len(k2); T = np.zeros((nk, 2, 2)); T[:, 0, 0] = 1; T[:, 1, 1] = 1
    def f(m, T):
        d = np.empty_like(T); d[:, 0, :] = T[:, 1, :]; d[:, 1, :] = -(k2 + m)[:, None] * T[:, 0, :]; return d
    for i in range(len(sg) - 1):
        h = sg[i + 1] - sg[i]; a0, a1 = m2[i], m2[i + 1]; am = 0.5 * (a0 + a1)
        a = f(a0, T); b = f(am, T + h / 2 * a); c = f(am, T + h / 2 * b); d = f(a1, T + h * c)
        T = T + h / 6 * (a + 2 * b + 2 * c + d)
    tr_ = T[:, 0, 0] + T[:, 1, 1]
    lam = np.abs(tr_) / 2 + np.sqrt(np.maximum(tr_ * tr_ / 4 - 1, 0))
    return np.log(np.maximum(lam, 1.0)) / (sg[-1] - sg[0])


def orbit_period_r2(E_, Q_, mH2, n=3000):
    T_, rmin, rmax = E.radial_period(E_, Q_, extra_m2=mH2)
    h = T_ / n; x = np.array([rmin, 0.0]); v = np.array([0.0, Q_ / rmin]); r = np.empty(n + 1); r[0] = x @ x
    def acc(x): return -(E.LX * (x @ x) + mH2) * x
    for i in range(n):
        k1v = acc(x); k1x = v; k2v = acc(x + h / 2 * k1x); k2x = v + h / 2 * k1v
        k3v = acc(x + h / 2 * k2x); k3x = v + h / 2 * k2v; k4v = acc(x + h * k3x); k4x = v + h * k3v
        x = x + h / 6 * (k1x + 2 * k2x + 2 * k3x + k4x); v = v + h / 6 * (k1v + 2 * k2v + 2 * k3v + k4v); r[i + 1] = x @ x
    return np.linspace(0, T_, n + 1), r, T_, rmin, rmax


def summ(mu):
    i = int(np.argmax(mu))
    pos = kgrid[mu > 0.02]
    return dict(mu_max=float(mu[i]), k_at_max=float(kgrid[i]), mu_matched_band=float(mu[(kgrid >= 6.74) & (kgrid <= 7.27)].max()),
                k_lowest_band=float(pos.min()) if len(pos) else None,
                mu_first_band=float(mu[(kgrid >= pos.min()) & (kgrid <= pos.min() + 1.0)].max()) if len(pos) else 0.0)


# pericentres: local minima of |chibar|^2
mins = np.flatnonzero((rb2[1:-1] < rb2[:-2]) & (rb2[1:-1] <= rb2[2:])) + 1
rows = []
G = {v: np.zeros_like(kgrid) for v in ('V1', 'V3', 'V4', 'V5')}
G_from84 = {v: np.zeros_like(kgrid) for v in ('V1', 'V3', 'V4', 'V5')}
for a, b in zip(mins[:-1], mins[1:]):
    sl = slice(a, b + 1); sg = s[sl]
    if sg[-1] - sg[0] < 0.2: continue
    Sm = float(np.mean(S[sl])); P2 = float(np.mean(tr['phi2'][sl]))
    mu2 = monodromy_series(sg, E.G2 * rb2[sl])
    mu3 = monodromy_series(sg, E.G2 * tr['X2'][sl] + 3 * E.LP * tr['phi2'][sl])
    mu4 = monodromy_series(sg, E.G2 * rb2[sl] + E.G2 * Sm + 3 * E.LP * P2)
    sg1, r1, T1, rmin1, rmax1 = orbit_period_r2(E0[a], Q0[a], 0.0)
    mu1 = monodromy_series(sg1, E.G2 * r1)
    mH2 = 2 * E.LX * Sm + E.G2 * P2
    EH = E0[a] + 0.5 * mH2 * rb2[a]
    sg5, r5, T5, rmin5, rmax5 = orbit_period_r2(EH, Q0[a], mH2)
    mu5 = monodromy_series(sg5, E.G2 * (r5 + Sm) + 3 * E.LP * P2)
    row = dict(tau_start=float(tau[a]), tau_end=float(tau[b]), T_measured=float(sg[-1] - sg[0]), T_free=float(T1), T_hartree=float(T5),
               E0_over_rho0=float(E0[a] / tr['rho'][0]), rmax_over_rmin_free=float(rmax1 / rmin1),
               rbar2_max_over_min_measured=float(rb2[sl].max() / rb2[sl].min()),
               X2_max_over_min_measured=float(tr['X2'][sl].max() / tr['X2'][sl].min()),
               S_mean=Sm, g2S=E.G2 * Sm, rbar2_mean=float(np.mean(rb2[sl])), g2_rbar2_min=float(E.G2 * rb2[sl].min()),
               V1_free=summ(mu1), V2_zero_mode=summ(mu2), V3_hartree_measured=summ(mu3), V4_const_shift=summ(mu4), V5_hartree_orbit=summ(mu5))
    rows.append(row)
    for v, mu, T_ in (('V1', mu1, T1), ('V3', mu3, sg[-1] - sg[0]), ('V4', mu4, sg[-1] - sg[0]), ('V5', mu5, T5)):
        G[v] += mu * (sg[-1] - sg[0])  # amplitude e-folds accumulated over the measured period length
        if tau[a] >= 84: G_from84[v] += mu * (sg[-1] - sg[0])
    print('%6.1f-%6.1f T %.3f/%.3f  E0 %.3f ecc %.2f  g2S %5.1f | mu1 %.3f@%.2f  mu2 %.3f@%.2f  mu3 %.3f@%.2f  mu4 %.3f@%.2f  mu5 %.3f@%.2f | band(6.74-7.27) %.3f %.3f %.3f'
          % (row['tau_start'], row['tau_end'], row['T_measured'], T1, row['E0_over_rho0'], row['rmax_over_rmin_free'], row['g2S'],
             row['V1_free']['mu_max'], row['V1_free']['k_at_max'], row['V2_zero_mode']['mu_max'], row['V2_zero_mode']['k_at_max'],
             row['V3_hartree_measured']['mu_max'], row['V3_hartree_measured']['k_at_max'], row['V4_const_shift']['mu_max'], row['V4_const_shift']['k_at_max'],
             row['V5_hartree_orbit']['mu_max'], row['V5_hartree_orbit']['k_at_max'],
             row['V2_zero_mode']['mu_matched_band'], row['V3_hartree_measured']['mu_matched_band'], row['V4_const_shift']['mu_matched_band']), flush=True)
cum = {v: dict(max_amplitude_efolds=float(G[v].max()), k_at_max=float(kgrid[int(np.argmax(G[v]))]), matched_band_efolds=float(G[v][(kgrid >= 6.74) & (kgrid <= 7.27)].max()),
             from_tau84_max=float(G_from84[v].max()), from_tau84_k=float(kgrid[int(np.argmax(G_from84[v]))])) for v in G}
print(json.dumps(cum, indent=0))
out.write_text(json.dumps(dict(label='OWN period-by-period Floquet analysis of the reproduction N48 history', kgrid=[float(kgrid[0]), float(kgrid[-1]), 0.02],
                               cumulative_amplitude_efolds_tau26_762=cum, periods=rows, seconds=time.time() - t0), indent=1) + '\n')
