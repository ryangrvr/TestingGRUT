"""en_neutral_skeptic item 1: 2x2 factorial separation of "pump depletion" from "fluctuation-mass detuning" for the
same-seed linear delta-phi, plus re-runs of en_neutral's two key controls (measured zero-mode pump L2, measured
Hartree mass L3) with own code and a k-resolved breakdown of where their tau-800 energy sits.

OWN code (sk_lib.py). Inputs: the en_neutral reproduction run of Codex's unmodified engine (N48, L12, seed 1):
state.tsv (engine trace, Delta tau = 0.05) and initial.bin (tau-20 field).

Pumps m2(s) for phi_k'' = -(k_lat^2 + m2(s)) phi_k :
  A   fixed free orbit from the tau-20 mean (Codex's LINEAR_PHI_CONTROL; must give 2.0917)
  B   A + g2 Sbar(s) + 3 lamP <phi^2>bar(s)       (undepleted pump, measured fluctuation mass: DETUNING ONLY)
  C   depleted free orbit: consecutive exact free quartic periods whose (E, Q) follow the measured zero mode's
      per-period mean (E0, Q0); each period at its own free frequency                    (DEPLETION ONLY)
  Cp  as C but (E0, Q0) taken at the measured pericentres (en_neutral's V1 convention)
  D   C + g2 Sbar + 3 lamP <phi^2>bar                                                   (BOTH)
  E   measured Hartree g2 <|chi_X|^2>(s) + 3 lamP <phi^2>(s)              (= en_neutral L3, own re-run)
  F   measured zero mode g2 |chibar(s)|^2                                   (= en_neutral L2, own re-run)
  G   measured zero mode + g2 Sbar + 3 lamP <phi^2>bar                      (= en_neutral L4 analogue)
Sbar = per-measured-period (pericentre to pericentre) mean of S = <|chi_X|^2> - |chibar|^2, linearly interpolated
between period mid-points. Energy convention as Codex: raw, interaction included with the |chi|^2 the pump uses.
usage: python3 -I factorial.py <rundir> <out_json>
"""
import sys, json, time
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sk_lib as K

rundir, out = Path(sys.argv[1]), Path(sys.argv[2])
t0 = time.time()
tr = K.read_trace(rundir / 'state.tsv')
seed = K.phi_seed_by_k2(rundir / 'initial.bin')
rho0 = float(tr['rho'][0])
s_tr, tau_tr = tr['eta'], tr['tau']
a0 = float(tr['a'][0]); Cc = float((tr['a'][-1] - tr['a'][0]) / tr['eta'][-1])
rho0_own = K.energy_density_realspace(seed['q'], seed['p'], seed['L'])

# fine grid: trace grid with midpoints inserted; RK4 stage points left/mid/right
s_f = np.empty(2 * len(s_tr) - 1); s_f[0::2] = s_tr; s_f[1::2] = 0.5 * (s_tr[1:] + s_tr[:-1])
sL, sR = s_f[:-1], s_f[1:]; sM = 0.5 * (sL + sR)
rec_tr = np.arange(0, len(s_tr), 10); rec_f = 2 * rec_tr
tau_rec = tau_tr[rec_tr]

# measured series
rb2 = tr['q0'] ** 2 + tr['q1'] ** 2
X2, P2 = tr['X2'], tr['phi2']
S = X2 - rb2
E0m = 0.5 * (tr['p0'] ** 2 + tr['p1'] ** 2) + K.LX / 4 * rb2 ** 2
Q0m = tr['q0'] * tr['p1'] - tr['q1'] * tr['p0']
peri = np.flatnonzero((rb2[1:-1] < rb2[:-2]) & (rb2[1:-1] <= rb2[2:])) + 1
spl = {nm: CubicSpline(s_tr, v) for nm, v in (('rb2', rb2), ('X2', X2), ('P2', P2))}

def period_means(v):
    mids, vals = [], []
    edges = np.r_[0, peri]
    for a, b in zip(edges[:-1], edges[1:]):
        sl = slice(a, b + 1)
        mids.append(0.5 * (s_tr[a] + s_tr[b])); vals.append(np.trapezoid(v[sl], s_tr[sl]) / (s_tr[b] - s_tr[a]))
    # last partial period
    a = peri[-1]; mids.append(0.5 * (s_tr[a] + s_tr[-1])); vals.append(np.trapezoid(v[a:], s_tr[a:]) / (s_tr[-1] - s_tr[a]))
    return np.array(mids), np.array(vals)

mS, vS = period_means(S); mP, vP = period_means(P2)
mE, vE = period_means(E0m); mQ, vQ = period_means(Q0m)
Sbar = lambda s: np.interp(s, mS, vS)
Pbar = lambda s: np.interp(s, mP, vP)
shift = lambda s: K.G2 * Sbar(s) + 3 * K.LP * Pbar(s)

# ---- A: fixed free orbit ----
x0 = np.array([tr['q0'][0], tr['q1'][0]]); v0 = np.array([tr['p0'][0], tr['p1'][0]])
solA = K.free_orbit_dense(x0, v0, s_tr[-1] + 1e-9)
def rA(s):
    y = solA.sol(s); return y[0] ** 2 + y[1] ** 2

# ---- C / Cp: depleted free orbit, concatenated exact free periods ----
def build_depleted(Et, Qt):
    # first segment: the original orbit up to its first pericentre
    yA = solA.sol(np.linspace(0, 1.2, 24001)); rv = yA[0] * yA[2] + yA[1] * yA[3]
    j = np.flatnonzero((rv[:-1] < 0) & (rv[1:] >= 0))[0]
    from scipy.optimize import brentq
    sp0 = brentq(lambda s: (lambda y: y[0] * y[2] + y[1] * y[3])(solA.sol(s)), 1.2 * j / 24000, 1.2 * (j + 1) / 24000, xtol=1e-14)
    segs = [(0.0, sp0, lambda s: rA(s), None)]
    s0 = sp0; info = []
    while s0 < s_tr[-1]:
        E_, Q_ = float(Et(s0)), float(Qt(s0))
        sol, T, rmin, rmax = K.one_free_period(E_, Q_)
        segs.append((s0, s0 + T, (lambda sol_, st: (lambda s: (lambda y: y[0] ** 2 + y[1] ** 2)(sol_.sol(s - st))))(sol, s0), None))
        info.append(dict(s_start=s0, T=T, E=E_, Q=Q_, ecc=rmax / rmin))
        s0 += T
    starts = np.array([sg[0] for sg in segs])
    def r(s):
        s = np.atleast_1d(s); out_ = np.empty_like(s)
        idx = np.clip(np.searchsorted(starts, s, side='right') - 1, 0, len(segs) - 1)
        for i in np.unique(idx):
            m = idx == i; out_[m] = segs[i][2](s[m])
        return out_
    return r, info

rC, infoC = build_depleted(lambda s: np.interp(s, mE, vE), lambda s: np.interp(s, mQ, vQ))
pE, pQ = s_tr[peri], None
rCp, infoCp = build_depleted(lambda s: np.interp(s, s_tr[np.r_[0, peri]], E0m[np.r_[0, peri]]),
                             lambda s: np.interp(s, s_tr[np.r_[0, peri]], Q0m[np.r_[0, peri]]))

pumps = {
    'A_fixed_orbit': (lambda s: K.G2 * rA(s), lambda s: K.G2 * rA(s)),
    'B_fixed_orbit_plus_Sbar': (lambda s: K.G2 * rA(s) + shift(s), lambda s: K.G2 * (rA(s) + Sbar(s))),
    'C_depleted_free_orbit_periodmean': (lambda s: K.G2 * rC(s), lambda s: K.G2 * rC(s)),
    'Cp_depleted_free_orbit_pericentre': (lambda s: K.G2 * rCp(s), lambda s: K.G2 * rCp(s)),
    'D_depleted_plus_Sbar': (lambda s: K.G2 * rC(s) + shift(s), lambda s: K.G2 * (rC(s) + Sbar(s))),
    'E_measured_hartree': (lambda s: K.G2 * spl['X2'](s) + 3 * K.LP * spl['P2'](s), lambda s: K.G2 * spl['X2'](s)),
    'F_measured_zero_mode': (lambda s: K.G2 * spl['rb2'](s), lambda s: K.G2 * spl['rb2'](s)),
    'G_measured_zero_mode_plus_Sbar': (lambda s: K.G2 * spl['rb2'](s) + shift(s), lambda s: K.G2 * (spl['rb2'](s) + Sbar(s))),
}
only = sys.argv[3].split(',') if len(sys.argv) > 3 else list(pumps)
kl = np.sqrt(seed['k2'])
kbins = [(0, 2.0), (2.0, 2.5), (2.5, 3.0), (3.0, 5.0), (5.0, 6.0), (6.0, 6.5), (6.5, 6.74), (6.74, 7.27), (7.27, 7.5), (7.5, 20)]
res = dict(label='OWN (sk_lib): same-seed linear delta-phi on factorial pumps; trace = en_neutral reproduction of Codex engine N48',
           rho0_trace=rho0, rho0_own_realspace=rho0_own, nmodes=seed['nmodes'], unique_k2=int(len(seed['k2'])),
           a0=a0, C=Cc, n_measured_periods=int(len(peri) - 1), controls={})
nl = (tr['Kphi'] + tr['Gphi'] + tr['Uphi'] + tr['Umix']) / tr['rho']
last_period = (s_tr[rec_tr] >= s_tr[peri[-2]]) & (s_tr[rec_tr] <= s_tr[peri[-1]])
for name in only:
    m2f, mintf = pumps[name]
    t1 = time.time()
    mL, mM, mR = m2f(sL), m2f(sM), m2f(sR)
    rec = K.propagate(seed['k2'], mL, mM, mR, s_f, rec_f)
    vals = np.empty(len(rec_f)); fin_e = None
    mint_rec = mintf(s_f[rec_f])
    for j, i in enumerate(rec_f):
        e = K.energy_from_T(rec[i], seed, mint_rec[j]); vals[j] = e.sum() / rho0
        if j == len(rec_f) - 1: fin_e = e
    kb = {f'{a}-{b}': float(fin_e[(kl >= a) & (kl < b)].sum() / fin_e.sum()) for a, b in kbins}
    detT = float(np.max(abs(np.linalg.det(rec[rec_f[-1]]) - 1)))
    sv = np.linalg.svd(rec[rec_f[-1]], compute_uv=False)[:, 0]          # largest singular value per unique k^2
    isv = int(np.argmax(sv)); band = (kl >= 6.74) & (kl <= 7.27)
    amp = dict(max_ln_sigma=float(np.log(sv[isv])), k_lat_at_max=float(kl[isv]), matched_band_max_ln_sigma=float(np.log(sv[band].max())))
    def at(x): return float(vals[int(np.argmin(abs(tau_rec - x)))])
    res['controls'][name] = dict(
        at={str(x): at(x) for x in (20, 40, 60, 70, 84, 100, 150, 200, 300, 400, 500, 600, 700, 800)},
        final=float(vals[-1]), max=float(vals.max()), last_measured_period_mean=float(vals[last_period].mean()),
        final_energy_share_by_klat=kb, k_lat_energy_weighted_final=float((kl * fin_e).sum() / fin_e.sum()),
        transfer_amplification_tau800=amp,
        det_err=detT, seconds=time.time() - t1)
    print(name, 'final %.5g lastperiod %.5g max %.5g det %.1e' % (vals[-1], vals[last_period].mean(), vals.max(), detT),
          {x: round(at(x), 5) for x in (60, 100, 200, 400, 600, 800)}, {a: round(b, 4) for a, b in kb.items()}, amp, '%.1fs' % (time.time() - t1), flush=True)
res['nonlinear_N48_trace'] = dict(final=float(nl[-1]), last_measured_period_mean=float(nl[rec_tr][last_period].mean()))
res['depleted_orbit_segments_periodmean'] = [dict(tau_start=float(np.interp(d['s_start'], s_tr, tau_tr)), T=d['T'], E_over_rho0=d['E'] / rho0,
                                                  Q=d['Q'], rmax_over_rmin=d['ecc']) for d in infoC]
res['seconds'] = time.time() - t0
prev = json.loads(out.read_text()) if out.exists() else None
if prev is not None and set(only) != set(pumps):
    prev['controls'].update(res['controls']); res['controls'] = prev['controls']
out.write_text(json.dumps(res, indent=1) + '\n')
print('rho0 trace %.10f own %.10f; nonlinear final %.6g, last-period mean %.6g' % (rho0, rho0_own, nl[-1], res['nonlinear_N48_trace']['last_measured_period_mean']))
for d in res['depleted_orbit_segments_periodmean'][::4]: print(d)
