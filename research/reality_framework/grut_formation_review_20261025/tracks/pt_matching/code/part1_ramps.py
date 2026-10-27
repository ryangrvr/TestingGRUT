"""Part 1b: free-matching ramps of the starting homogeneous orbit (independent code).
m(s) = 30 f(s/T), lambda(s) = lambda_X (1 - f)  [free matching], or lambda fixed [mass-only control],
f = (1 - cos pi u)/2 (Codex's ramp) or a C-infinity ramp. Final free antiparticle number
n_- = (E_free/m - q)/2 per unit charge; sudden limit T = 0; adiabatic error versus T."""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from pm_common import *

MOUT = 30.0


def phase_states(y4, Tr, phases, nper=8192):
    """States at t = phase*T_r after the tau=20 state (Codex's radial-phase convention)."""
    par = make_par(mode=0, lamramp=0, T=1e-9, mout=0.0, lamX0=LX, g2=0.0, lamP=0.0)
    out = []
    for ph in phases:
        y = y8(y4[0], y4[1], 0.0, y4[2], y4[3], 0.0, s=1.0)
        n = int(round(ph * nper))
        if n:
            c_fixed(par, y, Tr / nper, n, 8, store=False)
        out.append(np.array([y[0], y[1], y[3], y[4]]))
    return out


def ramp(y4, T, lamramp=1, ramp_type=0, h0=0.04, wb0=40.0):
    q = charge(y4)
    if T == 0:
        Ef = 0.5 * (y4[2] ** 2 + y4[3] ** 2 + MOUT ** 2 * (y4[0] ** 2 + y4[1] ** 2))
        return dict(nminus=(Ef / MOUT - q) / 2 / q, work=None, work_rel=None, charge_rel=0.0, quartic_frac=None)
    par = make_par(mode=0, lamramp=lamramp, ramp_type=ramp_type, T=T, mout=MOUT, lamX0=LX, g2=0.0, lamP=0.0)
    y = y8(y4[0], y4[1], 0.0, y4[2], y4[3], 0.0, s=0.0, ps=0.0)
    E0 = c_energy(par, y)
    c_run(par, y, T, h0, wb0, 0.0, 8)
    Ef = c_energy(par, y)          # at s=T: m=30, lambda = 0 (lamramp) or lambda_X
    work = -y[7]
    r2 = y[0] ** 2 + y[1] ** 2
    Efree = 0.5 * (y[3] ** 2 + y[4] ** 2 + MOUT ** 2 * r2)
    qf = y[0] * y[4] - y[1] * y[3]
    res = dict(nminus=(Efree / MOUT - q) / 2 / q, work=work, work_rel=abs(Ef - E0 - work) / Ef, charge_rel=abs(qf / q - 1),
               quartic_frac=(Ef - Efree) / Ef)
    if not lamramp:
        orb = orbit_quadrature(Ef, qf, LX, MOUT ** 2)
        res['final_J_r_per_q'] = orb['J_r'] / q
    return res


def main():
    t0 = time.time()
    y4 = tau20_state(); q = charge(y4)
    E = 0.5 * (y4[2] ** 2 + y4[3] ** 2) + 0.25 * LX * (y4[0] ** 2 + y4[1] ** 2) ** 2
    orb = orbit_quadrature(E, q, LX); Tr = orb['T_r']; Jq = orb['J_r'] / q
    print('J_r/q', Jq, 'T_r', Tr)
    phases8 = np.arange(8) / 8
    st8 = phase_states(y4, Tr, phases8)
    rows = []
    for T in [0.0, 0.2, 2.0, 20.0, 80.0]:
        for ph, s in zip(phases8, st8):
            r = ramp(s, T); r.update(duration=T, phase=ph); rows.append(r)
        a = np.array([r['nminus'] for r in rows if r['duration'] == T])
        print('T', T, 'min', a.min(), 'max', a.max(), 'mean', a.mean(), 'spread', np.ptp(a), time.time() - t0)
    summary = []
    for T in [0.0, 0.2, 2.0, 20.0, 80.0]:
        a = np.array([r['nminus'] for r in rows if r['duration'] == T])
        summary.append(dict(duration=T, min=a.min(), max=a.max(), mean=a.mean(), spread=np.ptp(a), max_abs_dev_from_Jr=np.max(np.abs(a - Jq))))
    # mass-only control (lambda fixed), T = 80, 4 phases
    mass_only = []
    for ph, s in zip(phases8[::2], st8[::2]):
        r = ramp(s, 80.0, lamramp=0); r.update(phase=ph); mass_only.append(r)
        print('mass-only', ph, r['final_J_r_per_q'], r['final_J_r_per_q'] - Jq, r['nminus'], r['quartic_frac'])
    # sudden switch, fine phase scan
    ph512 = np.arange(512) / 512
    st512 = phase_states(y4, Tr, ph512)
    sudden = np.array([ramp(s, 0.0)['nminus'] for s in st512])
    print('sudden fine: min', sudden.min(), 'at', ph512[sudden.argmin()], 'max', sudden.max(), 'at', ph512[sudden.argmax()], 'mean', sudden.mean())
    # analytic sudden extremes: n_- = ((p^2/2 + m^2 r^2/2)/m - q)/2 along the orbit
    # T dependence of the adiabatic error, cosine (Codex) and C-infinity ramps, 16 phases
    ph16 = np.arange(16) / 16; st16 = phase_states(y4, Tr, ph16)
    tdep = []
    for rt, Ts in [(0, [1, 2, 5, 10, 20, 40, 80, 160, 320]), (1, [1, 2, 5, 10, 20, 40, 80, 160])]:
        for T in Ts:
            a = np.array([ramp(s, float(T), ramp_type=rt)['nminus'] for s in st16])
            d = a - Jq
            tdep.append(dict(ramp_type=['cosine', 'Cinfinity'][rt], T=T, T_over_Tr=T / Tr, mean_dev=d.mean(), max_abs_dev=np.abs(d).max(),
                             rms_dev=np.sqrt(np.mean(d ** 2)), min=a.min(), max=a.max()))
            print(tdep[-1], time.time() - t0)
    # power-law fit for cosine ramp in the asymptotic range
    cs = [t for t in tdep if t['ramp_type'] == 'cosine' and t['T'] >= 20]
    slope = np.polyfit(np.log([t['T'] for t in cs]), np.log([t['max_abs_dev'] for t in cs]), 1)[0]
    # precision control: halve the step for T=80
    prec = [abs(ramp(s, 80.0, h0=0.02)['nminus'] - ramp(s, 80.0)['nminus']) for s in st8[:2]]
    save_json('matching_ramps.json', dict(
        J_r_per_q=Jq, T_r=Tr, rows=rows, summary=summary, mass_only=mass_only,
        sudden_fine=dict(n=512, min=sudden.min(), argmin_phase=ph512[sudden.argmin()], max=sudden.max(), argmax_phase=ph512[sudden.argmax()],
                         mean=sudden.mean(), mean_minus_Jr=sudden.mean() - Jq),
        T_dependence=tdep, cosine_powerlaw_slope_T_ge_20=slope, step_halving_change_T80=prec,
        codex=dict(sudden_range=[0.07693913569092221, 18.17569609823104], T20_spread=0.00035612459385570006, T80_spread=1.998900231803269e-05,
                   T80_mean=1.7637252784479496, mass_only_final_Jr=[1.7637247989859752, 1.7637253976993104, 1.7637257576229226, 1.7637251604443436])))
    print('slope', slope, 'prec', prec, 'elapsed', time.time() - t0)


if __name__ == '__main__':
    main()
