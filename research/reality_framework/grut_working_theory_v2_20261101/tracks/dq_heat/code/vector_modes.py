#!/usr/bin/env python3
"""dq_heat (Claude, own code), item 1: longitudinal hidden-vector production on the actual M1r background.

Mode equation: the longitudinal Proca mode of hw_gates/code/gate1_vector.py (copied evolve_mode, generalised to an
arbitrary background), S_k = 1/2 int deta [z^2 |B'|^2 - a^2 m^2 |B|^2], z^2 = a^2 m^2/(k^2 + a^2 m^2), m(eta) = e F(eta),
F^2 = 2 (|X|^2 + 9|Y|^2 + 36|Sigma|^2) (charges 1, 3, 6; single-Higgs-direction approximation, as gate 1).
Pair (B, Pi = z^2 B'): B' = Pi/z^2, Pi' = -a^2 m^2 B; adiabatic vacuum at k/aH = 1000; n_k from the Hamiltonian over the
frequency.  N_A = int dlnk k^3 n_k/(2 pi^2) per comoving H_i^3 (a = 1 at the preparation slice, the frame of N_X, E_X).

Background (conformal time eta, eta = 0 at the slice):
  inflation  eta < 0 : de Sitter H = H_i, a = 1/(1 - eta), F = F_inf of the M1r prepared state (275.93 H_i);
  reheating  0 <= eta <= eta_end : the recorded rptri_dq run (dense 0-200, coarse 200-2000), cubic splines of a and a^2 F^2;
  radiation  eta > eta_end : a = a_end + K (eta - eta_end) (K = a^2 H = 1.98581), a^2 F^2 = 2 <chi^2> + 72 |Sigma_vac|^2 a^2.
Variants: 'full' (the recorded F including the X oscillation), 'noX' (X part of F^2 removed: the Y + Sigma Higgs norm only),
'inst' (gate 1's instant reheating, a = 1 + eta, K = 1, with the P1 profile m^2 = m_vac^2 + (m_inf^2 - m_vac^2)/a^2; with
spot-point F_inf it reproduces gate1's N_A).
Late read-out: n_k at a = 20 a_* and 30 a_*, a_* = (K/m_vac)^{1/2} (H(a_*) = m_vac); UV modes also at eta_end.
Deposit: one vector -> X Xbar at <a> = (pi K/(2 Gamma_A))^{1/2}, Gamma_A = e^2 m_A beta^3/(48 pi), m_A = 6 e h = 240 e H_i;
R_vec = N_A m_A <a> / E_X.
Usage: python3 -B vector_modes.py run full:1e-6 ... ; python3 -B vector_modes.py assemble [out.json]
"""
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/dq_heat')
import os
BGTAG = os.environ.get('DQ_BG', 'm1r')          # 'm1r' (portals vanish at M_Y) or 'planck' (vanish at M_P)
BG = json.load(open(RES / ('background.json' if BGTAG == 'm1r' else 'background_planck.json')))
EX = BG['validation']['EX']
H_HAT = 40.0
M_VAC_PER_E = 6 * H_HAT
F_INF = BG['prep']['F_inf']
K_LATE = BG['late']['K_end']
MX_OVER_H = None    # set per H_i where needed (beta^3 = 1 to 1e-30)

_bg = {}


def load_bg():
    if _bg:
        return _bg
    d = np.load(SP / 'runs' / ('m1r_background.npz' if BGTAG == 'm1r' else 'planck_background.npz'))
    # dense 0-200 then coarse beyond 200
    ed, ad, F2d, zd = d['d_eta'], d['d_a'], d['d_F2'], d['d_z']
    ec, ac, F2c, zc = d['c_eta'], d['c_a'], d['c_F2'], d['c_z']
    tc, td = d['c_t'], d['d_t']
    sel = tc > td[-1] + 1e-9
    # coarse eta re-based on the dense eta at tau = 200 (the coarse trapezoid is less accurate early)
    e200c = np.interp(td[-1], tc, ec)
    eta = np.concatenate([ed, ec[sel] - e200c + ed[-1]])
    a = np.concatenate([ad, ac[sel]])
    F2 = np.concatenate([F2d, F2c[sel]])
    z = np.concatenate([zd, zc[sel]])
    X2 = 2 * np.abs(z[:, 0]) ** 2
    _bg.update(eta=eta, a=a, F2=F2, aF2_full=a * a * F2, aF2_noX=a * a * (F2 - X2))
    _bg['spl_a'] = CubicSpline(eta, a)
    _bg['spl_full'] = CubicSpline(eta, _bg['aF2_full'])
    _bg['spl_noX'] = CubicSpline(eta, _bg['aF2_noX'])
    _bg['eta_end'] = float(eta[-1]); _bg['a_end'] = float(a[-1])
    sel = d['c_t'] > 1600
    _bg['chi2'] = float(np.mean((d['c_a'][sel] * np.abs(d['c_z'][sel, 0])) ** 2))
    _bg['S2'] = float(np.mean(np.abs(d['c_z'][sel, 2]) ** 2))
    return _bg


def make_bg(variant, F_inf=F_INF, F_vac=M_VAC_PER_E):
    """returns functions a(eta), A2(eta) = a^2 F^2 (m = e F) and the late constant K."""
    if variant == 'inst':
        def a_of(eta):
            return 1.0 / (1.0 - eta) if eta < 0 else 1.0 + eta

        def A2(eta):
            if eta < 0:
                return a_of(eta) ** 2 * F_inf ** 2
            a = a_of(eta)
            return a * a * (F_vac ** 2 + (F_inf ** 2 - F_vac ** 2) / (a * a))
        return a_of, A2, 1.0
    B = load_bg()
    spl = B['spl_full'] if variant == 'full' else B['spl_noX']
    e_end, a_end = B['eta_end'], B['a_end']
    late_c = (2 * B['chi2'] if variant == 'full' else 0.0)
    if variant == 'noX':          # the same Higgs norm without X also during inflation (no artificial jump at the slice)
        F_inf = float(np.sqrt(B['aF2_noX'][0]))

    def a_of(eta):
        if eta < 0:
            return 1.0 / (1.0 - eta)
        if eta <= e_end:
            return float(B['spl_a'](eta))
        return a_end + K_LATE * (eta - e_end)

    def A2(eta):
        if eta < 0:
            return a_of(eta) ** 2 * F_inf ** 2
        if eta <= e_end:
            return float(spl(eta))
        a = a_of(eta)
        return late_c + 72 * B['S2'] * a * a
    return a_of, A2, K_LATE


def eta_of_a_late(a, variant):
    if variant == 'inst':
        return a - 1.0
    B = load_bg()
    if a <= B['a_end']:
        return float(np.interp(a, B['a'], B['eta']))
    return B['eta_end'] + (a - B['a_end']) / K_LATE


def evolve(k, e, variant, ends_eta):
    a_of, A2, K = make_bg(variant)
    a0 = k / 1000.0
    eta0 = 1.0 - 1.0 / a0
    e2 = e * e

    def rhs(eta, y):
        am2 = e2 * A2(eta)
        z2 = am2 / (k * k + am2)
        return [y[2] / z2, y[3] / z2, -am2 * y[0], -am2 * y[1]]

    am2 = e2 * A2(eta0); z2 = am2 / (k * k + am2); om = np.sqrt(k * k + am2)
    Bc = 1.0 / np.sqrt(2 * z2 * om); Pc = -1j * np.sqrt(z2 * om / 2)
    y = np.array([Bc.real, Bc.imag, Pc.real, Pc.imag])
    B = None if variant == 'inst' else load_bg()
    cuts = [eta0, 0.0] + ([B['eta_end']] if B else []) + list(ends_eta)
    segs = sorted(set(c for c in cuts if c >= eta0))
    out = {}
    for t1, t2 in zip(segs[:-1], segs[1:]):
        ms = 0.02 if (B is not None and 0 <= t1 < B['eta_end']) else np.inf
        sol = solve_ivp(rhs, (t1, t2), y, method='DOP853', rtol=1e-10, atol=1e-13 * float(np.max(np.abs(y))), max_step=ms)
        if not sol.success:
            raise RuntimeError(sol.message)
        y = sol.y[:, -1]
        am2 = e2 * A2(t2); z2 = am2 / (k * k + am2); om = np.sqrt(k * k + am2)
        Bc = y[0] + 1j * y[1]; Pc = y[2] + 1j * y[3]
        out[t2] = float((abs(Pc) ** 2 / (2 * z2) + am2 * abs(Bc) ** 2 / 2) / om - 0.5)
        out[('state', t2)] = (Bc, Pc, am2, z2)
    return out


def n_conformal(k, e, variant, eta, state):
    """number of the conformal variable v = z B (v' = Pi/z + z' B) with omega^2 = k^2 + a^2 m^2 (z''/z dropped: it vanishes
    in radiation for k >> a m).  For k >> a_end H_end this is the adiabatic invariant that the late (B, Pi) number tends to."""
    a_of, A2, K = make_bg(variant)
    Bc, Pc, am2, z2 = state
    h = 1e-4
    A2p = e * e * (A2(eta + h) - A2(eta - h)) / (2 * h)
    z = np.sqrt(z2); dz2 = A2p * k * k / (k * k + am2) ** 2; zp = dz2 / (2 * z)
    v = z * Bc; vp = Pc / z + zp * Bc
    om = np.sqrt(k * k + am2)
    return float((abs(vp) ** 2 + om * om * abs(v) ** 2) / (2 * om) - 0.5)


K_UV = 0.3      # modes with k >= K_UV are read one unit of eta after the recorded run (k (eta_end + 1) >= 13)


def job(args):
    """n_k from the conformal number of v = z B (conserved in radiation for k >> a m and equal to the (B, Pi) number once
    z' -> 0); the (B, Pi) number is kept for comparison (it carries a +-(z'/z)/omega cross term)."""
    k, e, variant = args[:3]
    force_late = len(args) > 3 and args[3]
    K = 1.0 if variant == 'inst' else K_LATE
    m = M_VAC_PER_E * e
    a_star = np.sqrt(K / m)
    a1, a2 = max(20 * a_star, 30.0), max(30 * a_star, 45.0)
    ends = [eta_of_a_late(a1, variant), eta_of_a_late(a2, variant)]
    uv = (variant != 'inst') and k >= K_UV and not force_late
    if uv:
        e1 = load_bg()['eta_end'] + 1.0
        r = evolve(k, e, variant, [e1])
        nv = n_conformal(k, e, variant, e1, r[('state', e1)])
        return dict(k=k, e=e, variant=variant, n_k=nv, n_k_BPi=r[e1], readout='conformal_eta_end+1', rel_change=0.0)
    r = evolve(k, e, variant, ends)
    n1 = n_conformal(k, e, variant, ends[0], r[('state', ends[0])])
    n2 = n_conformal(k, e, variant, ends[1], r[('state', ends[1])])
    row = dict(k=k, e=e, variant=variant, n_k=n2, n_k_20=n1, n_k_BPi=r[ends[1]], rel_change=abs(n2 - n1) / max(abs(n2), 1e-300),
               readout='conformal_late_20_30_astar')
    return row


def decay(e, K):
    m = M_VAC_PER_E * e
    G = e * e * m / (48 * np.pi)
    a_mean = np.sqrt(np.pi * K / (2 * G))
    return dict(m_A=m, Gamma_A=G, a_mean=a_mean, a_unit=np.sqrt(K / G), deposit_per_vector_EX=m * a_mean / EX)


def spectrum_rows(rows, kgrid):
    lk = np.log(kgrid)
    n = np.array([r['n_k'] for r in rows])
    integ = kgrid ** 3 * n / (2 * np.pi ** 2)
    return float(np.trapezoid(integ, lk)), integ


KG = np.geomspace(1e-5, 30.0, 61)
E_LIST = [1e-7, 1e-6, 1.5e-6, 2e-6, 3e-6, 5e-6, 1e-5]
SETS = ([(e, 'full') for e in E_LIST] + [(e, 'noX') for e in (1e-6, 1e-5)] + [(e, 'inst') for e in (1e-7, 1e-6)]) if BGTAG == 'm1r' else [(e, 'full') for e in E_LIST[1:]]


def setfile(e, var):
    return SP / 'vec_rows' / ('rows_%s%s_%.3e.json' % ('' if BGTAG == 'm1r' else 'planck_', var, e))


def run_sets(specs):
    (SP / 'vec_rows').mkdir(parents=True, exist_ok=True)
    load_bg()
    with Pool(2) as pool:
        for e, var in specs:
            t0 = time.time()
            rows = pool.map(job, [(k, e, var) for k in KG], chunksize=2)
            json.dump(rows, open(setfile(e, var), 'w'), default=float)
            print('set', var, e, 'done %.1fs' % (time.time() - t0), flush=True)


def assemble(out):
    load_bg()
    res = dict(method=__doc__, EX=EX, K=K_LATE, F_inf=F_INF, F_vac=M_VAC_PER_E, kgrid=KG.tolist(), cases=[])
    for e, var in SETS:
        if not setfile(e, var).exists():
            print('missing', var, e); continue
        rr = json.load(open(setfile(e, var)))
        rr.sort(key=lambda r: r['k'])
        NA, integ = spectrum_rows(rr, KG)
        K = 1.0 if var == 'inst' else K_LATE
        dc = decay(e, K)
        m = dc['m_A']; a_star = np.sqrt(K / m); k_star = a_star * m
        lo = KG <= 1.0
        NA_lo = float(np.trapezoid(integ[lo], np.log(KG[lo])))
        NA_hi = NA - NA_lo
        bpi = max(abs(r['n_k'] - r['n_k_BPi']) / max(abs(r['n_k']), 1e-300) for r in rr if 'n_k_BPi' in r and r['n_k'] > 1e-3)
        ipk = int(np.argmax(integ))
        # band split as kin_source: k/k* in [0.2, 5]
        band = (KG >= 0.2 * k_star) & (KG <= 5 * k_star)
        kb = KG[band]; NA_band = float(np.trapezoid(integ[band], np.log(kb)))
        case = dict(e=e, variant=var, N_A=NA, N_A_k_le_1=NA_lo, N_A_k_gt_1=NA_hi, N_A_band_0p2_5kstar=NA_band, peak_k=float(KG[ipk]),
                    k_star=k_star, a_star=a_star, max_rel_change_20_30=max(r['rel_change'] for r in rr), max_BPi_vs_conformal=bpi,
                    R_vec=NA * dc['deposit_per_vector_EX'], R_vec_k_le_1=NA_lo * dc['deposit_per_vector_EX'],
                    R_vec_k_gt_1=NA_hi * dc['deposit_per_vector_EX'], N_A_per_net_charge=NA / abs(BG['validation']['NX_400_1000']),
                    **dc, spectrum=[dict(k=float(k), n_k=r['n_k'], dN_dlnk=float(i)) for k, r, i in zip(KG, rr, integ)])
        res['cases'].append(case)
        print('e=%.2e %-4s N_A=%.5f (k<=1 %.5f, k>1 %.3e, band %.4f) peak k=%.3g k*=%.3g a*=%.1f <a>=%.4e dep/vec=%.4f R_vec=%.4f chg=%.1e' % (
            e, var, NA, NA_lo, NA_hi, NA_band, case['peak_k'], k_star, a_star, dc['a_mean'], dc['deposit_per_vector_EX'], case['R_vec'],
            case['max_rel_change_20_30']), flush=True)
    full = [c for c in res['cases'] if c['variant'] == 'full' and c['e'] >= 1e-6]
    le = np.log([c['e'] for c in full])
    res['scaling_1e-6_to_1e-5'] = dict(dlnNA_dlne=float(np.polyfit(le, np.log([c['N_A'] for c in full]), 1)[0]),
                                      dlnNA_kle1_dlne=float(np.polyfit(le, np.log([c['N_A_k_le_1'] for c in full]), 1)[0]),
                                      dlnNA_kgt1_dlne=float(np.polyfit(le, np.log([max(c['N_A_k_gt_1'], 1e-300) for c in full]), 1)[0]),
                                      dlnR_dlne=float(np.polyfit(le, np.log([c['R_vec'] for c in full]), 1)[0]),
                                      dlnDep_dlne=float(np.polyfit(le, np.log([c['deposit_per_vector_EX'] for c in full]), 1)[0]),
                                      R_vec_times_e_over_1e6=[c['R_vec'] * c['e'] / 1e-6 for c in full],
                                      NA_times_sqrt_e_over_1e6=[c['N_A'] * np.sqrt(c['e'] / 1e-6) for c in full])
    json.dump(res, open(out, 'w'), indent=1, default=float)
    print(json.dumps(res['scaling_1e-6_to_1e-5'], indent=1))


if __name__ == '__main__':
    # python3 -B vector_modes.py run full:1e-6 noX:1e-5 ...   |   python3 -B vector_modes.py assemble [out.json]
    if sys.argv[1] == 'run':
        run_sets([(float(x.split(':')[1]), x.split(':')[0]) for x in sys.argv[2:]])
    else:
        assemble(sys.argv[2] if len(sys.argv) > 2 else str(RES / ('vector_modes.json' if BGTAG == 'm1r' else 'vector_modes_planck.json')))
