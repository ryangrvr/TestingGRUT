"""Track `dissipation`, script 3: basin of the Minkowski end state in (bare constant lambda) x (initial data),
for the three KV dissipative models, with stated grids and convergence checks.

Grids (all stated in ../results/basin_maps.json):
  KV16a: lambda in [-1,1] step 0.05 (41) x f_i in [0.20,3.00] step 0.05 (57); r_r,i = 10 (radiation, units E_*^4);
         dust r_d,i in {0, 10, 100, 1000}; gamma in {1/(32 pi^2) (graviton value, eq 5.5), 0.01 (paper figures)}.
         lambda is in units E_*^4 = 9 M^4 (E_P^4 = (8 pi)^2 M^4 = 70.2 E_*^4).
         Classification (rigorous, using monotonicity of r_V, s1 check E):
           'crunch'  : r_V < 0 or h <= 0 reached (r_V can never increase; with diluting matter h^2 -> r_V < 0);
           'mink'    : r_V > 0, h > 0 and r_d < 1e-6 r_V (remaining dust drain < 1e-6 r_V; radiation is traceless);
           'invalid' : h_i^2 = r_V + r_r + r_d <= 0;  'undecided': none of the above by tau_max.
  KSV16: lambda in [-0.12, 1.0] (Minkowski equilibrium with f > 0 exists only for lambda >= -1/8) x f_i in
         [0.10, 1.50]; h_i = 0.55, hdot_i = 0, r_M,i = 0 (paper slice); gamma_f = gamma_h = 1.
         Outcomes: 'dS' (h > 0, |hdot| < 1e-3 h^2 at tau_max) / 'collapse' (h crosses 0) / 'undecided'.
         Separatrix f_sep(lambda) by bisection (tol 1e-10); asymptotic r_V near it -> fine-tuning slope.
  K22:   mu = 4; lambda in [-1, 5] step 0.25 (25) x xi_i in [0.1, 6.0] step 0.1 (60); h_i = 9; r_M from Friedmann;
         Gamma = g~ |h| h^4 (g~ = 1e-3, paper) and Gamma_pp = gamma |h| (hdot+2h^2)^2 (gamma = 1e-3).
         Outcomes: 'mink' / 'crunch' / 'singular' (xi reaches mu/2 with r_V > 0: xi-dot -> infinity) / 'invalid'.
Convergence: (i) rtol 1e-9 vs 1e-11 on every grid point; (ii) tau_max x10 on all undecided points;
             (iii) boundary located by bisection vs the analytic boundary where one exists.
Run: python3 -I s3_basin.py   (about 10-20 min on 4 cores)
"""
import json, os, sys, time
import numpy as np
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from diss_models import *  # noqa

GPHYS = 1.0 / (32.0 * np.pi**2)


# =========================================================================================== KV16a
def kv_classify(args):
    lam, f1, rr1, rd1, gam, tmax, rtol = args
    rV1 = kv16a_rV(f1, lam)
    h2 = rV1 + rr1 + rd1
    if h2 <= 0:
        return 'invalid', None
    if rd1 == 0.0:
        pass  # exact theorem: sign(r_V,i) decides; still integrate as a check
    y0 = [np.sqrt(h2), f1, rd1, 0.0]

    def e_rv(t, y, *a):
        return kv16a_rV(y[1], lam)
    e_rv.terminal = True; e_rv.direction = -1

    def e_h(t, y, *a):
        return y[0]
    e_h.terminal = True; e_h.direction = -1

    def e_cert(t, y, *a):
        return kv16a_rV(y[1], lam) - 1e6 * y[2]
    e_cert.terminal = True; e_cert.direction = 1

    if rV1 < 0:
        return 'crunch', 1.0
    if rV1 - 1e6 * rd1 > 0:   # already certified (radiation only): integrate to tmax as a check
        sol = integrate(kv16a_rhs, 1.0, tmax, y0, (lam, gam), rtol=rtol, atol=1e-14, events=[e_rv, e_h])
        if sol.status == 1:
            return 'crunch', float(sol.t[-1])
        return 'mink', float(kv16a_rV(sol.y[1, -1], lam))
    sol = integrate(kv16a_rhs, 1.0, tmax, y0, (lam, gam), rtol=rtol, atol=1e-14, events=[e_rv, e_h, e_cert])
    if sol.status == 1:
        if len(sol.t_events[2]) > 0:
            return 'mink', float(kv16a_rV(sol.y[1, -1], lam))
        return 'crunch', float(sol.t[-1])
    return 'undecided', float(kv16a_rV(sol.y[1, -1], lam))


def kv_boundary_bisect(args):
    """smallest r_V,i in the basin at fixed lambda (bisection in f_i)."""
    lam, rr1, rd1, gam, tmax = args
    lo, hi = 0.2, 3.0   # lo: inside (large r_V); hi: outside
    if kv_classify((lam, lo, rr1, rd1, gam, tmax, 1e-10))[0] != 'mink':
        return None
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        c = kv_classify((lam, mid, rr1, rd1, gam, tmax, 1e-10))[0]
        if c == 'mink':
            lo = mid
        elif c in ('crunch', 'invalid'):
            hi = mid
        else:
            return {'lambda': lam, 'status': 'undecided at f=%.6f' % mid}
        if hi - lo < 1e-9:
            break
    return {'lambda': lam, 'f_boundary': lo, 'r_V,i_at_boundary': float(kv16a_rV(lo, lam)),
            'analytic_f_star_rad_only': float(kv16a_fbar(lam))}


# =========================================================================================== KSV16
def ev_turn(tau, y, *a):
    return y[2]
ev_turn.terminal = True; ev_turn.direction = -1


def ksv_outcome(args):
    lam, f1, h1, tmax, rtol = args
    y0 = ksv_initial(h1, 0.0, f1, 0.0, lam)
    sol = integrate(ksv_rhs, 1.0, tmax, y0, (1.0, 1.0, 1 / 3, lam), rtol=rtol, atol=1e-14, events=[ev_turn])
    yy = sol.y[:, -1]
    if sol.status == 1:
        return 'collapse', float(sol.t[-1]), float(yy[3])
    hdot = ksv_epsp(yy[0]) - yy[3] - 2 * yy[2]**2
    if yy[2] > 0 and abs(hdot) < 1e-3 * yy[2]**2:
        return 'dS', float(ksv_rVtilde(yy, lam)), float(yy[3])
    return 'undecided', float(ksv_rVtilde(yy, lam)), float(yy[3])


def ksv_separatrix(args):
    """all collapse <-> dS transitions along f_i at fixed lambda; bisection; r_V on the dS side near each."""
    lam, h1, tmax = args
    fs = np.arange(0.10, 1.5001, 0.05)
    outs = [ksv_outcome((lam, f, h1, tmax, 1e-10))[0] for f in fs]
    trans = []
    for i in range(len(fs) - 1):
        a, b = outs[i], outs[i + 1]
        if {a, b} == {'collapse', 'dS'}:
            lo, hi = fs[i], fs[i + 1]          # lo side has outcome a
            for _ in range(36):
                mid = 0.5 * (lo + hi)
                o = ksv_outcome((lam, mid, h1, tmax, 1e-11))[0]
                if o == a:
                    lo = mid
                elif o == b:
                    hi = mid
                else:
                    break
            fsep = 0.5 * (lo + hi)
            side = +1 if b == 'dS' else -1     # direction into the dS side
            near = []
            for d in (1e-3, 1e-4, 1e-5):
                o = ksv_outcome((lam, fsep + side * d, h1, tmax, 1e-11))
                o4 = ksv_outcome((lam, fsep + side * d, h1, 4 * tmax, 1e-11))
                near.append({'delta_f': d, 'outcome': o[0], 'rVtilde(tau_max)': o[1],
                             'rVtilde(4 tau_max)': o4[1], 'u_eff(4 tau_max)': o4[2]})
            trans.append({'f_sep': fsep, 'from_to_increasing_f': [a, b], 'dS_side_near_separatrix': near,
                          'is_Minkowski_separatrix': bool(near[-1]['rVtilde(4 tau_max)'] is not None and
                                                          abs(near[-1]['rVtilde(4 tau_max)']) < 3e-3)})
    u0, f0 = ksv_u0(lam)
    return {'lambda': lam, 'u0(lambda)': u0, 'f0(lambda)': f0, 'transitions': trans, 'coarse_outcomes': outs}


# =========================================================================================== K22
def k22_classify(args):
    lam, xi1, mu, h1, gt, source, tmax, rtol = args
    y0 = k22_initial(h1, xi1, lam, mu)
    if y0[2] < 0:
        return 'invalid', None
    rV1 = k22_rV(xi1, lam, mu)
    if rV1 < 0:
        return 'crunch', 1.0
    if abs(2 * xi1 - mu) < 1e-12:
        return 'singular', 1.0

    def e_rv(t, y, *a):
        return k22_rV(y[1], lam, mu)
    e_rv.terminal = True; e_rv.direction = -1

    def e_h(t, y, *a):
        return y[0]
    e_h.terminal = True; e_h.direction = -1

    def e_sing(t, y, *a):   # approaching the turning point of r_V(xi): xi -> mu/2
        return abs(2 * y[1] - mu) - 1e-6
    e_sing.terminal = True

    def e_cert(t, y, *a):   # radiation negligible and its remaining drain < r_V/10 (see docstring)
        rV = k22_rV(y[1], lam, mu)
        hM2 = 8 * np.pi * y[2] / 3
        bound = gt * hM2**2 / 8 if source == 'adhoc' else 0.0
        return min(rV - 1e3 * y[2], rV - 10 * bound)
    e_cert.terminal = True; e_cert.direction = 1

    sol = integrate(k22_rhs, 1.0, tmax, y0, (lam, mu, 1 / 3, gt, source), rtol=rtol, atol=1e-16,
                    events=[e_rv, e_h, e_sing, e_cert])
    if sol.status == 1:
        if len(sol.t_events[3]):
            return 'mink', float(k22_rV(sol.y[1, -1], lam, mu))
        if len(sol.t_events[2]):
            return 'singular', float(sol.t[-1])
        return 'crunch', float(sol.t[-1])
    if sol.status == -1:
        return 'singular', float(sol.t[-1])
    return 'undecided', float(k22_rV(sol.y[1, -1], lam, mu))


def frac(grid, key):
    flat = [c for row in grid for c in row]
    valid = [c for c in flat if c != 'invalid']
    return float(sum(c == key for c in valid)) / max(len(valid), 1)


def ksv_only(pool, out, T0):
    lamsK = np.round(np.concatenate([[-0.12, -0.1, -0.05], np.arange(0.0, 1.0001, 0.1)]), 4)
    seps = pool.map(ksv_separatrix, [(l, 0.55, 1500.0) for l in lamsK])
    out['KSV16']['separatrix_by_lambda'] = seps
    return out


if __name__ == '__main__':
    T0 = time.time()
    pool = Pool(4)
    if len(sys.argv) > 1 and sys.argv[1] == 'ksv_separatrix_only':
        fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../results/basin_maps.json')
        out = ksv_only(pool, json.load(open(fn)), T0)
        json.dump(out, open(fn, 'w'), indent=0)
        print('ksv separatrix updated', time.time() - T0)
        sys.exit(0)
    out = {'labels': '[C] computed here unless marked [D]'}
    # ------------------------------------------------------------ KV16a
    lams = np.round(np.arange(-1.0, 1.0001, 0.05), 4)
    f1s = np.round(np.arange(0.20, 3.0001, 0.05), 4)
    kvres = {}
    for gam_name, gam in (('gamma_phys=1/(32pi^2)', GPHYS), ('gamma=0.01', 0.01)):
        for rd1 in (0.0, 10.0, 100.0, 1000.0):
            key = '%s, r_r,i=10, r_d,i=%g' % (gam_name, rd1)
            tmax = 1e3 if rd1 == 0 else 1e6
            jobs = [(lam, f1, 10.0, rd1, gam, tmax, 1e-9) for lam in lams for f1 in f1s]
            r9 = pool.map(kv_classify, jobs, chunksize=16)
            jobs11 = [(lam, f1, 10.0, rd1, gam, tmax, 1e-11) for lam in lams for f1 in f1s]
            r11 = pool.map(kv_classify, jobs11, chunksize=16)
            grid = [[r9[i * len(f1s) + j][0] for j in range(len(f1s))] for i in range(len(lams))]
            nmis = sum(a[0] != b[0] for a, b in zip(r9, r11))
            und = [(j[0], j[1]) for j, r in zip(jobs, r9) if r[0] == 'undecided']
            und10 = pool.map(kv_classify, [(l, f, 10.0, rd1, gam, tmax * 10, 1e-10) for (l, f) in und]) if und else []
            # analytic radiation-only prediction: mink iff r_V,i > 0
            if rd1 == 0:
                pred = [['mink' if kv16a_rV(f, l) > 0 else 'crunch' for f in f1s] for l in lams]
                agree = sum(grid[i][j] == pred[i][j] for i in range(len(lams)) for j in range(len(f1s)))
            else:
                agree = None
            bis = pool.map(kv_boundary_bisect, [(l, 10.0, rd1, gam, tmax) for l in (-1.0, -0.5, 0.0, 0.5, 1.0)])
            kvres[key] = {'grid_lambda': lams.tolist(), 'grid_f_i': f1s.tolist(), 'classes': grid,
                          'fraction_mink_of_valid': frac(grid, 'mink'), 'fraction_crunch_of_valid': frac(grid, 'crunch'),
                          'n_undecided': len(und), 'undecided_after_10x_tmax': [r[0] for r in und10],
                          'n_class_changes_rtol_1e-9_vs_1e-11': int(nmis),
                          'agreement_with_analytic_rV_i>0_(rad only)': agree,
                          'boundary_by_bisection': bis, 'tau_max': tmax}
            print(key, 'mink frac', kvres[key]['fraction_mink_of_valid'], 'undec', len(und), 'mis', nmis,
                  'agree', agree, time.time() - T0, flush=True)
    out['KV16a'] = kvres
    out['KV16a_note_[D]'] = ('Radiation only: r_V obeys dr_V/dxi = -4 gamma r_V^2 exactly, so the basin is exactly '
                             '{r_V,i >= 0, h_i > 0}: f_i <= f_*(lambda) = sqrt((lambda + sqrt(lambda^2+12))/2). With dust, '
                             'R = (rho_d + 4 rho_V)/K and the drain pushes r_V below 0 unless r_V,i exceeds ~ gamma r_d,i^2/24.')
    # ------------------------------------------------------------ KSV16
    lamsK = np.round(np.concatenate([[-0.12, -0.1, -0.05], np.arange(0.0, 1.0001, 0.1)]), 4)
    f1K = np.round(np.arange(0.10, 1.5001, 0.05), 4)
    jobs = [(l, f, 0.55, 1500.0, 1e-9) for l in lamsK for f in f1K]
    rK = pool.map(ksv_outcome, jobs, chunksize=8)
    rK11 = pool.map(ksv_outcome, [(l, f, 0.55, 1500.0, 1e-11) for l in lamsK for f in f1K], chunksize=8)
    gridK = [[rK[i * len(f1K) + j][0] for j in range(len(f1K))] for i in range(len(lamsK))]
    rVinf = [[rK[i * len(f1K) + j][1] if rK[i * len(f1K) + j][0] == 'dS' else None for j in range(len(f1K))]
             for i in range(len(lamsK))]
    seps = pool.map(ksv_separatrix, [(l, 0.55, 1500.0) for l in lamsK])
    dsv = [v for row in rVinf for v in row if v is not None]
    # (h_i, f_i) map at lambda = 0
    hsK = np.round(np.arange(0.45, 0.9001, 0.025), 4)
    f2K = np.round(np.arange(0.10, 0.6001, 0.025), 4)
    rHF = pool.map(ksv_outcome, [(0.0, f, h, 1500.0, 1e-9) for f in f2K for h in hsK], chunksize=8)
    gridHF = [[rHF[i * len(hsK) + j][0] for j in range(len(hsK))] for i in range(len(f2K))]
    out['KSV16'] = {'grid_lambda': lamsK.tolist(), 'grid_f_i': f1K.tolist(), 'h_i': 0.55, 'tau_max': 1500.0,
                    'classes': gridK, 'rVtilde_at_tau_max_for_dS': rVinf,
                    'fraction_dS': frac(gridK, 'dS'), 'fraction_collapse': frac(gridK, 'collapse'),
                    'fraction_undecided': frac(gridK, 'undecided'),
                    'n_class_changes_rtol_1e-9_vs_1e-11': int(sum(a[0] != b[0] for a, b in zip(rK, rK11))),
                    'min_rVtilde_among_dS_grid_points': min(dsv) if dsv else None,
                    'median_rVtilde_among_dS_grid_points': float(np.median(dsv)) if dsv else None,
                    'separatrix_by_lambda': seps,
                    'h_f_map_lambda0': {'grid_h_i': hsK.tolist(), 'grid_f_i': f2K.tolist(), 'classes': gridHF,
                                        'fraction_dS': frac(gridHF, 'dS'), 'fraction_collapse': frac(gridHF, 'collapse')},
                    'note_[D]': 'u_eff is non-decreasing (s >= 0); every static state (de Sitter family) has s = 0, '
                                'so u_eff freezes below u0 (dS) or overshoots (collapse); Minkowski = separatrix.'}
    print('KSV16 done', time.time() - T0, flush=True)
    # ------------------------------------------------------------ K22
    lamsB = np.round(np.arange(-1.0, 5.0001, 0.25), 4)
    xisB = np.round(np.arange(0.1, 6.0001, 0.1), 4)
    k22res = {}
    for source, gt in (('adhoc', 1e-3), ('pp', 1e-3)):
        jobs = [(l, x, 4.0, 9.0, gt, source, 1e5, 1e-9) for l in lamsB for x in xisB]
        r = pool.map(k22_classify, jobs, chunksize=16)
        r11 = pool.map(k22_classify, [(l, x, 4.0, 9.0, gt, source, 1e5, 1e-11) for l in lamsB for x in xisB],
                       chunksize=16)
        grid = [[r[i * len(xisB) + j][0] for j in range(len(xisB))] for i in range(len(lamsB))]
        und = [(j[0], j[1]) for j, rr in zip(jobs, r) if rr[0] == 'undecided']
        und10 = pool.map(k22_classify, [(l, x, 4.0, 9.0, gt, source, 1e6, 1e-10) for (l, x) in und]) if und else []
        lam_le4 = [i for i, l in enumerate(lamsB) if l <= 4.0]
        lam_gt4 = [i for i, l in enumerate(lamsB) if l > 4.0]
        k22res[source] = {'gamma': gt, 'classes': grid, 'grid_lambda': lamsB.tolist(), 'grid_xi_i': xisB.tolist(),
                          'fraction_mink_of_valid': frac(grid, 'mink'),
                          'fraction_mink_of_valid_lambda<=mu^2/4': frac([grid[i] for i in lam_le4], 'mink'),
                          'fraction_mink_of_valid_lambda>mu^2/4': frac([grid[i] for i in lam_gt4], 'mink'),
                          'classes_present_lambda>mu^2/4': sorted(set(c for i in lam_gt4 for c in grid[i])),
                          'n_undecided': len(und), 'undecided_after_10x_tmax': [x[0] for x in und10],
                          'n_class_changes_rtol_1e-9_vs_1e-11': int(sum(a[0] != b[0] for a, b in zip(r, r11)))}
        print('K22', source, k22res[source]['fraction_mink_of_valid'], len(und), time.time() - T0, flush=True)
    out['K22'] = k22res
    out['K22_note_[D]'] = ('mu = 4: r_V(xi) = (xi-2)^2 + lambda - 4. Minkowski needs a zero of r_V reachable '
                           'monotonically: lambda <= mu^2/4 = 4 and r_V(xi_i) >= 0. For lambda > 4 the flow reaches '
                           'xi = mu/2 with r_V > 0, where xi-dot = -Gamma/(2 xi - mu) diverges: the model has no '
                           'continuation (classified singular).')
    pool.close()
    out['runtime_s'] = time.time() - T0
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../results/basin_maps.json'), 'w') as fh:
        json.dump(out, fh, indent=0)
    print('done', time.time() - T0)
