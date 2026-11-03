#!/usr/bin/env python3
"""e2e_relic item 1 (Claude): the corrected pair-fraction map p(S/Q) = (Xbar per net charge after freeze-out) of the free-X
branch, as a callable function plus a table (../results/pmap.json).

What changed against the inherited maps (chain_mfi.p_of_SQ, inf_chain chain.p_of_SQ):
  * The gap 6.74 < S/Q < 17.1 is filled with the reviewed case-L conditional moment network (fx_chain/code/kin/m1_kin.py,
    read-only; driver code/pmap_runs.py), not bridged log-linearly.
  * ONE start convention throughout the central map: the 'relaxed' onset start (x0 = ln(a_on/a_ref) interpolated over
    fx_chain's heat_onset rows, R = 0 ... 49.6), which is the convention of the M1r gas map below the gap.  The inherited
    maps switched at S/Q = 6.74 from 'relaxed' (wt_critic/fx_chain gas map) to 'aref' (forced equilibrium at a_ref:
    inf_chain pmap_kin, mf_chain skeptic).  The two differ x1.85 at 6.74 (the 'junction' the mf_chain skeptic found),
    x1.41 at 8, x1.07 at 14, x1.03 at 21.8.  Above S/Q = 21.8 no relaxed rows exist (the onset table ends at R = 49.6);
    there the aref rows of inf_chain are used times 1 - 0.0280 (21.797/S)^2, a two-parameter fit of the computed
    relaxed/aref ratio at S/Q = 14-21.8 [E].  Every one of these is a model choice.
  * The visible-clock yield dependence is explicit: p(S/Q, f) = p(S/Q, 1) f^beta(S/Q), f = Y_Q / 5.07941423e-9, with
    beta computed at S/Q = 5, 7, 10, 21.8 (pmap_runs.py, f = 0.1 and 10) and at 155-958 from inf_chain's self-consistent
    runs.  beta = -0.38 (S/Q = 5) ... -0.02 (S/Q >= 155).  The abundance-calibrated yield f = 0.9374/(1 + 2p) is solved by
    p_selfconsistent().  This replaces inf_chain's sc_ratio interpolation (identical where both exist).
[C] on a conditional recipe (forced chemical equilibrium at the start, no objects, MB fugacities): inherited, not validated.
Usage: python3 -B pmap.py   -> writes ../results/pmap.json and prints the table.
       from pmap import p_of_SQ, p_selfconsistent, M_freeX
"""
import json, math, os
import numpy as np

YQ_M1 = 5.07941423e-9          # M1 calibration Y_Q = D/(0.9374 m_X), D = 4.4e-10 GeV [supplied]

# ---- relaxed-start rows at f = 1 (S/Q, p, source)
ROWS_RELAXED = [
    (3.6388366029832127, 7.544122020366876e-05, 'fx_chain gas map (wt_critic S/Q), relaxed'),
    (4.314953051486718, 0.0002698829640998089, 'fx_chain gas map, relaxed'),
    (4.936310238047122, 0.0009101892722290548, 'fx_chain gas map, relaxed'),
    (5.2682376922802785, 0.0017163247840089558, 'fx_chain gas map, relaxed'),
    (5.573114471699787, 0.002990076314594871, 'fx_chain gas map, relaxed'),
    (5.872535041261307, 0.004967989559836962, 'fx_chain gas map, relaxed'),
    (6.16695544373095, 0.007831057621247556, 'fx_chain gas map, relaxed'),
    (6.742, 0.016618986077024223, 'e2e_relic network, relaxed (= gas map 0.016617)'),
    (7.0, 0.022011356178076858, 'e2e_relic network, relaxed'),
    (7.5, 0.03486343500701146, 'e2e_relic network, relaxed'),
    (8.0, 0.050444186483864645, 'e2e_relic network, relaxed'),
    (9.0, 0.08759137134643431, 'e2e_relic network, relaxed'),
    (10.0, 0.129946303964089, 'e2e_relic network, relaxed'),
    (12.0, 0.22293447922712467, 'e2e_relic network, relaxed'),
    (14.0, 0.3214958560962796, 'e2e_relic network, relaxed'),
    (16.0, 0.4228011813397113, 'e2e_relic network, relaxed'),
    (17.131, 0.48080997258747643, 'e2e_relic network, relaxed'),
    (21.797, 0.7230608321005416, 'e2e_relic network, relaxed'),
]
# ---- aref-start rows (forced equilibrium at a_ref), f = 1
ROWS_AREF = [
    (6.742, 0.030805401666354627, 'mf_chain skeptic sk_pmap'), (8.0, 0.07109283966498119, 'sk_pmap'),
    (10.0, 0.1539440580819891, 'sk_pmap'), (12.0, 0.24721977237645992, 'sk_pmap'), (14.0, 0.3452512496109597, 'sk_pmap'),
    (15.0, 0.39529561359438364, 'sk_pmap'), (16.0, 0.44580846026132825, 'sk_pmap'),
    (17.13126030489702, 0.5033859970367478, 'inf_chain pmap_kin'), (21.79713734357875, 0.7438843557228119, 'inf_chain pmap_kin'),
    (41.608775427103964, 1.780789405424379, 'inf_chain pmap_kin'), (94.4312451625062, 4.543070189881549, 'inf_chain pmap_kin'),
    (232.62293436915292, 11.684549404379283, 'inf_chain pmap_kin'), (530.1016499067407, 26.845567625897147, 'inf_chain pmap_kin'),
    (1307.7296071865874, 65.83102022585068, 'inf_chain pmap_kin'),
]
# single aref point below the gap (fx_chain kin_m1 'm1_gas_aref_R497', S/Q = 4.936): start-convention systematic x3.5 there
AREF_LOW = (4.936310238047122, 0.003216211852310832)
# ---- yield exponent beta = d ln p / d ln Y_Q (relaxed network, f = 0.1 and 10; inf_chain self-consistent runs at >= 155)
BETA = [(5.0, -0.3832), (7.0, -0.1924), (10.0, -0.0832), (21.797, -0.0351),
        (155.0, -0.02127), (302.0, -0.02022), (447.0, -0.01985), (600.0, -0.01963), (958.0, -0.01943)]
RATIO_FIT = (21.797, 0.02799)   # relaxed/aref = 1 - 0.02799 (21.797/S)^2 above 21.797 [E]


def _central_rows():
    rows = [(s, p) for s, p, _ in ROWS_RELAXED]
    a, c = RATIO_FIT
    rows += [(s, p * (1 - c * (a / s) ** 2)) for s, p, _ in ROWS_AREF if s > a + 1e-3]
    rows.sort()
    return np.log([r[0] for r in rows]), np.log([r[1] for r in rows])


_LX, _LY = _central_rows()


def _loglin(lx, X, Y):
    if lx <= X[0]:
        return Y[0] + (Y[1] - Y[0]) / (X[1] - X[0]) * (lx - X[0])
    if lx >= X[-1]:
        return Y[-1] + (Y[-1] - Y[-2]) / (X[-1] - X[-2]) * (lx - X[-1])
    return float(np.interp(lx, X, Y))


def beta(SQ):
    """d ln p / d ln Y_Q, linear in ln S/Q between computed points; held constant outside [5, 958] (flagged [E])."""
    X = [math.log(a) for a, _ in BETA]; Y = [b for _, b in BETA]
    return float(np.interp(math.log(SQ), X, Y))


def p_of_SQ(SQ, YQ_factor=1.0):
    """Central corrected pair fraction (relaxed start), Xbar per net charge, at visible-clock yield Y_Q = YQ_factor * 5.0794e-9."""
    p1 = math.exp(_loglin(math.log(SQ), _LX, _LY))
    return p1 * YQ_factor ** beta(SQ)


def p_aref(SQ, YQ_factor=1.0):
    """Same network with the a_ref (forced-equilibrium) start: the upper edge of the start-convention band, where computed."""
    rows = sorted([(s, p) for s, p, _ in ROWS_AREF] + [AREF_LOW])
    X = np.log([r[0] for r in rows]); Y = np.log([r[1] for r in rows])
    return math.exp(_loglin(math.log(SQ), X, Y)) * YQ_factor ** beta(SQ)


def p_selfconsistent(SQ, YQ_ratio_extra=1.0, tol=1e-10):
    """p at the abundance-calibrated yield Y_Q = Y_Q,M1 * 0.9374/M, M = 1 + 2p (free X), times an optional extra factor
    (e.g. if the chain's own Y_Q differs from the calibrated one).  Fixed-point iteration (converges in < 10 steps)."""
    p = p_of_SQ(SQ)
    for _ in range(60):
        f = 0.9374 / (1 + 2 * p) * YQ_ratio_extra
        pn = p_of_SQ(SQ, f)
        if abs(pn - p) < tol * max(p, 1e-30):
            return pn, f
        p = pn
    return p, f


def M_freeX(SQ, YQ_factor=None):
    """Mass per net charge of the free-X relic in units of m_X: 1 + 2p.  YQ_factor None -> abundance-calibrated."""
    p = p_selfconsistent(SQ)[0] if YQ_factor is None else p_of_SQ(SQ, YQ_factor)
    return 1 + 2 * p


def table(grid=None):
    grid = grid if grid is not None else sorted(set([3.4716, 3.6, 4.0, 4.16, 4.5, 4.955, 5.0, 5.5, 6.0, 6.5, 6.742, 7.0, 7.5,
                                                     8, 9, 10, 11, 12, 12.2, 13, 14, 14.7, 15, 15.9, 16, 17.131, 19, 21.797,
                                                     25, 30, 41.6, 61.5, 94.4, 155, 232.6, 302, 530, 958, 1307.7]))
    old = None
    try:                                                     # inherited bridge, for comparison only (read-only import)
        import importlib.util, sys
        here = os.path.dirname(os.path.abspath(__file__))
        sys.path.insert(0, os.path.join(here, '..', '..', 'mf_ideal', 'code'))
        spec = importlib.util.spec_from_file_location('chain_mfi_ro', os.path.join(here, '..', '..', 'mf_ideal', 'code', 'chain_mfi.py'))
        old = importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
    except Exception:
        old = None
    out = []
    for s in grid:
        psc, f = p_selfconsistent(s)
        r = dict(S_per_Q=s, p_fixed_yield=p_of_SQ(s), p_aref_fixed_yield=p_aref(s) if s >= 4.936 else None,
                 beta=beta(s), p_selfconsistent=psc, YQ_factor_selfconsistent=f, M_selfconsistent=1 + 2 * psc,
                 asym_share=1 / (1 + 2 * psc), Xbar_over_X=psc / (1 + psc))
        if old is not None:
            r['p_inherited_chain_mfi'] = old.p_of_SQ(s)
        out.append(r)
    return out


if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    tab = table()
    res = dict(method=__doc__, rows_relaxed=ROWS_RELAXED, rows_aref=ROWS_AREF, aref_low=AREF_LOW, beta_points=BETA,
               ratio_fit=dict(form='p_relaxed/p_aref = 1 - c (S0/S)^2', S0=RATIO_FIT[0], c=RATIO_FIT[1]),
               YQ_M1=YQ_M1, table=tab)
    json.dump(res, open(os.path.join(here, '..', 'results', 'pmap.json'), 'w'), indent=1, default=float)
    print('%8s %11s %11s %11s %8s %11s %8s' % ('S/Q', 'p(f=1)', 'p_aref', 'p_inher', 'beta', 'p_sc', 'M_sc'))
    for r in tab:
        print('%8.3f %11.4e %11s %11s %8.4f %11.4e %8.4f' % (
            r['S_per_Q'], r['p_fixed_yield'], '%.4e' % r['p_aref_fixed_yield'] if r['p_aref_fixed_yield'] else '-',
            '%.4e' % r['p_inherited_chain_mfi'] if 'p_inherited_chain_mfi' in r else '-', r['beta'], r['p_selfconsistent'],
            r['M_selfconsistent']))
