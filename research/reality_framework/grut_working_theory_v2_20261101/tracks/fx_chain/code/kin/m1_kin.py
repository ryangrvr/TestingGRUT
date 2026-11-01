"""fx_chain item 4 (Claude): the case-L moment network at the GRUT-M1 object and heat, from the relaxed initial state.

Engine: the conditional fast-s moment network of REVIEW_KINETIC (Codex 34f23f8 moment_controls.py over the ported
kin_solve modules, copied unmodified into this directory; reviewed and reproduced in tracks/kin_review).  This driver
patches the inputs at run time; every change is a labelled model choice:
  * M1 object: Q = 2e5, omega = 0.8817204920 m_X, E/Q = 0.9374 m_X; binding object Q = 71263 (E/Q = m_X,
    omega = 0.94017); M1 E(Q) table from tracks/m1_vacuum (replaces the old symmetric-potential table);
  * M1 calibration Y_Q = 5.07941423e-9 (the 'old' 4.8294e-9 is kept for the reproduction runs);
  * M1 false-vacuum energy V(0) = lP v^4/4 + eps v^4 = 95 m_X^4 (old 75);
  * capture kernel (rules b, c) scaled geometrically by (r90_M1/r90_old)^2 = 1.041 [E];
  * relaxed start: the network starts at the relaxed onset a_on (heat_onset.py: circular condensate + thermal plasma
    at the energy per charge (1 + R) E_X/N_X), with rho = rho_c (1 + R)(a_ref/a_on)^4 + V(0), n_Q = n_Q0 (a_ref/a_on)^3,
    in forced chemical equilibrium there (the relaxed state IS chemically relaxed, REVIEW_TRANSPORT);
  * M1 scalar cubic V'''(v) = 6(lP - eps)v in the dominant pair-changing channel ss -> X Xbar: T-dependent factor
    0.904-0.998 from cubic_m1.py (applied to every M1 case unless the name ends in _nocubic).
Outputs p = P/n_Q at the stop (pairs per net charge), f = 2.02855 p (old convention) and the frozen-test term
2.13356 p.  [C] conditional recipes; not a full Boltzmann solution or a formation calculation.
Usage: python3 m1_kin.py <case-name> [<case-name> ...]   (cases below), or 'all'.
"""
import os, sys, json, time, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
SP = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad'
os.environ.setdefault('GRUT_KINETIC_WORKDIR', SP + '/wt2/fx_chain/kwork')
import numpy as np
import ks_common as K

TR = '/home/user/TestingGRUT/research/reality_framework/grut_formation_review_20261025/tracks'
OUT = os.path.join(HERE, '..', '..', 'results')
SOL = json.load(open(TR + '/m1_vacuum/results/soliton.json'))
HO = json.load(open(os.path.join(OUT, 'heat_onset.json')))

OLD = dict(OMEGA_REF=K.OMEGA_REF, E_OVER_QM_REF=K.E_OVER_QM_REF, Q_MIN=K.Q_MIN, YQ=K.YQ, V0=K.V0, Objects=K.Objects)
M1 = dict(OMEGA_REF=SOL['primary']['Omega'], E_OVER_QM_REF=SOL['primary']['E_over_Qm'], Q_MIN=SOL['Q_min']['Q'],
          YQ=5.07941423e-9, V0=K.V0 * 95.0 / 75.0)
Q_BIND, OM_BIND = SOL['E_over_Q_equals_mX']['Q'], SOL['E_over_Q_equals_mX']['Omega']
R90_M1 = SOL['primary_resolution_runs']['L60_h0.0125']['r90']
R90_OLD = K.THM['objects']['reference_Q200000']['r90_over_mX']
CAPTURE_SCALE = {'old': 1.0, 'm1': (R90_M1 / R90_OLD) ** 2}


class M1Objects:
    """M1 stable branch E(Q), Omega(Q) (tracks/m1_vacuum soliton.json EQ_scan + Q_min + binding point)."""
    def __init__(self):
        rows = [(r['Q'], r['Omega'], r['E_over_Qm'], r['r90']) for r in SOL['EQ_scan']]
        rows.append((SOL['Q_min']['Q'], SOL['Q_min']['Omega'], float('nan'), float('nan')))
        rows.append((Q_BIND, OM_BIND, 1.0, float('nan')))
        rows = sorted(r for r in rows if np.isfinite(r[2]))
        self.Q = np.array([r[0] for r in rows]); self.Om = np.array([r[1] for r in rows])
        self.EQ = np.array([r[2] for r in rows]); self.r90 = np.array([r[3] for r in rows])

    def omega(self, Q): return float(np.interp(np.log(Q), np.log(self.Q), self.Om))
    def e_over_q(self, Q): return float(np.interp(np.log(Q), np.log(self.Q), self.EQ))
    def r90_over_m(self, Q): return float(np.interp(np.log(Q), np.log(self.Q), np.nan_to_num(self.r90, nan=5.92)))


def configure(variant):
    P = M1 if variant == 'm1' else OLD
    for k in ('OMEGA_REF', 'E_OVER_QM_REF', 'Q_MIN', 'YQ', 'V0'):
        setattr(K, k, P[k])
    K.Objects = M1Objects if variant == 'm1' else OLD['Objects']


# import the network after the module-level defaults exist (OBJ/RATES are built at import with the old table;
# configure() + reset_obj() swap the object table per case)
import stage2_moments as S2
import moment_controls as MC
_cap0 = S2.capture_kernel
SCALE = dict(capture=1.0, cubic=False)


def capture_scaled(t, extrap):
    sv, ein = _cap0(t, extrap)
    return sv * SCALE['capture'], ein


S2.capture_kernel = capture_scaled
MC.capture_kernel = capture_scaled
_rates0 = S2.Solver.rates


_CUB = json.load(open(os.path.join(OUT, 'cubic_m1.json')))
_CUB_T = np.log([r['T_over_mX'] for r in _CUB['rows']]); _CUB_R = np.array([r['ratio'] for r in _CUB['rows']])


def cubic_ratio(t):
    """M1/untilted gamma(ss -> X Xbar) from cubic_m1.py; held at the end values outside 0.03-30 m_X (threshold
    value 0.9025 below)."""
    return float(np.interp(np.log(t), _CUB_T, _CUB_R))


def rates_scaled(self, st):
    R = _rates0(self, st)
    if SCALE['cubic']:
        g, dp, dn = R['ss_XXb']; R['ss_XXb'] = (g * cubic_ratio(st['t']), dp, dn)
    return R


S2.Solver.rates = rates_scaled


class Case(S2.Case):
    pass


class Control(MC.Control):
    """moment_controls.Control with (i) rules a/b/c and a chosen object charge, (ii) a start at x0 = ln(a0/a_ref)."""

    def __init__(self, R, formation='gas', Q_L=None, x0=0.0, rtol=2e-8, **kw):
        self.x0 = x0
        super().__init__(R, formation='gas' if formation == 'gas' else 'instant', rtol=rtol, **kw)
        if formation in ('instant', 'frag'):
            rule = 'b' if formation == 'instant' else 'c'
            self.case = S2.Case('control', R, rule=rule, e_cap_variant='finite', catalysis=False, Q_L=Q_L)
            self.helper = S2.Solver(self.case); self.B = self.helper.B; self.c = self.case
            self.inv = MC.Inverter(self.c.cap_p, self.c.e_cap, self.c.cap_m, self.c.e_cm)
        self.formation = formation

    def run(self):
        x0 = self.x0
        rho = (K.RHO_C * (1.0 + self.R) * math.exp(-4 * x0) + K.V0) / K.MX ** 4
        nq = K.N_Q0 * math.exp(-3 * x0) / K.MX ** 3
        initial = MC.initial_equilibrium(nq, rho, self.c.cap_p, self.c.e_cap)
        self.initial = initial
        t = initial['t']; up = initial['up']; um = initial['um']; nc = initial['n_cond']
        self.active = [nc > 1e-12, False]
        w = np.array([np.log(t), nc / nq if self.active[0] else up, um, 0.]); self.E0 = rho * math.exp(3 * x0)
        # the rest is moment_controls.Control.run with x = x0 instead of 0 (copied, unchanged otherwise)
        x = x0; trajectories = []; switches = []; segments = []; stopped = 'xend'; self.sols = []
        from scipy.integrate import solve_ivp
        for segment in range(20):
            def cold(x, w): return w[0] - np.log(self.tstop)
            cold.direction = -1; cold.terminal = True
            def tiny(x, w):
                st, _ = self.state(x, w); return st['P'] / max(st['Q'], 1e-300) - 1e-18
            tiny.direction = -1; tiny.terminal = True
            events = [cold, tiny]
            for j in range(2):
                active = self.active[j]; cap = self.c.cap_p if j == 0 else self.c.cap_m
                def transition(x, w, j=j, active=active, cap=cap): return w[j + 1] if active else w[j + 1] - cap
                transition.terminal = True; transition.direction = -1 if active else 1; events.append(transition)
            def closure(x, w): return self.closure_rate(x, w) - 1
            closure.direction = -1; closure.terminal = True
            if self.closure_stop: events.append(closure)
            start_active = list(self.active)
            sol = solve_ivp(self.rhs, (x, self.xend), w, method=self.method, jac=self.jac, rtol=self.rtol,
                            atol=[1e-11, 1e-15, 1e-15, 1e-11], max_step=.08, events=events, dense_output=True)
            if not sol.success: raise RuntimeError(sol.message)
            self.sols.append((list(self.active), sol))
            samples = np.linspace(x, sol.t[-1], max(3, int((sol.t[-1] - x) * 35)))
            for xx in samples:
                ww = sol.sol(xx); st, _ = self.state(xx, ww); self.current_x = xx
                fac, el, rev = self.recycling_rates(st) if self.formation != 'gas' else (1., None, None)
                p = st['P'] / st['Q']; energy = (np.exp(3 * xx) * st['rho_tot'] / self.E0 + ww[3] - 1)
                trajectories.append(dict(x=float(xx), T_over_m=st['t'], up=st['up'], um=st['um'], pairs_per_net=p,
                                         condensed_fraction=st['n_c'] / st['Q'], charge_error=st['Q'] / (nq * np.exp(-3 * (xx - x0))) - 1,
                                         energy_work_error=energy, energy_per_net_m=st['rho_tot'] / st['Q'],
                                         visible_T_GeV=self.B.T_of_x(xx), n_s_per_net=st['n_s'] / st['Q'],
                                         Xbar_gas_per_net=st['n_minus'] / st['Q']))
            segments.append(dict(start_x=x, end_x=float(sol.t[-1]), active=start_active, nfev=sol.nfev))
            x = float(sol.t[-1]); w = sol.y[:, -1]
            if sol.t_events[0].size: stopped = f'T={self.tstop}m'; break
            if sol.t_events[1].size: stopped = 'pair_ratio=1e-18'; break
            if self.closure_stop and sol.t_events[4].size: stopped = 'kinetic_shape_rate/H=1'; break
            event = next((j for j in range(2) if sol.t_events[j + 2].size), None)
            if event is None: break
            self.active[event] = not self.active[event]
            w[event + 1] = 0. if self.active[event] else (self.c.cap_p if event == 0 else self.c.cap_m)
            switches.append(dict(x=x, species=event, new_active=self.active[event]))
        else:
            raise RuntimeError('Too many active-set transitions')
        final = trajectories[-1]
        return dict(R=self.R, formation=self.formation, x0=x0, initial=dict(t=initial['t'], up=initial['up'],
                    cond_fraction=initial['n_cond'] / nq, rho_over_mX4=rho, nq_over_mX3=nq),
                    final=final, stop=stopped, segments=segments, switches=switches,
                    checks=dict(max_charge_error=max(abs(r['charge_error']) for r in trajectories),
                                max_energy_work_error=max(abs(r['energy_work_error']) for r in trajectories)),
                    trajectories=trajectories[::4] + [final])

    def channel_budget(self, npts=60):
        """Per-channel |dP/dln a| / n_Q along the solution (which reactions hold P near equilibrium)."""
        rows = []
        for active, sol in self.sols:
            self.active = list(active)
            for xx in np.linspace(sol.t[0], sol.t[-1], max(3, int(npts * (sol.t[-1] - sol.t[0]) / 10) + 3)):
                ww = sol.sol(xx); st, _ = self.state(xx, ww); self.current_x = xx
                H = self.B.H_of_x(xx); R = self.helper.rates(st)
                ch = {}
                for name, (g, dp, dn) in R.items():
                    base = name[:-3] if name.endswith('_cc') else name
                    ch[base] = ch.get(base, 0.0) + g * dp / (H * K.MX ** 3) / st['Q']
                rows.append(dict(x=float(xx), T_over_m=st['t'], p=st['P'] / st['Q'], dlnP_rates={k: float(v / max(st['P'] / st['Q'], 1e-300)) for k, v in ch.items()}))
        return rows


CASES = {}


def add(name, variant, R, formation='gas', Q_L=None, start='relaxed', cubic=None, Rtag=None):
    CASES[name] = dict(variant=variant, R=R, formation=formation, Q_L=Q_L, start=start, cubic=cubic, Rtag=Rtag)


R497 = HO['R']['tot_e1e6_band']; R114 = HO['R']['med_e1e7']; R091 = HO['R']['med_e1e6']; R417 = HO['R']['tot_e1e7_band']
R570 = HO['R']['tot_e1e6_broad']
# reproduction (old object, old Y_Q, old V0, a_ref start)
add('repro_gas_R497', 'old', R497, 'gas', start='aref')
add('repro_gas_R5', 'old', 5.0, 'gas', start='aref')
# M1, a_ref start (only the M1 inputs changed)
for nm, R, tag in (('R497', R497, 'tot_e1e6_band_4.97'), ('R114', R114, 'med_e1e7_1.14'), ('R091', R091, 'med_e1e6'),
                   ('R570', R570, 'tot_e1e6_broad'), ('R417', R417, 'tot_e1e7_band_41.7'), ('R0', 0.0, 'R0')):
    add('m1_gas_aref_' + nm, 'm1', R, 'gas', start='aref', Rtag=tag)
    add('m1_gas_relaxed_' + nm, 'm1', R, 'gas', start='relaxed', Rtag=tag)
    add('m1_gas_spin_' + nm, 'm1', R, 'gas', start='spinodal', Rtag=tag)
    add('m1_inst_relaxed_' + nm, 'm1', R, 'instant', start='relaxed', Rtag=tag)
    add('m1_bind_relaxed_' + nm, 'm1', R, 'instant', Q_L=Q_BIND * 1.0000001, start='relaxed', Rtag=tag)
    add('m1_frag_relaxed_' + nm, 'm1', R, 'frag', Q_L=2e5, start='relaxed', Rtag=tag)
add('m1_gas_relaxed_R497_oldYQ', 'm1', R497, 'gas', start='relaxed', Rtag='tot_e1e6_band_4.97')
add('m1_gas_relaxed_R497_nocubic', 'm1', R497, 'gas', start='relaxed', Rtag='tot_e1e6_band_4.97')
add('m1_gas_relaxed_R114_nocubic', 'm1', R114, 'gas', start='relaxed', Rtag='med_e1e7_1.14')
add('m1_gas_aref_R497_nocubic', 'm1', R497, 'gas', start='aref', Rtag='tot_e1e6_band_4.97')
add('m1_gas_relaxed_R497_oldV0', 'm1', R497, 'gas', start='relaxed', Rtag='tot_e1e6_band_4.97')


def x0_for(c):
    if c['start'] == 'aref':
        return 0.0
    if c['Rtag'] == 'interp':                       # R scan: interpolate ln(a_on/a_ref) in R over the computed rows
        rows = sorted([(r['R'], math.log(r['onset_untilted']['a_over_aref'])) for r in HO['chain']
                       if r['label'].startswith('heavy_release') and r['phi_in_equilibrium']])
        return float(np.interp(c['R'], [r[0] for r in rows], [r[1] for r in rows]))
    row = [r for r in HO['chain'] if r['label'] == 'heavy_release:' + c['Rtag'] and r['phi_in_equilibrium']][0]
    if c['start'] == 'relaxed':
        return math.log(row['onset_untilted']['a_over_aref'])
    hist = [h for h in HO['history_M1'] if h['label'] == 'heavy_release:' + c['Rtag']][0]
    return math.log(row['onset_untilted']['a_over_aref'] * hist['crossings']['spinodal_233']['a_over_aon'])


def run_case(name, budget=False):
    if name not in CASES and name.startswith('m1_gas_relaxed_Rx'):
        add(name, 'm1', float(name[len('m1_gas_relaxed_Rx'):]), 'gas', start='relaxed', Rtag='interp')
    capmul = 1.0
    if name not in CASES and '_capx' in name:          # capture-efficiency sensitivity: <name>_capx<multiplier>
        base, mul = name.split('_capx'); capmul = float(mul); CASES[name] = dict(CASES[base])
    c = CASES[name]
    configure(c['variant'])
    if name.endswith('_oldYQ'): K.YQ = OLD['YQ']
    if name.endswith('_oldV0'): K.V0 = OLD['V0']
    S2.OBJ = K.Objects()
    SCALE['capture'] = CAPTURE_SCALE['m1' if c['variant'] == 'm1' else 'old'] * capmul
    SCALE['cubic'] = c['variant'] == 'm1' and not name.endswith('_nocubic')
    x0 = x0_for(c)
    t0 = time.time()
    ctl = Control(c['R'], formation=c['formation'], Q_L=c['Q_L'], x0=x0,
                  rtol=float(os.environ.get('FX_RTOL', '2e-8')) if c['formation'] == 'gas' else float(os.environ.get('FX_RTOL_OBJ', '2e-6')))
    o = ctl.run()
    if budget:
        o['channel_budget'] = ctl.channel_budget()
    p = o['final']['pairs_per_net']
    o.update(name=name, case=c, x0=x0, a0_over_aref=math.exp(x0), p=p, f_old_convention=2.02855 * p,
             frozen_test_term=2.13356 * p, eta_min=1 - (0.01 - 2.13356 * p) / 0.0667805,
             YQ=K.YQ, V0_GeV4=K.V0, omega_cap=ctl.c.cap_p, e_cap=ctl.c.e_cap, capture_scale=SCALE['capture'],
             cubic_applied=SCALE['cubic'], rtol=ctl.rtol,
             seconds=time.time() - t0)
    return o


if __name__ == '__main__':
    names = sys.argv[1:]
    budget = '--budget' in names
    names = [n for n in names if not n.startswith('--')]
    if names == ['all']:
        names = list(CASES)
    for n in names:
        try:
            o = run_case(n, budget)
        except Exception as e:
            print('FAIL', n, repr(e), flush=True)
            json.dump(dict(name=n, error=repr(e)), open(os.path.join(SP, 'wt2/fx_chain/kin_%s.json' % n), 'w'))
            continue
        json.dump(o, open(os.path.join(SP, 'wt2/fx_chain/kin_%s.json' % n), 'w'), default=float)
        f = o['final']
        print('%-30s x0=%+.3f start t=%.3f cond=%.3f | stop=%s T/m=%.4f p=%.4e f=%.4e test=%.4e eta>=%.4f | dQ=%.1e dE=%.1e %.0fs'
              % (n, o['x0'], o['initial']['t'], o['initial']['cond_fraction'], o['stop'], f['T_over_m'], o['p'],
                 o['f_old_convention'], o['frozen_test_term'], o['eta_min'], o['checks']['max_charge_error'],
                 o['checks']['max_energy_work_error'], o['seconds']), flush=True)
