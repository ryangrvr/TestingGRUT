# Continued by Codex from Claude partial code at e71730c24fe6a0b28ddd72cd30982c5b39c9aedf.
# Modifications and limitations are recorded in this packet.
#!/usr/bin/env python3
"""kin_solve (Claude, own code): fixed inputs, interface loaders, visible background and the rate-table
interpolators the kinetic solver uses.

Inputs are read from the three builder interfaces of this review (tracks/kin_source, kin_rates, kin_thermo),
whose provenance is Codex 225028b INPUTS.json and 9877518 (read by those tracks via git show).  Nothing here
re-derives a kernel; this module turns the tabulated kernels into callables with stated interpolation and
extrapolation rules.  Units: GeV unless a name says otherwise; hidden thermodynamics in m_X = 1 (t = T/m_X).

Reuse of our earlier tracks, stated: the coarse SM g*(T) table is the one of tracks/ls_threshold/code/lsrelic.py
and kin_thermo/code/kt_common.py (copied, [Q own]); kin_thermo/code/thermo.py is imported by thermo_state.py for
the Bose state functions (its series was checked there against an independent quadrature to 2e-15).
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import kve
from scipy.interpolate import PchipInterpolator

HERE = Path(__file__).resolve().parent
from bootstrap import ROOT,WORK,INPUT_ROOT
TRACK = ROOT
TRACKS = INPUT_ROOT / 'tracks'
RES = WORK / 'results'
SCRATCH = WORK
CACHE = SCRATCH / 'cache'
for p in (RES, CACHE):
    p.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(TRACKS / 'kin_thermo' / 'code'))     # thermo.py (Bose state functions), kt_common constants

# ------------------------------------------------------------------ interfaces (read once)
SRC = json.loads((TRACKS / 'kin_source' / 'results' / 'source_interface.json').read_text())
RAT = json.loads((TRACKS / 'kin_rates' / 'results' / 'rates_interface.json').read_text())
THM = json.loads((TRACKS / 'kin_thermo' / 'results' / 'thermo_interface.json').read_text())

# ------------------------------------------------------------------ fixed point [Q INPUTS.json 225028b via kin_source]
FP = SRC['meta']['fixed_point']
MX = FP['mX_GeV']; MS = FP['ms_GeV']; RS = MS / MX
GX, LX, LP, V = FP['gX'], FP['lambda_X'], FP['lambda_phi'], FP['v_GeV']
YQ, DM_PER_S, F_ALLOW = FP['YQ'], FP['DM_per_entropy_GeV'], FP['allowance_f']
OMEGA_REF, E_OVER_QM_REF, Q_REF = FP['omega_over_m'], FP['E_over_Qm'], FP['Q_ref']
N_X, E_X, EPS = abs(FP['N_X']), FP['E_X'], abs(FP['eps'])
H_I = FP['H_GeV']
MRED = 2.435e18
LAM3 = RAT['meta']['parameters']['lam3_sXXb_GeV']
GAMMA_S_VAC = RAT['meta']['parameters']['Gamma_s_vac_GeV']

ONSET = SRC['release_population']['onset']
N_Q0 = ONSET['n_Q_GeV3']; RHO_C = ONSET['rho_c_GeV4']; V0 = ONSET['V0_GeV4']; N_C = ONSET['n_c_GeV3']
A_REF = ONSET['a_ref_from_charge']
E_REF_GEV = THM['objects']['reference_Q200000']['E_GeV']
E_MINUS_QOMEGA = THM['objects']['reference_Q200000']['E_minus_Q_omega_GeV']
Q_MIN = THM['objects']['Q_min']['Q']

# source energy ratios [C kin_source, CONFIRMED against Codex]
R_VEC = {k: SRC['vector_sources'][k]['per_sign_totals']['R'] for k in SRC['vector_sources']}
R_MED = {'e1e-7': SRC['mediator_debris']['totals']['e1e-7']['R_mediator_total'],
         'e1e-6': SRC['mediator_debris']['totals']['e1e-6']['R_mediator_total']}
PAIRS_PER_NET = {k: SRC['vector_sources'][k]['per_sign_totals']['pairs_per_net_charge'] for k in SRC['vector_sources']}
PAIRS_PER_NET_MED = SRC['mediator_debris']['totals']['e1e-7']['pairs_per_net_charge_total']


def rho_onset(R, with_V0=True):
    """Hidden energy density at a_ref for total debris ratio R (release + sources [+ V0])."""
    return RHO_C * (1.0 + R) + (V0 if with_V0 else 0.0)


def f_old(pairs_per_net):
    """The old convention: f = 2 m_X Y_Q pairs/(DM per entropy)."""
    return 2 * MX * YQ * pairs_per_net / DM_PER_S


def f_bench(pairs_per_net, e_over_qm=E_OVER_QM_REF):
    """Benchmark convention (Codex clarification 7): symmetric mass 2 m_X per pair over the net-charge object mass
    Y_Q E(Q_ref)/Q_ref, the same denominator in every run."""
    return 2 * pairs_per_net / e_over_qm


# ------------------------------------------------------------------ visible background
_TS = np.log([1e-6, 5e-4, 1e-3, .1, .15, .3, 1., 5., 100., 200., 1e19])
_G = [3.36, 3.36, 10.75, 10.75, 17.25, 61.75, 72.25, 86.25, 96.25, 106.75, 106.75]
_GS = [3.91, 3.91, 10.75, 10.75, 17.25, 61.75, 72.25, 86.25, 96.25, 106.75, 106.75]


class Background:
    """Visible radiation: s_vis(T) = 2pi^2/45 g_s T^3, H = sqrt(pi^2 g_rho/90) T^2/M_Pl.  a_ref is where
    Y_Q s_vis = n_Q0.  mode 'table': the coarse SM table (lsrelic/kin_thermo [Q own]); 'const': g* = 106.75."""

    def __init__(self, mode='table'):
        self.mode = mode
        # The instant-object diagnostic can outlive the old 10-keV table
        # floor. Never silently clamp its Hubble rate at that floor.
        lT = np.linspace(np.log(1e-12), np.log(1e6), 8001)
        T = np.exp(lT)
        gs = np.array([self.g_s(x) for x in T])
        lna = -np.log(gs ** (1 / 3) * T)                   # ln a + const (entropy conservation)
        # monotone decreasing in T; build the inverse ln a -> ln T
        order = np.argsort(lna)
        self._lna_grid = lna[order]; self._lT_grid = lT[order]
        self.T_ref = self.T_of_charge(N_Q0)
        self._lna_ref = -np.log(self.g_s(self.T_ref) ** (1 / 3) * self.T_ref)
        self.H_ref = self.hubble(self.T_ref)

    def g_rho(self, T):
        return float(np.interp(np.log(T), _TS, _G)) if self.mode == 'table' else 106.75

    def g_s(self, T):
        return float(np.interp(np.log(T), _TS, _GS)) if self.mode == 'table' else 106.75

    def entropy(self, T):
        return 2 * np.pi ** 2 / 45 * self.g_s(T) * T ** 3

    def hubble(self, T):
        return np.sqrt(np.pi ** 2 * self.g_rho(T) / 90) * T * T / MRED

    def T_of_charge(self, nQ):
        from scipy.optimize import brentq
        return brentq(lambda T: YQ * self.entropy(T) - nQ, 1e-6, 1e6, xtol=1e-14, rtol=1e-14)

    def T_of_x(self, x):
        """Visible temperature at ln(a/a_ref) = x."""
        where=self._lna_ref+x
        if not self._lna_grid[0]<=where<=self._lna_grid[-1]:
            raise ValueError('Visible temperature outside background grid')
        return float(np.exp(np.interp(where, self._lna_grid, self._lT_grid)))

    def H_of_x(self, x):
        return self.hubble(self.T_of_x(x))

    def nQ_of_x(self, x):
        return YQ * self.entropy(self.T_of_x(x))


# ------------------------------------------------------------------ kin_rates kernels as callables (m_X units for T)
def _loglog(T, y, name):
    T = np.asarray(T, float); y = np.asarray(y, float)
    ok = y > 0
    if ok.sum() < 2:
        raise ValueError('not enough positive points for %s' % name)
    return PchipInterpolator(np.log(T[ok]), np.log(y[ok]), extrapolate=False), float(T[ok].min()), float(T[ok].max())


def n_eq_mb(m_over_m, t):
    """MB equilibrium density at mu = 0 for mass r (m_X units), per bosonic degree, in m_X^3."""
    r = m_over_m
    return r * r * t * kve(2, r / t) * np.exp(-r / t) / (2 * np.pi ** 2)

def log_n_eq_mb(r,t):
    return 2*np.log(r)+np.log(t)+np.log(kve(2,r/t))-r/t-np.log(2*np.pi**2)


class Rates:
    """MB rate densities gamma(T) (GeV^4) at mu = 0 from kin_rates, with the Bose/MB factors, as functions of
    t = T/m_X.  Chemical potentials and detailed balance are applied by the solver (exact for MB)."""

    def __init__(self):
        c22 = RAT['channels_22']['ss_XXb']['thermal_MB']
        self.T22 = np.array(c22['T_over_mX'])
        self.g_ss = np.array(c22['gamma_sub'])                        # ss -> X Xb, no pole: sub = nores
        self._ss, self._ss_lo, self._ss_hi = _loglog(self.T22, self.g_ss, 'ss_XXb')
        b = RAT['channels_22']['ss_XXb']['bose_over_MB_mu0']
        self._bT = np.array(sorted(float(k) for k in b)); self._bss = np.array([b[str(k) if str(k) in b else k]['bose_over_MB'] for k in sorted(b, key=float)])
        bd = RAT['decay_s_XXb']['bose_over_MB_mu0']
        self._bdec = np.array([bd[k]['decay_BE_over_MB'] for k in sorted(bd, key=float)])
        bn = RAT['decay_s_XXb']['bose_factors']['net_width_bose_over_MB_mu0']
        self._bnet = np.array([bn[k] for k in sorted(bn, key=float)])
        el = RAT['channels_22']
        self.elastic = {}
        for ch, col in (('XXb_XXb', 'sigmav_nores'), ('XX_XX', 'sigmav_sub'), ('Xs_Xs', 'sigmav_nores'), ('ss_ss', 'sigmav_sub')):
            th = el[ch]['thermal_MB']
            self.elastic[ch] = _loglog(th['T_over_mX'], th[col], ch)
        # n <-> 2 subtracted remainders (P-changing channels); legs (incoming masses r_i) for the low-T extrapolation
        self.n2 = {}
        legs_in = {'A_XXXbXb_XXb': [1, 1, 1, 1], 'B_XXXXb_XX': [1, 1, 1, 1], 'E_XXbs_XXb': [1, 1, RS], 'E2_XXs_XX': [1, 1, RS], 'F_sss_ss': [RS] * 3,
                   'G2_Xss_Xs': [1, RS, RS], 'L_sss_XXb': [RS] * 3, 'N_ssss_ss': [RS] * 4, 'H_XXbss_XXb': [1, 1, RS, RS],
                   'J_XXXbs_Xs': [1, 1, 1, RS], 'K_XXbss_ss': [1, 1, RS, RS]}
        for ch, legs in legs_in.items():
            th = RAT['channels_n2'][ch]['thermal_MB']
            T = np.array(th['T_over_mX']); g = np.array(th['gamma_sub']); res = np.array(th['remainder_resolved'], bool)
            bound = np.array(th['remainder_bound']); err = np.array(th['gamma_sub_err'])
            use = res & (g > 0)
            entry = dict(T=T, g=g, resolved=res, bound=bound, err=err, legs=legs, Sin=RAT['channels_n2'][ch]['Sin'])
            if use.sum() >= 2:
                entry['interp'], entry['lo'], entry['hi'] = _loglog(T[use], g[use], ch)
                entry['Tmax_resolved'] = float(T[use].max())
            else:
                entry['interp'] = None
            self.n2[ch] = entry
        nr = RAT['checks']['nr_limits']
        self.A_old = nr['A_XXXbXb_XXb']['A_old']; self.B_old = nr['B_XXXXb_XX']['A_old']
        # decay: closed form checked against the table
        dec = RAT['decay_s_XXb']['MB']
        self._dec_T = np.array(dec['T_over_mX']); self._dec_g = np.array(dec['gamma_decay_MB'])
        # in-medium
        mt = RAT['masses_T']['oneloop']
        self.T_close = mt['T_close_s_decay_over_mX']; self.T_c = mt['T_c_over_mX']
        self._mT = np.array(mt['T_over_mX']); self._mX_T = np.array(mt['mX_GeV']) / MX; self._ms_T = np.array(mt['ms_GeV']) / MS

    # --- Bose factors (mu = 0 tables at 5 T; constant outside; 1 below 0.3 m where they are 1.00-1.05)
    def _bose(self, arr, t):
        lt = np.clip(np.log(t), np.log(self._bT[0]), np.log(self._bT[-1]))
        return float(np.exp(np.interp(lt, np.log(self._bT), np.log(arr))))

    def bose_ss(self, t): return self._bose(self._bss, t)
    def bose_dec_gross(self, t): return self._bose(self._bdec, t)
    def bose_dec_net(self, t): return self._bose(self._bnet, t)

    # --- ss -> X Xb MB rate density at mu = 0, GeV^4
    def gamma_ss(self, t):
        if t < self._ss_lo:            # threshold law: sigma v -> SV_SS const, gamma ∝ n_s^eq(t)^2
            return float(np.exp(self._ss(np.log(self._ss_lo)))) * (n_eq_mb(RS, t) / n_eq_mb(RS, self._ss_lo)) ** 2
        if t > self._ss_hi:
            return float(np.exp(self._ss(np.log(self._ss_hi))))
        return float(np.exp(self._ss(np.log(t))))

    def log_gamma_ss(self,t):
        if t<self._ss_lo:
            return float(self._ss(np.log(self._ss_lo)))+2*(log_n_eq_mb(RS,t)-log_n_eq_mb(RS,self._ss_lo))
        return float(self._ss(np.log(min(t,self._ss_hi))))

    def log_gamma_n2(self,ch,t):
        e=self.n2[ch]
        if e['interp'] is None or t>e['hi']:return None
        if t<e['lo']:
            return float(e['interp'](np.log(e['lo'])))+sum(log_n_eq_mb(r,t)-log_n_eq_mb(r,e['lo']) for r in e['legs'])
        v=float(e['interp'](np.log(t)));return v if np.isfinite(v) else None

    # --- n -> 2 subtracted remainders where resolved (0 where not; bound separately), GeV^4 at mu = 0
    def gamma_n2(self, ch, t):
        e = self.n2[ch]
        if e['interp'] is None:
            return 0.0
        if t < e['lo']:
            scale = 1.0
            for r in e['legs']:
                scale *= n_eq_mb(r, t) / n_eq_mb(r, e['lo'])
            return float(np.exp(e['interp'](np.log(e['lo'])))) * scale
        if t > e['hi']:
            return 0.0                     # above the last resolved point the remainder is not resolved -> 0 (bound carried)
        v = e['interp'](np.log(t))
        return float(np.exp(v)) if np.isfinite(v) else 0.0

    def bound_n2(self, ch, t):
        e = self.n2[ch]
        return float(np.interp(np.log(t), np.log(e['T']), e['bound']))

    def A_eff_ratio(self, ch, t):
        """A_eff(T)/A_old = S_in gamma_sub/(n_eq^4 A_old): the energy-dependence factor of the 4->2 kernels."""
        e = self.n2[ch]; A_old = self.A_old if ch.startswith('A') else self.B_old
        g = self.gamma_n2(ch, t)
        n4 = (n_eq_mb(1.0, t) * MX ** 3) ** 4
        return e['Sin'] * g / (n4 * A_old) if n4 > 0 else 0.0

    # --- s decay, MB (closed form [D], checked against the table)
    def gamma_dec_mb(self, t):
        T = t * MX
        return GAMMA_S_VAC * MS * MS * T * kve(1, MS / T) * np.exp(-MS / T) / (2 * np.pi ** 2)

    def mean_m_over_E(self, t):
        return float(kve(1, RS / t) / kve(2, RS / t))

    def masses_T(self, t):
        """In-medium m_X(T)/m_X, m_s(T)/m_s (one loop, kin_rates); vacuum below the table."""
        if t <= self._mT[0]:
            return 1.0, 1.0
        return float(np.interp(t, self._mT, self._mX_T)), float(np.interp(t, self._mT, self._ms_T))


# ------------------------------------------------------------------ kin_thermo object tables as callables
class Objects:
    def __init__(self):
        rows = [r for r in THM['objects']['table'] if r['branch'].startswith('stable')]
        Q = np.array([r['Q'] for r in rows]); o = np.argsort(Q)
        self.Q = Q[o]
        self.Om = np.array([r['Omega'] for r in rows])[o]
        self.EQ = np.array([r['E_over_QmX'] for r in rows])[o]
        self.EmQw = np.array([r['E_minus_Q_omega_GeV'] for r in rows])[o]
        self.r90 = np.array([r['r90_over_mX'] for r in rows])[o]

    def omega(self, Q):
        return float(np.interp(np.log(Q), np.log(self.Q), self.Om))

    def e_over_q(self, Q):
        return float(np.interp(np.log(Q), np.log(self.Q), self.EQ))

    def r90_over_m(self, Q):
        return float(np.interp(np.log(Q), np.log(self.Q), self.r90))


def conversion_table(extrap='hold', n=120):
    """MB thermal averages of the capture kernel: <sigma v>_f(t) (GeV^-2) over an MB Xbar distribution at
    temperature t (m_X units), and <E_in>(t)/m_X; cached.  Built from kin_thermo kernels.sv_forward (the dense
    sigma v(p) table of our conversion track, read as data)."""
    f = CACHE / ('conv_%s.json' % extrap)
    if f.exists():
        d = json.loads(f.read_text())
        return np.array(d['t']), np.array(d['sv']), np.array(d['Ein'])
    # Reuse the completed interface's certified cold forward table directly;
    # do not import the unavailable old conversion checkout. High-p choices
    # remain labelled recipes, not certified upper/lower bounds.
    raw=THM['kernels']['conversion_forward']
    ps=np.array(raw['p_minus_over_mX']);vs=np.array(raw['sigma_v_GeV2'])
    mask=ps<=.70800001;ps=ps[mask];vs=vs[mask]
    order=np.argsort(ps);ps=ps[order];vs=vs[order]
    interp=PchipInterpolator(np.log(ps),np.log(vs))
    radius=THM['objects']['reference_Q200000']['r90_over_mX']/MX
    def sv_forward(p,extrap):
        if p>ps[-1] and extrap=='geometric':return np.pi*radius**2*p/np.sqrt(1+p*p)
        return float(np.exp(interp(np.log(np.clip(p,ps[0],ps[-1])))))
    from scipy.integrate import quad
    ts = np.geomspace(1e-3, 30.0, n)
    sv = np.zeros(n); Ein = np.zeros(n)
    for i, t in enumerate(ts):
        pts = [0.0, 0.05, 0.2, 0.5, 0.708, 1.5, 3.0, 8.0, 30.0, 150.0, 600.0]
        pts = [p for p in pts if p < 60 * max(t, 0.05) + 1.0] + [60 * max(t, 0.05) + 1.0]
        w = lambda p: p * p * np.exp(-(np.sqrt(1 + p * p) - 1) / t)
        den = sum(quad(w, a, b, epsabs=0, epsrel=1e-9, limit=200)[0] for a, b in zip(pts[:-1], pts[1:]))
        num = sum(quad(lambda p: w(p) * sv_forward(p, extrap), a, b, epsabs=0, epsrel=1e-9, limit=200)[0] for a, b in zip(pts[:-1], pts[1:]))
        ein = sum(quad(lambda p: w(p) * sv_forward(p, extrap) * np.sqrt(1 + p * p), a, b, epsabs=0, epsrel=1e-9, limit=200)[0] for a, b in zip(pts[:-1], pts[1:]))
        sv[i] = num / den; Ein[i] = ein / num
    f.write_text(json.dumps(dict(t=ts.tolist(), sv=sv.tolist(), Ein=Ein.tolist())))
    return ts, sv, Ein


def saturated_table(n=241):
    """Bose state functions of one X degree at the saturation cap u = 1 - 1e-7 (direct quadrature, kin_thermo
    thermo.quad_species), tabulated in t for the gas_only / fragmentation rules; cached."""
    f = CACHE / 'saturated.json'
    if f.exists():
        d = json.loads(f.read_text())
        return {k: np.array(v) for k, v in d.items()}
    import thermo as TH
    ts = np.geomspace(1e-4, 40.0, n)
    out = dict(t=ts, n=np.zeros(n), P=np.zeros(n), rho=np.zeros(n), dn_dt=np.zeros(n), drho_dt=np.zeros(n))
    for i, t in enumerate(ts):
        q = TH.quad_species(t, 1.0, 1.0 - 1e-7, derivs=True)
        for k in ('n', 'P', 'rho', 'dn_dt', 'drho_dt'):
            out[k][i] = q[k]
    f.write_text(json.dumps({k: v.tolist() for k, v in out.items()}))
    return out


if __name__ == '__main__':
    B = Background('table'); Bc = Background('const')
    print('T_ref table %.4f GeV (const %.4f); H_ref %.4e (%.4e) GeV' % (B.T_ref, Bc.T_ref, B.H_ref, Bc.H_ref))
    print('kin_source a_ref', A_REF, 'T_vis_ref', ONSET['T_vis_ref_GeV'], 'kin_thermo H_ref', THM['expansion']['background']['H_ref_GeV'])
    for x in (0, 1, 2, 3, 5, 7, 9):
        print('x=%d T_vis=%.4e H=%.3e nQ=%.3e' % (x, B.T_of_x(x), B.H_of_x(x), B.nQ_of_x(x)))
    R = Rates()
    for t in (0.1, 0.3, 1.0, 3.0, 7.344, 12.0):
        print('t=%6.3f gamma_ss=%.3e (bose %.2f) gamma_dec closed=%.3e table=%.3e  F=%.2e L=%.2e A=%.2e B=%.2e E2=%.2e G2=%.2e' % (
            t, R.gamma_ss(t), R.bose_ss(t), R.gamma_dec_mb(t), float(np.exp(np.interp(np.log(t), np.log(R._dec_T), np.log(R._dec_g)))),
            R.gamma_n2('F_sss_ss', t), R.gamma_n2('L_sss_XXb', t), R.gamma_n2('A_XXXbXb_XXb', t), R.gamma_n2('B_XXXXb_XX', t),
            R.gamma_n2('E2_XXs_XX', t), R.gamma_n2('G2_Xss_Xs', t)))
    for t in (0.03, 0.1, 0.3, 1.0):
        print('A_eff/A_old at t=%.2f: %.3f  B: %.3f' % (t, R.A_eff_ratio('A_XXXbXb_XXb', t), R.A_eff_ratio('B_XXXXb_XX', t)))
    O = Objects()
    print('omega(Q_ref)=%.10f E/Q=%.10f  omega(1e6)=%.4f E/Q(1e6)=%.4f  Q_min=%.0f' % (O.omega(Q_REF), O.e_over_q(Q_REF), O.omega(1e6), O.e_over_q(1e6), Q_MIN))
    import time
    t0 = time.time(); ts, sv, Ein = conversion_table('hold'); print('conversion table hold built %.1fs' % (time.time() - t0))
    for t in (0.03, 0.1, 0.3, 1.0, 7.34):
        print('  t=%.2f <sv>_f(MB)=%.1f GeV^-2  <E_in>=%.3f' % (t, np.interp(np.log(t), np.log(ts), sv), np.interp(np.log(t), np.log(ts), Ein)))
    t0 = time.time(); ts, sv, Ein = conversion_table('geometric'); print('conversion table geometric built %.1fs' % (time.time() - t0))
    t0 = time.time(); S = saturated_table(); print('saturated table built %.1fs' % (time.time() - t0))
    for t in (0.05, 0.1, 0.3, 1.0):
        print('  t=%.2f n_sat=%.6e rho_sat=%.6e' % (t, np.interp(np.log(t), np.log(S['t']), S['n']), np.interp(np.log(t), np.log(S['t']), S['rho'])))
