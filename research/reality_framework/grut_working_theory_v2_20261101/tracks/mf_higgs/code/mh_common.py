#!/usr/bin/env python3
"""Track mf_higgs (Claude, multi-field inflaton pass, BRIEF_MULTIFIELD.md; GitHub issue #13 context).  Own code.

Higgs survival with xi_H = 0.1 on mf_hybrid's benchmarks HFa / HFb.  Every choice here is a model choice (the benchmark,
xi_H, the survival criterion, the plasma model, the fluid closure); none follows from closure.

Common module:
  * SM input: tracks/higgs_reheat/code/smpot.py through tracks/inf_higgs/code/ih_common.SMTab (read-only imports;
    the SM tables are cached in this track's scratchpad, IH_SCRATCH).
  * The benchmark histories of tracks/mf_hybrid/results/bench_*.json (read-only).
  * A post-waterfall fluid (own implementation of the mf_hybrid skeptic's Pauli-blocked fluid, sk_pauli.py, with one
    addition: the SM NUMBER density, so that the hard-gas induced Higgs mass m_ind^2 = c rho_SM / E_SM^2 can be formed).
    The production tables (tracks/mf_hybrid/results/axial_tables.json) and the Pauli saturation fractions
    (tracks/mf_hybrid/skeptic/results/sk_axial*.json) are read, not recomputed.
  * Background tables for the Higgs modes in H_i units (a = 1, H = H_i at sigma = sigma_c, mf_hybrid's slice):
      EXP     : the exported bench history (exact homogeneous stage to t = 0.6 with raw R; fluid mean after, plus the
                R oscillations 9 H^2 f_c cos(omega_c t) rebuilt from the exported fractions).  This is the unblocked
                (K ~ 2) corner the mf_hybrid skeptic flagged.
      UNBL    : bench to t = 0.19, then this track's fluid without Pauli blocking (control; reproduces K ~ 2).
      PAULI0/1/10, PAULIb : bench to t = 0.19, then the Pauli-blocked fluid with N-thermalisation rate
                Gamma_th = 0 / 1 / 10 H_i (perturbative tails free), or with the tails blocked too.
Units: fields in reduced Planck units, time in 1/H_i, R in H_i^2, rho in H_i^4 (rho = 3 mP^2 H^2, mP = M_P/H_ref).
"""
import sys
sys.dont_write_bytecode = True
import json
import math
import os
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
TRACK = HERE.parent
RES = TRACK / 'results'
TRACKS = TRACK.parent
SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/mf_higgs')
SP.mkdir(parents=True, exist_ok=True)
os.environ['IH_SCRATCH'] = str(SP)
sys.path.insert(0, str(TRACKS / 'inf_higgs' / 'code'))
import ih_common as IHC          # noqa: E402  read-only reuse (SMTab; smpot underneath)
S = IHC.S                        # higgs_reheat smpot (read-only)

MPL = 2.435e18
H_REF = 1.30e11                  # mf_hybrid design H (bench H_i = 1.30002e11)
mP_REF = MPL / H_REF
PDG24 = dict(Mt=172.57, sMt=0.29, Mh=125.20, sMh=0.11, als=0.1180, sals=0.0009)   # PDG 2024 [R]
H0C_LN = math.log(0.05 / 2.248e-4)        # ln(k_*/a0 H0), h = 0.674 [R]
MNU, VEW = 0.05e-9, 174.1                 # seesaw m-tilde (mf_hybrid's choice) and v_EW (GeV)
C_HARD = 0.11 * 2.70 ** 2 / (math.pi ** 2 * 106.75 / 30)   # 0.0228: inf_higgs/inf_reheat hard-gas coefficient [E]
C_HARD_SK = 0.032                         # inf_higgs skeptic N11 composition value (0.031-0.034) [E]
MB = TRACKS / 'mf_hybrid'


def jdump(obj, path):
    def conv(z):
        if isinstance(z, (np.floating, np.integer)):
            return float(z)
        if isinstance(z, np.ndarray):
            return z.tolist()
        return str(z)
    json.dump(obj, open(path, 'w'), indent=1, default=conv)


# --------------------------------------------------------------------------------------------- bench
_B = {}


def load_bench(name):
    if name not in _B:
        d = json.load(open(MB / 'results' / f'bench_{name}.json'))
        h = {k: np.array([np.nan if x is None else x for x in v], dtype=float) for k, v in d['history'].items()}
        d['h'] = h
        i = int(np.argmin(np.abs(h['t'] - d['waterfall']['t_eps1'])))
        d['t_eps1'] = float(d['waterfall']['t_eps1'])
        _B[name] = d
    return _B[name]


def Nstar_K(name, K):
    """N_* (pivot exit -> slice) for a radiation era with K (mf_hybrid's bookkeeping: N_* = const - (1/2) ln K)."""
    d = load_bench(name)
    return d['cmb']['Nstar'] - 0.5 * math.log(K / d['reheating']['K'])


def N_obs(name, K):
    """e-folds from the exit of today's horizon to the slice (a = 1 at sigma = sigma_c)."""
    return Nstar_K(name, K) + H0C_LN


# --------------------------------------------------------------------------------------------- fluid
TAB = json.load(open(MB / 'results' / 'axial_tables.json'))
FS = json.load(open(MB / 'skeptic' / 'results' / 'sk_axial_Fsat_vs_A.json'))
AXM = json.load(open(MB / 'skeptic' / 'results' / 'sk_axial.json'))['multi_period']


class Chan:
    """axial channel condensate (mass m) -> N (M, f): non-perturbative rate (m/pi) eps(A)/(S^2/2) from the track's
    exact half-period table, perturbative tail (M/f)^2 m beta/16pi, Pauli saturation F_sat(A) from the skeptic's
    multi-period runs (interpolated in ln A, slope 2 outside), seesaw rest width m_nu M^2/(8 pi v^2)."""

    def __init__(self, tag, m, M, f, grest_scale=1.0):
        self.tag, self.m, self.M, self.f = tag, m, M, f
        key = min(TAB.keys(), key=lambda s: abs(float(s) - M / m))
        self.lA = np.log(TAB[key]['A'])
        self.le = np.log(np.maximum(TAB[key]['eps'], 1e-300))
        pts = [(v['A'], v['Fsat']) for k, v in FS.items() if k.startswith(tag + '_')]
        akey = {'N1': 'sigma->N1', 'N2': 'chi->N2', 'N3': 'chi->N3', 'N4': 'chi->N4'}[tag]
        rows = AXM[akey]['rows']
        pts.append((AXM[akey]['A'], float(np.mean([r['F_tot'] for r in rows[6:]]))))
        pts.sort()
        self.flA = np.log([p[0] for p in pts])
        self.flF = np.log([p[1] for p in pts])
        b2 = 1 - 4 * M * M / (m * m)
        self.Gp = (M / f) ** 2 * m * math.sqrt(b2) / (16 * math.pi) if b2 > 0 else 0.0
        self.Grest = grest_scale * MNU * (M * H_REF) ** 2 / (8 * math.pi * VEW ** 2) / H_REF

    def amp(self, rho):
        Sig = math.sqrt(2 * max(rho, 0.0)) / self.m
        return Sig / (2 * self.f), Sig / self.m

    def rate(self, rho):
        if rho <= 0:
            return 0.0
        A, Sm = self.amp(rho)
        lA = math.log(A)
        if lA < self.lA[0]:
            return self.Gp
        le = float(np.interp(lA, self.lA, self.le)) if lA <= self.lA[-1] else self.le[-1] + 4 * (lA - self.lA[-1])
        return (self.m / math.pi) * math.exp(le) / (0.5 * Sm * Sm)

    def fsat(self, rho):
        A, _ = self.amp(rho)
        lA = math.log(max(A, 1e-30))
        if lA < self.flA[0]:
            lF = self.flF[0] + 2 * (lA - self.flA[0])
        elif lA > self.flA[-1]:
            lF = self.flF[-1] + 2 * (lA - self.flA[-1])
        else:
            lF = float(np.interp(lA, self.flA, self.flF))
        return min(math.exp(lF), 1.0)

    def mu5(self, rho):
        return math.sqrt(2 * max(rho, 0.0)) / (2 * self.f)


def bench_state(name, t0=0.19):
    d = load_bench(name)
    h = d['h']
    i = int(np.argmin(np.abs(h['t'] - t0)))
    return dict(t=float(h['t'][i]), a=float(h['a'][i]), H=float(h['H'][i]), fs=float(h['f_sig'][i]),
                fc=float(h['f_chi'][i]), fr=float(h['f_rad'][i]))


def fluid(name='HFa', blocking=True, Gth=0.0, pert_free=True, grest_scale=1.0, t_end=3e4, start=None, n_out=1500):
    """Post-waterfall fluid from the bench state at t ~ 0.19 (the skeptic's start).  Returns arrays (H_i units).
    Comoving state: ln a, rho_sig a^3, rho_chi a^3, E_j a^4, n_j a^3, B_j a^4 (j = N1..N4), rho_SM a^4, n_SM a^3."""
    mP = mP_REF
    d = load_bench(name)
    ax = d['params']
    cs = [Chan('N1', ax['m_sig_over_Hi'], ax['axial_sigma'][0]['M_over_Hi'], ax['axial_sigma'][0]['f_over_Hi'], grest_scale)]
    cc = [Chan(t, ax['m_chi_over_Hi'], c['M_over_Hi'], c['f_over_Hi'], grest_scale)
          for t, c in zip(('N2', 'N3', 'N4'), ax['axial_chi'])]
    allc = [(c, 0) for c in cs] + [(c, 1) for c in cc]
    nC = len(allc)
    st = start or bench_state(name)
    rho0 = 3 * mP ** 2 * st['H'] ** 2
    a0 = st['a']

    def unpack(y):
        a = math.exp(y[0])
        rs, rc = max(y[1], 0) / a ** 3, max(y[2], 0) / a ** 3
        E = [max(y[3 + j], 0) / a ** 4 for j in range(nC)]
        n = [max(y[3 + nC + j], 1e-300) / a ** 3 for j in range(nC)]
        B = [max(y[3 + 2 * nC + j], 0) / a ** 4 for j in range(nC)]
        rsm = max(y[3 + 3 * nC], 0) / a ** 4
        nsm = max(y[4 + 3 * nC], 0) / a ** 3
        return a, rs, rc, E, n, B, rsm, nsm

    def rhs(t, y):
        a, rs, rc, E, n, B, rsm, nsm = unpack(y)
        rho = rs + rc + sum(E) + rsm
        H = math.sqrt(rho / (3 * mP ** 2))
        P, Gl = [], []
        for j, (c, w) in enumerate(allc):
            rcond = rs if w == 0 else rc
            G = c.rate(rcond)
            if blocking:
                Es = c.fsat(rcond) * rcond
                bl = max(0.0, 1.0 - B[j] / Es) if Es > 0 else 0.0
            else:
                bl = 1.0
            Pnp = G * rcond * bl
            Pp = c.Gp * rcond * (1.0 if (pert_free or not blocking) else bl)
            P.append(max(Pnp, Pp))
            Emean = E[j] / n[j] if n[j] > 0 else c.M
            Gl.append(c.Grest * c.M / max(c.M, Emean))
        Gs = sum(P[j] for j, (c, w) in enumerate(allc) if w == 0) / max(rs, 1e-300)
        Gc = sum(P[j] for j, (c, w) in enumerate(allc) if w == 1) / max(rc, 1e-300)
        dy = [H, -Gs * rs * a ** 3, -Gc * rc * a ** 3]
        dE, dn, dB = [], [], []
        for j, (c, w) in enumerate(allc):
            rcond = rs if w == 0 else rc
            Einj = max(0.75 * c.mu5(rcond), c.M)
            Gband = Gl[j] if not blocking else c.Grest * c.M / max(c.M, 0.75 * c.mu5(rcond))
            dE.append((P[j] - Gl[j] * E[j]) * a ** 4)
            dn.append((P[j] / Einj - Gl[j] * n[j]) * a ** 3)
            Gcw = Gs if w == 0 else Gc
            dB.append((P[j] - (Gband + Gth + 2 * H + 2 * Gcw) * B[j]) * a ** 4)
        dsm = sum(Gl[j] * E[j] for j in range(nC)) * a ** 4
        dnsm = sum(2 * Gl[j] * n[j] for j in range(nC)) * a ** 3
        return dy + dE + dn + dB + [dsm, dnsm]

    y0 = [math.log(a0), st['fs'] * rho0 * a0 ** 3, st['fc'] * rho0 * a0 ** 3]
    Er0 = st['fr'] * rho0
    E0inj = max(allc[1][0].M, 0.75 * allc[1][0].mu5(st['fc'] * rho0))
    y0 += [(Er0 if j == 1 else 0.0) * a0 ** 4 for j in range(nC)]
    y0 += [(Er0 / E0inj if j == 1 else 1e-300) * a0 ** 3 for j in range(nC)]
    y0 += [(Er0 if j == 1 else 0.0) * a0 ** 4 for j in range(nC)]
    y0 += [0.0, 0.0]
    ts = np.unique(np.concatenate([np.linspace(st['t'], 0.4, 300), np.linspace(0.4, 3.0, 400),
                                   np.geomspace(3.0, t_end, n_out)]))
    sol = solve_ivp(rhs, (st['t'], t_end), y0, method='LSODA', rtol=1e-7, atol=1e-30, t_eval=ts, max_step=0.05,
                    first_step=1e-5)
    rec = {k: [] for k in ('t', 'a', 'H', 'f_sig', 'f_chi', 'f_rad', 'f_SM', 'rho_SM', 'n_SM', 'E_SM', 'K', 'rho',
                           'E_N')}
    for i, t in enumerate(sol.t):
        a, rs, rc, E, n, B, rsm, nsm = unpack(sol.y[:, i])
        rho = rs + rc + sum(E) + rsm
        H = math.sqrt(rho / (3 * mP ** 2))
        rec['t'].append(float(t)); rec['a'].append(a); rec['H'].append(H); rec['rho'].append(rho)
        rec['f_sig'].append(rs / rho); rec['f_chi'].append(rc / rho); rec['f_rad'].append((sum(E) + rsm) / rho)
        rec['f_SM'].append(rsm / rho); rec['rho_SM'].append(rsm); rec['n_SM'].append(nsm)
        rec['E_SM'].append(rsm / nsm if nsm > 0 else np.nan); rec['K'].append(a * a * H)
        nt = sum(n)
        rec['E_N'].append(sum(E) / nt if nt > 0 else np.nan)
    out = {k: np.array(v) for k, v in rec.items()}
    out['R_mean'] = 3 * out['H'] ** 2 * (out['f_sig'] + out['f_chi'])
    out['R_osc_chi'] = 9 * out['H'] ** 2 * out['f_chi']
    out['R_osc_sig'] = 9 * out['H'] ** 2 * out['f_sig']
    out['meta'] = dict(name=name, blocking=blocking, Gth=Gth, pert_free=pert_free, grest_scale=grest_scale,
                       start=st, K_final=float(out['K'][-1]), success=bool(sol.success))
    return out


def fluid_cached(name='HFa', blocking=True, Gth=0.0, pert_free=True, grest_scale=1.0, t_end=3e4):
    key = f'fl_{name}_b{int(blocking)}_G{Gth:g}_p{int(pert_free)}_g{grest_scale:.4g}_te{t_end:g}'
    p = SP / 'fluid' / (key + '.npz')
    p.parent.mkdir(exist_ok=True)
    if p.exists():
        z = np.load(p, allow_pickle=True)
        out = {k: z[k] for k in z.files if k != 'meta'}
        out['meta'] = json.loads(str(z['meta']))
        return out
    out = fluid(name, blocking=blocking, Gth=Gth, pert_free=pert_free, grest_scale=grest_scale, t_end=t_end)
    np.savez_compressed(p, meta=json.dumps(out['meta'], default=float), **{k: v for k, v in out.items() if k != 'meta'})
    return out


def first_time(t, x, thr):
    i = np.where(np.asarray(x) >= thr)[0]
    return float(t[i[0]]) if len(i) else None


# --------------------------------------------------------------------------------------------- backgrounds
HISTS = {'EXP': None, 'UNBL': dict(blocking=False), 'PAULI0': dict(blocking=True, Gth=0.0),
         'PAULI1': dict(blocking=True, Gth=1.0), 'PAULI10': dict(blocking=True, Gth=10.0),
         'PAULIb': dict(blocking=True, Gth=0.0, pert_free=False)}
OMEGA_CHI = 480.0      # R line of the chi oscillation (mf_hybrid: 470-490 H_i)
OMEGA_SIG = 1200.0     # 2 m_sigma


class Background:
    """Tabulated background for the Higgs modes.  t from the bench start (deep de Sitter) to t_end.
    R(t) = raw (exact) where the bench is homogeneous, otherwise mean + rebuilt oscillation (osc=True).
    Plasma: hard-gas m_ind^2(t) = c rho_SM / E_SM^2 from this track's fluid (EXP uses the UNBL fluid's SM)."""

    def __init__(self, name='HFa', hist='EXP', osc=True, grest_scale=1.0, Hi_GeV=None, t_end=200.0):
        d = load_bench(name)
        h = d['h']
        self.name, self.hist, self.osc = name, hist, osc
        self.Hi_GeV = d['H_i_GeV'] if Hi_GeV is None else Hi_GeV
        self.tau_endinf = d['t_eps1']
        self.t_end = t_end
        tb, Nb, Hb, Rb = h['t'], h['N'], h['H'], h['R']
        fl_kw = HISTS[hist] if hist != 'EXP' else HISTS['UNBL']
        self.fl = fluid_cached(name, grest_scale=grest_scale, **fl_kw)
        fl = self.fl
        if hist == 'EXP':
            self.t_join = 0.6
            m = tb <= t_end * 1.5
            self.tt = tb[m]; self.lna = Nb[m]; self.HH = Hb[m]
            self.Rm = Rb[m].copy()
            self.fchi = h['f_chi'][m]; self.fsig = h['f_sig'][m]
            hom = self.tt <= self.t_join + 1e-9
            self.Rm_is_raw = hom
            self.K = d['reheating']['K']
        else:
            self.t_join = fl['t'][0]
            m1 = tb < self.t_join - 1e-9
            m2 = fl['t'] <= t_end * 1.5
            self.tt = np.concatenate([tb[m1], fl['t'][m2]])
            self.lna = np.concatenate([Nb[m1], np.log(fl['a'][m2])])
            self.HH = np.concatenate([Hb[m1], fl['H'][m2]])
            self.Rm = np.concatenate([Rb[m1], fl['R_mean'][m2]])
            self.fchi = np.concatenate([h['f_chi'][m1], fl['f_chi'][m2]])
            self.fsig = np.concatenate([h['f_sig'][m1], fl['f_sig'][m2]])
            self.Rm_is_raw = np.concatenate([np.ones(m1.sum(), bool), np.zeros(m2.sum(), bool)])
            self.K = fl['meta']['K_final']
        # mean R for the criterion: running mean over 0.03 H^-1 (>= 2 chi periods) where raw; analytic mean elsewhere
        self.Rbar_tab = self.Rm.copy()
        raw = self.Rm_is_raw & (self.tt > self.tau_endinf - 0.03)
        if raw.any():
            tr = self.tt[self.Rm_is_raw]
            Rr = self.Rm[self.Rm_is_raw]
            cR = np.concatenate([[0.0], np.cumsum(0.5 * (Rr[1:] + Rr[:-1]) * np.diff(tr))])
            for i in np.where(raw)[0]:
                lo, hi = max(self.tt[i] - 0.015, tr[0]), min(self.tt[i] + 0.015, tr[-1])
                if hi > lo:
                    self.Rbar_tab[i] = (np.interp(hi, tr, cR) - np.interp(lo, tr, cR)) / (hi - lo)
        # plasma: hard-gas induced mass (H_i^2 units at this H_i)
        mPr = MPL / self.Hi_GeV
        fl_t = fl['t']
        rho_SM_Hi = fl['f_SM'] * 3 * fl['H'] ** 2 * mPr ** 2          # H_i^4 at this H_i (fixed dimensionless bg)
        E = fl['E_SM']
        self.pl_t = fl_t
        self.pl_m2 = np.where(np.isfinite(E) & (E > 0), C_HARD * rho_SM_Hi / np.where(E > 0, E, 1.0) ** 2, 0.0)
        self.pl_fSM = fl['f_SM']
        self.pl_E = E
        # N_* = const + (1/2) ln H_i at fixed K (radiation bookkeeping of mf_hybrid)
        self.N_obs = N_obs(name, self.K) + 0.5 * math.log(self.Hi_GeV / d['H_i_GeV'])

    def a_H_R(self, t):
        """a, H, R (raw or rebuilt) at time t (scalars)"""
        lna = np.interp(t, self.tt, self.lna)
        H = np.interp(t, self.tt, self.HH)
        R = np.interp(t, self.tt, self.Rm)
        if self.osc and self.t_join < t <= self.t_osc_max:
            fc = np.interp(t, self.tt, self.fchi)
            fs = np.interp(t, self.tt, self.fsig)
            R = R + 9 * H * H * (fc * math.cos(OMEGA_CHI * (t - self.t_join)) + fs * math.cos(OMEGA_SIG * (t - self.t_join)))
        return math.exp(lna), H, R

    def dt_limit(self, t):
        """time-step cap resolving the R lines (raw table spacing in the homogeneous stage; rebuilt lines after)."""
        if t <= self.t_join:
            return 2e-4 if t > self.tau_endinf - 0.05 else 1.0
        if not self.osc or t > self.t_osc_max:
            return 1.0
        fs = np.interp(t, self.tt, self.fsig)
        return (2 * math.pi / OMEGA_SIG) / 16 if fs > 1e-6 else (2 * math.pi / OMEGA_CHI) / 16

    t_osc_max = 6.0          # rebuilt R oscillations kept to t = 6 H_i^-1 (their parametric effect is ~1e-6, checked)

    def Rbar(self, t):
        return np.interp(t, self.tt, self.Rbar_tab)

    def plasma_m2(self, t):
        return np.interp(t, self.pl_t, self.pl_m2, left=0.0, right=self.pl_m2[-1])
