#!/usr/bin/env python3
"""Track e2e_higgs (Claude, end-to-end pass, BRIEF_E2E.md; GitHub issue #13 context).

Higgs survival (xi_H = 0.1, read as xi(M_t)) on the final histories.  Every choice here is a model choice (benchmark
HFa, xi_H, the argued danger criterion D, the closures of the histories, the medium-mass model); none follows from
closure.  Nothing here is a validation, an exclusion, a new theory or a GRUT prediction.

Reuse (read-only imports, nothing copied or edited):
  * tracks/mf_higgs/code/mh_modes.py   -- Higgs mode equations (RK4, Bunch-Davies, Hartree) and the danger exponent D
  * tracks/mf_higgs/code/mh_common.py  -- bench loader, N_obs bookkeeping (N_* = const - ln K/2 + ln H_i/2)
  * tracks/mf_higgs/code/mh_inflation.py -- Hawking-Moss bubble count (constant H and instantaneous)
  * tracks/inf_higgs/code/ih_common.SMTab on tracks/higgs_reheat/code/smpot.py -- SM input (3-loop running, flat
    2-loop lambda_eff, xi running)
What changes here is ONLY the background passed to mh_modes.run (class FinalBG, same interface as
mh_common.Background) and the medium mass column (both mf_higgs-skeptic corrections: the N-gas term
m^2 = sum_j y_j^2 rho_Nj / E_N^2, and the exact eps_H = 1 endpoint in the instantaneous count, see eh_infl.py).

Histories (dimensionless, H_i = 1 at sigma = sigma_c, a = 1 there; bench convention):
  'KIN'      : mf_hybrid bench (raw homogeneous R) to t = 0.19, then e2e_kinetic fluid_HFa.json (a, H, cycle-mean R).
  'LAT_P0'   : bench to 0.19, then e2e_lattice hist_HFa.json lattice block (primary P0 = local Pauli blocking) to
               t = 3 (box-mean a, H, cycle-mean R_cyc), then continued with e2e_kinetic's eps_H(t) [A] integrated
               from the lattice end state (dH/dt = -eps H^2, R = 6 H^2 (2 - eps)).
  'LAT_GM'   : same with the lattice alternate P0GM (skeptic's box-mean saturation).
  'LATF_P0'  : lattice P0 followed by its own HB.stage2 fluid continuation (K ~ 3.9), as exported.
Scratch: /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/e2e_higgs
"""
import sys
sys.dont_write_bytecode = True
import os
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
TRACK = HERE.parent
RES = TRACK / 'results'
TRACKS = TRACK.parent
SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/e2e_higgs')
SP.mkdir(parents=True, exist_ok=True)
# import ih_common FIRST with this track's scratch so that SM caches land here, not in mf_higgs's scratch
os.environ['IH_SCRATCH'] = str(SP)
sys.path.insert(0, str(TRACKS / 'inf_higgs' / 'code'))
import ih_common as IHC  # noqa: E402,F401
sys.path.insert(0, str(TRACKS / 'mf_higgs' / 'code'))
import mh_common as MC   # noqa: E402  (its own os.environ line runs after ih_common is already imported)
import mh_modes as MM    # noqa: E402

MPL = MC.MPL
PDG = MC.PDG24
H_BENCH = 1.30002426398e11
VEW = 174.1

KINF = TRACKS / 'e2e_kinetic' / 'results' / 'fluid_HFa.json'
LATF = TRACKS / 'e2e_lattice' / 'results' / 'hist_HFa.json'
SMF = TRACKS / 'e2e_sm' / 'results' / 'sm_card.json'

# e2e_kinetic K(H_i) (its H_i table, reference.json / digest) [Q]
K_KIN = {5e10: 9.19, 1.3e11: 16.35, 2.4e11: 17.49}


def jdump(obj, path):
    MC.jdump(obj, path)


def _arr(v):
    return np.array([np.nan if x is None else x for x in v], dtype=float)


_K = None


def kin():
    global _K
    if _K is None:
        d = json.load(open(KINF))
        _K = {k: _arr(v) for k, v in d['history'].items()}
        _K['_summary'] = d['summary']
    return _K


_L = None


def lat():
    global _L
    if _L is None:
        _L = json.load(open(LATF))
    return _L


def sm_y2(Hi):
    """e2e_sm S1 Yukawas y_j^2 = mtilde_j M_j / v^2 (N1..N4, M/H_i = 250, 400, 200, 120); y^2 ~ H_i at fixed mtilde."""
    d = json.load(open(SMF))['singlets']
    y2 = np.array(d['y2']['value'], dtype=float)
    return y2 * Hi / 1.3e11


def E_N_inj(Hi):
    """injection energy bound (3/4) mu5,max at eps_H = 1 (mf_higgs skeptic sk_plasma convention, f2 = 2000 H_i);
    held fixed, not redshifted (conservative: m^2 ~ 1/E^2).  H_i units."""
    b = MC.load_bench('HFa')
    h = b['h']
    mP = MPL / Hi
    ie = int(np.argmin(np.abs(h['t'] - b['t_eps1'])))
    rho = 3 * mP ** 2 * h['H'][ie] ** 2
    mu5 = math.sqrt(2 * h['f_chi'][ie] * rho) / (2 * 2000.0)
    return 0.75 * mu5


class FinalBG:
    """Background with the mh_common.Background interface (a_H_R, dt_limit, Rbar, plasma_m2, attributes)."""

    t_raw_end = 0.19

    def __init__(self, hist='KIN', Hi_GeV=1.3e11, K=None, t_end=200.0):
        b = MC.load_bench('HFa')
        hb = b['h']
        self.name, self.hist, self.Hi_GeV, self.t_end = 'HFa', hist, Hi_GeV, t_end
        self.tau_endinf = float(b['t_eps1'])          # exact homogeneous eps_H = 1 (0.19266)
        self.t_join = self.t_raw_end
        m1 = hb['t'] <= self.t_raw_end + 1e-12
        tb, lnab, Hb, Rb = hb['t'][m1], hb['N'][m1], hb['H'][m1], hb['R'][m1]
        k = kin()
        rms = None
        fN_hist = None      # (t, f_N total) of the history itself, for the N-gas medium
        fSM_hist = None
        if hist == 'KIN':
            m2 = k['t'] > self.t_raw_end + 1e-12
            t2, lna2, H2, R2 = k['t'][m2], np.log(k['a'][m2]), k['H'][m2], k['R'][m2]
            Kh = float(k['_summary']['K'])
        else:
            L = lat()
            if hist in ('LAT_P0', 'LATF_P0'):
                blk = L['lattice']
            elif hist == 'LAT_GM':
                blk = L['alternates']['P0GM']['lattice']
            else:
                raise ValueError(hist)
            lt = _arr(blk['t'])
            ml = lt > self.t_raw_end + 1e-12
            tl, lnal, Hl, Rl = lt[ml], np.log(_arr(blk['a'])[ml]), _arr(blk['H'])[ml], _arr(blk['R_cyc'])[ml]
            rmsl = _arr(blk['R_rms'])[ml]
            fNl = _arr(blk['f_rad'])[ml] - _arr(blk['f_SM'])[ml]
            fSMl = _arr(blk['f_SM'])[ml]
            if hist == 'LATF_P0':
                F = L['fluid']
                ft = _arr(F['t'])
                mf = ft > tl[-1] + 1e-9
                tf, lnaf, Hf, Rf = ft[mf], np.log(_arr(F['a'])[mf]), _arr(F['H'])[mf], _arr(F['R_cyc'])[mf]
                fNf = (_arr(F['f_rad']) - _arr(F['f_SM']))[mf]
                fSMf = _arr(F['f_SM'])[mf]
                Kh = float(np.array(F['K_running'], dtype=float)[-1])
            else:
                # continuation with e2e_kinetic's eps_H(t) from the lattice end state [A]
                kt, keps = k['t'], k['eps_H']
                t0, y0 = tl[-1], [lnal[-1], math.log(Hl[-1])]
                lkt = np.log(kt)

                def eps(t):
                    return float(np.interp(math.log(t), lkt, keps))

                def rhs(t, y):
                    return [math.exp(y[1]), -eps(t) * math.exp(y[1])]
                te_ = np.unique(np.concatenate([np.linspace(t0, 10, 400)[1:], np.exp(np.linspace(math.log(10), math.log(kt[-1]), 400))]))
                sol = solve_ivp(rhs, (t0, kt[-1]), y0, t_eval=te_, rtol=1e-9, atol=1e-12, method='LSODA')
                tf, lnaf, Hf = sol.t, sol.y[0], np.exp(sol.y[1])
                Rf = 6 * Hf ** 2 * (2 - np.array([eps(x) for x in tf]))
                # N and SM fractions after t_sw: e2e_kinetic's, rescaled to match the lattice at t_sw [A]
                fNk = np.interp(tf, kt, k['f_N'])
                fSMk = np.interp(tf, kt, k['f_SM'])
                rN = fNl[-1] / max(np.interp(t0, kt, k['f_N']), 1e-30)
                rS = fSMl[-1] / max(np.interp(t0, kt, k['f_SM']), 1e-30)
                fNf, fSMf = np.minimum(fNk * rN, 1.0), np.minimum(fSMk * rS, 1.0)
                Kh = float(np.exp(2 * lnaf[-1]) * Hf[-1])
            t2 = np.concatenate([tl, tf]); lna2 = np.concatenate([lnal, lnaf])
            H2 = np.concatenate([Hl, Hf]); R2 = np.concatenate([Rl, Rf])
            rms = (t2, np.concatenate([rmsl, np.zeros(len(tf))]))
            fN_hist = (t2, np.concatenate([fNl, fNf]))
            fSM_hist = (t2, np.concatenate([fSMl, fSMf]))
        self.K_hist = Kh
        self.K = Kh if K is None else K
        # ln a rebuilt from H after t_raw_end (trapezoid), so that the friction 3H and a^3 in the Wronskian agree;
        # the tabulated a of the histories agrees with this to <= 1% (difference logged as lna_mismatch)
        lna_int = lnab[-1] + np.concatenate([[0.0], np.cumsum(0.5 * (H2[1:] + H2[:-1]) * np.diff(t2))]) \
            + 0.5 * (Hb[-1] + H2[0]) * (t2[0] - tb[-1])
        self.lna_mismatch = float(np.max(np.abs(lna_int - lna2)))
        self.tt = np.concatenate([tb, t2]); self.lna = np.concatenate([lnab, lna_int])
        self.HH = np.concatenate([Hb, H2]); self.Rm = np.concatenate([Rb, R2])
        # criterion mean R: running mean (+-0.015) of the raw bench R up to t_raw_end, the history's cycle mean after
        hr = b['h']
        mr = hr['t'] <= 0.6 + 1e-9                     # bench homogeneous stage (raw) for the running mean
        tr, Rr = hr['t'][mr], hr['R'][mr]
        cR = np.concatenate([[0.0], np.cumsum(0.5 * (Rr[1:] + Rr[:-1]) * np.diff(tr))])
        self.Rbar_tab = self.Rm.copy()
        for i in np.where((self.tt <= self.t_raw_end) & (self.tt > self.tau_endinf - 0.05))[0]:
            lo, hi = max(self.tt[i] - 0.015, tr[0]), min(self.tt[i] + 0.015, tr[-1])
            self.Rbar_tab[i] = (np.interp(hi, tr, cR) - np.interp(lo, tr, cR)) / (hi - lo)
        self.Rrms_tab = np.zeros(len(self.tt))
        if rms is not None:
            self.Rrms_tab = np.interp(self.tt, rms[0], rms[1], left=0.0, right=0.0)
        # ---- medium masses (H_i^2 at this H_i), both [E]
        mP = MPL / Hi_GeV
        rho = 3 * mP ** 2 * self.HH ** 2
        E_N = E_N_inj(Hi_GeV)
        y2 = sm_y2(Hi_GeV)
        kt = k['t']
        fj = np.array([np.interp(self.tt, kt, k[f'f_N{j}'], left=0.0) for j in (1, 2, 3, 4)])
        fN_k = np.interp(self.tt, kt, k['f_N'], left=0.0)
        if fN_hist is None:
            fN = fN_k
        else:
            fN = np.interp(self.tt, fN_hist[0], fN_hist[1], left=0.0)
            fN = np.where(self.tt <= self.t_raw_end, 0.0, np.maximum(fN, 0.0))
        share = np.where(fN_k > 0, fj / np.where(fN_k > 0, fN_k, 1.0), 0.25)
        y2eff = np.sum(y2[:, None] * share, axis=0)
        y2eff = np.where(fN_k > 0, y2eff, y2.min())
        self.m2_N = y2eff * fN * rho / E_N ** 2
        self.m2_N_ymin = y2.min() * fN * rho / E_N ** 2
        if hist == 'KIN':
            # e2e_kinetic's own medium column kH Tcal^2 (computed at 1.3e11); Tcal^2/H_i^2 ~ 1/H_i at fixed dimensionless bg
            self.m2_SM = np.interp(self.tt, kt, k['m_H2_medium'], left=0.0) * (1.3e11 / Hi_GeV)
        else:
            # SM hard gas c rho_SM / (E_N/2)^2 (mf_higgs / skeptic formula, energies fixed) on the lattice's own f_SM
            fSM = np.interp(self.tt, fSM_hist[0], fSM_hist[1], left=0.0)
            fSM = np.where(self.tt <= self.t_raw_end, 0.0, np.maximum(fSM, 0.0))
            self.m2_SM = MC.C_HARD * fSM * rho / (E_N / 2) ** 2
        self.pl_m2 = self.m2_N + self.m2_SM
        self.E_N = E_N
        self.N_obs = MC.N_obs('HFa', self.K) + 0.5 * math.log(Hi_GeV / H_BENCH)

    def a_H_R(self, t):
        return (math.exp(float(np.interp(t, self.tt, self.lna))), float(np.interp(t, self.tt, self.HH)),
                float(np.interp(t, self.tt, self.Rm)))

    def dt_limit(self, t):
        if t <= self.t_raw_end:
            return 2e-4 if t > self.tau_endinf - 0.05 else 1.0
        return 1.0

    def Rbar(self, t):
        return float(np.interp(t, self.tt, self.Rbar_tab))

    def Rrms(self, t):
        return np.interp(t, self.tt, self.Rrms_tab)

    def plasma_m2(self, t):
        return float(np.interp(t, self.tt, self.pl_m2, left=0.0))

    def col(self, name, t):
        return np.interp(t, self.tt, getattr(self, name), left=0.0)


_SMT = {}


def smt(Mt, Mh=PDG['Mh'], als=PDG['als']):
    k = (round(Mt, 4), round(Mh, 3), round(als, 5))
    if k not in _SMT:
        _SMT[k] = IHC.SMTab(k[0], Mh=k[1], als=k[2])
    return _SMT[k]


def xi_at(Hi, xi_EW=0.1, Mt=172.57):
    return smt(Mt).xi(Hi, xi_EW)


def run_cached(hist, Hi, xi_EW=0.1, n=260, t_end=80.0, rerun=False):
    xi = xi_at(Hi, xi_EW)
    key = f'{hist}_Hi{Hi:.3g}_xiEW{xi_EW:.5g}_n{n}_te{t_end}'
    p = SP / 'runs' / (key + '.npz')
    p.parent.mkdir(exist_ok=True)
    if p.exists() and not rerun:
        z = np.load(p, allow_pickle=True)
        return dict(res=json.loads(str(z['res'])), rows=z['rows'], spec=z['spec'], ks=z['ks'], wq=z['wq'], key=key)
    bg = FinalBG(hist, Hi)
    R = MM.run(bg, xi, smt(172.57), t_end=t_end, n=n)
    R['res']['Hi_GeV'] = Hi
    R['res']['xi_EW'] = xi_EW
    R['res']['K_hist'] = bg.K_hist
    np.savez_compressed(p, res=json.dumps(R['res'], default=float), rows=R['rows'], spec=R['spec'], ks=R['ks'], wq=R['wq'])
    R['key'] = key
    return R
