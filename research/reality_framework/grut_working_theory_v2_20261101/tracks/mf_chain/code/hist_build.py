#!/usr/bin/env python3
"""mf_chain: skeptic-corrected post-slice histories of the mf_hybrid benchmarks (HFa, HFb), exported in the bench_*.json
format that mf_ideal's bench_run.py reads, at the design H_i or at a rescaled H_i (V0 free).

WHAT IS COPIED, AND WHAT CHANGED (stated, as the brief asks):
  * Stage 1 (homogeneous sigma, chi with axial decay friction, slice to t_sw): the equations are those of the mf_hybrid
    skeptic's sk_waterfall.bg (tracks/mf_hybrid/skeptic/code), re-typed here without the linear modes, with three
    changes: (i) the blocking variable B_j obeys the same equation as in the skeptic's fluid (unblocking by N decays,
    the N-thermalisation bracket Gamma_th, redshift and band-edge shrinkage) instead of the never-decaying cumulative
    energy; (ii) the perturbative tails (sigma -> N1 at 0.29 H_i, chi -> N4) act here too, free or blocked as in the
    fluid; (iii) every scale is a parameter, so H_i can be rescaled (rule below).  'track' friction (no blocking) is kept
    as the validation switch.
  * Stage 2 (fluid, t_sw to K convergence): a copy of the skeptic's sk_pauli.run (per-species N energy, number and
    in-band blocking energy; Pauli saturation F_sat(A) from the skeptic's exact multi-period mode runs; production tables
    eps_j(A) of mf_hybrid's axial_tables.json, both read-only), changed only in: the per-species initial state is taken
    from stage 1 (not 'all radiation in N2'); M_P/H_i and the seesaw rate scale are parameters; adaptive max_step for the
    late tail; and the record carries a, H, eps_H, R (cycle mean, w = 0 condensates) and the 2 m oscillation amplitudes.
  * Pre-slice rows (t < 0) are copied from the bench (exit-factor check only; the slice is unchanged).
  * The spectator s/s_* is integrated on the recorded background for the s check only (the engine evolves its own).
    m_s^2 model (mf_hybrid bench.py): HFa g_s^2 chi^2 (homogeneous) -> g_s^2 (v^2 + rho_chi/m_chi^2) (fluid);
    HFb 9 - (3/4) R (raw R in stage 1, cycle-mean R in stage 2).

RESCALING RULE (model choice, used only for the H_i iteration).  Hold fixed, in H_i units: m_chi, m_sig, the clock's
tree c1 and m2 (so eps_c ~ H^2 and A_s stays fixed to the precision the post-slice history needs), mu_R^2, the singlet
masses M_j and decay constants f_j, the spectator mass 3 H_i and |g_h| s_* = 12.30 H_i.  Hold fixed in M_P units: v,
sigma_c, s_*.  Bias: c/(lam v^3) fixed.  Then lam ~ H^2, g ~ H, the axial amplitude A = Sigma/2f ~ M_P/H_i, and the
seesaw N widths in H_i units ~ H_i (m_nu-tilde fixed at 0.05 eV).  The CMB is checked separately (cmb_rescale.py).

Units inside: H_ref = 1.30e11 GeV design, mP = M_P/H_ref' with H_ref' = H_GeV (the rescaled design scale).  The export is
in units of H_i = H(sigma = sigma_c) of the run (H_i/H_design = 1 + 2e-5), a = 1 there.
Usage: python3 -B hist_build.py NAME TAG [--H 1.3e11] [--fric pauli|track] [--Gth 0] [--tails free|blocked]
                                [--tsw 0.45] [--tend 1e9]
Labels: [C] on the [E] rates and blocking model of the mf_hybrid skeptic; every choice is a model choice.
"""
import argparse, json, math, os, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import mfc_common as M

sys.path.insert(0, str(M.SKC))
import sk_pauli as SPA                                   # noqa: E402  read-only: Chan (tables eps(A), F_sat(A))

MP_GEV, MNU, VEW = 2.435e18, 0.05e-9, 174.1
INTEG = []                                               # (t_lo, t_hi, integrator) of stage 2, for the record


# ------------------------------------------------------------------ model at a (possibly rescaled) design scale H_GeV
def model(name, H_GeV):
    d = M.load_bench(name)
    pr = d['params']
    mP0 = MP_GEV / M.H_REF
    mP = MP_GEV / H_GeV
    r = mP / mP0
    P = dict(name=name, H_GeV=H_GeV, mP=mP, lam=pr['lam'] / r ** 2, v=pr['v_MP'] * mP, g=pr['g'] / r, sc=pr['sigma_c_MP'] * mP,
             c1=pr['c1_Hi3'], m2=pr['m2_Hi2'], mu2=pr['mu_R2_Hi2'], cb=pr['c_bias_Hi3'] * r,
             m_chi=pr['m_chi_over_Hi'], m_sig=pr['m_sig_over_Hi'], smech=pr['smech'], s_star_MP=pr['s_star_MP'],
             sig_ch=[(c['M_over_Hi'], c['f_over_Hi']) for c in pr['axial_sigma']],
             chi_ch=[(c['M_over_Hi'], c['f_over_Hi']) for c in pr['axial_chi']])
    P['V0'] = 0.25 * P['lam'] * P['v'] ** 4
    P['g_s'] = 3.0 / P['v']
    # tree vacuum and the constant making V(vacuum) = 0
    s, x = 0.0, P['v']
    for _ in range(80):
        Vs, Vx = grad2(P, s, x, loop=False)
        a_, b_, c_ = hess2(P, s, x)
        det = a_ * b_ - c_ * c_
        s, x = s - (b_ * Vs - c_ * Vx) / det, x - (a_ * Vx - c_ * Vs) / det
    P['s0'], P['x0'] = s, x
    P['Vc'] = 0.0
    P['Vc'] = -V2(P, s, x, loop=False)
    return d, P


def M2c(P, s):
    return P['g'] ** 2 * s * s - P['lam'] * P['v'] ** 2


def V2(P, s, x, loop=True):
    d = s - P['sc']
    out = (0.25 * P['lam'] * (x * x - P['v'] ** 2) ** 2 + 0.5 * P['g'] ** 2 * s * s * x * x - P['cb'] * x
           + P['c1'] * d - 0.5 * P['m2'] * d * d + P['Vc'])
    if loop and s > P['sc']:
        m2 = M2c(P, s)
        out += m2 * m2 * (math.log(abs(m2) / P['mu2']) - 1.5) / (64 * math.pi ** 2)
    return out


def grad2(P, s, x, loop=True):
    d = s - P['sc']
    Vs = P['g'] ** 2 * s * x * x + P['c1'] - P['m2'] * d
    if loop and s > P['sc']:
        m2 = M2c(P, s)
        Vs += m2 * 2 * P['g'] ** 2 * s * (math.log(abs(m2) / P['mu2']) - 1.0) / (32 * math.pi ** 2)
    Vx = P['lam'] * (x * x - P['v'] ** 2) * x + P['g'] ** 2 * s * s * x - P['cb']
    return Vs, Vx


def hess2(P, s, x):
    return (P['g'] ** 2 * x * x - P['m2'], P['lam'] * (3 * x * x - P['v'] ** 2) + P['g'] ** 2 * s * s, 2 * P['g'] ** 2 * s * x)


def valley(P, s):
    """single-field valley potential above sigma_c (chi at its biased minimum c/M^2), and its slope."""
    m2 = M2c(P, s)
    d = s - P['sc']
    V = P['V0'] + P['c1'] * d - 0.5 * P['m2'] * d * d + m2 * m2 * (math.log(abs(m2) / P['mu2']) - 1.5) / (64 * math.pi ** 2) \
        - 0.5 * P['cb'] ** 2 / m2
    V1 = P['c1'] - P['m2'] * d + m2 * 2 * P['g'] ** 2 * s * (math.log(abs(m2) / P['mu2']) - 1.0) / (32 * math.pi ** 2) \
        + P['cb'] ** 2 * P['g'] ** 2 * s / m2 ** 2
    return V, V1


# ------------------------------------------------------------------ channels (skeptic's Chan, rescaled seesaw rate)
def channels(P):
    cs = [SPA.Chan('N1', P['m_sig'], M_, f_) for M_, f_ in P['sig_ch']]
    cc = [SPA.Chan(t, P['m_chi'], M_, f_) for t, (M_, f_) in zip(('N2', 'N3', 'N4'), P['chi_ch'])]
    for c in cs + cc:
        c.Grest = MNU * (c.M * P['H_GeV']) ** 2 / (8 * math.pi * VEW ** 2) / P['H_GeV']
    return [(c, 0) for c in cs] + [(c, 1) for c in cc]


def rates(P, allc, rs, rc, B, E, n, a, blocking, pert_free, Gth, H, Gs_prev=0.0, Gc_prev=0.0):
    """per-species production P_j (energy density per time), lab decay rates, band unblocking rates (fluid form)."""
    Pj, Gl, Gb = [], [], []
    for j, (c, w) in enumerate(allc):
        rcond = rs if w == 0 else rc
        G = c.rate(rcond)
        if blocking:
            Es = c.fsat(rcond) * rcond
            bl = max(0.0, 1.0 - B[j] / Es) if Es > 0 else 0.0
        else:
            bl = 1.0
        Pnp = G * rcond * bl
        if pert_free == 'band':
            # perturbative products sit at k = (m/2) beta; they are inside the Pauli band only while mu5,max > k
            kp = 0.5 * c.m * math.sqrt(max(1 - 4 * c.M ** 2 / c.m ** 2, 0.0))
            free = (not blocking) or c.mu5(rcond) <= kp          # skeptic's mu5(rho) = m Sigma/2f = mu5,max
        else:
            free = bool(pert_free) or not blocking
        Pp = c.Gp * rcond * (1.0 if free else bl)
        Pj.append(max(Pnp, Pp))
        Emean = E[j] / n[j] if n[j] > 0 else c.M
        Gl.append(c.Grest * c.M / max(c.M, Emean))
        Gb.append(c.Grest * c.M / max(c.M, 0.75 * c.mu5(rcond)) if blocking else Gl[-1])
    return Pj, Gl, Gb


# ------------------------------------------------------------------ stage 1: homogeneous
def stage1(P, fric='pauli', Gth=0.0, pert_free=True, t_sw=0.45, dt_rec=2e-4):
    allc = channels(P)
    nC = len(allc)
    mP = P['mP']
    blocking = (fric == 'pauli')

    def energies(s, sd, x, xd):
        Es = 0.5 * sd * sd + 0.5 * P['g'] ** 2 * x * x * (s - P['s0']) ** 2
        Et = 0.5 * (sd * sd + xd * xd) + V2(P, s, x)
        return max(Es, 0.0), max(Et - Es, 0.0)

    def unpack(y):
        a = math.exp(y[4])
        E = [max(v, 0) / a ** 4 for v in y[5:5 + nC]]
        n = [max(v, 1e-300) / a ** 3 for v in y[5 + nC:5 + 2 * nC]]
        B = [max(v, 0) / a ** 4 for v in y[5 + 2 * nC:5 + 3 * nC]]
        sm = max(y[5 + 3 * nC], 0) / a ** 4
        return a, E, n, B, sm

    def rhs(t, y, on):
        s, sd, x, xd = y[:4]
        a, E, n, B, sm = unpack(y)
        rr = sum(E) + sm
        Vv = V2(P, s, x)
        Kin = 0.5 * (sd * sd + xd * xd)
        H = math.sqrt(max(Kin + Vv + rr, 1e-300) / (3 * mP ** 2))
        Vs, Vx = grad2(P, s, x)
        dE, dn, dB = [0.0] * nC, [0.0] * nC, [0.0] * nC
        Gs = Gc = 0.0
        dsm = 0.0
        if on and fric != 'hubble':
            rs, rc = energies(s, sd, x, xd)
            Pj, Gl, Gb = rates(P, allc, rs, rc, B, E, n, a, blocking, pert_free, Gth, H)
            # friction coefficients so that the instantaneous transfer is Gamma phi_dot^2 (cycle average Gamma rho)
            for j, (c, w) in enumerate(allc):
                rcond = rs if w == 0 else rc
                G = Pj[j] / rcond if rcond > 0 else 0.0
                vel2 = sd * sd if w == 0 else xd * xd
                inj = G * vel2
                if w == 0:
                    Gs += G
                else:
                    Gc += G
                Einj = max(0.75 * c.mu5(rcond), c.M)
                dE[j] = (inj - Gl[j] * E[j]) * a ** 4
                dn[j] = (inj / Einj - Gl[j] * n[j]) * a ** 3
                dB[j] = (inj - (Gb[j] + Gth + 2 * H) * B[j]) * a ** 4
            dsm = sum(Gl[j] * E[j] for j in range(nC)) * a ** 4
        return [sd, -(3 * H + Gs) * sd - Vs, xd, -(3 * H + Gc) * xd - Vx, H] + dE + dn + dB + [dsm]

    # start 1 e-fold before sigma_c on the slow-roll attractor of the valley
    def Nsr(z):
        g = np.linspace(P['sc'] * (1 + 1e-6), z, 2001)
        f = [valley(P, q)[0] / (mP ** 2 * valley(P, q)[1]) for q in g]
        return float(np.trapezoid(f, g))
    hi = P['sc'] + 0.002 * mP
    s_start = brentq(lambda z: Nsr(z) - 1.0, P['sc'] * (1 + 1e-4), hi, xtol=1e-9 * mP)
    Vst, V1st = valley(P, s_start)
    sd0 = -V1st / math.sqrt(3 * Vst) * mP
    x0 = P['cb'] / M2c(P, s_start)
    y0 = [s_start, sd0, x0, 0.0, 0.0] + [0.0] * nC + [1e-300] * nC + [0.0] * nC + [0.0]
    ev = lambda t, y, on: y[0] - P['sc']
    ev.terminal, ev.direction = True, -1
    s1 = solve_ivp(rhs, (0, 50), y0, args=(False,), method='DOP853', rtol=1e-11, atol=1e-13, events=ev, max_step=0.01)
    yc = list(s1.y_events[0][0])
    yc[4] = 0.0
    tt = np.arange(0.0, t_sw + 0.5 * dt_rec, dt_rec)
    s2 = solve_ivp(rhs, (0, tt[-1]), yc, args=(True,), method='DOP853', rtol=1e-10, atol=1e-12, t_eval=tt,
                   max_step=0.02 / P['m_chi'], dense_output=False)
    Y = s2.y
    rec = dict(t=s2.t, s=Y[0], sd=Y[1], x=Y[2], xd=Y[3], la=Y[4])
    H, R, eps, fs, fc, fr, fsm, rsm = [], [], [], [], [], [], [], []
    for j in range(len(s2.t)):
        s, sd, x, xd = Y[:4, j]
        a, E, n, B, sm = unpack(Y[:, j])
        Vv = V2(P, s, x)
        Kin = 0.5 * (sd * sd + xd * xd)
        rr = sum(E) + sm
        rho = Kin + Vv + rr
        h = math.sqrt(rho / (3 * mP ** 2))
        Hd = -(sd * sd + xd * xd + (4.0 / 3.0) * rr) / (2 * mP ** 2)
        H.append(h); eps.append(-Hd / h ** 2)
        R.append((4 * Vv - (sd * sd + xd * xd)) / mP ** 2)
        es, ec = energies(s, sd, x, xd)
        fs.append(es / rho); fc.append(ec / rho); fr.append(rr / rho); fsm.append(sm / rho)
    rec.update(H=np.array(H), R=np.array(R), eps=np.array(eps), fs=np.array(fs), fc=np.array(fc), fr=np.array(fr),
               fsm=np.array(fsm))
    a, E, n, B, sm = unpack(Y[:, -1])
    s, sd, x, xd = Y[:4, -1]
    es, ec = energies(s, sd, x, xd)
    rec['end'] = dict(t=float(s2.t[-1]), a=a, rs=es, rc=ec, E=E, n=n, B=B, sm=sm, H=H[-1])
    rec['s_start'] = s_start
    return rec


# ------------------------------------------------------------------ stage 2: fluid (copy of sk_pauli.run, changes in header)
def stage2(P, st, blocking=True, Gth=0.0, pert_free=True, t_end=1e9, nrec=1500):
    allc = channels(P)
    nC = len(allc)
    mP = P['mP']

    def unpack(y):
        a = math.exp(y[0])
        rs, rc = max(y[1], 0) / a ** 3, max(y[2], 0) / a ** 3
        E = [max(y[3 + j], 0) / a ** 4 for j in range(nC)]
        n = [max(y[3 + nC + j], 1e-300) / a ** 3 for j in range(nC)]
        B = [max(y[3 + 2 * nC + j], 0) / a ** 4 for j in range(nC)]
        rsm = max(y[3 + 3 * nC], 0) / a ** 4
        return a, rs, rc, E, n, B, rsm

    def rhs(t, y):
        a, rs, rc, E, n, B, rsm = unpack(y)
        rho = rs + rc + sum(E) + rsm
        H = math.sqrt(rho / (3 * mP ** 2))
        Pj, Gl, Gb = rates(P, allc, rs, rc, B, E, n, a, blocking, pert_free, Gth, H)
        Gs = sum(Pj[j] for j, (c, w) in enumerate(allc) if w == 0) / max(rs, 1e-300)
        Gc = sum(Pj[j] for j, (c, w) in enumerate(allc) if w == 1) / max(rc, 1e-300)
        dy = [H, -Gs * rs * a ** 3, -Gc * rc * a ** 3]
        dE, dn, dB = [], [], []
        for j, (c, w) in enumerate(allc):
            rcond = rs if w == 0 else rc
            Einj = max(0.75 * c.mu5(rcond), c.M)
            dE.append((Pj[j] - Gl[j] * E[j]) * a ** 4)
            dn.append((Pj[j] / Einj - Gl[j] * n[j]) * a ** 3)
            Gcw = Gs if w == 0 else Gc
            dB.append((Pj[j] - (Gb[j] + Gth + 2 * H + 2 * Gcw) * B[j]) * a ** 4)
        dsm = sum(Gl[j] * E[j] for j in range(nC)) * a ** 4
        return dy + dE + dn + dB + [dsm]

    a0 = st['a']
    y0 = [math.log(a0), st['rs'] * a0 ** 3, st['rc'] * a0 ** 3]
    y0 += [st['E'][j] * a0 ** 4 for j in range(nC)]
    y0 += [max(st['n'][j], 1e-300) * a0 ** 3 for j in range(nC)]
    y0 += [st['B'][j] * a0 ** 4 if blocking else 0.0 for j in range(nC)]
    y0 += [st['sm'] * a0 ** 4]
    t0 = st['t']
    # segments with growing max_step; stop when the condensates hold < 1e-7 of rho
    grid = [1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8, 3e8, 1e9, 3e9, 1e10]
    edges = sorted(set([t0] + [e for e in grid if t0 < e < t_end] + [t_end]))
    T, Yall = [t0], [np.array(y0)]
    y = y0
    done = False
    for lo, hi in zip(edges[:-1], edges[1:]):
        ms = max(0.02, 0.02 * lo)
        te = np.linspace(lo, hi, 120) if hi <= 3.0 else np.geomspace(lo, hi, 70)
        sol = None
        for meth, mstep, rt in (('LSODA', ms, 1e-8), ('Radau', ms, 1e-7), ('RK45', min(ms, 0.002 * max(lo, 1.0)), 1e-7)):
            try:
                sol = solve_ivp(rhs, (lo, hi), y, method=meth, rtol=rt, atol=1e-30, t_eval=te, max_step=mstep)
            except Exception:
                sol = None
            if sol is not None and sol.status == 0 and isinstance(sol.y, np.ndarray) and sol.y.shape[1] == len(te):
                break
            sol = None
        if sol is None:
            raise RuntimeError('stage2 integration failed on [%g, %g]' % (lo, hi))
        INTEG.append((lo, hi, meth))
        for i in range(1, sol.y.shape[1]):
            T.append(float(sol.t[i])); Yall.append(sol.y[:, i])
        y = sol.y[:, -1]
        a, rs, rc, E, n, B, rsm = unpack(y)
        rho = rs + rc + sum(E) + rsm
        if (rs + rc) / rho < 1e-7 and hi >= 1e3:
            done = True
            break
    rec = dict(t=[], a=[], H=[], eps=[], R=[], Rs=[], Rc=[], frad=[], fsm=[], fchi=[], fsig=[], K=[], rc=[], rs=[])
    for t, yv in zip(T, Yall):
        a, rs, rc, E, n, B, rsm = unpack(yv)
        rho = rs + rc + sum(E) + rsm
        H = math.sqrt(rho / (3 * mP ** 2))
        p = (sum(E) + rsm) / 3.0
        rec['t'].append(float(t)); rec['a'].append(a); rec['H'].append(H)
        rec['eps'].append(1.5 * (rho + p) / rho)
        rec['R'].append((rs + rc) / mP ** 2); rec['Rs'].append(3 * rs / mP ** 2); rec['Rc'].append(3 * rc / mP ** 2)
        rec['frad'].append((sum(E) + rsm) / rho); rec['fsm'].append(rsm / rho)
        rec['fchi'].append(rc / rho); rec['fsig'].append(rs / rho); rec['K'].append(a * a * H)
        rec['rc'].append(rc); rec['rs'].append(rs)
    rec = {k: np.array(v) for k, v in rec.items()}
    rec['converged'] = done
    rec['integrators'] = list(INTEG)
    return rec


# ------------------------------------------------------------------ assemble the bench-format export
def build(name, tag, H_GeV=M.H_REF, fric='pauli', Gth=0.0, pert_free=True, t_sw=0.45, t_end=1e9, write=True):
    d, P = model(name, H_GeV)
    s1 = stage1(P, fric, Gth, pert_free, t_sw)
    s2 = stage2(P, s1['end'], blocking=(fric == 'pauli'), Gth=Gth, pert_free=pert_free, t_end=t_end)
    Hi = float(s1['H'][0])                                   # H at sigma = sigma_c (units of the design scale)
    mP = P['mP']
    # spectator mass^2 (design units)
    if P['smech'] == 'a':
        ms2_1 = (P['g_s'] * s1['x']) ** 2
        ms2_2 = P['g_s'] ** 2 * (P['v'] ** 2 + s2['rc'] / P['m_chi'] ** 2)
    else:
        ms2_1 = 9.0 - 0.75 * s1['R']
        ms2_2 = 9.0 - 0.75 * s2['R']
    t = np.concatenate([s1['t'], s2['t'][1:]])
    la = np.concatenate([s1['la'], np.log(s2['a'][1:])])
    H = np.concatenate([s1['H'], s2['H'][1:]])
    eps = np.concatenate([s1['eps'], s2['eps'][1:]])
    R = np.concatenate([s1['R'], s2['R'][1:]])
    ms2 = np.concatenate([ms2_1, ms2_2[1:]])
    fsig = np.concatenate([s1['fs'], s2['fsig'][1:]]); fchi = np.concatenate([s1['fc'], s2['fchi'][1:]])
    frad = np.concatenate([s1['fr'], s2['frad'][1:]]); fsm = np.concatenate([s1['fsm'], s2['fsm'][1:]])
    Ros = np.concatenate([np.full(len(s1['t']), np.nan), s2['Rs'][1:]])
    Roc = np.concatenate([np.full(len(s1['t']), np.nan), s2['Rc'][1:]])
    sig = np.concatenate([s1['s'] / mP, np.full(len(s2['t']) - 1, np.nan)])
    chi = np.concatenate([s1['x'] / mP, np.full(len(s2['t']) - 1, np.nan)])
    # convert to H_i units
    tH, HH, RH, msH = t * Hi, H / Hi, R / Hi ** 2, np.sqrt(np.maximum(ms2, 0)) / Hi
    RosH, RocH = Ros / Hi ** 2, Roc / Hi ** 2
    # spectator c = s/s_* on the recorded background (check only), exact to t = 60 H_i^-1, then flagged envelope
    from scipy.interpolate import interp1d
    fH = interp1d(tH, HH, kind='linear'); fm = interp1d(tH, msH ** 2, kind='linear')
    tmax = min(60.0, tH[-1])
    ss = solve_ivp(lambda q, y: [y[1], -3 * fH(q) * y[1] - fm(q) * y[0]], (0, tmax), [1.0, 0.0], method='DOP853', rtol=1e-10,
                   atol=1e-13, dense_output=True, max_step=0.01)
    c = np.where(tH <= tmax, ss.sol(np.minimum(tH, tmax))[0], np.nan)
    cd = np.where(tH <= tmax, ss.sol(np.minimum(tH, tmax))[1], np.nan)
    j60 = np.searchsorted(tH, tmax) - 1
    amp60 = math.hypot(c[j60], cd[j60] / msH[j60]) if j60 >= 0 else 1.0
    a = np.exp(la)
    env = amp60 * (a[j60] / a) ** 1.5
    s_out = np.where(tH <= tmax, c, env)
    env_flag = (tH > tmax).astype(float)
    sst = P['s_star_MP'] * mP                                # in design units
    rho_tot = 3 * mP ** 2 * H ** 2
    E_s = np.where(tH <= tmax, 0.5 * (cd * Hi) ** 2 + 0.5 * ms2 * c ** 2, 0.5 * ms2 * env ** 2) * sst ** 2
    f_s = E_s / rho_tot
    # pre-slice rows from the bench (identical inflationary history; exit factor only)
    h = d['history']
    tb = np.array(h['t'])
    pre = tb < 0
    hist = {}
    cols = ['t', 'N', 'H', 'eps_H', 'R', 'sigma_MP', 'chi_MP', 's_over_sstar', 'm_s', 'f_sig', 'f_chi', 'f_rad', 'f_SM', 'f_s', 'a',
            'K_running', 'R_osc_sig', 'R_osc_chi', 's_is_envelope']
    new = dict(t=tH, N=la, H=HH, eps_H=eps, R=RH, sigma_MP=sig, chi_MP=chi, s_over_sstar=s_out, m_s=msH, f_sig=fsig, f_chi=fchi,
               f_rad=frad, f_SM=fsm, f_s=f_s, a=a, K_running=a * a * HH, R_osc_sig=RosH, R_osc_chi=RocH, s_is_envelope=env_flag)
    for k in cols:
        # pre-slice inflation rows: H constant to 1e-7 in H_i units at any design scale (eps ~ 1e-8 (H/H_ref)^2)
        old = [h[k][i] for i in range(len(tb)) if pre[i]]
        vals = [None if (isinstance(x, float) and not np.isfinite(x)) else float('%.10g' % x) for x in new[k]]
        hist[k] = old + vals
    K = float(new['K_running'][-1])
    out = dict(name=tag, label='C on the [E] rates and the mf_hybrid skeptic\'s Pauli blocking model; model choices throughout',
               note='mf_chain corrected history of %s: friction %s, Gamma_th %g H_i, tails %s, t_sw %g; H design %.4g GeV' % (
                   name, fric, Gth, {True: 'free', False: 'blocked', 'band': 'band'}[pert_free], t_sw, H_GeV),
               slice='a = 1, H = H_i at sigma = sigma_c (as the bench)', H_i_GeV=Hi * H_GeV,
               params=dict(d['params'], H_i_GeV=Hi * H_GeV, s_star_GeV=P['s_star_MP'] * MP_GEV, lam=P['lam'], g=P['g'],
                           c_bias_Hi3=P['cb'], g_s=P['g_s'] if P['smech'] == 'a' else 0.0, design_H_GeV=H_GeV,
                           rescaling_rule='see hist_build.py header'),
               spectator=dict(Gamma_s_assumed=d['spectator'].get('Gamma_s_assumed', 1.3e-6)),
               reheating=dict(K=K, converged=bool(s2['converged']), integrators=[list(x) for x in s2['integrators']], t_end=float(tH[-1]), friction=fric, Gamma_th=Gth,
                              tails={True: 'free', False: 'blocked', 'band': 'blocked while in band'}[pert_free], t_sw=t_sw,
                              frad_t03=float(np.interp(0.3, tH, frad)), frad_t06=float(np.interp(0.6, tH, frad)),
                              frad_t3=float(np.interp(3.0, tH, frad)), frad_t100=float(np.interp(100.0, tH, frad)),
                              t_SM_1em3=float(tH[np.argmax(fsm >= 1e-3)]) if (fsm >= 1e-3).any() else None,
                              t_SM_half=float(tH[np.argmax(fsm >= 0.5)]) if (fsm >= 0.5).any() else None,
                              t_rad_99=float(tH[np.argmax(frad >= 0.99)]) if (frad >= 0.99).any() else None,
                              f_s_at_t60=float(np.interp(60.0, tH, f_s)), f_s_over_a_late=float(f_s[-200] / a[-200]) if len(a) > 200 else None),
               waterfall=dict(t_eps1=float(tH[np.argmax(eps >= 1.0)]), H_eps1=float(HH[np.argmax(eps >= 1.0)]), R_min=float(RH.min())),
               history_columns=d['history_columns'], history=hist, status='mf_chain corrected export (not mf_hybrid\'s)')
    if write:
        path = M.SP / 'hist' / ('bench_%s.json' % tag)
        M.dump(path, out)
        out['path'] = str(path)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('name'); ap.add_argument('tag')
    ap.add_argument('--H', type=float, default=M.H_REF); ap.add_argument('--fric', default='pauli')
    ap.add_argument('--Gth', type=float, default=0.0); ap.add_argument('--tails', default='free')
    ap.add_argument('--tsw', type=float, default=0.45); ap.add_argument('--tend', type=float, default=1e9)
    a = ap.parse_args()
    o = build(a.name, a.tag, a.H, a.fric, a.Gth, {'free': True, 'blocked': False, 'band': 'band'}[a.tails], a.tsw, a.tend)
    print(json.dumps(dict(path=o['path'], H_i_GeV=o['H_i_GeV'], reheating=o['reheating'], waterfall=o['waterfall']), indent=1))
