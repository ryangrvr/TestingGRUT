#!/usr/bin/env python3
"""inf_potential (Claude, own code): checks of the engine copy rpinf.c.
 (1) MODE 0 regression: rpinf and rpdq (dq_entropy, unmodified) on the M1r prepared state (lambda = 0.5, u = 10i, R = 12,
     (m_I, Gamma_I) = (3, 0.5), G from the slice constraint) must give byte-identical output.
 (2) MODE 4: the benchmark inflaton started at the slice S1 with the static preparation at u = 10i; the engine's H(tau)
     and phihat(tau) are compared with our Python background (inf_common) over the last e-fold of inflation and the first
     oscillations (the hidden sector's gravity is ~1e-10 of the total and is ignored in the Python reference).
Read-only import of dq_heat/skeptic/code/dq_model_copy.py for the stationary state.  Binaries are built in the scratch
directory.  Output: ../results/engine_test.json.
"""
import json, math, os, subprocess, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
TRK = HERE.parents[1]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(TRK / 'dq_heat' / 'skeptic' / 'code'))
import inf_common as ic
import dq_model_copy as fc
from exit_factor import model_from_bench

SP = Path(os.environ.get('INF_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/inf_potential'))
B = SP / 'build'
RES = HERE.parent / 'results'


def build():
    B.mkdir(parents=True, exist_ok=True)
    subprocess.run(['gcc', '-O2', '-o', str(B / 'rpinf'), str(HERE / 'rpinf.c'), '-lm'], check=True)
    subprocess.run(['gcc', '-O2', '-o', str(B / 'rpdq'), str(TRK / 'dq_entropy' / 'code' / 'rpdq.c'), '-lm'], check=True)


def spec_args(z):
    out = ['F']
    for c in z:
        out += ['%.17g' % c.real, '%.17g' % c.imag]
    out += ['0'] * 6
    return out


def common_args(mode, b, gam, eta2, p, G, end, nsamp):
    return [str(mode), '%.17g' % b, '%.17g' % gam, '%.17g' % p['xi'], '%.17g' % eta2, '0', '%.17g' % p['lx'], '%.17g' % p['ly'],
            '%.17g' % p['ls'], '%.17g' % p['kap'], '%.17g' % p['M'], '%.17g' % p['h'], '%.17g' % p['mu'], '0', '0',
            '%.17g' % G.real, '%.17g' % G.imag, '0', '%.17g' % end, '1e-10', '1e-13', str(nsamp), '0', '0']


def parse_S(txt):
    rows = [list(map(float, l.split()[1:])) for l in txt.splitlines() if l.startswith('S ')]
    return np.array(rows)


def main():
    build()
    p = fc.params(0.5)
    st = fc.stationary(fc.U_PREP, p, R=12.0)
    z = st['z']
    out = {}
    # (1) MODE 0 regression
    eta2 = fc.ETA2
    I0 = math.sqrt(6.0) / 3.0
    G = (fc.U_PREP - p['mu']) / (I0 / math.sqrt(eta2))
    args = common_args(0, 3.0, 0.5, eta2, p, G, 20.0, 200) + spec_args(z)
    r0 = subprocess.run([str(B / 'rpdq')] + args, capture_output=True, text=True).stdout
    r1 = subprocess.run([str(B / 'rpinf')] + args, capture_output=True, text=True).stdout
    S0 = parse_S(r0)
    out['mode0_identical'] = (r0 == r1)
    out['mode0_lines'] = len(r0.splitlines())
    out['mode0_final_NX_EX'] = [float(S0[-1][4]), float(S0[-1][8])]
    print('MODE 0 regression: identical =', r0 == r1, '| lines', len(r0.splitlines()), '| NX(20) =', S0[-1][4], flush=True)
    # (2) MODE 4 on the benchmarks
    out['mode4'] = {}
    for name in ('HT6t', 'T1', 'E1'):
        m, rec = model_from_bench(name)
        Hi = rec['H_i_GeV']; e2 = (Hi / ic.MPL_GEV) ** 2
        V0 = rec['V0_MP4']
        kind = {'T': 1, 'E': 2, 'HT': 3}[m.kind]
        env = dict(os.environ, RPINF_KIND=str(kind), RPINF_AMP='%.17g' % (V0 / e2), RPINF_PHIMIN='%.17g' % m.phi_min(),
                   RPINF_I0='%.17g' % rec['phiI_S1_MP'], RPINF_SIGN='%d' % (1 if rec['phidot_S1_over_HiMP'] > 0 else -1))
        if m.kind == 'HT':
            env.update(RPINF_P1='%.17g' % m.vev, RPINF_P2='%d' % m.pw, RPINF_P3='%.17g' % m.beta, RPINF_FTOP='%.17g' % m.ftop)
        else:
            env.update(RPINF_P1='%.17g' % m.M)
        Gb = (fc.U_PREP - p['mu']) / (rec['phiI_S1_MP'] / math.sqrt(e2))
        tend = 3.0
        args = common_args(4, rec['m_I_over_Hi'], 0.0, e2, p, Gb, tend, 300) + spec_args(z)
        res = subprocess.run([str(B / 'rpinf')] + args, capture_output=True, text=True, env=env)
        if res.returncode != 0:
            print(name, 'engine error', res.stderr[:300]); out['mode4'][name] = dict(error=res.stderr[:300]); continue
        S = parse_S(res.stdout)
        tau, a, H, ph, ure, uim = S[:, 0], S[:, 1], S[:, 2], S[:, 28], S[:, 26], S[:, 27]
        # Python reference from S1 (gamma = 0): post_inflation starting 1 e-fold before the end
        bg = ic.background(m, 75.0)
        hs = rec['h_slice'] if 'h_slice' in rec else Hi / ic.MPL_GEV / math.sqrt(V0)
        sol = ic.post_inflation(m, bg, 0.0, tend / hs, start_Nbefore=1.0, rtol=1e-11, atol=1e-15, dense=True)
        Y = sol.sol(tau / hs)
        rhoP = 0.5 * Y[2] ** 2 + np.array([m.v(x) for x in Y[1]])
        HP = np.sqrt(rhoP / 3) / hs
        phP = (Y[1] - m.phi_min()) / math.sqrt(e2)
        dH = float(np.max(np.abs(H / HP - 1)))
        dph = float(np.max(np.abs(ph - phP)) / np.max(np.abs(phP)))
        aend = None
        out['mode4'][name] = dict(max_rel_dH=dH, max_rel_dphihat=dph, u_at_0=[float(ure[0]), float(uim[0])],
                                  u_at_end=[float(ure[-1]), float(uim[-1])], NX_end=float(S[-1][4]), H_end=float(H[-1]),
                                  a_end=float(a[-1]), table=[dict(tau=float(t), a=float(x), H=float(h), phihat=float(q), u=[float(r), float(i)])
                                                             for t, x, h, q, r, i in zip(tau[::15], a[::15], H[::15], ph[::15], ure[::15], uim[::15])])
        print('%s MODE 4: max |H/H_py - 1| = %.2e, max |dphihat|/max = %.2e; u(0) = %.4f%+.4fi, u(3) = %.4f%+.4fi; N_X(3) = %.4g' % (
            name, dH, dph, ure[0], uim[0], ure[-1], uim[-1], S[-1][4]), flush=True)
    # (3) smoke test: Gamma_I = 0.5 H_i to tau = 300 from the S1 slice; N_X, E_X, eps (NOT a converged release analysis)
    out['smoke'] = {}
    for name in ('HT6t', 'T1', 'E1'):
        m, rec = model_from_bench(name)
        Hi = rec['H_i_GeV']; e2 = (Hi / ic.MPL_GEV) ** 2; V0 = rec['V0_MP4']
        kind = {'T': 1, 'E': 2, 'HT': 3}[m.kind]
        env = dict(os.environ, RPINF_KIND=str(kind), RPINF_AMP='%.17g' % (V0 / e2), RPINF_PHIMIN='%.17g' % m.phi_min(),
                   RPINF_I0='%.17g' % rec['phiI_S1_MP'], RPINF_SIGN='%d' % (1 if rec['phidot_S1_over_HiMP'] > 0 else -1))
        if m.kind == 'HT':
            env.update(RPINF_P1='%.17g' % m.vev, RPINF_P2='%d' % m.pw, RPINF_P3='%.17g' % m.beta, RPINF_FTOP='%.17g' % m.ftop)
        else:
            env.update(RPINF_P1='%.17g' % m.M)
        Gb = (fc.U_PREP - p['mu']) / (rec['phiI_S1_MP'] / math.sqrt(e2))
        args = common_args(4, rec['m_I_over_Hi'], 0.5, e2, p, Gb, 300.0, 300) + spec_args(z)
        res = subprocess.run([str(B / 'rpinf')] + args, capture_output=True, text=True, env=env)
        S = parse_S(res.stdout)
        # a at eps_H = 1 relative to S1: from our background
        bg = ic.background(m, 75.0)
        aend = math.exp(1.0)
        sel = S[:, 0] >= 200
        NX = float(np.mean(S[sel, 4])); EXv = float(np.mean(S[sel, 8])); epsv = float(np.mean(S[sel, 9]))
        out['smoke'][name] = dict(NX_S1norm=NX, EX_S1norm=EXv, eps=epsv, NX_endnorm=NX / aend ** 3, EX_endnorm=EXv / aend ** 4,
                                  NYL=float(np.mean(S[sel, 11])), NYH=float(np.mean(S[sel, 10])),
                                  NX_at=[dict(tau=float(t), NX=float(n)) for t, n in zip(S[::20, 0], S[::20, 4])])
        print('%s smoke (Gamma = 0.5, tau 200-300): N_X = %.5g (S1 norm), %.5g (end norm, x e^-3); E_X = %.5g (S1), %.5g (end, x e^-4); eps = %.4f; NYL = %.4g' % (
            name, NX, NX / aend ** 3, EXv, EXv / aend ** 4, epsv, out['smoke'][name]['NYL']), flush=True)
    json.dump(out, open(RES / 'engine_test.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
