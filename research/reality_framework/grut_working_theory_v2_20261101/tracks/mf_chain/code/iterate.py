#!/usr/bin/env python3
"""mf_chain: self-consistent H_i calibration of a hybrid benchmark (V0 is free in the hybrid).

Step n: build the history at design scale H_n (hist_build.py, rescaling rule in its header), run the whole chain
(run_bench.py -> mf_ideal's engine), read H_req,n = the H_i the dark matter needs at e = EREF (chain_mfi's calibration:
Y ~ H_i^{3/2} at fixed dimensionless history, relic M from the inherited p(S/Q) map or the Q-matter table).  Update
ln H by a damped secant on f(ln H) = ln(H_req/H_i) until |f| < tol.  The first step is the plain substitution
H_1 = H_req,0.  Every step is recorded (results/iterate_<tag>.json).

Usage: python3 -B iterate.py NAME TAG [--fric pauli|track] [--Gth 0] [--tails free|blocked|band] [--tsw 0.45 or 0.6]
                              [--H0 1.3e11] [--e 1e-6] [--tol 0.01] [--nmax 5]
Labels: [C] on the [E] rates; the rescaling rule and the calibration recipe are model choices.
"""
import argparse, json, math, subprocess, sys, os
import mfc_common as M
import hist_build as hb


def step(name, tag, H, a):
    pf = {'free': True, 'blocked': False, 'band': 'band'}[a.tails]
    o = hb.build(name, tag, H, a.fric, a.Gth, pf, a.tsw, a.tend)
    env = dict(os.environ, OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    r = subprocess.run([sys.executable, '-B', str(M.HERE / 'run_bench.py'), o['path'], tag], capture_output=True, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-2000:])
    ch = json.load(open(M.ENG / 'chain' / ('chain_%s.json' % tag)))
    lab = 'e1e-6' if abs(a.e - 1e-6) < 1e-12 else 'e3e-6'
    Hi = o['H_i_GeV']
    Hreq = ch[lab]['H_i_req_GeV']
    return dict(tag=tag, H_design_GeV=H, H_i_GeV=Hi, H_req_GeV=Hreq, ratio=Hreq / Hi, f=math.log(Hreq / Hi), K=ch['K'], eps=ch['eps'],
                N_X=ch['N_X'], N_YL=ch['N_YL'], SQ_1e6=ch['e1e-6']['S_per_Q'], SQ_3e6=ch['e3e-6']['S_per_Q'],
                branch_1e6=ch['e1e-6']['branch'], branch_3e6=ch['e3e-6']['branch'], M_1e6=ch['e1e-6']['M'], M_3e6=ch['e3e-6']['M'],
                SQ_inf=ch['S_per_Q_e_inf'], e_crit=ch['e_crit_bulk'], reheating=o['reheating'], waterfall=o['waterfall'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('name'); ap.add_argument('tag')
    ap.add_argument('--fric', default='pauli'); ap.add_argument('--Gth', type=float, default=0.0)
    ap.add_argument('--tails', default='free'); ap.add_argument('--tsw', type=float, default=0.45)
    ap.add_argument('--tend', type=float, default=1e9); ap.add_argument('--H0', type=float, default=M.H_REF)
    ap.add_argument('--e', type=float, default=1e-6); ap.add_argument('--tol', type=float, default=0.01)
    ap.add_argument('--nmax', type=int, default=5); ap.add_argument('--H1', type=float)
    a = ap.parse_args()
    out = dict(name=a.name, tag=a.tag, args=vars(a), steps=[])
    path = M.RES / ('iterate_%s.json' % a.tag)
    H = a.H0
    hist = []
    for n in range(a.nmax):
        s = step(a.name, '%s_it%d' % (a.tag, n), H, a)
        hist.append(s)
        out['steps'].append(s)
        M.dump(path, out)
        print('step %d: H_design %.4e -> H_i %.4e, K %.3f, S/Q(1e-6) %.2f %s, S/Q(3e-6) %.2f, M %.3f, H_req %.4e, ratio %.4f' % (
            n, H, s['H_i_GeV'], s['K'], s['SQ_1e6'], s['branch_1e6'], s['SQ_3e6'], s['M_1e6'], s['H_req_GeV'], s['ratio']), flush=True)
        if abs(s['f']) < a.tol:
            out['converged'] = True
            break
        if len(hist) == 1:
            lnH = math.log(a.H1) if a.H1 else math.log(s['H_req_GeV'])
        else:
            x0, x1 = math.log(hist[-2]['H_i_GeV']), math.log(hist[-1]['H_i_GeV'])
            f0, f1 = hist[-2]['f'], hist[-1]['f']
            slope = (f1 - f0) / (x1 - x0) if x1 != x0 else -1.0
            lnH = x1 - f1 / slope if slope < -0.05 else x1 + 0.7 * f1
        H = math.exp(lnH) * (H / hist[-1]['H_i_GeV'])          # design scale such that H_i lands at the target
    else:
        out['converged'] = False
    if len(hist) >= 2:
        x = [math.log(h['H_i_GeV']) for h in hist[-2:]]; f = [h['f'] for h in hist[-2:]]
        out['secant_fixed_point_GeV'] = math.exp(x[1] - f[1] * (x[1] - x[0]) / (f[1] - f[0])) if f[1] != f[0] else None
        out['contraction'] = (f[1] / f[0]) if f[0] else None
    out['final'] = hist[-1]
    M.dump(path, out)


if __name__ == '__main__':
    main()
