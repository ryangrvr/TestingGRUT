"""Skeptic runner for pt_m1tilt: one JSON per run in $SP/review/pt_m1tilt_skeptic/runs, restartable (checkpoint
after every grid chunk when the wall-clock budget is exceeded).  At most 2 worker processes.
Usage: python3 -B sk_run.py --budget 500 JOB [JOB ...]
  JOB = orbit:eps:phase[:key=val,...]   keys: growth (0 exp | 1 linear), epst, mend (M_end/Mapo), mstart (M0/Mperi),
        nper (steps per fastest period), method (8|6), rm (linear growth: M/Mperi at which dlnM/ds = eps*omega_r)
Measures (own formulas, sk_common.free_measures): windows of 60 X periods (2pi/(gM)) at 5, 10, 20, 40 Mapo;
a 400-point log grid in M with 4-radial-period short windows (n_-, N_phi, energies); H_ext and charge checks."""
import sys, os, json, time, argparse, numpy as np
from multiprocessing import Pool
from sk_common import *

T_START = time.time()
BUDGET = 500.0


def parse_job(spec):
    parts = spec.split(':')
    j = dict(orbit=parts[0], eps=float(parts[1]), phase=float(parts[2]), growth=0, epst=EPST, mend=40.0, mstart=0.5,
             nper=32, method=8, rm=1.0)
    for kv in parts[3:]:
        if not kv: continue
        for item in kv.split(','):
            k, v = item.split('=')
            j[k] = type(j[k])(float(v)) if k != 'orbit' else v
    return j


def job_name(j):
    n = f"{j['orbit']}_eps{j['eps']:.2e}_ph{j['phase']:.4f}_g{j['growth']}_et{j['epst']:g}_me{j['mend']:g}_ms{j['mstart']:g}_n{j['nper']}_m{j['method']}"
    if j['growth'] == 1:
        n += f"_rm{j['rm']:g}"
    return n


def run_job(j):
    name = job_name(j)
    fn = os.path.join(RUNS, name + '.json'); ck = os.path.join(RUNS, name + '.ckpt.json')
    if os.path.exists(fn):
        return name, 'exists'
    y4, info = orbit(j['orbit'])
    st = phase_state(y4, info['T_r'], j['phase'])
    q = charge(st)
    par = make_par(info, j['eps'], growth=j['growth'], epst=j['epst'], M_start_fac=j['mstart'], method=j['method'], rate_match_M=j['rm'])
    r2_0 = st[0] ** 2 + st[1] ** 2
    chi0 = chi_start(r2_0, par.M0, j['epst'])
    S = Sim(par, np.r_[st[:2], chi0], np.r_[st[2:], 0.0], nper=j['nper'])
    Mend = j['mend'] * info['Mapo']
    grid = np.geomspace(par.M0, Mend, 401)[1:]
    wins = [f * info['Mapo'] for f in (5.0, 10.0, 20.0, 40.0) if f <= j['mend'] + 1e-9]
    events = sorted([(Mg, 'g') for Mg in grid] + [(Mw, 'w') for Mw in wins])
    res = dict(job=j, name=name, q=q, M0=par.M0, kappa=par.kappa, om_fix=par.om_fix, Mperi=info['Mperi'], Mapo=info['Mapo'],
               J0_per_q=info['J_per_q'], omega_r0=info['omega_r'], T_r0=info['T_r'], chi0=chi0, chi0_over_Mperi=chi0 / info['Mperi'],
               trace=[], windows={}, i_event=0, runtime=0.0)
    if os.path.exists(ck):
        with open(ck) as f:
            c = json.load(f)
        res = c['res']; S.restore(c['state'])
    t0 = time.time()
    for i in range(res['i_event'], len(events)):
        Mg, kind = events[i]
        s_target = s_of_M(par, Mg)
        if s_target > S.s:
            S.to(s_target)
        if kind == 'g':
            # short window: 4 radial periods, radial period estimated from sqrt(omega_r0^2 + (2 g M)^2)
            Tr_est = 2 * np.pi / np.sqrt(info['omega_r'] ** 2 + (2 * GG * Mg) ** 2)
            ns = 200; ds = 4 * Tr_est / ns
            y0 = S.y.copy(); Y = S.sample(ds, ns); Y = np.vstack([y0, Y]); ts = Y[:, 6]
            m = free_measures(Y, ts, par, q, j['epst'])
            he, E = S.hext_err()
            qq = Y[:, 0] * Y[:, 4] - Y[:, 1] * Y[:, 3]
            res['trace'].append(dict(M=float(Mg), M_over_Mperi=float(Mg / info['Mperi']), chi_over_M=float(y0[2] / Mg),
                                     nminus=tav(m['nminus'], ts), Nphi=tav(m['Nphi'], ts), EX=tav(m['EX'], ts),
                                     Ephi_lin=tav(m['Ephi_lin'], ts), Ephi_exact=tav(m['Ephi_exact'], ts), Eint=tav(m['Eint'], ts),
                                     Equart=tav(m['Equart'], ts), dchi_rms=float(np.sqrt(np.mean((Y[:, 2] / Mg - 1) ** 2))),
                                     chi_min=float(Y[:, 2].min() / Mg), chi_max=float(Y[:, 2].max() / Mg),
                                     E=float(E), W=float(S.W()), hext_rel=float(he), qrel=float(np.max(np.abs(qq / q - 1))),
                                     s=float(S.s), steps=int(S.steps)))
        else:
            Tw = 60 * 2 * np.pi / (GG * Mg); ns = 2000; ds = Tw / ns
            y0 = S.y.copy(); Y = S.sample(ds, ns); Y = np.vstack([y0, Y]); ts = Y[:, 6]
            m = free_measures(Y, ts, par, q, j['epst'])
            he, E = S.hext_err()
            qq = Y[:, 0] * Y[:, 4] - Y[:, 1] * Y[:, 3]
            f = Mg / info['Mapo']
            res['windows'][f'{f:g}'] = dict(nminus=tav(m['nminus'], ts), nminus_actual=tav(m['nminus_actual'], ts),
                                            Nphi=tav(m['Nphi'], ts), nminus_std=float(np.std(m['nminus'])),
                                            EX_over_mXq=tav(m['EX'], ts) / (m['mX'][0] * q), Ephi_lin_over_mphiq=tav(m['Ephi_lin'], ts) / (m['mphi'][0] * q),
                                            Ephi_exact_over_mphiq=tav(m['Ephi_exact'], ts) / (m['mphi'][0] * q),
                                            Eint_over_mXq=tav(m['Eint'], ts) / (m['mX'][0] * q), Equart_over_EX=tav(m['Equart'], ts) / tav(m['EX'], ts),
                                            dchi_rms_over_M=float(np.sqrt(np.mean((Y[:, 2] / Mg - 1) ** 2))),
                                            E=float(E), W=float(S.W()), hext_rel=float(he), qrel=float(np.max(np.abs(qq / q - 1))),
                                            steps=int(S.steps), s=float(S.s))
        res['i_event'] = i + 1
        if time.time() - T_START > BUDGET:
            res['runtime'] += time.time() - t0
            with open(ck, 'w') as f:
                json.dump(dict(res=res, state=dict(y=S.y.tolist(), E0=S.E0, steps=S.steps)), f)
            return name, f'partial {i + 1}/{len(events)}'
    res['runtime'] += time.time() - t0
    fin = res['windows'][f"{j['mend']:g}"]
    res.update(nminus=fin['nminus'], Nphi=fin['Nphi'], steps=S.steps, hext_rel_final=fin['hext_rel'], qrel_final=fin['qrel'],
               max_abs_hext=max(abs(t['hext_rel']) for t in res['trace']), max_qrel=max(t['qrel'] for t in res['trace']))
    del res['i_event']
    with open(fn, 'w') as f:
        json.dump(res, f)
    if os.path.exists(ck):
        os.remove(ck)
    return name, f"done n-={fin['nminus']:.5f} Nphi={fin['Nphi']:.5f} steps={S.steps} t={res['runtime']:.0f}s"


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--budget', type=float, default=500.0); ap.add_argument('--nproc', type=int, default=2)
    ap.add_argument('jobs', nargs='+')
    A = ap.parse_args(); BUDGET = A.budget
    os.makedirs(RUNS, exist_ok=True)
    jobs = [parse_job(s) for s in A.jobs]
    with Pool(min(A.nproc, 2)) as p:
        for name, status in p.imap_unordered(run_job, jobs):
            print(status, name, flush=True)
