"""st_rerun sensitivity runs WITH CODEX'S OWN CODE (reproduction-type, not independent).

Imports Codex b73de64 timing_common (build/integrate/responses/crossings) from the extracted tree and
re-integrates the quiet histories:
  (1) on a dense time grid (dtau = 0.01 everywhere, plus Codex's early samples) at Codex's rtol 2e-10 and at 2e-12;
  (2) on a dense k grid (1.5..3.0 step 0.02) with the dense time grid, at rtol 2e-10.
Crossings use Codex's crossings() unchanged; the local-charge normalisation uses Codex's published
late_NX_normalization_through80 so only the sampling changes.
Nothing is written into the Codex tree: raw files go to GRUT_TIMING_WORKDIR (our scratch).

Usage: python3 -B dense_crossings.py <codex_code_dir> <published_RESPONSE_RESULTS.json> <out_dir> <job>
  job in: time_b3_r10 time_b5_r10 time_b3_r12 time_b5_r12 kgrid_b3 kgrid_b5
"""
import sys, json, time
from pathlib import Path
codex, published, outdir, job = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4]
sys.path.insert(0, codex)
from timing_common import build, integrate, responses, crossings, KS, M, WORK, np  # noqa: E402

outdir.mkdir(parents=True, exist_ok=True)
pub = {r['tag']: r for r in json.loads(published.read_text())['rows']}
early = np.load(M.ROOT / 'MODE_PROFILES.npz')['times']
_d = np.arange(80001) / 100.0
_d = _d[np.min(np.abs(_d[:, None] - early[None, :]), axis=1) > 1e-9]  # C integrator refuses dt<1e-14
times = np.unique(np.r_[early, _d])
SEEDS = [.03, 1e-6, 1e-9]


def crossing_table(times, a, amp, charge, ks):
    rows = []
    for s in SEEDS:
        rows.append(dict(phase_seed_radians=s,
                         amplitude_10percent=crossings(times, a, s * amp, ks, .1),
                         amplitude_unity=crossings(times, a, s * amp, ks, 1.),
                         local_charge_unity=crossings(times, a, s * charge, ks, 1.)))
    return rows


def run(tag, beta, decay, ks, rtol, label):
    exe = build(); t0 = time.monotonic()
    bg, vec, ctl = integrate(exe, label, beta, decay, ks, times, e=1e-7, rtol=rtol)
    amp, dq, nx = responses(beta, bg, vec)
    norm = pub[tag]['late_NX_normalization_through80']
    charge = dq / abs(norm)
    a = np.exp(bg[:, 15])
    (WORK / (label + '.bin')).unlink()  # raw 100+ MB; regenerable
    return dict(times=times, a=a, amp=amp, charge=charge, ks=np.asarray(ks), ctl=ctl,
                seconds=time.monotonic() - t0, bg_last=bg[-1], nx=nx)


def main():
    beta, decay, tag = (3., .5, 'quiet3') if 'b3' in job else (5., 1., 'quiet5')
    if job.startswith('time'):
        rtol = 2e-10 if job.endswith('r10') else 2e-12
        r = run(tag, beta, decay, KS, rtol, 'dense_' + job)
        table = crossing_table(r['times'], r['a'], r['amp'], r['charge'], r['ks'])
        # Same sampling as Codex for a like-for-like check of the integrator alone.
        coarse = np.isin(r['times'], np.unique(np.r_[early, np.arange(80., 800.00001, .5)]))
        tc = crossing_table(r['times'][coarse], r['a'][coarse], r['amp'][coarse], r['charge'][coarse], r['ks'])
        np.savez_compressed(outdir / (job + '.npz'), times=r['times'], a=r['a'], ks=r['ks'],
                            amp=r['amp'].astype(np.float64), charge=r['charge'], nx=r['nx'])
        res = dict(job=job, tag=tag, rtol=rtol, n_times=int(len(r['times'])), check=r['ctl']['check'],
                   seconds=r['seconds'], dense_crossings=table, codex_grid_crossings_this_run=tc,
                   final_a=float(r['a'][-1]))
    else:
        ks = np.round(np.arange(1.5, 3.00001, .02), 4)
        mx_amp = None; mx_chg = None; arg_amp = None; arg_chg = None; checks = []
        for start in range(0, len(ks), 19):
            kk = ks[start:start + 19]
            r = run(tag, beta, decay, kk, 2e-10, f'{job}_{start}')
            checks.append(r['ctl']['check'])
            A = np.abs(r['amp']); Q = np.abs(r['charge'])
            ia = A.argmax(1); iq = Q.argmax(1)
            ma = A[np.arange(len(A)), ia]; mq = Q[np.arange(len(Q)), iq]
            if mx_amp is None:
                mx_amp, mx_chg = ma, mq; arg_amp, arg_chg = kk[ia], kk[iq]
            else:
                ua = ma > mx_amp; uq = mq > mx_chg
                mx_amp = np.where(ua, ma, mx_amp); arg_amp = np.where(ua, kk[ia], arg_amp)
                mx_chg = np.where(uq, mq, mx_chg); arg_chg = np.where(uq, kk[iq], arg_chg)
            a = r['a']
        table = []
        for s in SEEDS:
            row = dict(phase_seed_radians=s)
            for name, series, arg, lev in [('amplitude_10percent', mx_amp, arg_amp, .1), ('amplitude_unity', mx_amp, arg_amp, 1.),
                                           ('local_charge_unity', mx_chg, arg_chg, 1.)]:
                w = np.flatnonzero(s * series >= lev)
                if not len(w):
                    row[name] = None; continue
                i = int(w[0])
                row[name] = dict(level=lev, k=float(arg[i]), tau_bracket=[float(times[i - 1]), float(times[i])],
                                 a_bracket=[float(a[i - 1]), float(a[i])], sampled_value=float(s * series[i]))
            table.append(row)
        np.savez_compressed(outdir / (job + '.npz'), times=times, a=a, ks=ks, max_amp=mx_amp, argk_amp=arg_amp,
                            max_charge=mx_chg, argk_charge=arg_chg)
        res = dict(job=job, tag=tag, rtol=2e-10, ks=ks, checks=checks, dense_k_dense_time_crossings=table)
    (outdir / (job + '.json')).write_text(json.dumps(res, indent=1, default=lambda x: x.tolist() if hasattr(x, 'tolist') else float(x)) + '\n')
    print(job, json.dumps(res.get('dense_crossings', res.get('dense_k_dense_time_crossings')), default=str)[:3000], flush=True)


if __name__ == '__main__':
    main()
