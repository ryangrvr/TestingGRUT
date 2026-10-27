"""st_rerun_skeptic: OWN-ENGINE background check (no Codex code).

Runs Claude's own reheating C engine (tracks/rp_reheating/code/rpreheat.c, mode 0, H from the Friedmann constraint,
Dormand-Prince 5(4)) for the two quiet histories (m_I,Gamma_I)/H_i = (3,0.5), (5,1) to tau = 800 at dtau = 0.01, with
our own preparation (rp.prepare).  Saves tau, a, H, N_X, Gauss, H_aux consistency to scratch npz, and compares with
Codex's published quiet3/quiet5.npz background (read-only) at the shared sample times.
Usage: python3 -B bg_own.py <codex_packet_dir> <scratch_dir> <out_json>
"""
import sys, json, os
from pathlib import Path
import numpy as np
packet, scratch, outp = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
os.environ['RP_SCRATCH'] = str(scratch / 'rp')
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'rp_reheating' / 'code'))
import rp  # noqa: E402  (our own engine)

out = {}
for tag, b, g in [('quiet3', 3., .5), ('quiet5', 5., 1.)]:
    res, S = rp.run('skeptic_' + tag, mode=0, b=b, gam=g, end=800., rtol=1e-11, atol=1e-14, nsamp=80000, save=False)
    C = dict(zip(rp.COLS, S.T))
    tau, a, H, NX = C['tau'], C['a'], C['H'], C['NX']
    np.savez_compressed(scratch / f'bg_{tag}.npz', tau=tau, a=a, H=H, NX=NX, R=C['R'], absX=C['absX'],
                        Gauss=C['Gauss'], Hcheck=C['Hcheck'])
    z = np.load(packet / f'{tag}.npz'); tc = z['times']; ac = np.exp(z['bg'][:, 15]); Hc = z['bg'][:, 16]
    sel = tc >= 80; idx = np.round(tc[sel] * 100).astype(int)
    assert np.allclose(tau[idx], tc[sel], atol=1e-9)
    ra = np.max(np.abs(a[idx] / ac[sel] - 1)); rH = np.max(np.abs(H[idx] / Hc[sel] - 1))
    Cown = a[-1] ** 2 * H[-1]; off = tau[-1] - a[-1] ** 2 / (2 * Cown)
    w = (tau >= 64) & (tau <= 80)
    nx_dense_mean = float(np.trapezoid(NX[w], tau[w]) / 16.)
    # Codex's normalization used a trapezoid on its sparse sample grid inside [64,80]
    sw = (tc >= 64) & (tc <= 80)
    nx_sparse = np.interp(tc[sw], tau, NX)
    nx_sparse_mean = float(np.trapezoid(nx_sparse, tc[sw]) / (tc[sw][-1] - tc[sw][0]))
    out[tag] = dict(seconds=res['seconds'], checks=res['checks'],
                    a_134=float(np.interp(134., tau, a)), a_147p5=float(np.interp(147.5, tau, a)),
                    a_772=float(np.interp(772., tau, a)), a_782=float(np.interp(782., tau, a)), a_800=float(a[-1]),
                    H_800=float(H[-1]), C_own=float(Cown), tau_star_own=float(off),
                    C_relrange_720_800=float(np.ptp((a * a * H)[tau >= 720]) / Cown),
                    max_rel_a_vs_codex_tau80_800=float(ra), max_rel_H_vs_codex_tau80_800=float(rH),
                    NX_mean_64_80_dense=nx_dense_mean, NX_mean_64_80_on_codex_sparse_grid=nx_sparse_mean,
                    n_codex_samples_in_64_80=int(sw.sum()), NX_rel_ptp_64_80=float(np.ptp(NX[w]) / abs(nx_dense_mean)),
                    NX_800_over_dense_mean_minus1=float(NX[-1] / nx_dense_mean - 1))
    print(tag, json.dumps(out[tag], default=float), flush=True)
outp.write_text(json.dumps(out, indent=1, default=float) + '\n')
