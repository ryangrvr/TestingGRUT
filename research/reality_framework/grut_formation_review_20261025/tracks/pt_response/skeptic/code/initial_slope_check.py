"""SKEPTIC (pt_response_skeptic), INDEPENDENT code, no Codex code.
pt_response reads the initial slopes Re G_N(Delta)/Delta = 0.986-0.997 as 'consistent with omega^2 = k^2 + Hartree m^2'.
Check what they measure exactly.  For q'' = -Omega^2(s) q, q(0) = 0, q'(0) = 1:
    q(Delta)/Delta = 1 - Omega0^2 Delta^2/6 - (Omega^2)'(s0) Delta^3/12 + O(Delta^4),
and for the shell-trace estimator with a white spatial source the expectation of Omega0^2 is
    kappa_bar^2(shell) + (1/2) tr Hbar_X(s0),
with kappa^2 the FD4 symbol averaged over the shell and Hbar_X the spatially averaged 2x2 X Hessian (the Hartree m^2
is (1/2) tr Hbar_X up to the g^2 phi^2 term).  So the slope is an O(Delta^3) identity about the INSTANTANEOUS mean
Hessian, not a quasiparticle or dispersion test.
Inputs: Hbar(s) record (pt_response, Codex engine + our diagnostics, Delta s = .01), the published RESPONSE json/npz
(Codex probes, s0 = 20 and 20.44, shells 1-15) and our low-noise traces (five phases, shells 1-4).
Usage: python3 -I initial_slope_check.py <hessian npz> <packet dir> <meas_traces.npz> <out json>
"""
import sys, json
import numpy as np
HZ, P, MEAS, OUT = sys.argv[1], sys.argv[2].rstrip('/') + '/', sys.argv[3], sys.argv[4]
N, L = 48, 12.; dx = L / N; D = .04
h = np.load(HZ); s = h['s']; Hb = h['Hbar']
trX = .5 * (Hb[:, 0, 0] + Hb[:, 1, 1])              # (1/2) tr of the X block
dtrX = np.gradient(trX, s)                            # its s-derivative (Delta s = .01 record)
# shells as in Codex's response.py: sid = floor(|k| / (2 pi / L) + .5), |k| the continuum wavenumber; all modes of the grid
m = np.fft.fftfreq(N, d=1. / N)                        # integers -24..23
mi, mj, mk = np.meshgrid(m, m, m, indexing='ij')
kc = 2 * np.pi / L * np.sqrt(mi ** 2 + mj ** 2 + mk ** 2)
u = 2 * np.pi * np.stack([mi, mj, mk]) / N
fd4 = np.sum(4 * np.sin(u / 2) ** 2 * (1 + np.sin(u / 2) ** 2 / 3), axis=0) / dx ** 2
sid = np.floor(kc / (2 * np.pi / L) + .5).astype(int)
sid[kc > 8. + 1e-12] = -1                             # Codex's mapping keeps |k| <= 8 only (full k-space)
kap2 = {n: float(fd4[sid == n].mean()) for n in range(1, 16)}
kap2_spread = {n: [float(fd4[sid == n].min()), float(fd4[sid == n].max())] for n in range(1, 16)}
nmodes_full = {n: int((sid == n).sum()) for n in range(1, 16)}

def pred(n, s0):
    i = int(round((s0 - s[0]) / .01)); assert abs(s[i] - s0) < 1e-9
    Om2 = kap2[n] + trX[i]
    p3 = 1 - Om2 * D ** 2 / 6 + Om2 ** 2 * D ** 4 / 120
    p4 = p3 - dtrX[i] * D ** 3 / 12
    return Om2, p3, p4, float(trX[i]), float(dtrX[i])

res = dict(kappa2_shell_mean=kap2, kappa2_shell_min_max=kap2_spread, modes_full_space=nmodes_full, published={}, lownoise={})
dev3 = []; dev4 = []
for tag, s0 in (('s20_seed7', 20.), ('s20.44_seed7', 20.44), ('s20_seed11', 20.)):
    J = json.load(open(P + f'RESPONSE_{tag}.json')); nr = np.load(P + f'RESPONSE_{tag}.npz')['normal']
    rows = []
    for r in J['rows']:
        n = r['shell']; meas = float(nr[1, n - 1].real / D)
        assert abs(meas - r['initial_slope']) < 1e-12
        Om2, p3, p4, tr0, dtr0 = pred(n, s0)
        rows.append(dict(shell=n, modes_codex_halfspace=r['modes'], modes_ours_fullspace=nmodes_full[n], measured=meas,
                         predicted_O3=p3, predicted_O4=p4, Omega2=Om2, half_tr_HbarX=tr0, d_half_tr_HbarX_ds=dtr0,
                         dev_O3=meas - p3, dev_O4=meas - p4))
        dev3.append(abs(meas - p3)); dev4.append(abs(meas - p4))
    res['published'][tag] = rows
mz = np.load(MEAS); TN = mz['tr_normal']; ph = mz['phases']
ln3 = []; ln4 = []
for i, s0 in enumerate(ph):
    rows = []
    for n in range(1, 5):
        meas = float(TN[i, n - 1, 1].real / D)
        Om2, p3, p4, tr0, dtr0 = pred(n, float(s0))
        rows.append(dict(shell=n, measured=meas, predicted_O3=p3, predicted_O4=p4, Omega2=Om2, half_tr_HbarX=tr0,
                         d_half_tr_HbarX_ds=dtr0, dev_O3=meas - p3, dev_O4=meas - p4))
        ln3.append(abs(meas - p3)); ln4.append(abs(meas - p4))
    res['lownoise'][f'{s0:g}'] = rows
res['summary'] = dict(published_max_abs_dev_O3=max(dev3), published_max_abs_dev_O4=max(dev4),
                      lownoise_max_abs_dev_O3=max(ln3), lownoise_max_abs_dev_O4=max(ln4),
                      half_tr_HbarX_at_s0={f'{s0:g}': float(trX[int(round((s0 - 20) / .01))]) for s0 in (20., 20.12, 20.24, 20.44, 20.68)},
                      codex_hartree_m2_at_s0={f'{s0:g}': float(h['hf'][int(round((s0 - 20) / .01))]) for s0 in (20., 20.12, 20.24, 20.44, 20.68)},
                      note='half tr HbarX = 2 lam <r^2> + g2 <phi^2> equals Codex Hartree_m2 exactly (same spatial mean); the slope tests the instantaneous mean Hessian at s0 to O(Delta^3), with the O(Delta^4) drift term visible.')
print(json.dumps(res['summary'], indent=1))
for tag, rows in res['published'].items():
    print(tag, 'modes half/full', [(r['modes_codex_halfspace'], r['modes_ours_fullspace']) for r in rows][:5])
    print(tag, 'dev O3 %s' % np.round([r['dev_O3'] for r in rows], 5).tolist())
    print(tag, 'dev O4 %s' % np.round([r['dev_O4'] for r in rows], 5).tolist())
for s0, rows in res['lownoise'].items():
    print('lownoise', s0, 'meas', np.round([r['measured'] for r in rows], 5).tolist(), 'dev O3', np.round([r['dev_O3'] for r in rows], 5).tolist(),
          'dev O4', np.round([r['dev_O4'] for r in rows], 5).tolist())
json.dump(res, open(OUT, 'w'), indent=1)
