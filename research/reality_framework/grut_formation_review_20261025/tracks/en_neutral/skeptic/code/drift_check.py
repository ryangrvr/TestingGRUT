"""en_neutral_skeptic item 1b: is it the band DRIFT of the depleting orbit (not the lower Floquet rate) that limits
growth? Same-seed linear delta-phi on a free quartic orbit held at FIXED depleted (E, Q) for the whole run:
(i) the late pericentre values (en_neutral V1 convention, last measured period), (ii) the late period-mean values.
If these give large growth while the slowly depleting orbit (factorial C/Cp) does not, the suppression comes from
the band moving through k-space faster than modes can grow. OWN code. usage: python3 -I drift_check.py <rundir> <out>"""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sk_lib as K
rd, out = Path(sys.argv[1]), Path(sys.argv[2])
tr = K.read_trace(rd / 'state.tsv'); rho0 = float(tr['rho'][0]); s_end = float(tr['eta'][-1])
seed = K.phi_seed_by_k2(rd / 'initial.bin'); kl = np.sqrt(seed['k2'])
rb2 = tr['q0'] ** 2 + tr['q1'] ** 2; s = tr['eta']
E0m = 0.5 * (tr['p0'] ** 2 + tr['p1'] ** 2) + K.LX / 4 * rb2 ** 2; Q0m = tr['q0'] * tr['p1'] - tr['q1'] * tr['p0']
per = np.flatnonzero((rb2[1:-1] < rb2[:-2]) & (rb2[1:-1] <= rb2[2:])) + 1
a, b = per[-2], per[-1]
cases = dict(late_pericentre=(float(E0m[a]), float(Q0m[a])),
             late_period_mean=(float(np.trapezoid(E0m[a:b + 1], s[a:b + 1]) / (s[b] - s[a])), float(np.trapezoid(Q0m[a:b + 1], s[a:b + 1]) / (s[b] - s[a]))))
res = dict(label='OWN fixed depleted orbit controls')
sg = np.linspace(0, s_end, 30001); sL, sR = sg[:-1], sg[1:]; sM = 0.5 * (sL + sR)
for name, (E, Q) in cases.items():
    sol, T, rmin, rmax = K.one_free_period(E, Q)
    r = lambda x: (lambda y: y[0] ** 2 + y[1] ** 2)(sol.sol(np.mod(x, T)))
    rec = K.propagate(seed['k2'], K.G2 * r(sL), K.G2 * r(sM), K.G2 * r(sR), sg, [len(sg) - 1])
    Tm = rec[len(sg) - 1]
    e = K.energy_from_T(Tm, seed, K.G2 * r(np.array([s_end]))[0])
    sv = np.linalg.svd(Tm, compute_uv=False)[:, 0]; i = int(np.argmax(sv))
    res[name] = dict(E_over_rho0=E / rho0, Q=Q, rmax_over_rmin=rmax / rmin, period=T, phi_plus_int_tau800=float(e.sum() / rho0),
                     max_ln_sigma=float(np.log(sv[i])), k_at_max=float(kl[i]))
    print(name, res[name], flush=True)
out.write_text(json.dumps(res, indent=1) + '\n')
