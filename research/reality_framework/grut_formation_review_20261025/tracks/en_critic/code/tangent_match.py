"""Independent integration of claim 5's charged-tangent (~1%) and neutral-transfer (1.5-2.7%) matching errors.

Own scipy DOP853 integration (tau 0->50) of the full background plus the 12-component preparation-phase tangent and the
2x2 phi transfer, using Codex's mode generator / neutral mass functions imported read-only from the sibling packets,
then own reduced (e=0 scalar-limit) tangent/transfer from tau 20. Same diagnostics as matching_check.py.

usage: GRUT_MODE_WORKDIR=... GRUT_REHEATING_WORKDIR=... python3 -B -I tangent_match.py <reality_framework_dir> <out_json>
"""
import sys, json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

RF, OUT = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(RF / 'grut_reheating_preparation_20261026' / 'code'))
sys.path.insert(0, str(RF / 'grut_reheating_modes_20261026' / 'code'))
import reheating as R  # noqa: E402
import modes as MO  # noqa: E402
import neutral as NE  # noqa: E402

LX, G2 = .001, .01


def main():
    beta, decay, ks = 3., .5, [2.13, 7.1]
    v0, _ = R.prepared(beta)
    plus, _ = R.prepared(beta, delta=1e-5); minus, _ = R.prepared(beta, delta=-1e-5)
    q = (plus[2:8] - minus[2:8]) / 2e-5 * np.sqrt(2); seed = np.r_[q, q]
    taus = np.r_[0., np.round(np.arange(200, 501) * .1, 10)]
    rows = []
    for k in ks:
        def rhs(t, y):
            v = y[:19]; out = np.empty_like(y)
            out[:19] = R.rhs(t, v, beta, decay, R.XI, R.ETA2)
            out[19:31] = MO.generator(v, beta, k, e=0.) @ y[19:31]
            a = np.exp(v[15]); m = NE.mass(v, beta); T = y[31:35].reshape(2, 2)
            out[31:35] = np.vstack([T[1] / a, -(k * k + m) * T[0] / a]).ravel()
            return out
        y0 = np.r_[v0, seed, np.eye(2).ravel()]
        sol = solve_ivp(rhs, (0, 50), y0, t_eval=taus, method='DOP853', rtol=1e-12, atol=1e-13)
        assert sol.success, sol.message
        Y = sol.y.T; v20 = Y[1, :19]
        a20 = np.exp(v20[15]); C = a20 * a20 * v20[16]
        z = v20[2:8:2] + 1j * v20[3:8:2]; vel = v20[8:14:2] + 1j * v20[9:14:2]
        chi = a20 * MO.realvec(z)[:2]; chip = a20 * a20 * MO.realvec(vel + v20[16] * z)[:2]
        tan20 = Y[1, 19:31]; T20 = Y[1, 31:35].reshape(2, 2)
        t = taus[1:]; s = 2 * (t - 20) / (a20 + np.sqrt(a20 * a20 + 2 * C * (t - 20)))

        def red(_, y):
            qq = y[:2]; pp = y[2:4]; r2 = qq @ qq; out = np.empty_like(y); out[:4] = np.r_[pp, -LX * r2 * qq]
            d = y[4:8]; M = (k * k + LX * r2) * np.eye(2) + 2 * LX * np.outer(qq, qq)
            out[4:8] = np.r_[d[2:], -M @ d[:2]]
            P = y[8:12].reshape(2, 2); out[8:12] = np.vstack([P[1], -(k * k + G2 * r2) * P[0]]).ravel()
            return out
        init = np.r_[chi, chip, tan20[:2], tan20[6:8], T20.ravel()]
        rs = solve_ivp(red, (0, s[-1]), init, t_eval=s, method='DOP853', rtol=1e-12, atol=1e-13).y.T
        target = Y[1:, 19:21]; approx = rs[:, 4:6]
        err = float(np.max(np.linalg.norm(approx - target, axis=1)) / np.max(np.linalg.norm(target, axis=1)))
        Pf = Y[1:, 31:35]; Pr = rs[:, 8:12]
        pherr = float(np.max(abs(Pr - Pf)) / np.max(abs(Pf)))
        rows.append(dict(k=k, phase_tangent_envelope_relative_error=err, neutral_transfer_relative_error=pherr))
        print(k, err, pherr, flush=True)
    codex = {2.13: (0.009774883339864232, 0.027171575954790524), 7.1: (0.010307871265675655, 0.014820004716620098)}
    for r in rows: r['codex'] = codex[r['k']]
    OUT.write_text(json.dumps(dict(method='own DOP853 on Codex full background + mode generator (e=0) + neutral mass, '
                                          'own reduced tangent from tau 20', rows=rows), indent=2) + '\n')


if __name__ == '__main__':
    main()
