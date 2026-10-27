"""Independent integration check of claim 5's background part (0.547%) and initial force mismatch (0.206%).

Own integrator (scipy DOP853 on Codex's full homogeneous reheating RHS, imported read-only from the sibling packet
grut_reheating_preparation_20261026) and own reduced quartic orbit, own apocentre normalisation. The full-model
equations themselves are Codex's (not re-derived here), so this checks the matching numerics and conventions, not
the full model.

usage: python3 -B -I background_match.py <reality_framework_dir> <packet_dir> <out_json>
"""
import sys, json, os
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

RF, PKT, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
sys.path.insert(0, str(RF / 'grut_reheating_preparation_20261026' / 'code'))
import reheating as R  # noqa: E402  (read-only import; GRUT_REHEATING_WORKDIR must point to scratch)

LX = .001


def rv(z):
    return np.sqrt(2) * np.array([z.real, z.imag]).T.ravel()


def light(v):
    a = np.exp(v[15]); z = v[2:8:2] + 1j * v[3:8:2]; vel = v[8:14:2] + 1j * v[9:14:2]
    chi = a * rv(z)[:2]; chip = a * a * rv(vel + v[16] * z)[:2]
    return a, chi, chip


def apocentre(q, p):
    E = .5 * p @ p + LX / 4 * (q @ q) ** 2; Q = q[0] * p[1] - q[1] * p[0]
    f = lambda r: Q * Q / (2 * r * r) + LX / 4 * r ** 4 - E
    rmax = (4 * E / LX) ** .25; rstar = (Q * Q / LX) ** (1 / 6)  # f is minimal (negative) at rstar
    return brentq(f, rstar, rmax)


def main():
    beta, decay = 3., .5
    v0, res = R.prepared(beta)
    taus = np.r_[0., np.round(np.arange(200, 501) * .1, 10)]
    sol = solve_ivp(R.rhs, (0, 50), v0, t_eval=taus, method='DOP853', rtol=1e-12, atol=1e-14,
                    args=(beta, decay, R.XI, R.ETA2))
    assert sol.success, sol.message
    V = sol.y.T
    a20, chi20, chip20 = light(V[1]); C = a20 * a20 * V[1, 16]
    m = np.load(PKT / 'MATCH_L12.npz')
    match_cmp = dict(a_own=a20, a_codex=float(m['a_match']), C_own=C, C_codex=float(m['C_match']),
                     chi_rel=float(np.linalg.norm(chi20 - m['light_mean_q'][:2]) / np.linalg.norm(m['light_mean_q'][:2])),
                     chip_rel=float(np.linalg.norm(chip20 - m['light_mean_p'][:2]) / np.linalg.norm(m['light_mean_p'][:2])))
    # reduced orbit in conformal time s from tau 20
    t = taus[1:]; s = 2 * (t - 20) / (a20 + np.sqrt(a20 * a20 + 2 * C * (t - 20)))
    rhs = lambda _, y: np.r_[y[2:], -LX * (y[:2] @ y[:2]) * y[:2]]
    red = solve_ivp(rhs, (0, s[-1]), np.r_[chi20, chip20], t_eval=s, method='DOP853', rtol=1e-13, atol=1e-12).y.T
    full = np.array([light(v)[1] for v in V[1:]])
    apo = apocentre(chi20, chip20)
    err = np.linalg.norm(red[:, :2] - full, axis=1) / apo
    # also own-clock check: full a(tau) vs a_m + C s
    afull = np.exp(V[1:, 15]); aclock = a20 + C * s
    # initial force mismatch, own evaluation of Codex's full potential force
    z, vel, pot, F, f, fd, Rc, hd = R.geometric(V[1], beta, R.XI, R.ETA2)
    fullforce = rv(F)[:2] * a20 ** 3 + (R.XI - 1 / 6) * a20 * a20 * Rc * chi20
    fbare = LX * (chi20 @ chi20) * chi20
    ferr = float(np.linalg.norm(fullforce - fbare) / np.linalg.norm(fbare))
    out = dict(method='own DOP853 on Codex full homogeneous RHS (reheating.rhs) + own reduced orbit; own apocentre',
               preparation_residual=res, match_vs_MATCH_L12=match_cmp,
               background_orbit_scaled_error_max=float(err.max()), at_tau=float(t[np.argmax(err)]),
               codex_value=0.005470871604787762, apocentre=apo,
               clock_max_relative_a_error_20_50=float(np.max(abs(aclock / afull - 1))),
               initial_force_relative_mismatch=ferr, codex_force=0.0020569577670522584, R_at_match=float(Rc))
    OUT.write_text(json.dumps(out, indent=2) + '\n'); print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
