"""Phase-sensitive linear to lattice match at the physical ALP hierarchy."""
from dataclasses import dataclass, replace
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "alp_conversion"))
from lattice import Box, gaussian
from model import MASS_EV, floquet, matching, physical_scales, rotation_energy, rotation_speed, write_json

DIRECTORY = HERE / "results"
START = .49
LEVELS = (1e-4, 1e-3)
K_MIN, K_MAX = 2., 3.
H = physical_scales()["H_star_over_m"]
VACUUM_SCALE = MASS_EV*1e-9/matching(2.)["f_a_gev"]
K_UNIT = 2*np.pi/96
WINDOW = 1500.


@dataclass(frozen=True)
class Grid:
    n: int = 48
    length: float = 24.
    dt: float = .05
    seam: int = 0
    seed: int = 314


BASE = Grid()
CASES = dict(base=BASE, spacing=replace(BASE, n=64),
             volume=replace(BASE, n=64, length=32.),
             time=replace(BASE, dt=.025), later=replace(BASE, seam=1))


def labels(box):
    scale = round(96/box.length)
    axes = tuple(scale*ix for ix in box.indices)
    mode_id = sum(ix*ix for ix in axes)
    keep = box.mask & (box.k >= K_MIN) & (box.k <= K_MAX)
    return axes, mode_id, keep


def measure(box, ids):
    _, index, keep = labels(box)
    multiplicity = np.broadcast_to(box.multiplicity, box.k.shape)[keep]
    weights = np.bincount(np.searchsorted(ids, index[keep]), weights=multiplicity,
                          minlength=len(ids))
    return weights


def modes():
    ids = []
    for grid in (CASES["spacing"], CASES["volume"]):
        _, index, keep = labels(Box(grid.n, grid.length))
        ids.extend(np.unique(index[keep]))
    return np.unique(ids)


def conformal_covariance(state, k2, a):
    x1, p1, x2, p2 = state[2:].reshape(4, -1)
    omega = np.sqrt(k2+a*a)
    u1, u2 = omega*x1, omega*x2
    v1, v2 = p1-H*x1/a, p2-H*x2/a
    return np.stack([u1*u1+u2*u2, v1*v1+v2*v2, u1*v1+u2*v2])


def linear(ids, rtol, sample=None):
    k2 = K_UNIT**2*ids
    w0 = np.sqrt(k2+START**2)
    initial = np.r_[0., START*rotation_speed(START), 1/np.sqrt(w0), np.zeros(len(ids)),
                    np.zeros(len(ids)), np.sqrt(w0)]
    box = Box(BASE.n, BASE.length)
    weights = measure(box, ids)
    prefactor = .25*VACUUM_SCALE**2/box.length**3

    def rhs(t, y):
        a = START+H*t
        angle, v = y[:2]
        x1, p1, x2, p2 = y[2:].reshape(4, -1)
        rate2 = k2+a*a*np.cos(angle)
        return np.r_[v, -2*H*v/a-a*a*np.sin(angle), p1, -rate2*x1, p2, -rate2*x2]

    def crossing(level):
        def event(t, y):
            a = START+H*t
            rho0 = .5*(y[1]/a)**2+2*np.sin(y[0]/2)**2
            c = conformal_covariance(y, k2, a)
            fraction = prefactor*np.dot(weights, c[0]+c[1])/(a**4*rho0)
            return fraction-level
        event.terminal = level == LEVELS[-1]
        event.direction = 1
        return event

    end = (.537-START)/H if sample is None else float(sample[-1])
    out = solve_ivp(rhs, (0., end), initial, t_eval=[end] if sample is None else sample,
                    events=[crossing(level) for level in LEVELS] if sample is None else None,
                    method="DOP853", rtol=rtol, atol=rtol/25, max_step=.5)
    if not out.success:
        raise RuntimeError(out.message)
    if sample is not None:
        return out.y.T, out.nfev
    if any(len(times) == 0 for times in out.t_events):
        raise RuntimeError("A registered linear handoff was not reached.")
    return np.array([t[0] for t in out.t_events]), np.array([y[0] for y in out.y_events]), out.nfev


def swept_tails(seams):
    a = np.unique(np.r_[np.linspace(.25, max(seams), 6145), START, seams])
    k = np.linspace(.30, 8., 1801)
    r = np.array([rotation_energy(v) for v in a])
    rate = floquet(k[None, :]/a[:, None], r[:, None])
    gain = cumulative_trapezoid(2*rate*a[:, None]/H, a, axis=0, initial=0.)
    rows = []
    for at in (START, *seams):
        j = np.searchsorted(a, at)
        exponent = gain[j]
        enhancement = .5*np.exp(exponent)*(-np.expm1(-exponent))**2
        spectrum = (VACUUM_SCALE**2/(4*np.pi**2*a[j]**3)*
                    k*k*np.sqrt(1+(k/a[j])**2)*enhancement)
        total = np.trapezoid(spectrum, k)
        outside = (k < K_MIN) | (k > K_MAX)
        omitted = np.trapezoid(spectrum*outside, k)
        rows.append(dict(a=float(a[j]), fluctuation_fraction=float(total/(2*r[j])),
                         omitted_fraction=float(omitted/total)))
    return rows


def initial_state(box, ids, state, a, seed):
    axes, mode_id, keep = labels(box)
    x1, p1, x2, p2 = state[2:].reshape(4, -1)
    indices = np.searchsorted(ids, mode_id[keep])
    draws = gaussian(*axes, seed)[keep], gaussian(*axes, seed+104729)[keep]
    norm = VACUUM_SCALE/np.sqrt(2*box.length**3)
    chi = np.zeros(box.k.shape, dtype=complex)
    speed = np.zeros_like(chi)
    chi[keep] = norm*(x1[indices]*draws[0]+x2[indices]*draws[1])
    speed[keep] = norm*(p1[indices]*draws[0]+p2[indices]*draws[1])
    angle = state[0]-2*np.pi*np.rint(state[0]/(2*np.pi))
    chi[0, 0, 0] = a*angle
    speed[0, 0, 0] = a*state[1]+H*angle
    return chi, speed


def evaluate():
    ids = modes()
    times, reference, evaluations = linear(ids, 3e-13)
    comparison, eval_other = linear(ids, 3e-11, times)
    scale = START+H*times
    k2 = K_UNIT**2*ids
    weights = measure(Box(BASE.n, BASE.length), ids)
    rows = []
    for a, first, second in zip(scale, reference, comparison):
        c1 = conformal_covariance(first, k2, a)
        c2 = conformal_covariance(second, k2, a)
        delta = c2-c1
        norm = np.sqrt(delta[0]**2+delta[1]**2+2*delta[2]**2)
        error = np.dot(weights, norm)/np.dot(weights, c1[0]+c1[1])
        grids = {}
        for name, config in CASES.items():
            if name in ("time", "later"):
                continue
            box = Box(config.n, config.length)
            w = measure(box, ids)
            e = c1[0]+c1[1]
            chi, speed = initial_state(box, ids, first, a, config.seed)
            theta = box.to_real(chi)/a
            _, density, _, fluct = box.energy(chi, speed, theta, a, H)
            grids[name] = dict(mode_participation=float(np.dot(w, e)**2/np.dot(w, e*e)),
                               expected_energy_over_base=float(
                                   np.dot(w, e)/box.length**3/(np.dot(weights, e)/BASE.length**3)),
                               realized_gradient_and_fluctuation_kinetic_fraction=float(fluct),
                               density=float(density))
        rows.append(dict(a=float(a), covariance_relative_error=float(error), grids=grids,
                         wrapped_phase_difference=float(
                             np.arctan2(np.sin(second[0]-first[0]), np.cos(second[0]-first[0])))))
    tails = swept_tails(scale)
    checks = dict(earlier=tails[0]["fluctuation_fraction"] < LEVELS[0]*1e-3,
                  band=max(x["omitted_fraction"] for x in tails[1:]) < .01,
                  covariance=max(row["covariance_relative_error"] for row in rows) < .02)
    result = dict(passed=all(checks.values()), checks=checks, h=H, vacuum_scale=VACUUM_SCALE,
                  mode_shells=len(ids), ode_evaluations=[evaluations, eval_other],
                  seams=rows, swept_tails=tails, levels=LEVELS)
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(DIRECTORY / "matched_states.npz", ids=ids, times=times, states=reference)
    write_json(DIRECTORY / "match.json", result)
    return result


if __name__ == "__main__":
    result = evaluate()
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit("The registered matched-state gate missed.")
