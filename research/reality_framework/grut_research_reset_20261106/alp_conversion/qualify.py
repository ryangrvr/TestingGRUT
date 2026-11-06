"""Homogeneous capture, exact Lamé bands, and a linear 3D wave check."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipj, ellipk

from lattice import Box
from model import HERE, Q, floquet, libration_action, rotation_speed, write_json


def homogeneous(h, phase=0., final=4.):
    a0 = .60
    start = rotation_speed(a0, phase)
    initial_energy = a0**3*(.5*start**2+2*np.sin(phase/2)**2)

    def derivative(t, y):
        theta, speed, _ = y
        a = np.sqrt(a0*a0+2*h*t)
        pressure = .5*speed*speed-2*np.sin(theta/2)**2
        return [speed, -3*h/a**2*speed-np.sin(theta), -3*h*a*pressure]

    a = np.linspace(3.5, final, 501)
    times = (a*a-a0*a0)/(2*h)
    sol = solve_ivp(derivative, (0., times[-1]), [phase, start, 0.],
                    t_eval=times, method="DOP853", rtol=3e-11, atol=1e-12,
                    max_step=.35)
    if not sol.success:
        raise RuntimeError(sol.message)
    theta, speed, work = sol.y
    theta = theta-2*np.pi*np.rint(np.mean(theta)/(2*np.pi))
    rho = .5*speed**2+2*np.sin(theta/2)**2
    exact = np.array([libration_action(b, r/2)/Q for b, r in zip(a, rho)])
    hubble = h/a**2
    frequency = np.sqrt(1+.75*hubble**2)
    free = .5*a**3*((speed+1.5*hubble*theta)**2+frequency**2*theta**2)/frequency/Q
    return dict(h=h, phase=phase, action_C=float(np.mean(exact)),
                particle_C=float(np.mean(free)), action_drift=float(np.ptp(exact)/np.mean(exact)),
                particle_drift=float(np.ptp(free)/np.mean(free)),
                peak_angle=float(max(abs(theta))),
                late_budget_error=float(np.max(abs(a**3*rho-initial_energy-work))/(a[-1]**3*rho[-1])))


def monodromy(p, r, rotating):
    parameter = 1/r if rotating else r
    frequency = np.sqrt(r) if rotating else 1.
    period = 2*ellipk(parameter)/frequency

    def derivative(t, flat):
        s = ellipj(frequency*t, parameter)[0]
        cosine = 1-2*s*s*(1 if rotating else r)
        matrix = np.array([[0., 1.], [-p*p-cosine, 0.]])
        return (matrix @ flat.reshape(2, 2)).ravel()

    sol = solve_ivp(derivative, (0., period), np.eye(2).ravel(), method="DOP853",
                    rtol=1e-12, atol=1e-14)
    matrix = sol.y[:, -1].reshape(2, 2)
    rate = np.arccosh(max(1., abs(np.trace(matrix))/2))/period
    return float(rate), period, matrix


def wave_test(p, r, period, matrix):
    box = Box(12, 2*np.pi/p)
    chi = np.zeros(box.k.shape, dtype=complex)
    speed = chi.copy()
    eig, basis = np.linalg.eig(matrix)
    vector = basis[:, np.argmax(abs(eig))].real
    initial = 1e-7*vector
    chi[0, 0, 1], speed[0, 0, 1] = .5*initial
    speed[0, 0, 0] = 2*np.sqrt(r)
    steps_per_period, periods = 320, 8
    dt = period/steps_per_period
    force, _ = box.acceleration(chi, 1.)
    for _ in range(steps_per_period*periods):
        force, _ = box.integrate_step(chi, speed, force, 1., dt)
    final = np.array([chi[0, 0, 1].real, speed[0, 0, 1].real])*2
    return float(np.log(np.linalg.norm(final)/np.linalg.norm(initial))/(periods*period))


def qualify(path):
    rows = [homogeneous(h) for h in (.02, .01, .005, .001)]
    rows.extend(homogeneous(.001, phase) for phase in (1., 2.))
    tests = []
    for r, rotating, momenta in ((2., True, (.7, np.sqrt(1.5), 1.7)),
                                 (.5, False, (.4, .85, 1.2))):
        for p in momenta:
            rate, period, matrix = monodromy(p, r, rotating)
            exact = float(floquet(p, r, rotating))
            tests.append(dict(energy=r, momentum=p, rotating=rotating,
                              monodromy=rate, analytic=exact,
                              determinant_error=abs(np.linalg.det(matrix)-1)))
    p = np.sqrt(1.5)
    rate, period, matrix = monodromy(p, 2., True)
    measured = wave_test(p, 2., period, matrix)
    wave_error = abs(measured/rate-1)
    passed = bool(max(abs(r["action_C"]/2-1) for r in rows if r["h"] == .001) < .03
                  and max(abs(t["monodromy"]-t["analytic"]) for t in tests) < 1e-8
                  and max(t["determinant_error"] for t in tests) < 1e-9
                  and wave_error < .02)
    result = dict(passed=passed, homogeneous=rows, floquet=tests,
                  wave=dict(expected=rate, measured=measured, error=wave_error))
    write_json(path, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "results")
    args = parser.parse_args()
    result = qualify(args.output / "qualification.json")
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit("Qualification missed.")


if __name__ == "__main__":
    main()
