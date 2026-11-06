"""Match the force conjugate to one propagating Majorana mass-phase mode."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "strong_drive"))
from propagation import (Config, Grid, Pulse, coherence, free_step,
                         initial_frame, kick, outgoing, phase_nodes,
                         spatial_basis, spin_sum, transform, write_json)

from response import (continuum_force, kernel, parameter_integral, profile,
                      spectral_work, split_tangent, static_reference)

BASE = Config(bands=5, bloch_nodes=8, radial_nodes=4, radial_cut=5., dt=.04)
FINE = Config(bands=8, bloch_nodes=12, radial_nodes=5, radial_cut=8., dt=.02)
CASES = {
    "base": BASE,
    "time": replace(BASE, dt=.02),
    "spacing": replace(BASE, bloch_nodes=12, radial_nodes=5),
    "longitudinal": replace(BASE, bands=8),
    "transverse": replace(BASE, radial_cut=8.),
    "fine": FINE,
    "tapered": replace(FINE, regulator="tapered"),
}


class System:
    def __init__(self, config, *, fold=False):
        self.config = config
        self.grid = Grid.quadrature(config)
        if fold:
            if config.bloch_nodes % 2:
                raise ValueError("Bloch folding requires an even quadrature order.")
            positive = self.grid.z[:, config.bands] > 0
            self.grid = Grid(**{key: value[positive] for key, value in vars(self.grid).items()})
            self.grid.weight *= 2
        self.weight = self.grid.weight[:, None]
        n = self.grid.z.shape[1]
        links = np.full(n-1, .5)
        self.tapered = config.regulator == "tapered"
        if self.tapered:
            links[[0, -1]] *= .25
            links[[1, -2]] *= .75
            c = np.diag(links, -1)+np.diag(links, 1)
            self.wave, self.basis = np.linalg.eigh(c)
        else:
            self.wave = phase_nodes(n, "open")
            self.basis = spatial_basis(n)
        self.free_coherence = (-1/self.grid.energy) @ self.basis**2
        self.pair_frequency, self.pair_strength, self.contact = self.linear_data(links)
        self.a, self.b, self.match = self.counterterms()

    def linear_data(self, links):
        g = self.grid
        nu = g.energy[:, 1:]+g.energy[:, :-1]
        amplitude2 = (spin_sum(g.positive[:, 1:], g.negative[:, :-1])
                      + spin_sum(g.positive[:, :-1], g.negative[:, 1:]))
        strength = self.weight*links**2*amplitude2
        c2_diagonal = np.pad(links**2, (1, 0))+np.pad(links**2, (0, 1))
        contact = float(np.sum(self.weight*c2_diagonal/g.energy))
        return nu.ravel(), strength.ravel(), contact

    def counterterms(self):
        nu, strength = self.pair_frequency, self.pair_strength
        bare_value = self.contact-np.sum(strength/nu)
        bare_slope = -np.sum(strength/nu**3)
        value, slope = static_reference(self.config.wave)
        a, b = value-bare_value, slope-bare_slope
        return a, b, dict(a=a, b=b, contact=self.contact,
                          cutoff_static=bare_value, cutoff_omega2_slope=bare_slope,
                          continuum_static=value, continuum_omega2_slope=slope)

    def rotate(self, frame, inverse=False):
        if not self.tapered:
            return transform(frame, "open", inverse=inverse)
        basis = self.basis if inverse else self.basis.T
        shape = frame.shape
        return (basis @ frame.reshape(*shape[:2], -1)).reshape(shape)

    def mass(self, amplitude):
        residual = np.expm1(1j*amplitude*self.wave)
        derivative = 1j*self.wave*(1+residual)
        return residual, derivative

    def trace(self, derivative, moment):
        # The two chiral blocks cancel the Majorana half in this off-diagonal trace.
        return float(np.sum(self.weight*(derivative*moment).real))

    def step_mass(self, frame, amplitude, dt):
        residual, derivative = self.mass(amplitude)
        frame = self.rotate(frame)
        frame = kick(frame, residual, dt/2)
        force = self.trace(derivative, coherence(frame))
        frame = kick(frame, residual, dt/2)
        return self.rotate(frame, inverse=True), force

    def loop_energy(self, frame, amplitude):
        free = outgoing(self.grid, frame, False)["energy"]
        residual, _ = self.mass(amplitude)
        return free+self.trace(residual, coherence(self.rotate(frame)))

    def tangent(self, amplitude):
        return split_tangent(amplitude, self.config.dt, self.pair_frequency,
                             self.pair_strength, self.contact)


def evolve_driven(system, pulse):
    grid, dt = system.grid, system.config.dt
    steps = round(2*pulse.duration/dt)
    time = -pulse.duration + (np.arange(steps)+.5)*dt
    amplitude, speed, acceleration = profile(time, pulse)
    raw, vacuum = np.empty(steps), np.empty(steps)
    half, whole = free_step(grid, dt/2), free_step(grid, dt)
    frame = half @ initial_frame(grid)
    for i, value in enumerate(amplitude):
        frame, raw[i] = system.step_mass(frame, value, dt)
        _, derivative = system.mass(value)
        vacuum[i] = system.trace(derivative, system.free_coherence)
        frame = (half if i == steps-1 else whole) @ frame
    energy = outgoing(grid, frame, False)
    data = dict(time=time, amplitude=amplitude, speed=speed,
                acceleration=acceleration, raw=raw, vacuum_contact=vacuum)
    return data, energy


def driven(case, directory):
    system = System(CASES[case])
    trace, out = evolve_driven(system, Pulse())
    finite_linear = system.tangent(trace["amplitude"])
    return record_driven(system, trace, out, finite_linear, directory, case,
                         lambda speed, force: system.config.dt*np.dot(speed, force))


def record_driven(system, trace, out, finite_linear, directory, case, integrate):
    a, v, acc, raw = (trace[key] for key in ("amplitude", "speed", "acceleration", "raw"))
    dt = system.config.dt
    continuum, before = continuum_force(a, dt, system.config.wave)
    local_match = raw+system.a*a-system.b*acc
    full_linear_match = raw-finite_linear+continuum
    trace.update(old=raw-trace["vacuum_contact"], local_match=local_match,
                 linear_match=full_linear_match, finite_linear=finite_linear,
                 continuum_linear=continuum, nonlinear=raw-finite_linear)
    terms = ("old", "raw", "local_match", "linear_match", "nonlinear")
    work = {term: float(integrate(v, trace[term])) for term in terms}
    peak = {term: float(np.max(np.abs(trace[term]))) for term in terms}
    row = dict(case=case, config=vars(system.config), match=system.match,
               work=work, peak_force=peak, outgoing=out,
               work_relative_error={key: abs(work[key]/out["energy"]-1)
                                    for key in ("raw", "old", "local_match", "linear_match")},
               continuum_pre_pulse=float(np.max(np.abs(before))))
    save(directory, case, trace, row)
    return row


def save(directory, name, trace, row):
    directory.mkdir(parents=True, exist_ok=True)
    temporary = directory / (name + ".tmp.npz")
    np.savez_compressed(temporary, **trace)
    temporary.replace(directory / (name + ".npz"))
    write_json(directory / (name + ".json"), row)


def qualify(directory):
    kernel_errors = {}
    for s in (-12., -3., -.001, .001, .4, 3.9, 4.1, 9., 60.):
        expected = parameter_integral(s)
        kernel_errors[str(s)] = abs(complex(kernel(s))-expected)/abs(expected)
    pulse = Pulse()
    want = spectral_work(pulse, BASE.wave)
    rows = []
    for dt, padding in ((.04, 128), (.04, 512), (.02, 512), (.02, 1024)):
        time = -pulse.duration + (np.arange(round(2*pulse.duration/dt))+.5)*dt
        a, v, _ = profile(time, pulse)
        f, pre = continuum_force(a, dt, BASE.wave, padding=padding)
        work = float(dt*np.dot(v, f))
        rows.append(dict(dt=dt, padding=padding, work=work,
                         error=abs(work/want-1),
                         pre_fraction=float(np.max(np.abs(pre))/np.max(np.abs(f)))))
    weak = replace(pulse, amplitude=.002, duration=6., frequency=2.5)
    small = System(Config(bands=3, radial_cut=3., dt=.02))
    trace, _ = evolve_driven(small, weak)
    tangent = small.tangent(trace["amplitude"])
    weak_error = float(np.max(np.abs(trace["raw"]-tangent))/np.max(np.abs(tangent)))
    grid = small.grid
    contact = np.sum(small.weight/grid.energy)
    tangent_static = np.sum(small.weight*spin_sum(grid.positive, grid.negative)/(2*grid.energy))
    constant_error = float(abs(contact-tangent_static)/contact)
    checks = dict(kernel=max(kernel_errors.values()) < 1e-9,
                  continuum=all(row["error"] < .002 and row["pre_fraction"] < .002
                                for row in rows[-2:]),
                  weak=weak_error < .002, constant=constant_error < 1e-12)
    result = dict(passed=all(checks.values()), checks=checks,
                  kernel_relative_errors=kernel_errors,
                  continuum=rows, spectral_work=want,
                  weak_force_error=weak_error, constant_static_error=constant_error)
    write_json(directory / "qualification.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("qualify", *CASES))
    parser.add_argument("--output", type=Path, default=HERE / "results")
    args = parser.parse_args()
    if args.case == "qualify":
        result = qualify(args.output)
    else:
        gate = args.output / "qualification.json"
        if not gate.exists() or not json.loads(gate.read_text())["passed"]:
            raise SystemExit("Pass qualification before running the source comparison.")
        result = driven(args.case, args.output)
    print(json.dumps(result, indent=2))
    if args.case == "qualify" and not result["passed"]:
        raise SystemExit("Qualification failed")


if __name__ == "__main__":
    main()
