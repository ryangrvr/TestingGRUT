"""Fourth-order reciprocal-source refinement on the registered grids."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.signal import fftconvolve

from gate import curve_difference, prescribed
from source import (Config, HERE, Pulse, System, coherence, free_step,
                    initial_frame, kick, outgoing, record_driven, write_json)
from response import profile

BASE4 = Config(bands=7, bloch_nodes=16, radial_nodes=8, radial_cut=8., dt=.04)
FINE4 = replace(BASE4, bands=9, bloch_nodes=20, radial_nodes=10)
CASES4 = {
    "base": BASE4,
    "time": replace(BASE4, dt=.02),
    "spacing": replace(BASE4, bloch_nodes=20, radial_nodes=10),
    "longitudinal": replace(BASE4, bands=9),
    "transverse": replace(BASE4, radial_cut=6.),
    "fine": FINE4,
    "tapered": replace(FINE4, regulator="tapered"),
}
FIRST = 1/(2-np.cbrt(2))
STAGES = np.array([FIRST, 1-2*FIRST, FIRST])
CENTERS = np.cumsum(STAGES)-STAGES/2


def evolve(system, pulse):
    g, dt = system.grid, system.config.dt
    steps = round(2*pulse.duration/dt)
    time = np.linspace(-pulse.duration, pulse.duration, steps+1)
    amplitude, speed, acceleration = profile(time, pulse)
    stages = profile(time[:-1, None]+dt*CENTERS, pulse)[0]
    free_outer = free_step(g, FIRST*dt/2)
    free_inner = free_step(g, (1-FIRST)*dt/2)
    raw, contact = np.zeros(steps+1), np.zeros(steps+1)
    frame = initial_frame(g)
    for i, stage_values in enumerate(stages):
        frame = free_outer @ frame
        for substep, (a, length) in enumerate(zip(stage_values, STAGES)):
            residual, _ = system.mass(a)
            frame = system.rotate(kick(system.rotate(frame), residual, dt*length),
                                  inverse=True)
            frame = (free_inner if substep < 2 else free_outer) @ frame
        _, derivative = system.mass(amplitude[i+1])
        raw[i+1] = system.trace(derivative, coherence(system.rotate(frame)))
        contact[i+1] = system.trace(derivative, system.free_coherence)
    trace = dict(time=time, amplitude=amplitude, speed=speed,
                 acceleration=acceleration, raw=raw, vacuum_contact=contact)
    return trace, outgoing(g, frame, False), stages


def tangent(system, amplitude, stage_amplitudes):
    steps = len(amplitude)-1
    dt = system.config.dt
    nu, strength = system.pair_frequency, system.pair_strength
    result = system.contact*amplitude
    for j, (center, length) in enumerate(zip(CENTERS, STAGES)):
        lag = (np.arange(1, steps+1)-center)*dt
        sine = np.zeros(steps)
        for start in range(0, len(nu), 256):
            batch = slice(start, start+256)
            sine += strength[batch] @ np.sin(nu[batch, None]*lag)
        result[1:] -= dt*length*fftconvolve(sine, stage_amplitudes[:, j])[:steps]
    return result


def run(case, directory):
    system = System(CASES4[case], fold=True)
    trace, out, stages = evolve(system, Pulse())
    linear = tangent(system, trace["amplitude"], stages)
    return record_driven(system, trace, out, linear, directory, case,
                         lambda speed, force: simpson(speed*force, x=trace["time"]))


def qualify(directory):
    pulse = Pulse(amplitude=.7, duration=4., frequency=2.3)
    config = Config(wave=1.4, bands=3, radial_nodes=2, radial_cut=3., dt=.08)
    coarse, out_coarse, _ = evolve(System(config, fold=True), pulse)
    full, out_full, _ = evolve(System(config), pulse)
    refinements = [coarse]
    for dt in (.04, .02, .01):
        refinements.append(evolve(System(replace(config, dt=dt), fold=True), pulse)[0])
    step_errors = [curve_difference(a, b, "raw") for a, b in
                   zip(refinements, refinements[1:])]
    ratios = [a/b for a, b in zip(step_errors, step_errors[1:])]
    peak = np.max(np.abs(full["raw"]))
    folding = max(np.max(np.abs(coarse["raw"]-full["raw"]))/peak,
                  abs(out_coarse["energy"]/out_full["energy"]-1))
    weak = replace(pulse, amplitude=.002)
    system = System(replace(config, dt=.04), fold=True)
    trace, _, stages = evolve(system, weak)
    lin = tangent(system, trace["amplitude"], stages)
    weak_error = np.max(np.abs(trace["raw"]-lin))/np.max(np.abs(lin))
    result = dict(passed=bool(folding < 1e-9 and min(ratios) > 8 and weak_error < .002),
                  folding_error=folding, timestep_differences=step_errors,
                  timestep_ratios=ratios, weak_tangent_error=weak_error)
    write_json(directory / "qualification.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("qualify", "gate", *CASES4))
    parser.add_argument("--output", type=Path, default=HERE / "results" / "refinement")
    args = parser.parse_args()
    if args.case == "qualify":
        result = qualify(args.output)
    else:
        if not json.loads((args.output / "qualification.json").read_text())["passed"]:
            raise SystemExit("Pass the fourth-order qualification before continuing.")
        result = (prescribed(args.output, cases=CASES4) if args.case == "gate"
                  else run(args.case, args.output))
    print(json.dumps(result, indent=2))
    if "passed" in result and not result["passed"]:
        raise SystemExit("Registered fourth-order gate missed.")


if __name__ == "__main__":
    main()
