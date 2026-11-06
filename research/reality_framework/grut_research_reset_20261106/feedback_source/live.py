"""Reciprocal feedback with the same matched loop force and Hamiltonian."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

from gate import LIVE_CASES
from source import CASES, HERE, Pulse, System, free_step, initial_frame, save
from response import profile

WAVE = 1.3
INERTIA = 1.
FREQUENCY2 = WAVE**2 + 1.1**2
REFERENCE = Pulse(amplitude=.85, duration=8., frequency=2.7)
END_TIME = 12.


def run(case, directory):
    config = replace(CASES[case], wave=WAVE)
    system = System(config)
    grid, dt = system.grid, config.dt
    inertia = INERTIA-system.b
    spring = INERTIA*FREQUENCY2+system.a
    if min(inertia, spring) <= 0:
        raise ValueError("The matched bare scalar quadratic form is unstable.")

    steps = round((END_TIME+REFERENCE.duration)/dt)
    time = -REFERENCE.duration + np.arange(steps+1)*dt
    midpoint = time[:-1]+dt/2
    a_ref, _, acc_ref = profile(midpoint, REFERENCE)
    drive = INERTIA*(acc_ref+FREQUENCY2*a_ref)
    reference, _, _ = profile(time, REFERENCE)
    frame = initial_frame(grid)
    free_half = free_step(grid, dt/2)

    amplitude, speed, work = np.zeros(steps+1), np.zeros(steps+1), np.zeros(steps+1)
    raw_force, matched_force = np.empty(steps), np.empty(steps)
    energy_t, e_loop, e_scalar, e_total, e_work = [time[0]], [0.], [0.], [0.], [0.]
    a, p = 0., 0.
    stride = max(1, round(.08/dt))
    for i, force_drive in enumerate(drive):
        a_previous = a
        a += dt*p/(2*inertia)
        frame, force = system.step_mass(free_half @ frame, a, dt)
        acceleration = (force_drive-spring*a-force)/inertia
        raw_force[i] = force
        matched_force[i] = force+system.a*a-system.b*acceleration
        p += dt*inertia*acceleration
        a += dt*p/(2*inertia)
        frame = free_half @ frame
        amplitude[i+1], speed[i+1] = a, p/inertia
        work[i+1] = work[i] + force_drive*(a-a_previous)

        if (i+1) % stride == 0 or i == steps-1:
            scalar = (p*p/inertia+spring*a*a)/2
            loop = system.loop_energy(frame, a)
            energy_t.append(time[i+1])
            e_loop.append(loop)
            e_scalar.append(scalar)
            e_total.append(scalar+loop)
            e_work.append(work[i+1])

    total, injected = np.array(e_total), np.array(e_work)
    peak = float(np.max(np.abs(total)))
    departure = float(np.max(np.abs(amplitude-reference)))
    defect = float(np.max(np.abs(total-injected))/peak)
    trace = dict(time=time, amplitude=amplitude, speed=speed, reference=reference,
                 drive_time=midpoint, external_force=drive, work=work,
                 raw_loop_force=raw_force, matched_loop_force=matched_force,
                 ledger_time=energy_t, loop_energy=e_loop, scalar_energy=e_scalar,
                 total_energy=total, injected_sampled=injected)
    row = dict(case=case, config=vars(config), match=system.match,
               bare_inertia=inertia, bare_spring=spring,
               relative_energy_defect=defect, peak_total_energy=peak,
               amplitude_departure=departure,
               final=dict(time=END_TIME, amplitude=a, speed=p/inertia,
                          scalar_bare_energy=e_scalar[-1], loop_energy=e_loop[-1],
                          total_energy=total[-1], drive_work=injected[-1]))
    save(directory, "live_"+case, trace, row)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=LIVE_CASES)
    parser.add_argument("--output", type=Path, default=HERE / "results")
    args = parser.parse_args()
    gate = args.output / "prescribed_gate.json"
    if not gate.exists() or not json.loads(gate.read_text())["passed"]:
        raise SystemExit("Pass the prescribed force gate before evolving feedback.")
    print(json.dumps(run(args.case, args.output), indent=2))


if __name__ == "__main__":
    main()
