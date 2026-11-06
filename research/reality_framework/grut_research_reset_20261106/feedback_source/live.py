"""Reciprocal feedback with the same matched loop force and Hamiltonian."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

from gate import LIVE_CASES, live
from higher_order import CASES4, CENTERS, STAGES
from source import CASES, HERE, Pulse, System, free_step, initial_frame, save
from response import profile

WAVE = 1.3
INERTIA = 1.
FREQUENCY2 = WAVE**2 + 1.1**2
REFERENCE = Pulse(amplitude=.85, duration=8., frequency=2.7)
END_TIME = 12.


def run(case, directory, fourth_order=False):
    cases = CASES4 if fourth_order else CASES
    config = replace(cases[case], wave=WAVE)
    system = System(config, fold=fourth_order)
    grid, dt = system.grid, config.dt
    inertia = INERTIA-system.b
    spring = INERTIA*FREQUENCY2+system.a
    if min(inertia, spring) <= 0:
        raise ValueError("The matched bare scalar quadratic form is unstable.")

    steps = round((END_TIME+REFERENCE.duration)/dt)
    time = -REFERENCE.duration + np.arange(steps+1)*dt
    lengths = dt*(STAGES if fourth_order else np.ones(1))
    centers = CENTERS if fourth_order else np.array([.5])
    stage_time = time[:-1, None]+dt*centers
    a_ref, _, acc_ref = profile(stage_time, REFERENCE)
    drive = INERTIA*(acc_ref+FREQUENCY2*a_ref)
    reference, _, _ = profile(time, REFERENCE)
    frame = initial_frame(grid)
    free_outer = free_step(grid, lengths[0]/2)
    bridges = [free_step(grid, (left+right)/2) for left, right in
               zip(lengths, lengths[1:])] + [free_outer]

    amplitude, speed, work = np.zeros(steps+1), np.zeros(steps+1), np.zeros(steps+1)
    no_fermions = np.zeros(steps+1)
    raw_force, matched_force = np.empty_like(drive), np.empty_like(drive)
    energy_t, e_loop, e_scalar, e_total, e_work = [time[0]], [0.], [0.], [0.], [0.]
    a, p, injected, a_test, v_test = 0., 0., 0., 0., 0.
    stride = max(1, round(.08/dt))
    for i, forces in enumerate(drive):
        frame = free_outer @ frame
        for j, (force_drive, length, free_after) in enumerate(zip(forces, lengths, bridges)):
            a_previous = a
            a += length*p/(2*inertia)
            frame, force = system.step_mass(frame, a, length)
            acceleration = (force_drive-spring*a-force)/inertia
            raw_force[i, j] = force
            matched_force[i, j] = force+system.a*a-system.b*acceleration
            p += length*inertia*acceleration
            a += length*p/(2*inertia)
            frame = free_after @ frame
            injected += force_drive*(a-a_previous)

            a_test += length*v_test/2
            v_test += length*(force_drive/INERTIA-FREQUENCY2*a_test)
            a_test += length*v_test/2
        amplitude[i+1], speed[i+1], work[i+1] = a, p/inertia, injected
        no_fermions[i+1] = a_test

        if (i+1) % stride == 0 or i == steps-1:
            scalar = (p*p/inertia+spring*a*a)/2
            loop = system.loop_energy(frame, a)
            energy_t.append(time[i+1])
            e_loop.append(loop)
            e_scalar.append(scalar)
            e_total.append(scalar+loop)
            e_work.append(injected)

    total, injected_sampled = np.array(e_total), np.array(e_work)
    peak = float(np.max(np.abs(total)))
    departure = float(np.max(np.abs(amplitude-reference)))
    defect = float(np.max(np.abs(total-injected_sampled))/peak)
    trace = dict(time=time, amplitude=amplitude, speed=speed, reference=reference,
                 no_fermions=no_fermions, drive_time=stage_time, external_force=drive,
                 work=work, raw_loop_force=raw_force, matched_loop_force=matched_force,
                 ledger_time=energy_t, loop_energy=e_loop, scalar_energy=e_scalar,
                 total_energy=total, injected_sampled=injected_sampled)
    row = dict(case=case, config=vars(config), match=system.match,
               bare_inertia=inertia, bare_spring=spring,
               relative_energy_defect=defect, peak_total_energy=peak,
               amplitude_departure=departure,
               zero_loop_numerical_error=float(np.max(np.abs(no_fermions-reference))),
               final=dict(time=END_TIME, amplitude=a, speed=p/inertia,
                          scalar_bare_energy=e_scalar[-1], loop_energy=e_loop[-1],
                          total_energy=total[-1], drive_work=injected))
    save(directory, "live_"+case, trace, row)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=(*LIVE_CASES, "gate"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--fourth-order", action="store_true")
    args = parser.parse_args()
    default = HERE / "results"
    if args.fourth_order:
        default = default / "refinement"
    directory = args.output or default
    gate = directory / "prescribed_gate.json"
    if not gate.exists() or not json.loads(gate.read_text())["passed"]:
        raise SystemExit("Pass the prescribed force gate before evolving feedback.")
    if args.case == "gate":
        cases = CASES4 if args.fourth_order else CASES
        result = live(directory, cases)
    else:
        result = run(args.case, directory, fourth_order=args.fourth_order)
    print(json.dumps(result, indent=2))
    if "passed" in result and not result["passed"]:
        raise SystemExit("Registered live gate missed.")


if __name__ == "__main__":
    main()
