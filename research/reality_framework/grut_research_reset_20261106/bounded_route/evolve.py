"""Fixed physical ALP window using the matched conformal states."""
import argparse
from dataclasses import asdict
import json
import signal

import numpy as np
from scipy.integrate import cumulative_simpson

from match import Box, CASES, DIRECTORY, H, START, WINDOW, initial_state
from model import write_json


def evolve(name, stop_requested):
    config = CASES[name]
    with np.load(DIRECTORY / "matched_states.npz") as f:
        ids, times, states = f["ids"], f["times"], f["states"]
    box = Box(config.n, config.length)
    start = times[config.seam]
    end = times[0]+WINDOW
    steps = int(np.ceil((end-start)/config.dt))
    dt = (end-start)/steps
    a = START+H*start
    chi, momentum = initial_state(box, ids, states[config.seam], a, config.seed)
    acceleration, theta = box.acceleration(chi, a)
    ledger, late, spectra = [], [], []
    stride = max(1, round(.20/dt))
    late_stride = max(1, round(.50/dt))
    for step in range(steps+1):
        t = start+step*dt
        a = START+H*t
        if step:
            acceleration, theta = box.integrate_step(chi, momentum, acceleration, a, dt)
        capped = stop_requested()
        final_window = t >= end-100 and (step % late_stride == 0 or step == steps)
        if step % stride == 0 or final_window or step == steps or capped:
            velocity, rho, pressure, fluctuations = box.energy(chi, momentum, theta, a, H)
            phase = chi[0, 0, 0].real/a
            coherent_kinetic = .5*abs(velocity[0, 0, 0])**2/rho
            ledger.append([t-times[0], a, a**3*rho, -3*a*a*pressure,
                           fluctuations, coherent_kinetic, phase])
            if final_window:
                row, spectrum = box.late_diagnostic(chi, velocity, theta, rho, a, H)
                late.append(row)
                spectra.append(spectrum)
        if capped:
            break
    ledger = np.asarray(ledger)
    late = np.asarray(late).reshape(-1, 9)
    spectra = np.asarray(spectra).reshape(len(late), -1) if len(late) else np.empty((0, 0))
    work = cumulative_simpson(ledger[:, 3], x=ledger[:, 1], initial=0.)
    defect = ledger[:, 2]-ledger[0, 2]-work
    row = dict(case=name, config=asdict(config), steps_completed=step,
               steps_planned=steps, completed=step == steps, dt=dt,
               initial_a=float(ledger[0, 1]), final_a=float(ledger[-1, 1]),
               initial_fluctuation_fraction=float(ledger[0, 4]),
               max_fluctuation_fraction=float(max(ledger[:, 4])),
               terminal_fluctuation_fraction=float(ledger[-1, 4]),
               terminal_coherent_kinetic_fraction=float(ledger[-1, 5]),
               energy_residual=float(max(abs(defect))/ledger[-1, 2]))
    if len(late):
        window = ledger[:, 0] >= WINDOW-100
        phase = ledger[window, 6]
        fixed = bool(np.ptp(np.rint(phase/(2*np.pi))) == 0)
        proxy = float(late[:, 1].mean())
        checks = dict(fixed_branch=fixed, drift=float(np.ptp(late[:, 1])/proxy),
                      dictionary=float(max(abs(late[:, 2]/late[:, 1]-1))),
                      anharmonic=float(max(late[:, 5])),
                      outside_branch=float(max(late[:, 6])),
                      cutoff_tail=float(max(late[:, 7])))
        row.update(quadratic_proxy=proxy, terminal_window=checks,
                   winding_span=float(np.ptp(phase)/(2*np.pi)),
                   coherent_kinetic_final_window_mean=float(np.mean(ledger[window, 5])))
    trace = dict(ledger=ledger, pressure_work=work, energy_defect=defect,
                 late=late, late_quadratic_spectrum=spectra)
    return row, trace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=CASES)
    args = parser.parse_args()
    matched = json.loads((DIRECTORY / "match.json").read_text())
    if not matched["passed"]:
        raise SystemExit("Pass the linear match before running the lattice window.")
    stopped = False

    def cap(_signum, _frame):
        nonlocal stopped
        stopped = True

    signal.signal(signal.SIGTERM, cap)
    row, trace = evolve(args.case, lambda: stopped)
    path = DIRECTORY / (args.case+".tmp.npz")
    np.savez_compressed(path, **trace)
    path.replace(DIRECTORY / (args.case+".npz"))
    write_json(DIRECTORY / (args.case+".json"), row)
    print(json.dumps(row, indent=2))


if __name__ == "__main__":
    main()
