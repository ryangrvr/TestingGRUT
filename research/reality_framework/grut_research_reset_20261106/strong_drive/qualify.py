"""Run the reference gates without evaluating the reserved strong pulse."""
import argparse
from dataclasses import replace
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag

from propagation import (Config, Grid, HERE, Pulse, cosine_sum, evolve, tangent,
                         write_json)

ARCHIVE = HERE.parent / "results" / "checks.json"


def weak_rows():
    archive = json.loads(ARCHIVE.read_text())["spatial"]
    rows = {}
    for label, reference in archive.items():
        tau, wave = reference["tau"], reference["k"]
        coarse = Config(bands=5, bloch_nodes=32, radial_nodes=8,
                        radial_cut=6., wave=wave)
        fine = replace(coarse, bloch_nodes=80, radial_nodes=24)
        dt = .06
        times = dt * (np.arange(round(8*tau/dt)) + .5)
        gaussian = np.exp(-.5*(times/tau)**2)*np.cos(2.4*times)

        def numerical(frequencies):
            return cosine_sum(frequencies, times, gaussian, dt)

        def analytic(frequencies):
            return tau*np.sqrt(2*np.pi)/2 * (
                np.exp(-.5*(tau*(frequencies-2.4))**2)
                + np.exp(-.5*(tau*(frequencies+2.4))**2))

        coarse_grid, fine_grid = Grid.quadrature(coarse), Grid.quadrature(fine)
        first, second = tangent(coarse_grid, analytic), tangent(fine_grid, numerical)
        exact = tangent(fine_grid, analytic)
        local_first = tangent(coarse_grid, analytic, local=True)
        local_second = tangent(fine_grid, numerical, local=True)
        want, local_want = reference["invariant_energy"], reference["local_energy"]
        rows[label] = dict(
            coarse=first["energy"], fine=second["energy"], reference=want,
            grid_error=abs(first["energy"]/second["energy"]-1),
            transform_error=abs(exact["energy"]/second["energy"]-1),
            relative_error=abs(second["energy"]/want-1),
            local_coarse=local_first["energy"], local_fine=local_second["energy"],
            local_reference=local_want,
            local_relative_error=abs(local_second["energy"]/local_want-1))
    return rows


def weak_compact():
    pulse = Pulse(amplitude=.002, duration=6., frequency=2.5)
    config = Config(bands=3, wave=1., dt=.02, radial_cut=3.)
    grid = Grid.quadrature(config)
    values, _ = pulse.samples(config.dt)
    positive_times = (np.arange(len(values)//2)+.5)*config.dt

    def midpoint(frequencies):
        return cosine_sum(frequencies, positive_times,
                          values[len(values)//2:]/pulse.amplitude, config.dt)

    rows = {}
    for arm in ("propagating", "local"):
        local = arm == "local"
        finite = evolve(grid, config, pulse, local=local)
        linear = tangent(grid, midpoint, local)
        normalized = finite["energy"]/pulse.amplitude**2
        rows[arm] = dict(finite_amplitude_energy_coefficient=normalized,
                         tangent=linear["energy"],
                         relative_error=abs(normalized/linear["energy"]-1))
    return rows


def one_sector_reference(grid, pulse, constant=False):
    """Independent dense Hamiltonian and unsplit DOP853 mode evolution."""
    n = grid.z.size
    h_free = block_diag(*grid.h0[0])
    energies, eigenvectors = np.linalg.eigh(h_free)
    initial = eigenvectors[:, :2*n].astype(complex)
    shape = initial.shape
    c = np.eye(n) if constant else .5*(np.eye(n, k=-1)+np.eye(n, k=1))
    cosine, spatial = np.linalg.eigh(c)
    m_x = np.kron([[0, 1], [1, 0]], np.eye(2))
    minus_m_y = np.kron([[0, 1j], [-1j, 0]], np.eye(2))

    def rhs(t, vector):
        phase = (pulse.amplitude*np.cos(np.pi*t/(2*pulse.duration))**4
                 * np.cos(pulse.frequency*t))
        angle = phase*cosine
        real = (spatial*(np.cos(angle)-1)) @ spatial.T
        imag = (spatial*np.sin(angle)) @ spatial.T
        h = h_free + np.kron(real, m_x) + np.kron(imag, minus_m_y)
        return (-1j * h @ vector.reshape(shape)).ravel()

    result = solve_ivp(rhs, (-pulse.duration, pulse.duration), initial.ravel(),
                       method="DOP853", rtol=1e-11, atol=1e-13,
                       t_eval=[pulse.duration])
    if not result.success:
        raise RuntimeError(result.message)
    amplitudes = eigenvectors.conj().T @ result.y[:, -1].reshape(shape)
    particles = np.sum(abs(amplitudes[2*n:])**2, axis=1)
    holes = 1 - np.sum(abs(amplitudes[:2*n])**2, axis=1)
    return dict(number=float(.5*(particles.sum()+holes.sum())),
                energy=float(.5*(energies[2*n:] @ particles-energies[:2*n] @ holes)),
                evaluations=result.nfev)


def sector_controls():
    config = Config(bands=2, wave=1.4, dt=.01, phases=8)
    z = (.19 + config.wave*np.arange(-2, 3))[None, :]
    grid = Grid.from_points(z, np.array([[.77]]), np.ones(1))
    pulse = Pulse(.7, 5., 2.3)
    reference = one_sector_reference(grid, pulse)
    split = evolve(grid, config, pulse)
    homogeneous = {
        "reference": one_sector_reference(grid, pulse, constant=True),
        "propagating": evolve(grid, config, pulse, constant=True),
        "local": evolve(grid, config, pulse, local=True, constant=True),
    }
    static = evolve(grid, config, replace(pulse, amplitude=0.))
    names = ("number", "energy", "work")

    def compact(row):
        return {name: row[name] for name in names if name in row}

    return dict(reference=reference, split=compact(split),
                split_error={m: abs(split[m]/reference[m]-1) for m in names[:2]},
                homogeneous={a: compact(row) for a, row in homogeneous.items()},
                static=compact(static))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=type(HERE), default=HERE / "results")
    args = parser.parse_args()
    gaussian, compact, controls = weak_rows(), weak_compact(), sector_controls()
    homogeneous = controls["homogeneous"]
    checks = dict(
        weak_gaussian=all(row["relative_error"] < .002 and
                          row["local_relative_error"] < .002 for row in gaussian.values()),
        weak_compact=all(row["relative_error"] < .005 for row in compact.values()),
        unsplit_sector=max(controls["split_error"].values()) < .002,
        static_vacuum=max(abs(v) for v in controls["static"].values()) < 1e-11,
        homogeneous_match=all(abs(homogeneous["propagating"][m] -
                                  homogeneous["local"][m]) < 1e-11
                              for m in ("number", "energy", "work")),
        homogeneous_reference=all(abs(homogeneous["propagating"][m] /
                                      homogeneous["reference"][m]-1) < .002
                                  for m in ("number", "energy")))
    result = dict(gaussian=gaussian, compact=compact, controls=controls,
                  checks=checks, passed=all(checks.values()))
    write_json(args.output / "qualification.json", result)
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit("Qualification failed. Reserved pulse stays unrun.")


if __name__ == "__main__":
    main()
