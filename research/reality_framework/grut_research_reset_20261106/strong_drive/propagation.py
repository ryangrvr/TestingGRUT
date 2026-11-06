"""Majorana vacuum production by a periodic mass phase, with a local control.

Natural units set M=1. Each Bloch sector includes the two transverse spins.
The local and propagating arms share the quadrature and the free in/out basis.
"""
import argparse
from dataclasses import dataclass, replace
import json
from pathlib import Path

import numpy as np
from scipy.fft import dst, fft, ifft
from scipy.special import roots_legendre

HERE = Path(__file__).resolve().parent
PANELS = (0., .5, 1., 2., 3., 4., 5., 6., 8.)
BINS = (0., .5, 1., 1.5, 2., 3., 4., 6., 9., 14., 25., np.inf)
IDENTITY = np.eye(4)
TANGENT = np.block([[np.zeros((2, 2)), 1j * np.eye(2)],
                    [-1j * np.eye(2), np.zeros((2, 2))]])


@dataclass(frozen=True)
class Pulse:
    amplitude: float = 1.2
    duration: float = 12.
    frequency: float = 2.4

    def samples(self, dt):
        steps = round(2 * self.duration / dt)
        if not np.isclose(steps * dt, 2 * self.duration):
            raise ValueError("timestep must tile the pulse")
        t = (np.arange(steps) + .5) * dt - self.duration
        alpha = np.pi / (2 * self.duration)
        c, s = np.cos(alpha * t), np.sin(alpha * t)
        carrier, quadrature = np.cos(self.frequency * t), np.sin(self.frequency * t)
        value = self.amplitude * c**4 * carrier
        derivative = -self.amplitude * (4 * alpha * c**3 * s * carrier
                                        + self.frequency * c**4 * quadrature)
        return value, derivative


@dataclass(frozen=True)
class Config:
    bands: int = 4
    bloch_nodes: int = 6
    radial_nodes: int = 3
    radial_cut: float = 4.
    phases: int = 16
    dt: float = .04
    regulator: str = "open"
    wave: float = 1.6


BASE = Config()
FINE = Config(bands=7, bloch_nodes=12, radial_nodes=5, radial_cut=6.,
              phases=32, dt=.02)
CASES = {
    "base": BASE,
    "time": replace(BASE, dt=.02),
    "spacing": replace(BASE, bloch_nodes=12, radial_nodes=5, phases=32),
    "longitudinal": replace(BASE, bands=7),
    "transverse": replace(BASE, radial_cut=6.),
    "fine": FINE,
    "regulator": replace(FINE, regulator="cyclic"),
}


@dataclass
class Grid:
    z: np.ndarray
    r: np.ndarray
    weight: np.ndarray
    h0: np.ndarray
    energy: np.ndarray
    negative: np.ndarray
    positive: np.ndarray

    @classmethod
    def from_points(cls, z, r, weight):
        r = np.broadcast_to(r, z.shape)
        e = np.sqrt(1 + r*r + z*z)
        h = np.zeros((*z.shape, 4, 4))
        h[..., 0, 0] = h[..., 3, 3] = z
        h[..., 1, 1] = h[..., 2, 2] = -z
        h[..., 0, 1] = h[..., 1, 0] = r
        h[..., 2, 3] = h[..., 3, 2] = -r
        h[..., 0, 2] = h[..., 2, 0] = h[..., 1, 3] = h[..., 3, 1] = 1.
        _, basis = np.linalg.eigh(h)
        return cls(z, r, weight, h, e, basis[..., :2], basis[..., 2:])

    @classmethod
    def quadrature(cls, config):
        x, w = roots_legendre(config.bloch_nodes)
        bloch, wb = x * config.wave/2, w * config.wave/2
        if config.radial_cut not in PANELS[1:]:
            raise ValueError("transverse cutoff must be a panel edge")
        x, w = roots_legendre(config.radial_nodes)
        radial, wr = [], []
        for a, b in zip(PANELS, PANELS[1:]):
            if b > config.radial_cut:
                break
            radius = (a+b)/2 + (b-a)/2*x
            radial.extend(radius)
            wr.extend((b-a)/2 * w * radius)
        z = np.repeat(bloch, len(radial))[:, None]
        z = z + config.wave * np.arange(-config.bands, config.bands+1)
        r = np.tile(radial, len(bloch))[:, None]
        weight = np.multiply.outer(wb, wr).ravel() / (2*np.pi)**2
        return cls.from_points(z, r, weight)


def spatial_basis(n):
    return dst(np.eye(n), type=1, axis=0, norm="ortho")


def phase_nodes(n, regulator):
    if regulator == "open":
        return np.cos(np.pi * np.arange(1, n+1) / (n+1))
    if regulator == "cyclic":
        return np.cos(2*np.pi * np.arange(n) / n)
    raise ValueError(regulator)


def transform(frame, regulator, inverse=False):
    if regulator == "open":
        return dst(frame, type=1, axis=1, norm="ortho", workers=1)
    operation = ifft if inverse else fft
    return operation(frame, axis=1, norm="ortho", workers=1)


def initial_frame(grid):
    sectors, n = grid.z.shape
    frame = np.zeros((sectors, n, 4, 2*n), dtype=complex)
    for band in range(n):
        frame[:, band, :, 2*band:2*band+2] = grid.negative[:, band]
    return frame


def free_step(grid, dt):
    angle = dt * grid.energy[..., None, None]
    return np.cos(angle) * IDENTITY - 1j * dt * np.sinc(angle/np.pi) * grid.h0


def kick(frame, residual, dt):
    angle = dt * np.abs(residual)
    c = np.cos(angle)[..., None, None]
    s = (dt * np.sinc(angle/np.pi) * residual)[..., None, None]
    top, bottom = frame[..., :2, :], frame[..., 2:, :]
    return np.concatenate((c*top - 1j*s*bottom, c*bottom - 1j*s.conj()*top), axis=-2)


def coherence(frame):
    return np.einsum("...ic,...ic->...", frame[..., :2, :].conj(), frame[..., 2:, :])


def projection(basis, frame):
    amplitudes = basis.swapaxes(-1, -2).conj() @ frame
    return np.sum(np.abs(amplitudes)**2, axis=(-2, -1))


def outgoing(grid, frame, local, details=False):
    particle = projection(grid.positive, frame)
    hole = 2 - projection(grid.negative, frame)
    norm_defect = float(np.max(np.abs(particle-hole)))
    if local:
        particle, hole = particle.mean(axis=0), hole.mean(axis=0)
    average = .5 * (particle+hole)
    weight = grid.weight[:, None]
    p = np.hypot(grid.r, grid.z)
    number = float(np.sum(weight*average))
    energy = float(np.sum(weight*grid.energy*average))
    row = dict(number=number, energy=energy,
               particle_number=float(np.sum(weight*particle)),
               hole_number=float(np.sum(weight*hole)),
               sector_number_defect=norm_defect if local else
               float(np.max(np.abs((particle-hole).sum(axis=1)))))
    if details:
        # Free outgoing moments. The pressure here is a projected diagnostic.
        row["pressure_free"] = float(np.sum(weight*average*p*p/(3*grid.energy)))
        row["negative_hole_roundoff"] = float(min(0., np.min(hole)))
        row["spectrum_number"] = np.histogram(p, bins=BINS, weights=weight*average)[0].tolist()
        row["spectrum_energy"] = np.histogram(p, bins=BINS,
                                               weights=weight*average*grid.energy)[0].tolist()
    return row


def evolve(grid, config, pulse, *, local=False, constant=False):
    n = grid.z.shape[1]
    frame = initial_frame(grid)
    vacuum = -1/grid.energy
    weight = grid.weight[:, None]
    if local:
        phase = np.pi * (np.arange(config.phases) + .5) / config.phases
        frame = np.broadcast_to(grid.negative, (config.phases, *grid.negative.shape))
        wave = np.ones_like(phase) if constant else np.cos(phase)
        wave = wave[:, None, None]
        weight = weight / config.phases
    else:
        wave = np.ones(n) if constant else phase_nodes(n, config.regulator)
        if config.regulator == "open":
            vacuum = vacuum @ (spatial_basis(n)**2)
        else:
            vacuum = vacuum.mean(axis=1, keepdims=True)
    dt = config.dt
    half_free, whole_free = free_step(grid, dt/2), free_step(grid, dt)
    values, derivatives = pulse.samples(dt)
    powers = np.empty(len(values))
    frame = half_free @ frame
    for step, (value, derivative) in enumerate(zip(values, derivatives)):
        if not local:
            frame = transform(frame, config.regulator)
        residual = np.expm1(1j * value * wave)
        frame = kick(frame, residual, dt/2)
        d_mass = 1j * wave * (1+residual) * derivative
        # Two conjugate chiral blocks cancel the factor 1/2 in the Majorana trace.
        powers[step] = np.sum(weight * (d_mass * (coherence(frame)-vacuum)).real)
        frame = kick(frame, residual, dt/2)
        if not local:
            frame = transform(frame, config.regulator, inverse=True)
        step_free = half_free if step == len(values)-1 else whole_free
        frame = step_free @ frame
    row = outgoing(grid, frame, local, details=True)
    row.update(work=float(dt*np.sum(powers)),
               positive_work=float(dt*np.maximum(powers, 0.).sum()),
               negative_work=float(dt*np.minimum(powers, 0.).sum()),
               peak_abs_power=float(np.max(np.abs(powers))))
    stride = max(1, len(values)//96)
    row["work_trace_dt"] = dt
    row["work_trace"] = [
        [step, float(powers[step]), float(total)]
        for step, total in enumerate(dt*np.cumsum(powers))
        if step % stride == 0 or step == len(values)-1
    ]
    return row


def spin_sum(positive, negative):
    overlap = positive.swapaxes(-1, -2).conj() @ TANGENT @ negative
    return np.sum(np.abs(overlap)**2, axis=(-2, -1))


def tangent(grid, transform_function, local=False):
    """First-order pair amplitudes of the split kick and exact free rotations."""
    if local:
        nu = 2*grid.energy
        squared = spin_sum(grid.positive, grid.negative) * .5
    else:
        nu = grid.energy[:, 1:] + grid.energy[:, :-1]
        squared = (spin_sum(grid.positive[:, 1:], grid.negative[:, :-1])
                   + spin_sum(grid.positive[:, :-1], grid.negative[:, 1:])) * .25
    weight = grid.weight[:, None] * squared * .5
    contribution = weight * transform_function(nu)**2
    return dict(number=float(np.sum(contribution)), energy=float(np.sum(nu*contribution)))


def cosine_sum(frequency, positive_times, signal, dt):
    shape = frequency.shape
    flat = frequency.ravel()
    out = np.empty_like(flat)
    for start in range(0, len(flat), 256):
        batch = flat[start:start+256]
        out[start:start+256] = 2*dt*(np.cos(batch[:, None]*positive_times) @ signal)
    return out.reshape(shape)


def write_json(path, row):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(row, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def run_case(case, output):
    config = CASES[case]
    grid = Grid.quadrature(config)
    result = dict(case=case, config=vars(config), sectors=len(grid.weight))
    for arm in ("propagating", "local"):
        result[arm] = evolve(grid, config, Pulse(), local=arm == "local")
        print(case, arm, {k: result[arm][k] for k in ("number", "work", "energy")}, flush=True)
    write_json(output / (case + ".json"), result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=CASES)
    parser.add_argument("--output", type=Path, default=HERE / "results")
    args = parser.parse_args()
    gate_path = args.output / "qualification.json"
    if not gate_path.is_file() or not json.loads(gate_path.read_text())["passed"]:
        raise SystemExit("Run and pass qualify.py before the reserved pulse.")
    run_case(args.case, args.output)


if __name__ == "__main__":
    main()
