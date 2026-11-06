"""Band-limited 3D cosine field in a prescribed radiation background."""
import argparse
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scipy.fft import irfftn, rfftn
from scipy.integrate import cumulative_simpson

from model import A_S, CASES, HERE, Q, rotation_speed, window, write_json

BINS = np.array([0., .25, .5, .75, 1., 1.25, 1.5, 2., 2.5, 3., 4., 5., 6., 8., 12.])


def mix(bits):
    """SplitMix64 permutation on unsigned labels."""
    bits = (bits ^ (bits >> np.uint64(30))) * np.uint64(0xbf58476d1ce4e5b9)
    bits = (bits ^ (bits >> np.uint64(27))) * np.uint64(0x94d049bb133111eb)
    return bits ^ (bits >> np.uint64(31))


def gaussian(ix, iy, iz, seed):
    """Hermitian complex Gaussian labelled by the physical wavevector."""
    positive = (ix > 0) | ((ix == 0) & ((iy > 0) | ((iy == 0) & (iz >= 0))))
    sign = np.where(positive, 1, -1)
    x = (sign*ix).astype(np.uint64)
    y = (sign*iy).astype(np.uint64)
    z = (sign*iz).astype(np.uint64)
    label = ((x*np.uint64(0x9e3779b185ebca87)) ^ (y*np.uint64(0xc2b2ae3d27d4eb4f)) ^
             (z*np.uint64(0x165667b19e3779f9)) ^ mix(np.array([seed], dtype=np.uint64))[0])
    scale = 2.**-53
    u = ((mix(label) >> np.uint64(11)).astype(float)+.5)*scale
    v = ((mix(label ^ np.uint64(0x4b33a62ed433d4a3)) >> np.uint64(11)).astype(float)+.5)*scale
    return np.sqrt(-np.log(u))*np.exp(2j*np.pi*sign*v)


class Box:
    def __init__(self, n, length):
        self.n, self.length = n, length
        self.shape = (n, n, n)
        line = np.rint(np.fft.fftfreq(n)*n).astype(int)
        self.indices = line[:, None, None], line[None, :, None], np.arange(n//2+1)[None, None, :]
        x, y, z = self.indices
        edge = np.maximum(np.maximum(abs(x), abs(y)), abs(z))
        self.mask = 3*edge < n
        self.tail = (edge >= .8*((n-1)//3)) & self.mask
        self.k2 = (2*np.pi/length)**2*(x*x+y*y+z*z)
        self.k = np.sqrt(self.k2)
        self.multiplicity = np.full(n//2+1, 2.)
        self.multiplicity[[0, -1]] = 1.
        self.bin_index = np.clip(np.searchsorted(BINS, self.k, side="right")-1,
                                 0, len(BINS)-2).ravel()

    def to_real(self, field):
        return irfftn(field, s=self.shape, axes=(0, 1, 2), norm="forward", workers=1)

    def to_modes(self, field):
        return rfftn(field, axes=(0, 1, 2), norm="forward", workers=1)

    def seed_modes(self, config):
        a = config.a_start
        speed = rotation_speed(a)
        omega = np.sqrt(1+self.k2/a**2)
        power = np.zeros_like(self.k)
        keep = self.mask & (self.k > 0)
        power[keep] = (2*np.pi**2*A_S*speed**2*config.seed_scale**2*window(self.k[keep])
                       / (self.k[keep]**3*self.length**3))
        spacing = round(96/self.length)
        if not np.isclose(spacing*self.length, 96):
            raise ValueError("Seed labels require a shared 96-unit momentum cell.")
        indices = [spacing*label for label in self.indices]
        theta = np.sqrt(power)/omega*gaussian(*indices, config.seed)
        velocity = np.sqrt(power)*gaussian(*indices, config.seed+104729)
        velocity[0, 0, 0] = speed
        return a*theta, a*a*velocity + config.h*theta

    def acceleration(self, chi, a):
        theta = self.to_real(chi)/a
        force = -self.k2*chi - a**3*self.to_modes(np.sin(theta))*self.mask
        return force, theta

    def integrate_step(self, chi, canonical_speed, force, a, dt):
        canonical_speed += .5*dt*force
        chi += dt*canonical_speed
        force, theta = self.acceleration(chi, a)
        canonical_speed += .5*dt*force
        return force, theta

    def energy(self, chi, canonical_speed, theta, a, h):
        angular_velocity = (canonical_speed-h/a*chi)/a**2
        kinetic = .5*np.sum(abs(angular_velocity)**2*self.multiplicity)
        gradient = .5*np.sum(self.k2*abs(chi)**2*self.multiplicity)/a**4
        potential = np.mean(2*np.sin(theta/2)**2)
        density = kinetic+gradient+potential
        pressure = kinetic-gradient/3-potential
        mean_kinetic = .5*abs(angular_velocity[0, 0, 0])**2
        inhomogeneous = (kinetic-mean_kinetic+gradient)/density
        return angular_velocity, density, pressure, inhomogeneous

    def late_diagnostic(self, chi, angular_velocity, theta, density, a, h):
        branch = 2*np.pi*np.rint((chi[0, 0, 0].real/a)/(2*np.pi))
        centered = theta-branch
        modes = chi/a
        modes[0, 0, 0] -= branch
        omega = np.sqrt(1+self.k2/a**2)
        hubble = h/a**2
        freq = np.sqrt(omega**2+.75*hubble**2)
        u_velocity = angular_velocity+1.5*hubble*modes
        action = .5*a**3*(abs(u_velocity)**2+freq**2*abs(modes)**2)/freq*self.multiplicity
        inst = .5*a**3*np.sum((abs(angular_velocity)**2+omega**2*abs(modes)**2)
                              /omega*self.multiplicity)
        count = float(action.sum())
        spec = np.bincount(self.bin_index, weights=action.ravel(), minlength=len(BINS)-1)
        exact = 2*np.sin(centered/2)**2
        anharmonic = np.mean(abs(exact-.5*centered**2))/density
        return (np.array([a, count/Q, inst/Q, action[0, 0, 0]/count,
                          np.sqrt(np.sum(self.k2*action)/count)/a,
                          anharmonic, np.mean(abs(centered) > np.pi),
                          np.sum(action[self.tail])/count, hubble]),
                spec/Q)


def evolve(config):
    box = Box(config.n, config.length)
    chi, speed = box.seed_modes(config)
    a0, h, dt = config.a_start, config.h, config.dt
    steps = round((config.a_end-a0)/(h*dt))
    if not np.isclose(a0+h*steps*dt, config.a_end):
        raise ValueError("The conformal steps must tile the requested interval.")
    a = a0
    force, theta = box.acceleration(chi, a)
    ledger, late, spectra = [], [], []
    ledger_stride = max(1, round(.10/dt))
    late_stride = max(1, round(.01/(h*dt)))
    for step in range(steps+1):
        if step:
            a = a0+h*step*dt
            force, theta = box.integrate_step(chi, speed, force, a, dt)
        sample_late = a >= 3.5-1e-12 and (step % late_stride == 0 or step == steps)
        if step % ledger_stride == 0 or step == steps or sample_late:
            velocity, rho, pressure, inhomogeneous = box.energy(chi, speed, theta, a, h)
            ledger.append([a, a**3*rho, -3*a*a*pressure, inhomogeneous,
                           a**3*velocity[0, 0, 0].real])
            if sample_late:
                row, spectrum = box.late_diagnostic(chi, velocity, theta, rho, a, h)
                late.append(row)
                spectra.append(spectrum)
    ledger, late, spectra = np.asarray(ledger), np.asarray(late), np.asarray(spectra)
    integrated = cumulative_simpson(ledger[:, 2], x=ledger[:, 0], initial=0.)
    defect = ledger[:, 1]-ledger[0, 1]-integrated
    action = late[:, 1]
    mean = np.mean(action)
    checks = dict(
        drift=float(np.ptp(action)/mean),
        dictionaries=float(np.max(abs(late[:, 2]/action-1))),
        anharmonic=float(np.max(late[:, 5])),
        out_of_branch=float(np.max(late[:, 6])),
        uv_tail=float(np.max(late[:, 7])),
        adiabatic=float(np.max(late[:, 8])),
        energy=float(np.max(abs(defect))/ledger[-1, 1]),
    )
    tolerances = dict(drift=.02, dictionaries=.02, anharmonic=.02,
                      out_of_branch=1e-4, uv_tail=.01, adiabatic=.01, energy=.01)
    failures = [name for name, value in checks.items() if value >= tolerances[name]]
    row = dict(config=asdict(config), steps=steps, particle_gate=not failures,
               misses=failures, checks=checks,
               C=float(mean) if not failures else None,
               late_quadratic_proxy=float(mean), late_proxy_min=float(min(action)),
               late_proxy_max=float(max(action)),
               coherent_fraction=float(np.mean(late[:, 3])) if not failures else None,
               rms_p_over_m_final=float(late[-1, 4]) if not failures else None,
               max_gradient_and_kinetic_fraction=float(np.max(ledger[:, 3])),
               initial_inhomogeneous_fraction=float(ledger[0, 3]),
               final_comoving_energy=float(ledger[-1, 1]),
               spectrum_over_charge=np.mean(spectra, axis=0).tolist() if not failures else None)
    trace = dict(ledger=ledger, pressure_work=integrated, energy_defect=defect,
                 late=late, late_quadratic_spectrum=spectra, comoving_edges=BINS)
    return row, trace


def run(name, directory):
    row, trace = evolve(CASES[name])
    path = directory / (name+".tmp.npz")
    directory.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **trace)
    path.replace(directory / (name+".npz"))
    write_json(directory / (name+".json"), row)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=CASES)
    parser.add_argument("--output", type=Path, default=HERE / "results")
    args = parser.parse_args()
    gate = args.output / "qualification.json"
    if not gate.exists() or not json.loads(gate.read_text())["passed"]:
        raise SystemExit("Qualify the homogeneous and linear ladder before the 3D runs.")
    print(json.dumps(run(args.case, args.output), indent=2))


if __name__ == "__main__":
    main()
