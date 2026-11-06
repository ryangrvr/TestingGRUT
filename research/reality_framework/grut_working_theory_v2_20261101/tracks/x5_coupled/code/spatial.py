#!/usr/bin/env python3
"""Self-consistent scalar work and all occupied Dirac orbitals in a flat 1D box."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

LENGTH = 2 * math.pi
MASS = 0.7
PHASE_SCALE = 1.15
SCALAR_MASS = 2.6
DURATION = 2.5


class Box:
    def __init__(self, points):
        self.points = points
        self.dx = LENGTH / points
        self.x = np.arange(points) * self.dx
        self.k = 2 * math.pi * np.fft.fftfreq(points, self.dx)
        self.freq = np.hypot(self.k, SCALAR_MASS)
        self.vac_energy = -float(np.hypot(self.k, MASS).sum())
        f = np.fft.fft(np.eye(points), axis=0, norm='ortho')
        self.momentum = f.conj().T @ (self.k[:, None] * f)

    def vacuum(self):
        negative = np.zeros((self.points, 2, self.points), complex)
        for i, k in enumerate(self.k):
            _, vec = np.linalg.eigh([[k, MASS], [MASS, -k]])
            negative[i, :, i] = vec[:, 0]
        return np.fft.ifft(negative, axis=0, norm='ortho')

    def matrix(self, phi):
        out = np.zeros((2*self.points, 2*self.points), complex)
        out[0::2, 0::2] = self.momentum
        out[1::2, 1::2] = -self.momentum
        even = 2*np.arange(self.points)
        mass = MASS * np.exp(1j*phi/PHASE_SCALE)
        out[even, even+1] = mass
        out[even+1, even] = mass.conj()
        return out

    def energies(self, phi, speed, orbitals):
        scalar_k = np.fft.fft(phi, norm='ortho')
        wave_energy = 0.5 * self.dx * (np.sum(speed*speed) + np.sum(abs(self.freq*scalar_k)**2))
        p = np.fft.fft(orbitals, axis=0, norm='ortho')
        kinetic = np.sum(self.k[:, None]*(abs(p[:,0,:])**2-abs(p[:,1,:])**2))
        mass = 2*MASS*np.sum((np.exp(1j*phi/PHASE_SCALE)[:, None] *
                              orbitals[:,0,:].conj()*orbitals[:,1,:]).real)
        return float(wave_energy), float(kinetic+mass-self.vac_energy)

    def projected_energy(self, phi, orbitals):
        energy, axes = np.linalg.eigh(self.matrix(phi))
        assert energy[self.points-1] < 0 < energy[self.points]
        weights = np.sum(abs(axes.conj().T @ orbitals.reshape(2*self.points,self.points))**2, axis=1)
        holes = 1 - weights[:self.points]
        particles = weights[self.points:]
        excitation = energy[self.points:] @ particles - energy[:self.points] @ holes
        return dict(energy=float(excitation), positive_number=float(sum(particles)),
                    hole_number=float(sum(holes)),
                    instantaneous_vacuum_shift=float(sum(energy[:self.points])-self.vac_energy),
                    max_occupation=float(max(weights)))


def axial(orbitals):
    return np.sum(abs(orbitals[:,0,:])**2-abs(orbitals[:,1,:])**2, axis=1)


def evolve(box, phi, speed, orbitals, dt, steps, record):
    c, s = np.cos(box.freq*dt/2), np.sin(box.freq*dt/2)
    fermion_freq = np.hypot(box.k, MASS)[:,None]
    fermion_c = np.cos(fermion_freq*dt/2)
    fermion_s = -1j*np.sin(fermion_freq*dt/2)/fermion_freq
    k = box.k[:,None]

    def free(field, velocity, states):
        q = np.fft.fft(field)
        p = np.fft.fft(velocity)
        field = np.fft.ifft(c*q+s*p/box.freq).real
        velocity = np.fft.ifft(c*p-s*box.freq*q).real
        wave = np.fft.fft(states,axis=0)
        top, bottom = wave[:,0,:], wave[:,1,:]
        wave = np.stack([fermion_c*top+fermion_s*(k*top+MASS*bottom),
                         fermion_c*bottom+fermion_s*(MASS*top-k*bottom)],axis=1)
        states = np.fft.ifft(wave,axis=0)
        return field, velocity, states

    for step in range(steps):
        phi, speed, orbitals = free(phi,speed,orbitals)
        phase = np.exp(1j*phi/PHASE_SCALE)
        # Put the full uniform vacuum Hamiltonian in the free step. Here
        # the residual mass vanishes at phi=0, preserving the exact vacuum.
        mass_delta = MASS*(phase-1)
        b = np.column_stack([mass_delta.real,-mass_delta.imag,np.zeros(box.points)])
        omega = np.linalg.norm(b,axis=1)
        axis = np.divide(b,omega[:,None],out=np.zeros_like(b),where=omega[:,None]>0)
        cross = np.sum(orbitals[:,0,:].conj()*orbitals[:,1,:],axis=1)
        bloch = np.column_stack([2*cross.real,2*cross.imag,axial(orbitals)])
        longitudinal = np.sum(bloch*axis,axis=1)[:,None]*axis
        sin_int = dt*np.sinc(2*omega*dt/np.pi)
        cos_int = omega*dt**2*np.sinc(omega*dt/np.pi)**2
        integrated = sin_int[:,None]*bloch+(dt-sin_int)[:,None]*longitudinal
        integrated += cos_int[:,None]*np.cross(axis,bloch)
        speed += (MASS/(PHASE_SCALE*box.dx)) * (
            phase.imag*integrated[:,0]+phase.real*integrated[:,1])

        mc = np.cos(omega*dt)[:,None]
        ms = -1j*dt*np.sinc(omega*dt/np.pi)[:,None]*mass_delta[:,None]
        top = mc*orbitals[:,0,:]+ms*orbitals[:,1,:]
        bottom = mc*orbitals[:,1,:]-ms.conj()*orbitals[:,0,:]
        orbitals = np.stack([top,bottom],axis=1)
        phi, speed, orbitals = free(phi,speed,orbitals)
        if record is not None:
            record(step+1,phi,speed,orbitals)
    return phi, speed, orbitals


def ode_control():
    box = Box(9)
    phi0 = 0.04*np.sin(box.x)
    speed0 = 0.6*np.cos(box.x)+0.17
    orb0 = box.vacuum()
    points = box.points
    duration = 0.5
    y0 = np.r_[phi0, speed0, orb0.ravel()]

    def rhs(_, y):
        phi, speed = y[:points].real, y[points:2*points].real
        orbitals = y[2*points:].reshape(points,2,points)
        mat = box.matrix(phi)
        state_gradient = 2*MASS/PHASE_SCALE * np.sum(
            (np.exp(1j*phi[:,None]/PHASE_SCALE)*orbitals[:,0,:].conj()*orbitals[:,1,:]).imag,
            axis=1)
        acceleration = -np.fft.ifft(box.freq**2*np.fft.fft(phi)).real+state_gradient/box.dx
        dpsi = -1j * (mat @ orbitals.reshape(2*points,points))
        return np.r_[speed, acceleration, dpsi.ravel()]

    ref = solve_ivp(rhs,(0,duration),y0,rtol=3e-13,atol=1e-14,method='DOP853')
    assert ref.success
    records = []
    for steps in (100,200,400):
        fields, speed, orb = evolve(box,phi0,speed0,orb0,duration/steps,steps,None)
        error = float(np.max(abs(np.r_[fields,speed,orb.ravel()]-ref.y[:,-1])))
        records.append(dict(steps=steps, max_absolute_error=error))
    assert records[-1]['max_absolute_error'] < 2e-5
    assert records[-1]['max_absolute_error'] < records[0]['max_absolute_error']/10
    return records


def run_case(points, steps, amplitude, preparation, reverse=False):
    box = Box(points)
    phi = np.zeros(points)
    speed = amplitude * (np.cos(box.x) if preparation=='wave' else np.ones(points))
    orbitals = box.vacuum()
    start = phi,speed,orbitals
    E0 = sum(box.energies(*start))
    samples = []

    def record(step, field, velocity, states):
        if step % (steps//50) != 0 and step != steps:
            return
        scalar, fermion = box.energies(field,velocity,states)
        samples.append(dict(t=step*DURATION/steps, scalar_energy=scalar, fermion_energy=fermion,
                            budget_error=(scalar+fermion-E0)/max(1,abs(E0)),
                            mean_phase=float(np.mean(field)/PHASE_SCALE),
                            phase_std=float(np.std(field)/PHASE_SCALE)))

    record(0,*start)
    final = evolve(box,*start,DURATION/steps,steps,record)
    projected = box.projected_energy(final[0],final[2])
    expected_e = projected['energy']+projected['instantaneous_vacuum_shift']
    assert abs(expected_e-samples[-1]['fermion_energy']) < 1e-9
    matrix = final[2].reshape(2*points,points)
    ortho = float(np.max(abs(matrix.conj().T @ matrix-np.eye(points))))
    momentum = np.fft.fft(final[2],axis=0,norm='ortho')
    off_diag = ~np.eye(points,dtype=bool)
    mixing = float(np.sum(abs(momentum.transpose(0,2,1)[off_diag])**2)/points)
    inversion_error = None
    if reverse:
        back = evolve(box,*final,-DURATION/steps,steps,None)
        inversion_error = max(float(np.max(abs(a-b))) for a,b in zip(back,start,strict=True))
        assert inversion_error < 2e-10
    return dict(points=points, steps=steps, initial_amplitude=amplitude, preparation=preparation,
                initial_energy=E0, final=projected, max_budget_error=max(abs(s['budget_error']) for s in samples),
                max_mean_phase=max(abs(s['mean_phase']) for s in samples),
                max_phase_std=max(s['phase_std'] for s in samples),
                orthonormality_error=ortho, momentum_mixing=mixing, inversion_error=inversion_error,
                samples=samples)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    outcomes = dict(
        wave17=run_case(17,500,1.2,'wave'),
        wave33=run_case(33,500,1.2,'wave',True),
        wave65=run_case(65,500,1.2,'wave'),
        wave33_half_step=run_case(33,1000,1.2,'wave'),
        mean33=run_case(33,500,0,'mean'),
        uniform17=run_case(17,500,1.2,'mean'),
        weak33=run_case(33,500,0.15,'wave'),
        weaker33=run_case(33,500,0.075,'wave'),
    )
    baseline = outcomes['wave33']
    fine = outcomes['wave65']
    half = outcomes['wave33_half_step']
    assert baseline['final']['energy'] > 1e-3
    assert abs(outcomes['mean33']['final']['energy']) < 1e-9
    assert outcomes['uniform17']['momentum_mixing'] < 1e-20
    assert max(row['orthonormality_error'] for row in outcomes.values()) < 1e-10
    assert max(row['max_budget_error'] for row in outcomes.values()) < 2e-4
    assert max(abs(case['final']['energy']/baseline['final']['energy']-1) for case in (fine,half)) < 0.02
    ode = ode_control()
    output = dict(
        scope='closed finite flat 1D Dirac/scalar control',
        parameters=dict(length=LENGTH, fermion_mass=MASS, phase_scale=PHASE_SCALE,
                        scalar_mass=SCALAR_MASS, duration=DURATION),
        unrotated_ODE=ode, cases=outcomes, scientific_acceptance=False,
    )
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    for name,row in outcomes.items():
        print(name, 'energy', row['final']['energy'], 'conservation', row['max_budget_error'],flush=True)


if __name__ == '__main__':
    main()
