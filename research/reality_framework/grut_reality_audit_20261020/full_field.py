"""Separate full-field test of Candidate 1 bulk waves and growth.

Periodic Fourier spatial differentiation and classical RK4 time integration.
No rotating-frame forces, imposed background, tangent solver or damping.
Fields depend on one coordinate: plane perturbations of uniform 3D bulk.
This is not an unrestricted 3D droplet/collision stability calculation.
"""
from pathlib import Path
import hashlib
import json
import platform
import time
import numpy as np

ROOT=Path(__file__).resolve().parent


def execute(spec):
    n=spec['N']; length=2*np.pi/spec['k']; dt=spec['dt']
    x=np.arange(n)*length/n
    cosine=np.cos(spec['k']*x)
    k2=(2*np.pi*np.fft.rfftfreq(n,d=length/n))**2
    a,b=1.,.1
    if spec['background']=='dense':
        f=(a/b)**.25; s=0.; omega=(a*b)**.25
    else:
        f=.5; s=np.sqrt(1-f*f/a); omega=np.sqrt(s*s+b*f*f)
    state=np.zeros((6,n))
    state[0]=f; state[2]=s
    if spec['seed']=='charged':state[0]+=spec['epsilon']*cosine
    else:state[2]+=spec['epsilon']*cosine
    state[4]=-omega*state[0]

    def lap(fields):
        return np.fft.irfft(-k2*np.fft.rfft(fields,axis=-1),n=n,axis=-1)

    def rhs(y):
        u,w,organizer=y[:3]
        f2=u*u+w*w
        force=lap(y[:3])
        force[:2]-=(organizer*organizer+b*f2)*y[:2]
        force[2]-=(a*(organizer*organizer-1)+f2)*organizer
        return np.concatenate((y[3:],force),axis=0)

    def measure(y,t):
        u,w,organizer,pu,pw,ps=y
        f2=u*u+w*w
        potential=.25*a*(organizer*organizer-1)**2+.5*organizer*organizer*f2+.25*b*f2*f2
        energy=length*np.mean(.5*np.sum(y[3:]**2,axis=0)+potential-.5*np.sum(y[:3]*lap(y[:3]),axis=0))
        charge=length*np.mean(w*pu-u*pw)
        z=u+1j*w; p=pu+1j*pw; mean=z.mean()
        phase=mean.conjugate()/abs(mean)
        rotating=z*phase
        measured_omega=-float((p.mean()/mean).imag)
        rotating_velocity=phase*(p+1j*measured_omega*z)
        values=[2*np.mean(v*cosine) for v in [rotating.real,rotating.imag,organizer,
                    rotating_velocity.real,rotating_velocity.imag,ps]]
        return [float(t),float(energy),float(charge),*map(float,values),
                float(np.linalg.norm(values)),float(abs(mean)),measured_omega]

    steps=round(spec['T']/dt); period=round(spec['sample_period']/dt)
    assert abs(steps*dt-spec['T'])<1e-9 and abs(period*dt-spec['sample_period'])<1e-9
    trace=[]; start=time.monotonic()
    for step in range(steps+1):
        if step%period==0:
            row=measure(state,step*dt)
            if not np.isfinite(row).all():raise RuntimeError('Nonfinite full-field state')
            trace.append(row)
        if step<steps:
            k1=rhs(state)
            k2v=rhs(state+.5*dt*k1)
            k3=rhs(state+.5*dt*k2v)
            k4=rhs(state+dt*k3)
            state+=(dt/6)*(k1+2*k2v+2*k3+k4)
    array=np.array(trace)
    output=dict(specification=spec,geometry='Periodic plane-dependent solution of the full classical scalar PDE; energy and charge per unit transverse area.',
                integrator='Classical explicit RK4',spatial_operator='Fourier pseudospectral; no filtering or imposed tangent evolution',
                columns=['time','energy','charge','eta_mode','xi_mode','S_mode','eta_velocity_mode','xi_velocity_mode','S_velocity_mode',
                         'mode_norm','mean_charged_amplitude','mean_rotation_frequency'],
                trace=trace,maximum_relative_energy_drift=float(np.max(abs(array[:,1]/array[0,1]-1))),
                maximum_relative_charge_drift=float(np.max(abs(array[:,2]/array[0,2]-1))),
                elapsed_seconds=time.monotonic()-start,
                runtime=dict(python=platform.python_version(),numpy=np.__version__),
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    target=ROOT/'results'/(spec['id']+'.json'); target.parent.mkdir(exist_ok=True)
    target.write_text(json.dumps(output,separators=(',',':'))+'\n')
    return dict(id=spec['id'],seconds=output['elapsed_seconds'],energy=output['maximum_relative_energy_drift'],charge=output['maximum_relative_charge_drift'])


if __name__=='__main__':
    import sys
    protocol=json.loads((ROOT/'protocol.json').read_text())
    chosen=next(s for s in protocol['runs'] if s['id']==sys.argv[1])
    print(json.dumps(execute(chosen)))
