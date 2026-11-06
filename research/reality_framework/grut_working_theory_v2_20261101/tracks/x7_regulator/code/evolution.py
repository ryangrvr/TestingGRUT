#!/usr/bin/env python3
"""Evolve the physical and signed vacuum-response sectors with a common source."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass, replace
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid

from regulator import Configuration, HEAVIER, MASSES, RadialSlab, coefficients
from slab import BETA, State as SlabState, adjoint

DURATION = 3.
CASES = {
    'cylinder_9': (9,12,0.5,49,3,1200,MASSES),
    'cylinder_12': (12,16,0.5,65,3,1200,MASSES),
    'cylinder_16': (16,22,0.5,81,3,1200,MASSES),
    'longitudinal': (16,16,0.5,81,3,1200,MASSES),
    'transverse_cutoff': (12,22,0.5,65,3,1200,MASSES),
    'transverse_bins': (12,16,0.25,65,3,1200,MASSES),
    'half_step': (12,16,0.5,65,3,2400,MASSES),
    'spatial_quadrature': (12,16,0.5,129,3,1200,MASSES),
    'scalar_bandwidth': (12,16,0.5,65,5,1200,MASSES),
    'heavier': (12,16,0.5,65,3,1200,HEAVIER),
}
RADIAL_CASES = {
    'uniform_fine': ((12,16,0.125,65,3,1200,MASSES),None),
    'gauss_reference': ((12,16,0.5,65,3,1200,MASSES),8),
    'gauss_refined': ((12,16,0.5,65,3,1200,MASSES),12),
    'gauss_cutoff': ((16,22,0.5,81,3,1200,MASSES),8),
    'gauss_cutoff_refined': ((16,22,0.5,81,3,1200,MASSES),12),
    'gauss_heavier': ((16,22,0.5,81,3,1200,HEAVIER),8),
}


@dataclass(frozen=True)
class State:
    a: float
    phi: np.ndarray
    pi: np.ndarray
    orbitals: tuple


@dataclass(frozen=True)
class Expectation:
    kinetic: np.ndarray
    mass: float
    comoving_momentum: float
    force: np.ndarray


class RegulatedSlab:
    def __init__(self,cfg,auxiliaries,gauss_order=None):
        self.cfg = cfg
        self.masses,self.signed = coefficients(cfg.mass,auxiliaries)
        self.sectors = [RadialSlab(replace(cfg,mass=m,subtract=False,expand=False),gauss_order)
                        for m in self.masses]
        self.spatial = self.sectors[0]
        self.weight = self.spatial.weight
        self.qx2 = self.spatial.k**2
        self.qperp2 = self.spatial.perp[:,None]**2
        self.q2 = self.spatial.q2
        self.mass2 = self.masses[:,None,None]**2
        self.signed_weight = self.signed[:,None,None]*self.weight[None,:,None]

    def matching(self,a):
        p2 = self.q2/a**2
        e = np.sqrt(self.mass2+p2)
        w = self.signed_weight*self.mass2/(4*a**3)
        z = np.sum(w/e**3)/self.cfg.phase_scale**2
        slope = 3*np.sum(w*self.mass2/e**5)/self.cfg.phase_scale**2
        zg = np.sum(w*p2/e**5)
        zg_slope = -5*np.sum(w*self.mass2*p2/e**7)
        if z>=1:
            raise ValueError('Nonpositive matched kinetic coefficient')
        return float(1-z),float(slope),float(zg),float(zg_slope),e

    def observe(self,state):
        phase = state.phi/self.cfg.phase_scale
        kinetic = np.zeros(2)
        mass = momentum = 0.
        force = np.zeros_like(phase)
        for c,box,occ in zip(self.signed,self.sectors,state.orbitals,strict=True):
            if c==0:
                continue
            blocks = box.blocks(occ)
            diagonal = blocks @ adjoint(blocks)

            def trace(matrix):
                return box.weight @ np.einsum('...ij,...ji->...',matrix,diagonal).real.sum(axis=1)

            kinetic += c*np.array([trace(box.kx),trace(box.kt)])
            cov = box.cov(occ)
            mass += c*(trace(box.cfg.mass*BETA)+np.einsum('ij,ji->',box.residual(phase),cov).real)
            count = np.trace(diagonal,axis1=-2,axis2=-1).real-2
            momentum += c*(box.weight @ (count @ box.k))/2
            force += c*box.mass_derivative(phase,cov)
        return Expectation(kinetic,mass,momentum,self.spatial.project_scalar(force))

    def measure(self,state,expectation,work=False):
        cfg,a,phi = self.cfg,state.a,state.phi
        A,slope,zg,zg_slope,e = self.matching(a)
        v = state.pi/(a**3*A)
        dx = self.spatial.derivative(phi)
        speed2 = np.mean(v*v)
        grad2 = np.mean(dx*dx)/a**2
        potential = cfg.vacuum_density+0.5*cfg.scalar_mass**2*np.mean(phi*phi)
        density_scalar = potential+0.5*A*(speed2+grad2)
        pressure_scalar = 0.5*A*speed2-potential+slope*speed2/6-(A+slope)*grad2/6
        Kx,Kt = expectation.kinetic/a
        rho_F = (Kx+Kt+expectation.mass+2*np.sum(self.signed_weight*e))/(2*a**3)
        px = (Kx+2*np.sum(self.signed_weight*self.qx2/(a**2*e)))/(2*a**3)
        pt = (Kt+2*np.sum(self.signed_weight*self.qperp2/(a**2*e)))/(4*a**3)
        P_F = (px+2*pt)/3
        rho = density_scalar+rho_F
        P = pressure_scalar+P_F
        if rho<=0 or 3*cfg.planck**2+zg/2<=0:
            raise ValueError('Invalid density or gravitational matching')
        H = math.sqrt(rho/(3*cfg.planck**2+zg/2)) if cfg.expand else 0.
        z = 1-A
        rho_ct = 0.5*z*(speed2+grad2)
        P_ct = (z/2-slope/6)*speed2-(z-slope)*grad2/6
        result = dict(a=a,H=H,rho=rho,P=P,rho_F=rho_F,P_F=P_F,rho_F_matched=rho_F-rho_ct,
                      P_F_matched=P_F-P_ct,P_F_x=px,P_F_perp=pt,
                      P_x=pressure_scalar+2*A*grad2/3+px,
                      P_perp=pressure_scalar-A*grad2/3+pt,
                      rho_scalar=density_scalar,scalar_gradient=0.5*A*grad2,
                      momentum=expectation.comoving_momentum/a**4-A*np.mean(v*dx)/a,
                      A=A,Z_g=zg)
        if work:
            force = expectation.force/a**3
            laplace = self.spatial.derivative(dx)/a**2
            acceleration = (A*laplace-cfg.scalar_mass**2*phi-force-H*(3*A+slope)*v)/A
            matched = force-z*acceleration-H*(3*z-slope)*v+z*laplace
            hdot = -(3*(rho+P)+zg_slope*H*H/2)/(6*cfg.planck**2+zg) if cfg.expand else 0.
            result.update(Jv=np.mean(v*force),Jv_matched=np.mean(v*matched),
                          acceleration=acceleration,Hdot=hdot)
        return result

    def initial(self,phi,speed):
        cfg = self.cfg
        phi = self.spatial.project_scalar(phi)
        speed = self.spatial.project_scalar(speed)
        A,*_ = self.matching(1.)
        pi = A*speed
        gx = self.spatial.derivative(phi)
        classical = cfg.vacuum_density+0.5*np.mean(speed**2+cfg.scalar_mass**2*phi**2+gx**2)
        hubble = math.sqrt(classical/(3*cfg.planck**2)) if cfg.expand else 0.
        spectra = [np.linalg.eigh(box.full_hamiltonian(1.,phi/cfg.phase_scale)) for box in self.sectors]
        for box,(e,_) in zip(self.sectors,spectra,strict=True):
            if np.any(e[:,box.d-1]>=0) or np.any(e[:,box.d]<=0):
                raise ValueError('Prepared state has no occupied/unoccupied gap')
        first = []
        for box,(e,basis) in zip(self.sectors,spectra,strict=True):
            _,mix = box.first_graph(1.,phi,speed,hubble,e,basis)
            first.append(box.graph_orbitals(basis,mix))
        preliminary = State(1.,phi,pi,tuple(first))
        diagnostics = self.measure(preliminary,self.observe(preliminary),work=True)
        H,Hdot,acceleration = (diagnostics[key] for key in ('H','Hdot','acceleration'))

        def second_graph(box,e,basis):
            gap,mix = box.first_graph(1.,phi,speed,H,e,basis)

            def embedded(t):
                a = math.exp(H*t+0.5*Hdot*t*t)
                field = phi+speed*t+0.5*acceleration*t*t
                velocity = speed+acceleration*t
                ev,u = np.linalg.eigh(box.full_hamiltonian(a,field/cfg.phase_scale))
                _,x = box.first_graph(a,field,velocity,H+Hdot*t,ev,u)
                block = u[:,:,box.d:] @ x @ adjoint(u[:,:,:box.d])
                return block+adjoint(block)

            eps = 2e-5
            tangent = (embedded(eps)-embedded(-eps))/(2*eps)
            mix += 1j*(adjoint(basis[:,:,box.d:]) @ tangent @ basis[:,:,:box.d])/gap
            return box.graph_orbitals(basis,mix)

        return replace(preliminary,orbitals=tuple(second_graph(box,*spectrum)
                       for box,spectrum in zip(self.sectors,spectra,strict=True)))

    def linear(self,state,dt):
        cfg,a = self.cfg,state.a
        A,*_ = self.matching(a)
        inertia = a**3*A
        omega = np.sqrt(cfg.scalar_mass**2/A+(self.spatial.scalar_k/a)**2)
        cos,sin = np.cos(omega*dt),np.sin(omega*dt)
        phi,v = np.fft.fft(state.phi),np.fft.fft(state.pi)/inertia
        occupied = tuple(box.free_step(a,occ,dt)
                         for box,occ in zip(self.sectors,state.orbitals,strict=True))
        return State(a,np.fft.ifft(cos*phi+sin*v/omega).real,
                     inertia*np.fft.ifft(cos*v-sin*omega*phi).real,occupied)

    def matter(self,state,dt):
        mid = self.linear(state,dt/2)
        phase = mid.phi/self.cfg.phase_scale
        occupied,impulses = zip(*(box.local(phase,occ,dt,True)
                               for box,occ in zip(self.sectors,mid.orbitals,strict=True)),strict=True)
        return self.linear(replace(mid,pi=mid.pi+self.signed @ np.array(impulses),orbitals=occupied),dt/2)

    def expand(self,state,expectation,dt):
        if not self.cfg.expand:
            return state,np.zeros(3)
        a = state.a
        H = self.measure(state,expectation)['H']
        middle = self.measure(replace(state,a=a*math.exp(H*dt/2)),expectation)
        dvol = a**3*math.expm1(3*middle['H']*dt)
        pressure = np.array([middle[key] for key in ('P','P_F','P_F_matched')])
        return replace(state,a=a*math.exp(middle['H']*dt)),dvol*pressure

    def step(self,state,expectation,dt):
        before,first = self.expand(state,expectation,dt/2)
        after = self.matter(before,dt)
        expectation = self.observe(after)
        state,second = self.expand(after,expectation,dt/2)
        return state,expectation,first+second


def evaluate(cfg,steps,masses,duration=DURATION,gauss_order=None):
    system = RegulatedSlab(cfg,masses,gauss_order)
    x = 2*np.pi*system.spatial.x/cfg.length
    state = system.initial(0.15*np.sin(x),0.9*np.cos(x))
    expectation = system.observe(state)
    values = []
    pressure = np.zeros(3)
    dt = duration/steps
    for i in range(steps+1):
        result = system.measure(state,expectation,work=True)
        result.pop('acceleration')
        values.append(dict(t=i*dt,work_P=pressure[0],work_P_F=pressure[1],
                           work_P_F_matched=pressure[2],**result))
        if i<steps:
            state,expectation,increment = system.step(state,expectation,dt)
            pressure += increment
    d = {k:np.array([v[k] for v in values]) for k in values[0]}
    vol = d['a']**3
    scale = vol[0]*d['rho'][0]
    energy = vol*d['rho']+d['work_P']
    momentum = d['a']**4*d['momentum']
    orthogonal = []
    for box,occ in zip(system.sectors,state.orbitals,strict=True):
        amplitudes = occ.reshape(len(box.perp),2*box.d,box.d)
        orthogonal.append(float(np.max(abs(adjoint(amplitudes) @ amplitudes-np.eye(box.d)))))
    summary = dict(
        configuration=asdict(cfg), auxiliary_masses=list(masses), dt=dt, elapsed=duration,
        transverse_rule=dict(gauss_order_per_panel=gauss_order,nodes=len(system.spatial.perp)),
        final={k:float(v[-1]) for k,v in d.items()},
        max_pressure_work_error=float(np.max(abs(energy/scale-1))),
        max_comoving_momentum_error=float(np.max(abs(momentum-momentum[0]))/scale),
        orthonormality_error_by_sector=orthogonal,
        max_fermion_anisotropic_stress=float(np.max(abs(d['P_F_x']-d['P_F_perp']))),
        max_scalar_gradient_fraction=float(np.max(d['scalar_gradient']/d['rho'])),
        min_scalar_inertia=float(np.min(d['A'])),
    )
    for suffix in ('','_matched'):
        instantaneous = vol*d['Jv'+suffix]
        work = cumulative_trapezoid(instantaneous,d['t'],initial=0)
        stress = vol*d['rho_F'+suffix]
        stress = stress-stress[0]+d['work_P_F'+suffix]
        d['comoving_work'+suffix] = instantaneous
        d['integrated_transfer'+suffix] = work
        summary.update({
            'integrated_transfer'+suffix:float(work[-1]),
            'integrated_absolute_transfer'+suffix:float(trapezoid(abs(instantaneous),d['t'])),
            'transfer_from_stress'+suffix:float(stress[-1]),
            'max_fermion_work_error'+suffix:float(np.max(abs(work-stress))/scale),
            'peak_instantaneous_work'+suffix:float(np.max(abs(instantaneous))),
        })
    return summary,d


def pulse(masses,dt):
    cfg = replace(Configuration(),longitudinal=9,transverse_cut=12.,transverse_bins=24,
                  points=49,expand=False)
    system = RegulatedSlab(cfg,masses)
    nsteps = round(8/dt)
    dt = 8/nsteps
    end_energies = []
    for box in system.sectors:
        zero = np.zeros_like(box.x)
        _,basis = np.linalg.eigh(box.full_hamiltonian(1.,zero))
        occ = basis[:,:,:box.d].reshape(len(box.perp),box.d,2,box.d)
        for i in range(nsteps):
            t = (i+0.5)*dt
            phase = 0.22*np.sin(np.pi*t/8)**4*np.cos(box.x-2.6*t)
            occ = box.free_step(1.,occ,dt/2)
            occ,_ = box.local(phase,occ,dt,False)
            occ = box.free_step(1.,occ,dt/2)
        end_energies.append(box.measure(SlabState(1.,zero,zero,occ))['rho_F'])
    weighted = system.signed*np.array(end_energies)
    ratio = sum(abs(weighted[1:]))/weighted[0]
    return dict(dt=dt,auxiliary_masses=list(masses),end_energy=end_energies,
                weighted_end_energy=weighted.tolist(),auxiliary_ratio=float(ratio),gate_pass=bool(ratio<0.01))


def check_expectation():
    cfg = replace(Configuration(),longitudinal=3,transverse_cut=4.,transverse_bins=8,expand=False)
    system = RegulatedSlab(cfg,MASSES)
    box = system.spatial
    state = system.initial(0.15*np.sin(box.x),0.9*np.cos(box.x))
    obs = system.observe(state)
    result = {}
    for a in (1.,1.2):
        own = system.measure(replace(state,a=a),obs)
        reference = [s.measure(SlabState(a,state.phi,state.pi,o),work=True)
                     for s,o in zip(system.sectors,state.orbitals,strict=True)]
        result[str(a)]={k:float(own[k]-system.signed @ np.array([r[k] for r in reference]))
                        for k in ('rho_F','P_F','P_F_x','P_F_perp')}
        result[str(a)]['force_error'] = float(np.max(abs(obs.force-
            sum(c*s.project_scalar(s.mass_derivative(state.phi/cfg.phase_scale,s.cov(o)))
                for c,s,o in zip(system.signed,system.sectors,state.orbitals,strict=True)))))
    return result


def run_case(args):
    name,work = args
    settings,gauss_order = RADIAL_CASES[name] if name in RADIAL_CASES else (CASES[name],None)
    n,qt,spacing,points,scalar_modes,steps,auxiliaries = settings
    cfg = replace(Configuration(),longitudinal=n,transverse_cut=qt,transverse_bins=round(qt/spacing),
                  points=points,scalar_modes=scalar_modes)
    summary,history = evaluate(cfg,steps,auxiliaries,gauss_order=gauss_order)
    (work/f'{name}.json').write_text(json.dumps(summary,allow_nan=False))
    np.savez_compressed(work/f'{name}.npz',**history)
    print(name,'work',summary['integrated_transfer_matched'],'budget',summary['max_pressure_work_error'],flush=True)
    return name


def compare(a,b,suffix='_matched'):
    r,d = a
    s,q = b
    work = 'comoving_work'+suffix
    transfer = 'integrated_transfer'+suffix
    absolute = 'integrated_absolute_transfer'+suffix
    instant = np.max(abs(np.interp(d['t'],q['t'],q[work])-d[work]))/r['peak_instantaneous_work'+suffix]
    cum = np.max(abs(np.interp(d['t'],q['t'],q[transfer])-d[transfer]))/np.max(abs(d[transfer]))
    final = abs(s[transfer]/r[transfer]-1)
    abs_work = abs(s[absolute]/r[absolute]-1)
    return dict(instantaneous_peak_scaled_error=float(instant),cumulative_error=float(cum),
                terminal_transfer_error=float(final),absolute_transfer_error=float(abs_work),
                gate_pass=bool(max(instant,cum,final,abs_work)<0.01))


def gather(work,radial=False):
    history = {}
    for name in (['transverse_bins',*RADIAL_CASES] if radial else CASES):
        with np.load(work/f'{name}.npz') as source:
            history[name] = (json.loads((work/f'{name}.json').read_text()),dict(source))
    if radial:
        pairs = dict(uniform_refinement=('uniform_fine','transverse_bins'),
                     cross_quadrature=('gauss_refined','uniform_fine'),
                     gauss_order=('gauss_refined','gauss_reference'),
                     gauss_order_wide=('gauss_cutoff_refined','gauss_cutoff'),
                     gauss_cutoff=('gauss_cutoff_refined','gauss_refined'),
                     gauss_heavier=('gauss_cutoff','gauss_heavier'))
    else:
        pairs = dict(cutoff_9_to_12=('cylinder_12','cylinder_9'),cutoff_12_to_16=('cylinder_16','cylinder_12'))
        pairs.update({name:('cylinder_12',name) for name in list(CASES)[3:]})
    result = dict(
        scope='synthetic expanding 3+1 slab with one resolved spatial axis and a homogeneous metric',
        cases={name:row[0] for name,row in history.items()},
        comparisons={name:compare(history[a],history[b]) for name,(a,b) in pairs.items()},
        unallocated_comparisons={name:compare(history[a],history[b],'') for name,(a,b) in pairs.items()},
        stress_crosscheck=check_expectation(),
        scientific_acceptance=False,
    )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--workers',type=int,default=5)
    parser.add_argument('--only',nargs='+',choices=[*CASES,*RADIAL_CASES])
    parser.add_argument('--gather',action='store_true')
    parser.add_argument('--pulse',action='store_true')
    parser.add_argument('--radial',action='store_true')
    args = parser.parse_args()
    args.work.mkdir(parents=True,exist_ok=True)
    if not args.gather and not args.pulse:
        names = args.only or (RADIAL_CASES if args.radial else CASES)
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(run_case,[(name,args.work) for name in names]))
        if args.only:
            return
    if args.pulse:
        result = dict(cases=[pulse(masses,dt) for masses,dt in ((MASSES,0.01),(MASSES,0.005),(HEAVIER,0.01))])
    else:
        result = gather(args.work,args.radial)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Completed',args.output,flush=True)


if __name__=='__main__':
    main()
