#!/usr/bin/env python3
"""Signed auxiliary response and independent static finite-volume references."""
import argparse
from dataclasses import replace
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.special import roots_legendre

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'x6_slab/code'))
from slab import AP, Configuration, Slab, adjoint

MASSES = (3.2,4.8)
HEAVIER = (4.8,7.2)
HARMONICS = (1,2,3)


def transverse_rule(cfg,gauss_order=None):
    if gauss_order is None:
        edges = np.linspace(0,cfg.transverse_cut,cfg.transverse_bins+1)
        return (edges[:-1]+edges[1:])/2,np.diff(edges**2)/(4*np.pi*cfg.length)
    edges = np.r_[0.,[x for x in (1,2,4,8,12,16) if x<cfg.transverse_cut],cfg.transverse_cut]
    x,w = roots_legendre(gauss_order)
    half = np.diff(edges)/2
    p = (edges[:-1,None]+edges[1:,None])/2+half[:,None]*x
    return p.ravel(),(half[:,None]*w*p/(2*np.pi*cfg.length)).ravel()


class RadialSlab(Slab):
    def __init__(self,cfg,gauss_order):
        super().__init__(cfg)
        if gauss_order is not None:
            self.perp,self.weight = transverse_rule(cfg,gauss_order)
            self.q2 = self.k[None,:]**2+self.perp[:,None]**2
            self.kt = self.perp[:,None,None,None]*AP


def coefficients(mass, auxiliaries=MASSES):
    low,high = auxiliaries
    if not 0<=mass<low<high:
        raise ValueError('Require a physical mass below two distinct auxiliary masses')
    squared = np.square([mass,low,high])
    w = squared[0]*np.array([1.,-(squared[2]-squared[0])/(squared[2]-squared[1]),
                             (squared[1]-squared[0])/(squared[2]-squared[1])])
    return np.sqrt(squared),np.r_[1.,w[1:]/squared[1:]]


def response_parts(cfg,masses,signed,a,order,gauss_order=None,frequency_squared=0.):
    q = 2*np.pi/cfg.length*np.arange(-cfg.longitudinal,cfg.longitudinal+1)
    p,weight = transverse_rule(cfg,gauss_order)
    shape = (p[:,None]**2+q**2)/a**2
    energy = np.sqrt(masses[:,None,None]**2+shape)
    mass_weight = signed*masses**2/(cfg.phase_scale**2*a**3)
    lo,hi = energy[:,:,:-order],energy[:,:,order:]
    temporal = np.sum(mass_weight[:,None,None]*weight[None,:,None]/(4*energy**3))
    pair = lo+hi
    if np.any(pair**2<=frequency_squared):
        raise ValueError('Response check is restricted to frequencies below the pair threshold')
    interior = np.sum(mass_weight[:,None,None]*weight[None,:,None]*pair/(2*lo*hi*(pair**2-frequency_squared)))
    edge = (0.5*(1/energy[:,:,:order]+1/energy[:,:,-order:])).sum(axis=2)
    edge = np.sum(mass_weight[:,None]*weight[None,:]*edge)
    return edge,interior,temporal


def cylinder(cfg,masses,signed,a,order,gauss_order=None):
    edge,interior,temporal = response_parts(cfg,masses,signed,a,order,gauss_order)
    edge *= (a*cfg.length/(2*np.pi*order))**2
    return dict(edge=float(edge),interior=float(interior),local=float(temporal),
                total=float(edge+interior),matched=float(edge+interior-temporal))


def continuum_periodic(cfg,masses,signed,a,order,modes=1200,nodes=48,frequency_squared=0.):
    x,w = roots_legendre(nodes)
    x = (1+x)/2
    q = 2*np.pi/(a*cfg.length)*np.arange(-modes,modes+1)
    k = 2*np.pi*order/(a*cfg.length)
    shift = x[:,None]*(2*q*k+k*k)-(x*(1-x))[:,None]*frequency_squared
    mass_weight = signed*masses**2
    active = mass_weight != 0
    base = masses[active,None,None]**2+q[None,None,:]**2
    # Form the difference before adding over the signed regulator species.
    before = np.sqrt(base)
    shifted = base+shift[None,:,:]
    if np.any(shifted<=0):
        raise ValueError('Response check is restricted to frequencies below the pair threshold')
    after = np.sqrt(shifted)
    difference = -shift[None,:,:]/(before*after*(before+after))
    integral = np.sum(mass_weight[active,None,None]*w[None,:,None]*difference)/2
    return float(integral/(8*np.pi*a*cfg.length*cfg.phase_scale**2))


def continuum_bulk(cfg, masses, signed, a, k):
    x,w = roots_legendre(96)
    mass_weight = signed*masses**2
    active = mass_weight != 0
    terms = np.log1p((1-x*x)[None,:]*(k/(2*a*masses[active,None]))**2)
    return float(-np.sum(mass_weight[active,None]*w*terms)/(32*np.pi**2*cfg.phase_scale**2))


def eigen_curvature(cfg, masses, signed, a, order):
    result = 0.
    eps = 0.004
    for mass,coefficient in zip(masses,signed,strict=True):
        box = Slab(replace(cfg,mass=mass))
        k = order*2*math.pi/cfg.length
        wave = np.cos(k*box.x)/cfg.phase_scale
        h0 = box.full_hamiltonian(a,0*wave)
        reference = np.linalg.eigvalsh(h0)[:,:box.d]

        def energy_delta(amplitude):
            h = box.full_hamiltonian(a,amplitude*wave)
            return 0.5*(box.weight @ (np.linalg.eigvalsh(h)[:,:box.d]-reference).sum(axis=1))/a**3

        result += coefficient*(energy_delta(eps)+energy_delta(-eps))/eps**2*2*a*a/k**2
    return float(result)


def spectral_curvature(cfg,masses,signed,a,order,frequency_squared):
    curvature=0.
    for mass,c in zip(masses,signed,strict=True):
        box=Slab(replace(cfg,mass=mass))
        wave=np.cos(2*np.pi*order*box.x/cfg.length)/cfg.phase_scale
        energy,vectors=np.linalg.eigh(box.full_hamiltonian(a,0*wave))
        occupied,empty=vectors[:,:,:box.d],vectors[:,:,box.d:]
        vertex=np.kron(box.mass_matrix(1j*mass*wave),np.eye(2))
        seagull=np.kron(box.mass_matrix(-mass*wave**2),np.eye(2))
        diagonal=np.trace(adjoint(occupied) @ seagull @ occupied,axis1=-2,axis2=-1).real
        transitions=abs(adjoint(empty) @ vertex @ occupied)**2
        gap=energy[:,box.d:,None]-energy[:,None,:box.d]
        virtual=np.sum(2*transitions*gap/(gap**2-frequency_squared),axis=(1,2))
        curvature+=c*(box.weight @ (diagonal-virtual))/a**3
    return float(curvature)


def frequency_check():
    """Separate the finite virtual response from the truncation of that response."""
    cfg = Configuration()
    light,signed_light = coefficients(cfg.mass)
    heavy,signed_heavy = coefficients(cfg.mass,HEAVIER)
    rows=[]
    for a in (1.,1.15):
        for r in HARMONICS:
            for frequency_squared in (-2.6**2,0.,2.6**2):
                s=(2*np.pi*r/(cfg.length*a))**2-frequency_squared
                # The physical sector is identical in both regulators.
                spectra = ((light[1:],signed_light[1:]),(heavy[1:],signed_heavy[1:]))
                reference = [
                    continuum_periodic(cfg,m,c,a,r,modes=2400,nodes=96,
                                       frequency_squared=frequency_squared)*s
                    for m,c in spectra
                ]
                infinite_difference = reference[1]-reference[0]
                cylinders = []
                for n,qt in ((12,16),(16,22),(48,64),(96,128),(192,256)):
                    config=replace(cfg,longitudinal=n,transverse_cut=qt,transverse_bins=8*qt)
                    parts=[response_parts(config,m,c,a,r,frequency_squared=frequency_squared)
                           for m,c in spectra]
                    edge,interior,local=np.subtract(parts[1],parts[0])
                    remainder=s*(interior-local)
                    cylinders.append(dict(longitudinal=n,transverse_cut=qt,edge=float(edge),
                                          interior_difference=float(remainder),
                                          total_difference=float(edge+remainder),
                                          interior_cutoff_error=float(remainder-infinite_difference)))
                rows.append(dict(a=a,harmonic=r,frequency_squared=frequency_squared,
                                 unbounded_auxiliary_difference=infinite_difference,
                                 cylinders=cylinders))
    return rows


def static_check():
    cfg = Configuration()
    masses,signed = coefficients(cfg.mass)
    m_heavy,c_heavy = coefficients(cfg.mass,HEAVIER)
    m_high,c_high = coefficients(cfg.mass,(9.6,14.4))
    geometries = ((6,8),(9,12),(12,16),(16,22),(24,32),(48,64),(96,128),(192,256))
    references,cylinders = [],[]
    for a in (1.,1.2):
        for order in HARMONICS:
            reference = continuum_periodic(cfg,masses,signed,a,order)
            refined = continuum_periodic(cfg,masses,signed,a,order,modes=2400,nodes=96)
            heavier = continuum_periodic(cfg,m_heavy,c_heavy,a,order)
            high = continuum_periodic(cfg,m_high,c_high,a,order)
            physical = continuum_periodic(cfg,masses[:1],signed[:1],a,order)
            bulk = continuum_bulk(cfg,masses,signed,a,float(order))
            large = continuum_periodic(replace(cfg,length=8*cfg.length),masses,signed,a,8*order,
                                       modes=2400,nodes=96)
            references.append(dict(a=a,harmonic=order,periodic=reference,refined=refined,
                                   heavier=heavier,highest_auxiliary_mass=high,physical=physical,
                                   bulk=bulk,large_periodic=large))
        for n,qt in geometries:
            for spacing in (0.5,0.125):
                config = replace(cfg,longitudinal=n,transverse_cut=qt,transverse_bins=round(qt/spacing))
                rows=[]
                for order in HARMONICS:
                    row = cylinder(config,masses,signed,a,order)
                    raw = cylinder(config,masses[:1],signed[:1],a,order)
                    reference = next(r['refined'] for r in references if r['a']==a and r['harmonic']==order)
                    rows.append(dict(harmonic=order,**row,unregulated_edge=raw['edge'],
                                     continuum_relative_error=abs(row['matched']/reference-1)))
                cylinders.append(dict(a=a,longitudinal=n,transverse_cut=qt,spacing=spacing,harmonics=rows))
    eigen=[]
    for n,qt in geometries[:2]:
        config = replace(cfg,longitudinal=n,transverse_cut=qt,transverse_bins=2*qt)
        for order in HARMONICS:
            measured = eigen_curvature(config,masses,signed,1.,order)
            formula = cylinder(config,masses,signed,1.,order)['total']
            eigen.append(dict(longitudinal=n,harmonic=order,formula=formula,eigen_curvature=measured,
                              relative_error=abs(measured/formula-1)))

    independent_cuts=[]
    for n,qt in ((12,16),(16,16),(12,22),(16,22),(48,22),(16,64)):
        config=replace(cfg,longitudinal=n,transverse_cut=qt,transverse_bins=8*qt)
        independent_cuts.append(dict(longitudinal=n,transverse_cut=qt,
                                     harmonics=[cylinder(config,masses,signed,1.,r) for r in HARMONICS]))

    zero_mass,zero_signed=coefficients(0.)
    spectral=[]
    for freq2 in (-2.6**2,0.6**2,2.6**2):
        m,c=(masses[1:],signed[1:]) if freq2>1 else (masses,signed)
        for order in HARMONICS:
            edge,interior,_=response_parts(cfg,m,c,1.,order,frequency_squared=freq2)
            formula=edge+(order**2-freq2)*interior
            measured=spectral_curvature(cfg,m,c,1.,order,freq2)
            spectral.append(dict(harmonic=order,frequency_squared=freq2,formula=float(formula),
                                 spectral_curvature=measured,absolute_error=abs(formula-measured)))
    edge_pairs={}
    for a in (1.,1.2):
        lo,hi=(next(row['harmonics'] for row in cylinders
                    if row['a']==a and row['longitudinal']==n and row['spacing']==0.125)
               for n in (96,192))
        edge_pairs[str(a)]=[dict(harmonic=x['harmonic'],signed_order=math.log2(abs(x['edge']/y['edge'])),
                                 unsigned_order=math.log2(abs(x['unregulated_edge']/y['unregulated_edge'])))
                            for x,y in zip(lo,hi,strict=True)]
    return dict(masses=masses.tolist(),signed_coefficients=signed.tolist(),
                cancellation_residuals=[float(np.dot(signed,masses**n)) for n in (2,4)],
                references=references,cylinders=cylinders,independent_cutoffs=independent_cuts,
                curvature_checks=eigen,edge_orders_96_to_192=edge_pairs,
                spectral_curvature_checks=spectral,
                massless=cylinder(cfg,zero_mass,zero_signed,1.,1),
                massless_continuum=continuum_periodic(cfg,zero_mass,zero_signed,1.,1),
                frequency_response=frequency_check())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    result=static_check()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Cutoff edge orders',result['edge_orders_96_to_192'])


if __name__=='__main__':
    main()
