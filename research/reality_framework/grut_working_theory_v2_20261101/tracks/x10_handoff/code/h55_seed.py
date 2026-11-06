#!/usr/bin/env python3
"""Reconstruct the archived H55 scalar seed and its exact planar projections."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

TRACKS=Path(__file__).resolve().parents[2]
X2=TRACKS/'x2_strip/results'
sys.path.insert(0,str(TRACKS/'lf_fermion/code'))
from lfcommon import HB, vacuum_field


def archive():
    p=X2/'lower_H_lattice.json'
    card=json.loads(p.read_text())
    file=X2/'lower_H_diag.txt.gz'
    sha=hashlib.sha256(file.read_bytes()).hexdigest()
    if sha!=card['diag_sha256']:
        raise ValueError('Archived H55 diagnostic digest differs')
    with gzip.open(file,'rt') as f:
        names=f.readline().lstrip('# ').split()
        first=dict(zip(names,map(float,f.readline().split()),strict=True))
    return card['run'],first,sha


def reconstruct(run):
    P=run['P']
    start=run['start']
    a=start['a']
    n,L=run['N'],run['L']
    cutoff=run['kcut_mchi']*P['m_chi']
    rng=np.random.default_rng(run['seed'])
    Vss,Vxx,_=HB.hess2(P,start['s'],start['x'])
    columns={}
    for key,velocity,m2 in (('sigma','dsigma',Vss),('chi','dchi',Vxx),
                            ('spectator_a','dspectator_a',0.),('spectator_b','dspectator_b',0.)):
        columns[key],columns[velocity]=vacuum_field(n,L,m2,a,cutoff,rng)
    return columns


def axial(data):
    return [np.mean(data,axis=tuple(j for j in range(3) if j!=k)) for k in range(3)]


def derivative_fractions(field,L):
    n=len(field)
    k=2*np.pi*np.fft.fftfreq(n,d=L/n)
    q2=k[:,None,None]**2+k[None,:,None]**2+k[None,None,:]**2
    full=np.sum(abs(np.fft.fftn(field)/n**3)**2*q2)
    full_field=np.mean(field**2)
    out=[]
    for strip in axial(field):
        spectrum=np.fft.fft(strip)/n
        out.append(dict(variance_fraction=float(np.mean(strip**2)/full_field),
                        gradient_fraction=float(np.sum(abs(spectrum)**2*k*k)/full),
                        first_harmonic_gradient_fraction=float(2*abs(spectrum[1]*k[1])**2/full),
                        first_harmonic_coefficient=[float(spectrum[1].real),float(spectrum[1].imag)]))
    return out


def summarize(run,first,sha,columns):
    P=run['P']
    a=run['start']['a']
    L=run['L']
    n=run['N']
    dx=L/n
    derived=dict(sig_std=np.std(columns['sigma'])/P['v'],chi_std=np.std(columns['chi'])/P['v'],
                 xdot_std=np.std(columns['dchi'])/P['v'],
                 sa_std=np.std(columns['spectator_a'])/P['s_star'])

    def grad(field):
        power=sum(np.mean((np.roll(field,-1,axis=j)-field)**2)/dx**2 for j in range(3))
        return power/(2*a*a*P['v']**2)

    # This is the nearest-neighbor gradient of the reconstructed grid.
    derived['Gs']=grad(columns['sigma'])
    derived['Gx']=grad(columns['chi'])
    diagnostics={k:dict(reconstructed=float(v),archive=float(first[k]),
                        relative_difference=float(v/first[k]-1)) for k,v in derived.items()}
    axes={name:derivative_fractions(value,L) for name,value in columns.items()}
    return dict(scope='reconstruction of X2 H55 Gaussian scalar seed at its declared cutoff',
                diagnostic_sha256=sha,time=run['t0'],scale=a,box_length=L,grid=n,seed=run['seed'],
                k_first=2*np.pi/L,diagnostics=diagnostics,axis_projections=axes,
                scientific_acceptance=False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    run,first,sha=archive()
    report=summarize(run,first,sha,reconstruct(run))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('H55 seed reproduced',report['diagnostics'],flush=True)


if __name__=='__main__':
    main()
