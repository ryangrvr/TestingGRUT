#!/usr/bin/env python3
"""Exact free-Fermi moment feasibility at the declared X3 kinetic interface.

Units are the source's H_i units. Two helicities; n is total particle number,
rho includes rest energy. This does not test an interacting vacuum stress tensor.
"""
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, linprog

HERE = Path(__file__).resolve().parent
TRACK = HERE.parent
TRACKS = TRACK.parent
RES = TRACK / 'results'


def dump(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False, default=float) + '\n')


def source():
    path = TRACKS / 'x2_strip/results/lower_H_diag.txt.gz'
    meta = json.loads((TRACKS / 'x2_strip/results/lower_H_lattice.json').read_text())
    assert hashlib.sha256(path.read_bytes()).hexdigest() == meta['diag_sha256']
    with gzip.open(path, 'rt') as f:
        cols = f.readline().lstrip('# ').split()
        data = np.loadtxt(f)
    return meta, cols, data


def rho_min(n, mass, cap=1., g=2.):
    """Bathtub principle: occupy the lowest energies first, f=cap below kF."""
    k = (6 * math.pi**2 * n / (g * cap))**(1/3)
    w = math.hypot(k, mass)
    integral = (k*w*(2*k*k + mass*mass) - mass**4*math.asinh(k/mass))/8
    return g*cap*integral/(2*math.pi**2)


def main():
    meta, cols, data = source()
    v = meta['run']['P']['v']
    rows = []
    for time in [0.3, 0.5, 1., 1.5, 2., 2.5, data[-1, 0]]:
        d = dict(zip(cols, data[np.argmin(abs(data[:,0]-time))]))
        species = []
        for j, mass in enumerate([250.,400.,200.,120.],1):
            energy, number = d[f'rN{j}']*v*v, d[f'nN{j}']*v*v
            entry = {'species':f'N{j}', 'mass':mass, 'energy':energy, 'number':number,
                     'mean_energy':energy/number, 'caps':{}}
            for cap in [1., .5]:
                kf = (3*math.pi**2*number/cap)**(1/3)
                emin = rho_min(number,mass,cap)
                independent = cap/math.pi**2*quad(lambda k:k*k*math.hypot(k,mass),0,kf,
                                                   epsabs=1e-5,epsrel=1e-12)[0]
                nmax = brentq(lambda n:rho_min(n,mass,cap)-energy,number*1e-6,number*100,xtol=1e-4)
                entry['caps'][str(cap)] = {'E_min':emin, 'E_over_E_min':energy/emin,
                    'feasible':bool(energy>=emin), 'n_max_at_fixed_energy':nmax,
                    'minimum_fraction_number_removed':max(0.,1-nmax/number),
                    'quadrature_relative_error':abs(independent/emin-1)}
            species.append(entry)
        rows.append({'t':float(d['t']),'species':species})
    # Independent finite linear program: no Fermi-sphere formula in its objective.
    edges = np.linspace(0,12,4001); k=(edges[:-1]+edges[1:])/2
    wt=(edges[1:]**3-edges[:-1]**3)/(3*math.pi**2)
    n=.7; mass=1.3
    lp=linprog(wt*np.hypot(k,mass),A_eq=wt[None,:],b_eq=[n],bounds=(0,1),method='highs')
    assert lp.success
    check={'linear_program_rho':lp.fun,'formula_rho':rho_min(n,mass),
           'relative_error':abs(lp.fun/rho_min(n,mass)-1),
           'massless_normalization':'rho_min=g fcap kF^4/(8 pi^2); n=g fcap kF^3/(6 pi^2)',
           'g':2,'lp_bins':4000}
    assert check['relative_error']<2e-6
    assert max(s['caps'][str(c)]['quadrature_relative_error'] for r in rows for s in r['species'] for c in [1.,.5])<1e-10
    out={'source_sha256':meta['diag_sha256'],'H_GeV':meta['run']['H_GeV'],
         'scope':'free two-helicity kinetic particle moments, E=sqrt(p^2+M^2)',
         'source_scalar_budget_absmax':float(abs(data[:,cols.index('budget')]).max()),
         'checks':check,'rows':rows,'verdict':'native-moment splice fails: N1 violates physical Pauli bound'}
    dump(RES/'moment_feasibility.json',out)
    print(json.dumps({'last':rows[-1],'checks':check},indent=2),flush=True)


if __name__=='__main__': main()
