"""Re-solve the stationary M1 branch at reported frozen-profile snapshots."""
from pathlib import Path
import argparse
import importlib.util
import json
import os
import sys

ROOT=Path(__file__).resolve().parents[1]
RF=ROOT.parent
os.environ['GRUT_M1_KINETICS_WORK']=str(ROOT/'.work'/'inherited')
sys.path.insert(0,str(RF/'grut_m1_kinetics_20261028'/'code'))
import numpy as np
from scipy.integrate import simpson
from setup import A,B,C,EREF,G2,OMEGA,QREF
from cold_completion import Seed

spec=importlib.util.spec_from_file_location('m1_deformation',RF/'grut_matched_candidate_20261027'/'code'/'deformation.py')
D=importlib.util.module_from_spec(spec)
spec.loader.exec_module(D)


def stationary(q):
    sol=D.solve(q,.002,Seed(),length=160.,tol=2e-9)
    r=np.linspace(0,160,24001)
    f,fp,s,sp,j=sol.sol(r)
    omega=float(sol.p[0])
    potential=B/4*(s*s-1)**2-C*(s-1)**3+.5*s*s*f*f+A/4*f**4
    energy=float(4*np.pi/G2*simpson(r*r*(.5*(fp*fp+sp*sp)+.5*omega*omega*f*f+potential),x=r))
    return {'Q':q,'Omega':omega,'energy_mX':energy,
            'BVP_residual':float(sol.rms_residuals.max()),
            'tangent_minus_stationary_mX':float(EREF+OMEGA*(q-QREF)-energy)}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('tags',nargs='+')
    args=p.parse_args()
    result=[]
    for tag in args.tags:
        d=json.loads((ROOT/(tag+'.json')).read_text())
        a=d['snapshots']
        selected=[a[0],next(s for s in a if abs(s['tau']-1000)<1e-4),
                  next(s for s in a if abs(s['a_over_a0']-2)<1e-6),a[-1]]
        rows=[]
        for r in selected:
            b=stationary(r['Qcore'])
            b.update({'a_over_a0':r['a_over_a0'],'tau':r['tau'],
                      'Nplus':r['Nplus'],'Nminus':r['Nminus'],
                      'kinetic_mX':r['gas_kinetic_mX'],
                      'stationary_cold_excess_fraction':(b['energy_mX']+r['Nplus']+r['Nminus'])/EREF-1})
            rows.append(b)
        source=d['configuration'].get('source','parent')
        if source=='parent':
            pulse=json.loads((RF/'grut_m1_kinetics_20261028'/'PHYSICS_CHECKS.json').read_text())['source_energy_ledger']['pulse_initial_lattice_energy_mX']
        else:
            pulse=json.loads((ROOT/'sources'/f'ringdown_{source.removesuffix("early")}.json').read_text())['initial_pulse_energy_mX']
        remainder=pulse-rows[0]['energy_mX']-a[0]['gas_energy_mX']
        result.append({'tag':tag,'source':source,'pulse_initial_lattice_energy_mX':pulse,
                       'unassigned_initial_energy_mX':remainder,
                       'remainder_definition':'pulse lattice energy minus continuum stationary core energy at inferred charge minus kinetic gas source; includes excitation, finite-window and matching effects',
                       'snapshots':rows})
    output={'reference_mass_mX':EREF,'results':result,
            'scope':'Cold diagnostic only: stationary BVP evaluated at the charge inferred from the frozen-profile trajectory; no backreacted scattering or excitation release.'}
    (ROOT/'MASS_ANALYSIS.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps(output,indent=2))

if __name__=='__main__':
    main()
