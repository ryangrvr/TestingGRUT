#!/usr/bin/env python3
"""Conditional unflavoured baryon yield with the axial law's fixed CP invariants."""
import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np

from production import TRACKS, read_columns


def run(work, kinetic, law, choice):
    work.mkdir(parents=True,exist_ok=True)
    os.environ['EK_SCRATCH']=str(kinetic)
    sys.path.insert(0,str(TRACKS/'e2e_kinetic/code'))
    import ek_common
    sys.path.insert(0,str(TRACKS/'e2e_card/code'))
    from lepto_fix import rho_eq
    import sm_lepto

    record=ek_common.read_rec(kinetic/'native')
    ledger=read_columns(kinetic/'native_ledger.txt')
    kin=json.loads((kinetic/'summary.json').read_text())
    source=json.loads(law.read_text())
    H=source['H55_scale_GeV']
    M=np.array(source['sterile_masses_GeV'])/H
    entry=source[choice]
    if max(abs(np.array(kin['mtilde_eV'])/entry['mtilde_eV']-1))>1e-10:
        raise ValueError('Kinetic and CP rate normalizations differ')
    temp=record['T_th']
    history=dict(t=ledger['t'].tolist(),a=ledger['a'].tolist(),H=record['H'].tolist(),
                 T_thermal=temp.tolist())
    for j in range(1,5):
        eq=np.array([rho_eq(M[j-1],T) for T in temp])/ledger['rho']
        history[f'f_N{j}']=np.maximum(record[f'fN{j}']-eq,0.).tolist()
        parent_mass=600. if j==1 else 250.
        mean_energy=record[f'fN{j}']*ledger['rho']/np.maximum(record[f'nN{j}'],1e-100)
        history[f'A_N{j}']=(mean_energy/(.75*parent_mass)).tolist()
    hist=work/'history.json'
    hist.write_text(json.dumps(dict(history=history),allow_nan=False)+'\n')
    # Delta is an input-row sentinel for the inherited driver. All
    # epsilon values are taken from the axial law at the given H55 scale.
    see=work/'seesaw.json'
    see.write_text(json.dumps(dict(H_i_GeV=H,M_over_Hi=M.tolist(),parent=['sigma','chi','chi','chi'],
                                  scenarios=dict(axial=dict(rows=[dict(delta=1.,eps=entry['epsilon'],
                                     Gamma_rest_over_Hi=(np.array(entry['Gamma_rest_GeV'])/H).tolist())])))))
    values={}
    for wash in (.5,1.,2.):
        result=sm_lepto.run(str(hist),str(see),'axial',H,delta_ref=1.,washout_scale=wash)
        eta=result['etaB_per_unit_delta']
        values[str(wash)]=dict(eta_B=float(eta),ratio_to_observed=float(eta/source['benchmark_eta_B']),
                              thermal_only_eta_B=float(result['etaB_thermal_only_per_unit_delta']),
                              integrated_washout=float(result['W_total']))
    out=dict(scope='conditional unflavoured Boltzmann response using excess over equilibrium',
             law_case=choice,kinetic_boundary=kin['source_time'],CP_epsilon=entry['epsilon'],
             mtilde_eV=entry['mtilde_eV'],washout_cases=values,accepted_GRUT_law=False)
    (work/'summary.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(choice,values,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work',type=Path)
    parser.add_argument('kinetic',type=Path)
    parser.add_argument('law',type=Path)
    parser.add_argument('choice',choices=('literal_positive_parity','inherited_parity_native_widths',
                                         'inherited_parity_later_widths'))
    args=parser.parse_args()
    run(args.work.resolve(),args.kinetic.resolve(),args.law.resolve(),args.choice)


if __name__=='__main__':
    main()
