"""Exploratory post-result moments of the archived interacting preparations.

This does not change a registered verdict. It distinguishes the free reference
displacement in the Weyl formula from the interacting expectation of F.
"""
import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np

from continuum import interval_spectrum

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    card=json.loads((HERE/"protocol_m1_01.json").read_text())
    omega,h2,fn,_,meta=interval_spectrum(card,card["reference_modes"],card["quadrature_epsabs"])
    hs=fn/np.sqrt(2*np.pi*omega)
    rows=[]
    for cutoff in card["energy_cutoffs"]:
        data=np.load(HERE/f"raw_m1_01/boson_m05_E{cutoff}_q128.npz")
        states=data["basis"];ground=data["ground"]
        field=defaultdict(complex)
        for state,amplitude in zip(states,ground):
            field[tuple(state)]+=meta["mean"]*amplitude
            for j,n in enumerate(state):
                if hs[j]==0:continue
                raised=state.copy();raised[j]+=1
                field[tuple(raised)]+=hs[j]*np.sqrt(n+1)*amplitude
                if n:
                    lowered=state.copy();lowered[j]-=1
                    field[tuple(lowered)]+=hs[j]*np.sqrt(n)*amplitude
        norm=float(np.vdot(ground,ground).real)
        mean=float(np.vdot(ground,np.array([field[tuple(state)] for state in states])).real)
        second=float(sum(abs(z)**2 for z in field.values())+h2[states.shape[1]:].sum()*norm)
        rows.append({"cutoff":cutoff,"mean":mean,"variance":second-mean*mean,"state_norm":norm})
    result={"role":"exploratory diagnostic after M1-01, not a new registered acceptance test",
            "free_reference_displacement":meta["mean"],"interacting_moments":rows,
            "variance_mode_tail_bound":meta["tail_variance_bound"],
            "scope":"Moments in the archived truncated ground states; no continuum preparation certificate."}
    with args.output.open("x") as f:json.dump(result,f,indent=2);f.write("\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
