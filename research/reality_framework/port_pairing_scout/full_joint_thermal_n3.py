#!/usr/bin/env python3
"""Complete finite-T quantum trace of full port-fermion joint N=3 model.

Independent of the Born-Oppenheimer Gibbs approximation; exactly all
300 joint eigenstates are included. This does not resolve large N.
Run python full_joint_thermal_n3.py with numpy/scipy.
"""
import json
from pathlib import Path
import numpy as np
from full_joint import assemble,cycle_lengths

def run(sign):
    n=3
    H,match,fb,index=assemble(n,sign)
    E,V=np.linalg.eigh(H.toarray())
    nfermion=len(fb)
    mark=np.array([cycle_lengths(p,n)==(n,) for p in match])
    p=np.sum(abs(V.reshape(len(match),nfermion,-1)[mark])**2,
             axis=(0,1))
    assert abs(p.mean()-sum(mark)/len(match))<1e-12
    rows=[]
    for beta in (0.,1.,5.,20.,100.):
        weights=np.exp(-beta*(E-E.min()))
        prob=float(np.dot(weights,p)/weights.sum())
        rows.append(dict(beta=beta,connected=prob))
    return dict(N=n,sign=sign,states=H.shape[0],probabilities=rows)

def main():
    rows=[run(-1.),run(1.)]
    assert abs(rows[0]['probabilities'][1]['connected']-
               .6921014938439972)<1e-10
    assert abs(rows[0]['probabilities'][2]['connected']-
               .997021708621257)<1e-10
    assert rows[1]['probabilities'][2]['connected']<.001
    Path(__file__).with_name('full_joint_thermal_n3_results.json').write_text(
       json.dumps(rows,indent=2)+'\n')
    print(json.dumps(rows,indent=2))
    print('SCOUT_RG02_FULL_JOINT_THERMAL_N3_PASS')
if __name__=='__main__':main()
