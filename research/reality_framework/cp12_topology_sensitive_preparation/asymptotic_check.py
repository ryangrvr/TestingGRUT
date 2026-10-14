"""Closed-form cubic differences and their analytic asymptotic limits.
Numbers at large N evaluate formulas, not Hamiltonian time evolutions.
"""
from fractions import Fraction as F
from math import comb,sqrt
from pathlib import Path
import json,os,sys
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
from topology import cubic_leader,falling,hypergeometric_average


def moment_check():
    # Independently verify the finite factorial identity for arbitrary pools.
    for M in range(8,31):
        for a in [3,M//2,M-3]:
            b=M-a
            for Q in [3,M//2,M-3]:
                def h(k):
                    return F(falling(k,3),falling(a,3))-3*F(falling(k,2)*(Q-k),falling(a,2)*b)+3*F(k*falling(Q-k,2),a*falling(b,2))-F(falling(Q-k,3),falling(b,3))
                got=hypergeometric_average(lambda k:k**3*h(k),6,M,Q,a)
                expected=F(6*falling(Q,3)*falling(M-Q,3),falling(M,6))
                assert got==expected


def main():
    moment_check();rows=[];gamma=F(1,10)
    tri_limit=36864*sqrt(2)*float(gamma**3);square_limit=-18*tri_limit
    for j in range(9):
        n=24*2**j;blocks=n//6;M=comb(n,2);Q=3*n
        coarse=int(F(M,2)-7*blocks+18)//36;a=7*blocks+36*coarse
        assert abs(F(a)-F(M,2))<=18
        c=1-F(M*(M-1),2*a*(M-a));ss=1-c*c;assert ss>0
        tri_over_sin=cubic_leader(M,c)*2*blocks*F(6*falling(Q,3)*falling(M-Q,3),falling(M,6))
        tri=float(tri_over_sin)*sqrt(float(ss));square=-float(F(3*(Q-3),M-6))*tri
        rows.append(dict(N=n,triangle_difference_T3=tri,square_difference_T3=square,
                         N5_scaled_triangle=tri*n**5,N6_scaled_square=square*n**6,
                         triangle_scaled_over_limit=tri*n**5/tri_limit,square_scaled_over_limit=square*n**6/square_limit))
    assert abs(rows[-1]['triangle_scaled_over_limit']-1)<.01
    assert abs(rows[-1]['square_scaled_over_limit']-1)<.01
    out=dict(status='CP12_TOPOLOGY_ASYMPTOTIC_FORMULA_PASS',triangle_N5_limit=tri_limit,square_N6_limit=square_limit,rows=rows,
             scope='Closed-form derivative asymptotics. No claim about response at times growing with N or about finite-time assembly.')
    (Path(os.environ.get('CP12_OUTPUT_DIR',str(P)))/'asymptotic_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
