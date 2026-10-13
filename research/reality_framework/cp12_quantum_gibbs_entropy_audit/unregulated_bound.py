#!/usr/bin/env python3
"""CP12 unregulated rotor thermodynamic projector bound formula controls.

Checks the quantum relative entropy / electric Gaussian tail / low-rank
bound's finite-N algebra, not a substitute for auditing the proof.
"""
import json, math
from pathlib import Path
EPS=0.25; HOP=1.; GAMMA=.1; BETA=1.

def logchoose(n,k):
    return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)

def theta(s):
    assert s>0
    if s<.1:
        return math.sqrt(math.pi/s)*(1+2*sum(
             math.exp(-math.pi**2*j*j/s) for j in range(1,9)))
    return 1+2*sum(math.exp(-s*j*j) for j in range(1, max(20,int(15/math.sqrt(s)))))

def bound(N,nu,scale=8):
    assert nu>1 and N>=10 and N*(N-1)//2>nu*N
    Q=nu*N; M=N*(N-1)//2; K=scale*math.log(N); beta=BETA
    # Theorem: D(rho_H||rho_A)<=2 beta ||B|| <=4 beta(t+Gamma)Q.
    c0=logchoose(N+Q,Q)+2*Q*math.log(theta(beta*EPS/4))
    info=beta*K*N/2-c0
    q=1 if info<=0 else min(1,(4*beta*(abs(HOP)+GAMMA)*Q+math.log(2))/info)
    # Low electric configurations bounded by exp(2N/eps)*Theta(1/K)^(2Q).
    log_low_fiber=logchoose(N+Q,Q)+2*N/EPS+2*Q*math.log(theta(1/K))
    # |S|<=N! for one specified spatial graph class.
    log_low=math.lgamma(N+1)-logchoose(M,Q)+log_low_fiber+beta*(
           2*abs(HOP)*Q+GAMMA*Q*(M-Q)/M)
    return dict(N=N,nu=nu,K=K,high_energy_prob_upper=q,
         log_low_energy_prob_upper=log_low,
         low_rate=log_low/(N*math.log(N)))

def main():
    rows=[bound(N,nu) for nu in (2,3)
          for N in (100,1000,10000,100000,1000000,10**9)]
    for nu in (2,3):
        seq=[x for x in rows if x['nu']==nu]
        assert seq[-1]['high_energy_prob_upper']<1
        assert seq[-1]['low_rate']<seq[-2]['low_rate']
    out=dict(status="CP12_UNREGULATED_RELATIVE_ENTROPY_BOUND_ALGEBRA_CHECK_PASS",
        limitation="Exact formulas checked; proof and broad spatial-phase counting require partner review",rows=rows)
    Path(__file__).with_name("unregulated_results.json").write_text(json.dumps(out,indent=2)+"\n")
    print(out['status'])
    for row in rows:
        print(row['nu'],row['N'],round(row['high_energy_prob_upper'],5),round(row['low_rate'],5))
if __name__=="__main__":main()
