#!/usr/bin/env python3
"""Independent CP12-influenced classical degree/fragmentation tests.

These are EXACT statements about declared classical statistical surrogates,
NOT full CP12 quantum Gibbs or controlled t_hop=1 perturbation theory.
Standard library only; run: python check.py
"""
import collections,itertools,json,math
from pathlib import Path

def finite_injection(n):
    edges=list(itertools.combinations(range(n),2))
    reg=[];defects=collections.Counter();hist=collections.Counter()
    for graph in itertools.combinations(edges,n):
        d=[0]*n
        for a,b in graph:d[a]+=1;d[b]+=1
        W=sum(x*(x-1)//2 for x in d)
        hist[W]+=1
        if any(x!=2 for x in d):continue
        reg.append(graph)
        gs=set(graph)
        for x,y in graph:
            for a,b in ((x,y),(y,x)):
                for c in range(n):
                    e=tuple(sorted((a,c)))
                    if c in (a,b) or e in gs:continue
                    new=frozenset((gs-{tuple(sorted((a,b)))})|{e})
                    defects[new]+=1
    assert sum(defects.values())==2*n*(n-3)*len(reg)
    assert max(defects.values())<=3
    for graph in defects:
        d=[0]*n
        for a,b in graph:d[a]+=1;d[b]+=1
        assert sorted(d)==[1]+[2]*(n-2)+[3]
    out=[]
    for a in (0.,1.,2.,4.):
        z=sum(count*math.exp(-a*(w-n)) for w,count in hist.items())
        p=len(reg)/z
        bound=1/(1+(2*n*(n-3)/3)*math.exp(-a))
        assert p<=bound+1e-13
        out.append(dict(beta_lambda=a,exact_regular_probability=p,
                        theorem_upper_bound=bound))
    return dict(N=n,all_graphs=sum(hist.values()),regular=len(reg),
         defect_graphs=len(defects),max_preimage=max(defects.values()),
         forward_operations=sum(defects.values()),probabilities=out)

def single_cycle_prob(n,theta):
    # exp[theta sum_(l>=3) z^l/l] = sum c_n z^n
    # n*c_n = theta*sum_(l=3)^n c_(n-l).
    c=[0.]*(n+1);prefix=[0.]*(n+1)
    c[0]=1.;prefix[0]=1.
    for j in range(1,n+1):
        c[j]=theta*prefix[j-3]/j if j>=3 else 0.
        prefix[j]=prefix[j-1]+c[j]
    return theta/(n*c[n])

def saddle_distribution(a):
    # Supplementary standard leading sparse-graph saddle:
    # p_d proportional to z^d exp[-a*d(d-1)/2]/d!, with mean d=2.
    lg=[-a*d*(d-1)/2-math.lgamma(d+1) for d in range(40)]
    def calc(logz):
        v=[x+d*logz for d,x in enumerate(lg)]
        c=max(v);w=[math.exp(x-c) for x in v];z=sum(w)
        return sum(d*u for d,u in enumerate(w))/z,[u/z for u in w]
    lo,hi=-40.,80.
    for _ in range(130):
        mid=(lo+hi)/2
        if calc(mid)[0]>2:hi=mid
        else:lo=mid
    mean,p=calc((lo+hi)/2)
    assert abs(mean-2)<1e-12
    return dict(beta_lambda=a,p0=p[0],p1=p[1],p2=p[2],
       p3=p[3],p4=p[4],non_degree_2_fraction=1-p[2])

def disconnected_same_degree(n):
    assert n>=6
    path={(i,i+1) for i in range(n-1)}
    triangle_path={(0,1),(1,2),(0,2)}|{(i,i+1) for i in range(3,n-1)}
    def degrees(es):
        d=[0]*n
        for i,j in es:d[i]+=1;d[j]+=1
        return sorted(d)
    assert len(path)==len(triangle_path)==n-1
    assert degrees(path)==degrees(triangle_path)==[1,1]+[2]*(n-2)
    return dict(N=n,Q=n-1,min_W=n-2,
                spanning_path_and_disconnected_triangle_path_same_W=True)

def main():
    injection=[finite_injection(n) for n in (4,5,6,7)]
    cycles=[]
    for bd in (0.,1.,2.,4.):
        theta=math.exp(-bd)/2
        values=[dict(N=n,p=single_cycle_prob(n,theta))
                for n in (6,8,12,20,40,80,160,320,1000)]
        constant=math.gamma(theta+1)*math.exp(1.5*theta)
        assert abs(values[-1]["p"]*1000**theta/constant-1)<.02
        if bd==0.:assert abs(values[0]["p"]-6/7)<1e-13
        cycles.append(dict(beta_component_penalty=bd,
          cycle_fugacity_theta=theta,asymptotic_amplitude=constant,
          values=values))
    out=dict(theorems={
      "local_degree_injection":"P_2reg <= 1/[1+(2/3)N(N-3)exp(-beta lambda)]",
      "cycle_exact":"P_one_cycle=theta/[N*[z^N]exp(theta sum_(l>=3)z^l/l)]",
      "cycle_asymptotic":"P_one_cycle ~ Gamma(1+theta) exp(3theta/2) N^(-theta)",
      "theta":"exp(-beta Delta_component)/2"},
      conditions=["Theorem 1: labelled simple graphs, Q=N, Boltzmann e^(-beta lambda W).",
                  "Theorem 2: exactly degree-2 labelled simple graphs, fixed additive energy Delta per cycle.",
                  "Neither is the exact full CP12 interacting quantum Gibbs state."],
      finite_injection_checks=injection,cycle_results=cycles,
      degree_entropy_saddle=[saddle_distribution(a) for a in (0.,1.,2.,4.,6.,10.)],
      CP12_formal_beta_lambda_at_beta_one_t_equals_one=
           2/(2+.25/2)**3,
      formal_coefficient_warning="Only t_hop->0; NOT controlled at frozen t_hop=1.",
      q_equals_N_minus_one_degeneracy=
           [disconnected_same_degree(n) for n in (6,7,8,12)])
    Path(__file__).with_name("results.json").write_text(
       json.dumps(out,indent=2)+"\n")
    print(json.dumps({"injection":[(x["N"],x["regular"],x["defect_graphs"])
       for x in injection],
       "cycle_at_N1000":[(x["beta_component_penalty"],x["values"][-1]["p"])
         for x in cycles]},indent=2))
    print("CP12_REGULARITY_FRAGMENTATION_INDEPENDENT_PASS")
if __name__=="__main__":main()
