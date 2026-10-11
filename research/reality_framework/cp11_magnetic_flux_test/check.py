"""Independent static CP-11 magnetic-incidence comparator, NOT frozen CP-11 law.
Adds supplied Z2 phases on incidence columns to test fermionic loop energy.
Conditional filled-fermion sector only; no joint graph+gauge quantum dynamics.
Requires numpy, scipy. Run: python check.py
"""
import itertools, json, math
from pathlib import Path
import numpy as np
from scipy.integrate import quad

KAPPA=0.5

def energy(n,edges,negative_mask=0):
    B=np.zeros((n,len(edges)),complex)
    for j,(u,v) in enumerate(edges):
        B[u,j]=1.
        B[v,j]=-(-1 if (negative_mask>>j)&1 else 1)
    eig=np.linalg.eigvalsh(B@B.conj().T)
    assert eig.min()>-1e-10
    return KAPPA*len(edges)+n-float(np.sqrt(1+np.maximum(eig,0)).sum())

def all_signs(n,edges):
    return [energy(n,edges,bits) for bits in range(1<<len(edges))]

def square_bulk(k=512,pi_flux=True):
    kx=np.arange(k)*2*np.pi/k-np.pi
    c=np.cos(kx)
    if pi_flux:
        a=2*np.sqrt(c[:,None]**2+c[None,:]**2)
        avg=np.mean((np.sqrt(5+a)+np.sqrt(5-a))/2)
    else:
        a=2*(c[:,None]+c[None,:])
        avg=np.sqrt(5-a).mean()
    return 2-float(avg)

def regular_tree_bulk(d):
    r=2*math.sqrt(d-1)
    def rho(a):
        return d*math.sqrt(max(0,4*(d-1)-a*a))/(2*math.pi*(d*d-a*a))
    avg=quad(lambda a:rho(a)*math.sqrt(1+d-a),-r,r,
             epsabs=2e-12,epsrel=2e-12)[0]
    return KAPPA*d/2+1-avg

def main():
    path=[(0,1),(1,2),(2,3)]
    square=[(0,1),(1,2),(2,3),(0,3)]
    complete=list(itertools.combinations(range(4),2))
    forest=all_signs(4,path)
    assert max(forest)-min(forest)<1e-12
    bare_ring=energy(4,square,0)
    best_ring=min(all_signs(4,square))
    assert abs(bare_ring-(-.7001695926375442))<1e-11
    assert abs(best_ring-(-.7205662327304481))<1e-11
    complete_zero=energy(4,complete,0)
    complete_best=min(all_signs(4,complete))
    assert complete_best<best_ring
    pi=[square_bulk(k,True) for k in (128,256,512)]
    zero=[square_bulk(k,False) for k in (128,256,512)]
    assert max(pi)-min(pi)<1e-10
    assert max(zero)-min(zero)<1e-10
    assert abs(zero[-1]-(-.184539293475))<3e-9
    assert abs(pi[-1]-(-.188160689705))<3e-9
    tree3,tree4=regular_tree_bulk(3),regular_tree_bulk(4)
    assert abs(tree3-(-.197366056730))<3e-9
    assert abs(tree4-(-.186535041049))<3e-9
    assert zero[-1]>tree4>pi[-1]>tree3
    pairs=list(itertools.combinations(range(4),2))
    competitors=[]
    for mask in range(1<<len(pairs)):
        edges=[e for i,e in enumerate(pairs) if mask>>i&1]
        vals=all_signs(4,edges)
        competitors.append(dict(edges=len(edges),minimum=min(vals),zero=vals[0],mask=mask))
    winner=min(competitors,key=lambda x:x["minimum"])
    assert winner["edges"]==6  # flux-optimized complete graph at N=4
    result=dict(
        controls=dict(kappa=.5,m=1,t=1,finite_graph_n=4,
            phase_choice="binary Z2 incidence phases, SUPPLIED",
            full_graph_gauge_quantum_dynamics="NOT constructed"),
        four_vertex=dict(path_all_phases_identical=forest[0],
            square_zero=bare_ring,square_best_phase=best_ring,
            complete_zero=complete_zero,complete_best_phase=complete_best,
            full_binary_phase_graph_winner=winner),
        thermodynamic_energy_per_vertex=dict(
            square_zero_flux=zero[-1],square_pi_flux=pi[-1],
            four_regular_tree=tree4,three_regular_tree=tree3,
            square_flux_improvement=zero[-1]-pi[-1],
            square_pi_advantage_over_same_degree_tree=tree4-pi[-1],
            square_pi_disadvantage_to_three_degree_tree=pi[-1]-tree3),
        independent_uniform_phase_thermal_density={
            str(beta):1/(1+math.exp(beta*.5)) for beta in (1.,2.,5.)},
        verdict="Flux reverses same-degree loop disadvantage; does not select all graph classes or cure finite-temperature entropy.")
    print(json.dumps(result,indent=2))
    Path(__file__).with_name("results.json").write_text(json.dumps(result,indent=2)+"\n")
    print("CP11_MAGNETIC_INCIDENCE_COMPARATOR_PASS")
if __name__=="__main__":
    main()
