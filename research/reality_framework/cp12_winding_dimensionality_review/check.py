#!/usr/bin/env python3
"""Independent CP12 induced Wilson-channel and static torus comparator audit.

Does NOT calculate CP12's complete Feshbach eigenenergy, quantum geometry
probability, continuum spacetime or interacting many-body free energy.
Uses owner-verified K_9(a) numbers pinned in commit 8e866d0.

Dependencies: NumPy (spectral heat kernel). Stdlib for exact walk counting.
"""
import json, math
from collections import defaultdict
from pathlib import Path
import numpy as np

SOURCE = '8e866d045f5c7419afc203211fc97282b0e8b90b'
GAMMA = .1
EPSILON_E = .25
MASS = 1.

def return_walks_exact(a,b,steps):
    """Closed walks per start vertex, integer universal-cover enumeration."""
    v = {(0, 0): 1}
    for _ in range(steps):
        new = defaultdict(int)
        for (x,y), weight in v.items():
            for nx,ny in ((x+1,y),(x-1,y),(x,y+1),(x,y-1)):
                new[nx,ny] += weight
        v = new
    return sum(w for (x,y), w in v.items() if x % a == 0 and y % b == 0)

def heat_dims(a,b,s):
    """The SAME graph Laplacian defines heat propagation on both tori."""
    def dim_one(n):
        k = np.arange(n)
        ev = 2 - 2*np.cos(2*np.pi*k/n)
        terms = np.exp(-s*ev)
        p = terms.mean()
        slope = -(ev*terms).mean()/p
        return float(p),float(slope)
    p1,slope1 = dim_one(a)
    p2,slope2 = dim_one(b)
    return {'heat_return_probability':p1*p2,
            'spectral_dimension':-2*s*(slope1+slope2)}

def flux_channel(n,short_cycle_length,k9):
    """Isolated + and - one-flux channel, not complete Feshbach energy."""
    q=2*n; m=n*(n-1)//2
    a=GAMMA*(m-q)/m
    dflux=short_cycle_length*(EPSILON_E-a)
    ncycles=n//short_cycle_length
    total_coefficient=-2*ncycles*k9*k9/dflux
    gap_bound=2*MASS+EPSILON_E/2-a
    convergence_sufficient=gap_bound/(2*q)
    return dict(N=n,Q=q,M=m,a=a,K9=k9,
       short_noncontractible_cycles=ncycles,
       flux_state_diagonal_cost=dflux,
       one_flux_mixing_matrix_element=f'{k9:.14g} t_hop^18',
       isolated_plusminus_flux_energy_shift_coefficient_t36=total_coefficient,
       isolated_shift_per_vertex_coefficient_t36=total_coefficient/n,
       finite_graph_sufficient_Neumann_threshold=convergence_sufficient,
       log10_abs_per_vertex_shift_at_sufficient_threshold=(
           math.log10(abs(total_coefficient/n))+36*math.log10(convergence_sufficient)),
       caveat='Channel-only Schur complement second-order shift, not entire vacuum/free energy.' )

def main():
    x=[]
    for n,thin,balanced,K9 in (
       (324,(9,36),(18,18),.013089333246327634),
       (729,(9,81),(27,27),.013249963407633625),
       (2916,(9,324),(54,54),.01334710643865799)):
        assert thin[0]*thin[1] == balanced[0]*balanced[1] == n
        assert thin[0] == 9
        assert 2*n == 2*thin[0]*thin[1]
        item={'N':n,'thin':thin,'balanced':balanced,'fixed_degree':4,
              'square_count_each':n,'short_winding_length_thin':9,
              'short_winding_length_balanced':balanced[0],
              'fraction_short_winding_channels_thin':(n//9)/n,
              'flux':flux_channel(n,9,K9)}
        item['spectral_dimensions']={str(s):{'thin':heat_dims(*thin,s),
                                             'balanced':heat_dims(*balanced,s)}
                    for s in (5.,10.,20.,30.,50.)}
        x.append(item)
    assert abs(x[1]['spectral_dimensions']['20.0']['thin']['spectral_dimension']-1.009641)<1e-5
    assert abs(x[1]['spectral_dimensions']['20.0']['balanced']['spectral_dimension']-2.003657)<1e-5
    assert abs(x[2]['spectral_dimensions']['20.0']['balanced']['spectral_dimension']-2.012825)<1e-5
    assert abs(x[1]['flux']['isolated_plusminus_flux_energy_shift_coefficient_t36']+0.020990495375)<1e-10
    moments=[]
    for k in range(0,14):
        t=return_walks_exact(9,81,k)
        b=return_walks_exact(27,27,k)
        moments.append(dict(length=k,thin_returns_per_vertex=t,
                            balanced_returns_per_vertex=b,difference=t-b))
        if k<9:assert t==b
    assert moments[9]['difference']==2
    assert moments[10]['difference']==0
    assert moments[11]['difference']==242
    n=11664
    late={'N':n,'thin':(9,1296),'balanced':(108,108),'s_200':{
        'thin':heat_dims(9,1296,200), 'balanced':heat_dims(108,108,200)}}
    assert .999<late['s_200']['thin']['spectral_dimension']<1.002
    assert 2.<late['s_200']['balanced']['spectral_dimension']<2.005
    result={'base_commit':SOURCE, 'scaling_families':x,
      'exact_closed_walk_moments_729':moments,'late_dimension_crosscheck':late,
      'derivation':{'isolated_flux_level':'d=ell*(epsilon_E-a)',
       'minimal_vertex':'K_ell(a)*t_hop^(2ell)',
       'isolated_channel_ground_shift':'-2*(number_of_cycles)*K_ell(a)^2/d * t_hop^(4ell)',
       'walk_topology':'Tr(A^k)/N agrees through k=8, first differs by 2 at k=9',
       'limit_thin':'C9 x C_(N/9) is a 9-width quasi-one-dimensional strip',
       'limit_balanced':'C_sqrtN x C_sqrtN has diverging width in both directions'},
      'nonclaims':['Not complete CP12 Gibbs or fixed-graph ground energy.',
       'Not a proof of dimension selection by quantum dynamics.',
       'Not a justification for frozen t_hop=1 from small-t convergent expansion.',
       'Flux-channel second-order piece need not equal entire coefficient at order 36.',
       'Initial graph topology is supplied in both competitors.']}
    Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'N729':x[1], 'late_spectral_200':late,
             'walk_first_difference':moments[9]},indent=2))
    print('CP12_WINDING_DIMENSIONALITY_CHECK_PASS')

if __name__=='__main__':main()
