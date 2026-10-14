"""Unregulated exact amplitudes, beyond the upstream reachability test."""
from fractions import Fraction as F
from math import factorial,lcm,comb
from pathlib import Path
import json,sys,time
from loops import P,masks,incidence_masks,directed_support,cycle_table
from exact_destination import sea,powers,apply,gauss


def coefficient(v,scale,order,predicate):
    return F(sum(a*a for s,a in v.items() if predicate(s)),scale**(2*order)*factorial(order)**2)


def run():
    started=time.monotonic();n=3;g=7;s0=sea(n,g)
    sc,vs=powers(n,{s0:1},7)
    def prot(s):return masks(n,s)[0]==7
    def flux(s):return masks(n,s)[1]==7
    def oneflux(s):return masks(n,s)[1]!=0
    for j in range(3):assert not any(prot(s) for s in vs[j])
    for j in range(6):assert not any(flux(s) for s in vs[j])
    for j in range(3):assert not any(oneflux(s) for s in vs[j])
    cp=coefficient(vs[3],sc,3,prot);cf=coefficient(vs[6],sc,6,flux)
    ck=coefficient(vs[3],sc,3,oneflux)
    assert cp==2 and cf==F(2,225) and ck==F(2,3)
    flux_states=[s for s in vs[6] if flux(s)];assert len(flux_states)==2
    # Protected-geometry loss from a vacuum circulation requires three hops.
    initial={flux_states[0]:1};sp,loss=powers(n,initial,3)
    for j in range(3):assert not any(not prot(s) for s in loss[j])
    loss_coef=coefficient(loss[3],sp,3,lambda s:not prot(s))
    assert loss_coef>0
    # Occupied-flux identity can change at one hop, while all edges stay protected.
    flux_loss=coefficient(loss[1],sp,1,lambda s:not flux(s))
    assert flux_loss==6
    circulation_loss=coefficient(loss[1],sp,1,lambda s:not directed_support(3,*incidence_masks(3,s)))
    assert circulation_loss==3
    def hole_circulation(s):
        return bool(masks(3,s)[2]) and bool(directed_support(3,*incidence_masks(3,s)))
    for j in range(7):assert not any(hole_circulation(s) for s in vs[j])
    hc=coefficient(vs[7],sc,7,hole_circulation)
    assert hc==F(32,33075)
    # The leading all-hole state has ordinary quadratic geometry loss.
    hole_state=next(s for s in vs[3] if prot(s))
    sh,holepowers=powers(n,{hole_state:1},1)
    hole_loss=coefficient(holepowers[1],sh,1,lambda s:not prot(s))
    assert hole_loss==3
    result=dict(N=3,Q=3,protected_triangle_first_amplitude_order=3,
        protected_triangle_probability_leading_coefficient=str(cp),
        occupied_flux_triangle_first_amplitude_order=6,
        occupied_flux_triangle_probability_leading_coefficient=str(cf),
        any_occupied_flux_first_amplitude_order=3,
        occupied_flux_fraction_probability_T6_coefficient=str(ck/3),
        leading_flux_triangle_target_count=len(flux_states),
        vacuum_flux_triangle_to_nonflux_probability_T2_coefficient=str(flux_loss),
        unit_vacuum_flux_triangle_circulation_loss_probability_T2_coefficient=str(circulation_loss),
        vacuum_flux_triangle_geometric_loss_probability_T6_coefficient=str(loss_coef),
        hole_bearing_circulation_first_amplitude_order=7,
        hole_bearing_circulation_probability_T14_coefficient=str(hc),
        leading_hole_triangle_geometric_loss_probability_T2_coefficient=str(hole_loss),
        exact_support_sizes=[len(v) for v in vs],
        integer_H6_flux_amplitudes=[str(F(vs[6][s],sc**6)) for s in flux_states],
        scope='Exact rational amplitudes of the full unbounded frozen law, no electric cutoff. The flux target and hole target used for the local loss calculation are supplied basis controls, not autonomous preparations. Leading coefficients, not long-time predictions.',
        seconds=time.monotonic()-started)
    (P/'ONSET_EXACT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True);return result


if __name__=='__main__':run()
