"""Induced gauge-loop dynamics in CP-12: exact small checks, no law changes.

Run from any directory with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python run.py.
Frozen-law results use m=t=1, electric=.25, the complete N=3,Q=3 graph,
and only an explicitly converged electric regulator. Additional m or t values
are labelled perturbative parameter controls, never frozen-law predictions.
"""
from __future__ import annotations
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import json
import itertools
import math
import sys
import time

import numpy as np
import scipy.linalg as la
from scipy.sparse import diags

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / 'cp12_ground_connectivity'
sys.path.insert(0, str(SOURCE))
import exact_destination as ed
import frozen_model as fm


def nv(n, state):
    return (state[1] & ((1 << n)-1)).bit_count()


def hopping(n, state):
    g, bits, E = state
    for e, (i,j) in enumerate(ed.edges(n)):
        if not (g >> e & 1):
            continue
        for leg, vertex, orient in [(0,i,1),(1,j,-1)]:
            for dst, src, shift in [(vertex,n+e,1),(n+e,vertex,-1)]:
                z = ed.hop(bits,dst,src)
                if z is None:
                    continue
                b, sign = z
                f = list(E)
                f[2*e+leg] += shift
                yield (g,b,tuple(f)), orient*sign


def minimum_loop_paths(length, mass=F(1), electric=F(1,4), z=F(0)):
    """Exact signed minimal paths, retaining Q=1-P intermediate states.

    Coefficient of t**(2ell) in <W sea|H_eff(z)|sea>. Electric term is
    electric/2 * E^2 and z excludes the common kappa*Q constant.
    No electric cutoff or amplitude pruning is used.
    """
    n = length
    cycle = list(range(n))
    g = ed.mask(n,list(zip(cycle,cycle[1:]+cycle[:1])))
    start = ed.sea(n,g)
    target = fm.wilson_vector(n,cycle,start)
    target_E = target[2]
    vec = {start:F(1)}
    raw = {start:1}
    rows = []
    for h in range(1,2*length+1):
        out, ro = defaultdict(F), defaultdict(int)
        for s, amp in vec.items():
            for dest, sign in hopping(n,s):
                # Geodesic criterion: each altered half-leg moves directly
                # from zero to its final +/-1 value, once.
                if any(x != 0 and x != y for x,y in zip(dest[2],target_E)):
                    continue
                if sum(abs(x) for x in dest[2]) != h:
                    continue
                if h < 2*length and nv(n,dest)==0:
                    continue
                den = z - 2*mass*nv(n,dest) - electric*F(sum(x*x for x in dest[2]),2)
                out[dest] += amp*sign if h==2*length else amp*sign/den
        for s, amp in raw.items():
            for dest, sign in hopping(n,s):
                if any(x != 0 and x != y for x,y in zip(dest[2],target_E)):
                    continue
                if sum(abs(x) for x in dest[2]) != h:
                    continue
                if h < 2*length and nv(n,dest)==0:
                    continue
                ro[dest] += amp*sign
        vec = {s:a for s,a in out.items() if a}
        raw = {s:a for s,a in ro.items() if a}
        rows.append({'hops':h,'support':len(vec)})
    coefficient = vec.get(target,F(0))
    return {'ell':length,'mass':str(mass),'electric':str(electric),'z':str(z),
            'scope':('Full frozen-law N3Q3 coefficient when m=t=1, epsilon_E=.25; R is identically zero' if length==3 else
                     'Fixed cycle host with rewiring omitted; diagonal-denominator hopping-skeleton control'),
            'resolvent_minimum_order_coefficient':str(coefficient),
            'coefficient_float':float(coefficient),
            'signed_H_power_amplitude':raw.get(target,0),'supports':rows}


def reduced(n=3, q=3, L=3, mass=1., t=1., electric=.25):
    H, states, lookup, _ = fm.build(n,q,L)
    number = np.array([nv(n,s) for s in states],float)
    e2 = np.array([sum(x*x for x in s[2]) for s in states],float)
    diag0 = .5*q+2*number+.125*e2
    K = H-diags(diag0)
    A = t*K+diags(2*mass*number+electric/2*e2)
    p = np.flatnonzero(number==0)
    qq = np.flatnonzero(number>0)
    return A.toarray(), states, lookup, p, qq


def feshbach(A,p,qq,z):
    return A[np.ix_(p,p)] + A[np.ix_(p,qq)] @ la.solve(
        z*np.eye(len(qq))-A[np.ix_(qq,qq)], A[np.ix_(qq,p)],assume_a='sym')


def boundary(states,lookup):
    out=defaultdict(dict)
    for i,s in enumerate(states):
        for dest,sgn in hopping(3,s):
            if dest not in lookup:
                out[dest][i]=out[dest].get(i,0)+sgn
    rowmax=max(sum(abs(c) for c in v.values()) for v in out.values())
    cols=defaultdict(int)
    for v in out.values():
        for i,c in v.items():
            cols[i]+=abs(c)
    colmax=max(cols.values())
    return out, rowmax*colmax


def boundary_action(out,x):
    return math.sqrt(sum(sum(c*x[i] for i,c in row.items())**2 for row in out.values()))


def unbounded_controls(L=16):
    A,states,lookup,p,q=reduced(L=L)
    vals,V=la.eigh(A)
    qvals=la.eigvalsh(A[np.ix_(q,q)])
    bout,b2=boundary(states,lookup)
    # Fourier gauge-angle fibers give ||K|| <= sum singular values(B)
    # <= sqrt(N Tr BB*) = sqrt(2*N*Q)=sqrt(18).
    knorm=math.sqrt(18)
    tail=.125*(L+1)**2-knorm
    qtail=2+tail
    # Exact 2x2 block lower comparison, using ||B||^2<=||B||_1||B||_inf.
    delta=qtail-qvals[0]
    qlower=float(qvals[0]-2*b2/(math.sqrt(delta*delta+4*b2)+delta))
    assert qlower>0
    # Count/isolate the full ground level, then Temple bound its energy.
    separator=-1.5
    beta=b2/(tail-separator)
    assert vals[1]-beta>separator
    outside_ground_res=boundary_action(bout,V[:,0])
    finite_ground_res=float(la.norm(A@V[:,0]-vals[0]*V[:,0]))
    ground_res=math.hypot(outside_ground_res,finite_ground_res)
    temple_width=ground_res**2/(separator-vals[0])
    s=ed.sea(3,7);w=fm.wilson_vector(3,[0,1,2],s)
    qs=A[q,lookup[s]];qw=A[q,lookup[w]]
    rr=[]
    for z in [0.,float(vals[0])]:
        Bz=z*np.eye(len(q))-A[np.ix_(q,q)]
        x=la.solve(Bz,qs,assume_a='sym')
        full=np.zeros(len(states));full[q]=x
        leak=boundary_action(bout,full)
        finite=float(la.norm(Bz@x-qs))
        bound=float(la.norm(qw))*(leak+finite)/(qlower-z)
        rr.append({'z':z,'offdiagonal':float(qw@x),
                   'outside_resolvent_residual':leak,
                   'finite_solve_residual':finite,
                   'absolute_full_unbounded_offdiagonal_error_bound':bound})
    return {'L':L,'K_global_norm_upper':knorm,'boundary_norm_squared_upper_integer':b2,
            'omitted_space_H_lower':tail,'omitted_matter_space_H_lower':qtail,
            'full_unbounded_matter_block_lower':qlower,
            'second_full_eigenvalue_lower':separator,
            'ground_outside_residual':outside_ground_res,
            'ground_finite_residual':finite_ground_res,
            'ground_combined_residual':ground_res,'ground_Temple_interval_width_from_computed_residual':temple_width,
            'ground_energy_upper_numerical':float(vals[0]),
            'full_unbounded_feshbach_entries':rr,
            'arithmetic_status':'Analytic error formulas evaluated in ordinary float64, not interval arithmetic; roundoff is not formally certified.'}


def frozen_cutoffs():
    rows=[]
    for L in [1,2,3,4,6,8,12,16]:
        A, states, lookup, p, qq = reduced(L=L)
        vals,V = la.eigh(A)
        z=float(vals[0])
        qvals=la.eigvalsh(A[np.ix_(qq,qq)])
        Heff=feshbach(A,p,qq,z)
        projected=V[p,0]
        residual=float(la.norm(Heff@projected-z*projected))
        s=ed.sea(3,7)
        w=fm.wilson_vector(3,[0,1,2],s)
        pi={int(a):i for i,a in enumerate(p)}
        a,b=pi[lookup[s]],pi[lookup[w]]
        effective_offdiag=float(Heff[b,a])
        Wilson=float(sum(V[i,0]*V[lookup[w],0] for i,s in enumerate(states)
                         if (w:=fm.wilson_vector(3,[0,1,2],s)) in lookup))
        zerofesh=feshbach(A,p,qq,0.)
        first=zerofesh[b,a]
        # Unregulated ground residual: only hopping that leaves the cube.
        full_out=defaultdict(float)
        for state, amp in zip(states,V[:,0]):
            for dest,sgn in hopping(3,state):
                if dest not in lookup:
                    full_out[dest] += amp*sgn
        leak=float(math.sqrt(sum(x*x for x in full_out.values())))
        row={'L':L,'dimension':len(states),'vacuum_dimension':len(p),
             'ground_energy_excluding_kappaQ':z,'first_five_energies':vals[:5].tolist(),
             'eliminated_block_minimum':float(qvals[0]),
             'distance_ground_to_eliminated_spectrum':float(np.min(abs(qvals-z))),
             'distance_zero_to_eliminated_spectrum':float(np.min(abs(qvals))),
             'ground_vacuum_weight':float(projected@projected),
             'ground_mean_nv':float(sum(nv(3,s)*v*v for s,v in zip(states,V[:,0]))),
             'ground_real_Wilson_expectation':Wilson,
             'ground_vacuum_unitflux_to_neutral_amplitude_ratio':float(V[p[b],0]/V[p[a],0]),
             'exact_feshbach_ground_residual':residual,
             'feshbach_neutral_to_unitflux_at_ground_energy':effective_offdiag,
             'feshbach_neutral_to_unitflux_at_zero':float(first),
             'unregulated_ground_residual_norm':leak}
        rows.append(row)
    return rows


def perturbative_controls():
    # t controls are labelled law-family checks of the exact t^6 coefficient.
    coeff=F(minimum_loop_paths(3)['resolvent_minimum_order_coefficient'])
    rows=[]
    for t in [.5,.25,.125,.0625]:
        A,states,lookup,p,qq=reduced(L=3,t=t)
        H=feshbach(A,p,qq,0.)
        pi={int(a):i for i,a in enumerate(p)}
        s=ed.sea(3,7);w=fm.wilson_vector(3,[0,1,2],s)
        actual=float(H[pi[lookup[w]],pi[lookup[s]]])
        pred=float(coeff)*t**6
        rows.append({'t':t,'exact_feshbach_offdiag':actual,
                     'sixth_order':pred,'relative_error':actual/pred-1})
    return rows


def magnetic_laplacian(cycle,flux):
    n=len(cycle)
    B=np.zeros((n,n),complex)
    for e in range(n):
        i,j=cycle[e],cycle[(e+1)%n]
        phase=np.exp(1j*flux) if e==n-1 else 1.
        B[i,e]=1.
        B[j,e]=-phase
    return B@B.conj().T


def static_angle_controls():
    rows=[]
    for ell in [3,4,5,6]:
        catalan=math.comb(2*(ell-1),ell-1)//ell
        J=F(ell*catalan,2**(2*ell-1))
        path=minimum_loop_paths(ell,mass=F(1),electric=F(0))
        assert F(path['resolvent_minimum_order_coefficient'])==J
        # Direct single-particle spectrum: isolate flux dependence with
        # phase difference pi vs zero; leading difference is -4J*t^(2ell).
        # Use m=1 and t=.1-.25 only as controlled parameter comparison.
        rr=[]
        for t in [.25,.125]:
            f0=ell-np.sqrt(1+t*t*la.eigvalsh(magnetic_laplacian(list(range(ell)),0.))).sum()
            fp=ell-np.sqrt(1+t*t*la.eigvalsh(magnetic_laplacian(list(range(ell)),np.pi))).sum()
            delta=float(fp-f0)
            pred=-4*float(J)*t**(2*ell)
            rr.append({'t':t,'E_pi_minus_E_zero':delta,'leading_prediction':pred,
                       'ratio':delta/pred})
        rows.append({'ell':ell,'J_ell_at_m1':str(J),'geodesic_check':path,'direct_band_checks':rr})
    return rows


def coordination_controls():
    """Every N6Q12 graph; a static zero-angle fermion-band calculation only."""
    n,q=6,12
    edges=list(itertools.combinations(range(n),2))
    records=[]
    for present in itertools.combinations(range(len(edges)),q):
        B=np.zeros((n,q))
        degree=np.zeros(n,dtype=int)
        for a,e in enumerate(present):
            i,j=edges[e];B[i,a]=1;B[j,a]=-1;degree[i]+=1;degree[j]+=1
        L=B@B.T
        lam=la.eigvalsh(L)
        trace2=int(sum(degree*degree)+2*q)
        assert abs(float(lam@lam)-trace2)<1e-10
        records.append((tuple(degree.tolist()),int(sum(degree*degree)),lam,trace2))
    rows=[]
    for mass in [4.,8.]:
        exact=np.array([n*mass-np.sqrt(mass*mass+r[2]).sum() for r in records])
        sumd2=np.array([r[1] for r in records],int)
        fourth=np.array([-q/mass+r[3]/(8*mass**3) for r in records])
        sixth=np.array([fourth[i]-float(sum(r[2]**3))/(16*mass**5) for i,r in enumerate(records)])
        assert np.all(exact<=fourth+1e-12)
        assert np.all(exact>=sixth-1e-12)
        groups=[]
        for d2 in sorted(set(sumd2)):
            e=exact[sumd2==d2]
            groups.append({'sum_degree_squared':int(d2),'number_graphs':len(e),
                           'exact_energy_min':float(e.min()),'exact_energy_max':float(e.max()),
                           'fourth_order_energy':float(-q/mass+(d2+2*q)/(8*mass**3))})
        minimizers=np.flatnonzero(exact<=exact.min()+1e-12)
        regular_formula=5*mass-3*math.sqrt(mass*mass+4)-2*math.sqrt(mass*mass+6)
        assert abs(float(exact.min())-regular_formula)<1e-12
        assert all(records[i][0]==(4,4,4,4,4,4) for i in minimizers)
        rows.append({'mass':mass,'t':1.,'graphs':len(records),
                     'static_angle_choice':'zero incidence angles; not a physical rotor-angle eigenstate preparation',
                     'max_expansion_ratio_t2lambda_over_m2':float(max(max(r[2]) for r in records)/mass**2),
                     'minimizer_degree_sequences':sorted({records[i][0] for i in minimizers}),
                     'number_minimizers':len(minimizers),
                     'independent_regular_graph_closed_form_energy':regular_formula,
                     'max_absolute_E4_error':float(max(abs(exact-fourth))),
                     'max_absolute_E6_error':float(max(abs(exact-sixth))),
                     'all_exact_energies_between_E6_and_E4':True,'groups':groups})
    return {'N':n,'Q':q,'number_labelled_graphs':len(records),'controls':rows,
            'scope':'Exact static free-fermion band energy at supplied gauge angles. No graph dynamics or finite-electric adiabatic approximation is executed.'}


def main():
    started=time.time()
    exact=[minimum_loop_paths(ell) for ell in [3,4,5,6]]
    assert F(exact[0]['resolvent_minimum_order_coefficient'])==F(65536,653429)
    assert exact[0]['signed_H_power_amplitude']==-48
    # Independent, unpruned unregulated H powers from the pre-existing
    # exact_destination implementation; no geodesic restriction is supplied.
    sea=ed.sea(3,7)
    w=fm.wilson_vector(3,[0,1,2],sea)
    scale,powers=ed.powers(3,{sea:1},6)
    assert powers[-1][w]==-48*scale**6
    independent={'unpruned_H6_target_amplitude':str(F(powers[-1][w],scale**6)),
                 'full_H_power_supports':[len(v) for v in powers]}
    zeros=static_angle_controls()
    cutoffs=frozen_cutoffs()
    perturb=perturbative_controls()
    unbounded=unbounded_controls()
    coordination=coordination_controls()
    assert all(r['coefficient_float']>0 for r in exact)
    assert cutoffs[-1]['unregulated_ground_residual_norm']<1e-10
    assert abs(cutoffs[-1]['ground_energy_excluding_kappaQ']-cutoffs[-2]['ground_energy_excluding_kappaQ'])<1e-10
    assert max(r['exact_feshbach_ground_residual'] for r in cutoffs)<1e-11
    assert abs(perturb[-1]['relative_error'])<abs(perturb[0]['relative_error'])
    out={'status':'PASS','frozen_parameters':{'m':1,'t':1,'epsilon_E':.25,'kappa':.5,'Gamma':.1},
         'exact_geodesic_resolvent_coefficients':exact,
         'electric_zero_static_angle_controls':zeros,
         'frozen_triangle_cutoff_convergence':cutoffs,
         'perturbative_t_controls_not_frozen_predictions':perturb,
         'unbounded_operator_controls':unbounded,
         'independent_unregulated_matrix_power_check':independent,
         'static_coordination_controls':coordination,
         'seconds':time.time()-started}
    (HERE/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'exact':exact,'cutoff_last':cutoffs[-1],
                      't_controls':perturb,'unbounded_controls':unbounded,'seconds':out['seconds']},indent=2))


if __name__=='__main__':
    main()
