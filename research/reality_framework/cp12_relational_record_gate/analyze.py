"""Relational diagnostics of the SAME states published at 8ece3e2.

No new size, Hamiltonian, preparation, source, measurement or evolution time.
Reconstruct projected states with sparse exponentials; compare every replayed
count against the saved trajectory and inherit its full-law leakage bound.
All geometric diagnostics are final-time, label-invariant diagonal operators.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations, permutations
from math import comb
import argparse, hashlib, json, os, sys, time
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import expm_multiply

P = Path(__file__).resolve().parent
UP = P.parent / 'cp12_autonomous_flux_gate'
sys.path.insert(0, str(UP))
from dynamics import resources, rewire_matrix
from run import model, graphstats

METRICS = ['connected', 'has_cycle', 'largest_component_fraction',
           'nonisolated_vertex_fraction', 'cycle_vertex_fraction',
           'triangle_density', 'square_density', 'short_cycle_packing_fraction',
           'cycle_rank_fraction', 'degree_square_fraction', 'reachable_pair_fraction']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@lru_cache(None)
def distances(n, mask):
    es = tuple(combinations(range(n), 2))
    adj = [[] for _ in range(n)]
    for e, (a, b) in enumerate(es):
        if mask >> e & 1:
            adj[a].append(b); adj[b].append(a)
    out = np.full((n, n), n, dtype=np.int8)
    for a in range(n):
        out[a, a] = 0; todo = [a]
        for b in todo:
            for c in adj[b]:
                if out[a, c] == n:
                    out[a, c] = out[a, b] + 1; todo.append(c)
    return out


def features(n, mask):
    es = tuple(combinations(range(n), 2)); M = len(es)
    D = distances(n, mask)
    components = []; available = set(range(n))
    while available:
        a = min(available); C = set(np.flatnonzero(D[a] < n).tolist())
        components.append(C); available -= C
    cycles = set()
    for e, (a, b) in enumerate(es):
        if mask >> e & 1 and distances(n, mask ^ (1 << e))[a, b] < n:
            cycles.update([a, b])
    gs = graphstats(n, mask, es)
    p = mask.bit_count(); rank = p - n + len(components)
    vals = [gs[0], float(rank > 0), max(map(len, components))/n,
            sum(np.count_nonzero(D[a] < n) > 1 for a in range(n))/n,
            len(cycles)/n, gs[1]/comb(n, 3), gs[2]/(3*comb(n, 4)),
            gs[3]/(n//3), rank/(M-n+1), gs[4]/(n*(n-1)**2),
            (np.count_nonzero(D < n)-n)/(n*(n-1))]
    vals += [gs[5+r]/n for r in range(n)]
    return np.asarray(vals, float)


@lru_cache(None)
def catalog(n):
    """Partition all labelled simple graphs into exact permutation orbits."""
    es = tuple(combinations(range(n), 2)); M = len(es)
    maps = [tuple(es.index(tuple(sorted((p[a], p[b])))) for a, b in es)
            for p in permutations(range(n))]
    canonical = np.full(1 << M, -1, dtype=np.int32)
    shapes = {}
    for mask in range(1 << M):
        if canonical[mask] >= 0:
            continue
        occ = [e for e in range(M) if mask >> e & 1]
        orbit = {sum(1 << m[e] for e in occ) for m in maps}
        key = min(orbit)
        canonical[list(orbit)] = key
        shapes[key] = features(n, key)
    keys = sorted(shapes); ids = {g:i for i, g in enumerate(keys)}
    class_ids = np.array([ids[int(g)] for g in canonical])
    scores = np.array([shapes[g] for g in keys])
    pop = np.array([g.bit_count() for g in range(1 << M)])
    hist = np.zeros((M+1, len(keys)))
    for k in range(M+1):
        hist[k] = np.bincount(class_ids[pop == k], minlength=len(keys))/comb(M, k)
    return dict(n=n, M=M, keys=keys, canonical=canonical,
                ids=class_ids, scores=scores, pop=pop,
                uniform_hist=hist, uniform_scores=hist@scores)


def core_mask(n, s):
    g, bits, fields = s
    mask = sum(1 << e for e in range(len(fields)//2)
               if (g >> e & 1) and
               (not (bits >> (n+e) & 1) or fields[2*e] or fields[2*e+1]))
    # Edge Gauss law forces a hole to have nonzero electric field.
    assert mask == sum(1 << e for e in range(len(fields)//2)
                       if fields[2*e] or fields[2*e+1])
    assert mask & ~g == 0
    return mask


def conditional_kernel(values, weights, groups):
    """Uniform physical-configuration null, then condition on each group.

    Quotient orbit sizes count labelled configurations. Only configurations
    represented in the existing reachable signed-invariant basis are used.
    This is a Gauss-compatible diagnostic distribution, not a claim that a
    diagonal physical mixture has the original full Hamiltonian energy.
    """
    _, inverse = np.unique(groups, axis=0, return_inverse=True)
    z = np.bincount(inverse, weights=weights)
    hist = np.zeros((len(z), values.max()+1))
    np.add.at(hist, (inverse, values), weights)
    hist /= z[:, None]
    return hist, inverse


def setup(n, q, states, sizes, r, nv, e2):
    C = catalog(n); M = C['M']; K = len(C['keys'])
    names = METRICS + [f'ball_{j}_fraction' for j in range(n)]
    cores = np.array([core_mask(n, s) for s in states], dtype=np.int32)
    counts = C['pop'][cores]
    assert np.array_equal(counts, q-r)
    classes = C['ids'][cores]
    actual = C['scores'][classes]
    uniform = C['uniform_scores'][counts]
    graphs = np.array([s[0] for s in states], dtype=np.int32)
    # Uniform protection of p among the edges of the SAME parent graph.
    thinning = {}
    for g, k in sorted(set(zip(graphs.tolist(), counts.tolist()))):
        hist = np.zeros(K)
        edges = [e for e in range(M) if g >> e & 1]
        for subset in combinations(edges, k):
            mask = sum(1 << e for e in subset)
            hist[C['ids'][mask]] += 1
        hist /= comb(q, k)
        thinning[g, k] = hist
    pairs = sorted(thinning); pair_id = {x:j for j,x in enumerate(pairs)}
    thin_inverse = np.array([pair_id[g, k] for g, k in zip(graphs, counts)])
    thin_hist = np.array([thinning[x] for x in pairs])
    resource_keys = np.column_stack([nv, e2, r]).astype(int)
    gauss_hist, gauss_inverse = conditional_kernel(classes, sizes, resource_keys)
    gp_hist, gp_inverse = conditional_kernel(classes, sizes,
        np.column_stack([resource_keys, graphs]))
    # Pad empty final classes for consistent matrix products.
    def pad(x):
        return np.pad(x, ((0,0),(0,K-x.shape[1])))
    gauss_hist = pad(gauss_hist); gp_hist = pad(gp_hist)
    # Candidate bins are evaluated relative to the protected core before
    # adding a neutral connection. Empty/full parent endpoints need not be
    # in the same core component. Disconnection is its own explicit bin.
    bins = ['distance_2', 'distance_3_or_more', 'disconnected']
    pool = {}; shortcut = np.zeros((len(states), 3)); completion = shortcut.copy()
    absent = shortcut.copy(); core_possible = shortcut.copy()
    for g, mask in sorted(set(zip(graphs.tolist(), cores.tolist()))):
        D = distances(n, mask); aa = np.zeros(3); nn = aa.copy(); allc = aa.copy()
        for e, (a,b) in enumerate(combinations(range(n),2)):
            if mask >> e & 1:
                continue
            d = D[a,b]
            category = 2 if d == n else (0 if d == 2 else 1)
            assert d >= 2
            allc[category] += 1
            if g >> e & 1:
                nn[category] += 1
            else:
                aa[category] += 1
        denom = M-mask.bit_count()
        expected = (q-mask.bit_count())*allc/denom if denom else allc*0
        pool[g, mask] = nn/q, expected/q, aa/(M-q), allc/denom if denom else allc*0
    for j, (g, mask) in enumerate(zip(graphs, cores)):
        shortcut[j], completion[j], absent[j], core_possible[j] = pool[g, mask]
    return dict(catalog=C, names=names, cores=cores, counts=counts, classes=classes,
        actual=actual, uniform=uniform, graphs=graphs,
        controls=dict(uniform_core_size=(C['uniform_hist'], counts),
                      parent_thinning=(thin_hist, thin_inverse),
                      gauss_resources=(gauss_hist, gauss_inverse),
                      gauss_parent_resources=(gp_hist, gp_inverse)),
        shortcut_bins=bins, shortcut=shortcut, completion=completion,
        absent=absent, core_possible=core_possible)


def analyze(n, q, budget, selected, snapshot_dir=None):
    started = time.monotonic()
    oldfile = UP / f'N{n}_Q{q}_B{budget}_T{20 if n==4 and q==4 else 1 if q==12 else 10}.json'
    old = json.loads(oldfile.read_text())
    sub,H,B,states,sizes,psi = model(n,q,budget)
    r,k,nv,e2 = resources(n,q,sub,states)
    S = setup(n,q,states,sizes,r,nv,e2)
    print(json.dumps(dict(stage='setup',N=n,Q=q,dimension=len(states),
        core_shapes=len(S['catalog']['keys']),seconds=time.monotonic()-started)),flush=True)
    # Replay already published times, not an additional larger/time-extended run.
    step = old['time_step']; end = max(selected); grid_steps = round(end/step)
    A = H-q/2*eye(H.shape[0],format='csr')
    waves = expm_multiply(-1j*A,psi,start=0.,stop=end,num=grid_steps+1,
                           traceA=-1j*float(A.diagonal().sum()))
    R = rewire_matrix(n,q,sub,states,sizes,H=H)
    rows=[]; maximum_replay_error=0.; maximum_norm_drift=0.; maximum_energy_drift=0.
    energy_diag = q/2+2*nv+e2/8+.1/sub.M*r*(sub.M-q)
    for j,v in enumerate(waves):
        t=j*step; ref=old['rows'][j]; assert abs(ref['T']-t)<1e-10
        p=abs(v)**2; e=ref['state_error_upper']; delta=2*e+e*e
        errors=[abs(float(p@u/q)-ref[key]) for u,key in
                [(r,'removable_fraction'),(k,'flux_protected_fraction'),(nv,'hole_protected_fraction')]]
        maximum_replay_error=max(maximum_replay_error,*errors)
        maximum_norm_drift=max(maximum_norm_drift,abs(float(p.sum())-1))
        maximum_energy_drift=max(maximum_energy_drift,abs(float(np.vdot(v,H@v).real)-q/2))
        assert max(errors)<1e-9
        if min(abs(t-s) for s in selected)>1e-8:
            continue
        hist=np.bincount(S['classes'],weights=p,minlength=len(S['catalog']['keys']))
        means=p@S['actual']; controls={}
        for label,(kernel,inverse) in S['controls'].items():
            mass=np.bincount(inverse,weights=p,minlength=len(kernel))
            hh=mass@kernel; baseline=hh@S['catalog']['scores']
            assert abs(hh.sum()-p.sum())<1e-10
            controls[label]=dict(means=dict(zip(S['names'],map(float,baseline))),
                contrasts=dict(zip(S['names'],map(float,means-baseline))),
                core_shape_TV=float(abs(hist-hh).sum()/2),core_shape_distribution=hh.tolist())
        # These diagnostic mixtures match diagonal resource components, but
        # discard hopping/rewiring coherence; explicitly report their energy.
        incoherent_energy=float(p@energy_diag)
        Hv=H@v;Rv=R@v
        shortcut_actual=p@S['shortcut']; shortcut_uniform=p@S['completion']
        currents=np.array([2*np.vdot(S['shortcut'][:,c]*v,Rv).imag for c in range(3)])
        countdist=np.bincount(S['counts'],weights=p,minlength=sub.M+1)
        row=dict(T=t,state_error_upper=e,bounded_operator_error_upper=delta,
            means=dict(zip(S['names'],map(float,means))),
            core_shape_distribution=hist.tolist(),protected_edge_count_distribution=countdist.tolist(),
            controls=controls,full_graph_means=dict(zip(S['names'],map(float,p@S['catalog']['scores'][S['catalog']['ids'][S['graphs']]]))),
            shortcut=dict(bins=S['shortcut_bins'],neutral_edges_per_Q=shortcut_actual.tolist(),
                uniform_completion_per_Q=shortcut_uniform.tolist(),
                contrasts_per_Q=(shortcut_actual-shortcut_uniform).tolist(),
                absent_candidate_pool_fractions=(p@S['absent']).tolist(),
                available_core_complement_fractions=(p@S['core_possible']).tolist(),
                rewiring_only_neutral_edge_count_current_per_Q=currents.tolist()),
            projected_mean_energy=float(np.vdot(v,Hv).real),
            matched_resource_diagonal_mixture_energy=incoherent_energy,
            missing_coherence_energy=incoherent_energy-q/2,
            flux_protected_fraction=float(p@k/q),removable_fraction=float(p@r/q))
        rows.append(row)
        if snapshot_dir:
            d=Path(snapshot_dir);d.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(d/f'N{n}_Q{q}_B{budget}_T{t:g}.npz',wave=v)
        print(json.dumps(dict(stage='analysis',N=n,Q=q,T=t,
            cycle_probability=row['means']['has_cycle'],connected=row['means']['connected'],
            gauss_parent_TV=controls['gauss_parent_resources']['core_shape_TV'],
            seconds=time.monotonic()-started)),flush=True)
    # Full-graph fixed-Q null has the same mean energy as the initial pure
    # preparation ONLY when its neutral amplitudes and coherences are kept.
    C=S['catalog'];full_null=C['uniform_hist'][q]@C['scores']
    out=dict(N=n,Q=q,budget=budget,dimension=len(states),source_commit='8ece3e20a331b8fe0c745d6ab1d778185e56dade',
        input_trajectory=oldfile.name,input_trajectory_sha256=digest(oldfile),
        maximum_replay_count_difference=maximum_replay_error,
        maximum_replay_norm_drift=maximum_norm_drift,maximum_replay_energy_drift=maximum_energy_drift,
        metrics=S['names'],core_shape_canonical_masks=C['keys'],
        full_graph_uniform_fixed_Q_means=dict(zip(S['names'],map(float,full_null))),
        rows=rows,seconds=time.monotonic()-started,
        error_scope='All scores and conditional-null score functions are diagonal contractions. Analytic full-law cutoff errors inherited at identical times from the verified source trajectory. Expectation and contrast errors <=2e+e^2; core-shape-TV error <=2e+e^2. Ordinary floating point; no roundoff interval certificate. Rewiring-current bounds are separate, see derivation.',
        control_scope='Uniform core-size and parent-thinning nulls are combinatorial and need not obey flux Gauss constraints. Resource nulls use represented reachable physical configurations, weighted by labelled orbit size, conditional on n_v,E^2,r (and parent unlabelled graph). Extend these nulls by the identity on physical invariant configuration orbits outside the represented retained support. Their displayed distributions match diagonal resources, not full Hamiltonian energy; the projected diagonal-mixture energy mismatch is reported. The initial pure neutral state is an exactly mean-energy-matched fixed-Q reference.')
    path=P/f'RELATIONAL_N{n}_Q{q}_B{budget}.json'
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(stage='saved',file=path.name,seconds=out['seconds'])),flush=True)
    return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True)
    ap.add_argument('--q',type=int,required=True);ap.add_argument('--budget',type=int,required=True)
    ap.add_argument('--times',type=float,nargs='+',default=[0,1,3,5,10]);ap.add_argument('--snapshot-dir')
    a=ap.parse_args();analyze(a.n,a.q,a.budget,a.times,a.snapshot_dir)
