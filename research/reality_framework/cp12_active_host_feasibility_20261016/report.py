"""Generate a claim-scoped report directly from completed evidence."""
from pathlib import Path
import json
P=Path(__file__).resolve().parent


def read(name):return json.loads((P/name).read_text())
def interval(I):return f'[{I[0]:.8g}, {I[1]:.8g}]'


def main():
    dynamics={n:read(f'HOST_N{n}_Q{n+1}_B{192 if n==6 else 256}.json') for n in [4,5,6]}
    responses={}
    for n in [4,5,6]:
        candidates=[json.loads(p.read_text()) for p in P.glob(f'HOST_RESPONSE_N{n}_Q{n+1}_B*_L*.json')]
        responses[n]=max(candidates,key=lambda d:d['level'])
    geometry_rows=[];source_rows=[];comparison_rows=[];decisions=[]
    for n in [5,6]:
        data=dynamics[n];first,last=data['rows'][0],data['rows'][-1];D=n-2
        geometry_rows.append(f"| {n},{n+1} | {100*first['geometry']['connectivity_probability']:.6f}% | {100*last['geometry']['connectivity_probability']:.6f}% | {100*first['geometry']['diameter_distribution'][str(D)]:.6f}% | {100*last['geometry']['diameter_distribution'][str(D)]:.6f}% | {100*last['total_variation_from_uniform']:.6f}% |")
        response=responses[n];initial=response['initial'];qrow=next(r for r in response['rows'] if r['lag']==.25)
        for d,shell in qrow['shells'].items():source_rows.append(f"| {n} | {d} | {interval(shell['response']['interval'])} | {'yes' if shell['response']['zero_excluded'] else 'unresolved about zero'} |")
        near=qrow['shells']['1']['response']['absolute_interval'];far=qrow['shells'][str(max(response['initial_host_distances']))]['response']['absolute_interval']
        hierarchy=near[0]>far[1]
        decisions.append(dict(N=n,lag=.25,near_absolute_interval=near,farthest_absolute_interval=far,
            far_below_near_resolved=hierarchy,
            far_to_near_upper_bound=far[1]/near[0] if near[0]>0 else None,
            all_adjacent_distance_orders_resolved=all(qrow['shells'][str(d)]['response']['absolute_interval'][0]>qrow['shells'][str(d+1)]['response']['absolute_interval'][1] for d in range(1,max(response['initial_host_distances'])))))
        comparison_rows.append(f"| {n} | {100*initial['class_probability']:.6f}% | {100*initial['fixed_source_peripheral_class_probability']:.6f}% | {initial['one_labelled_host_probability']:.8g} | {interval(response['source_variance_interval'])} | {interval(qrow['host_retention_interval'])} |")
    cost_rows=[]
    for n,data in dynamics.items():
        r=data['rows'][-1];cost_rows.append(f"| {n} | {interval(r['complete_record_work_interval'])} | {interval(r['complete_record_entropy_interval_bits'])} |")
    six=dynamics[6]['rows'][-1];r6=responses[6];checks=read('RESPONSE_CHECK.json')
    source_summary='\n'.join(f"- N={r['N']}: farthest/nearest absolute response is bounded above by {r['far_to_near_upper_bound']:.7g}; all adjacent distance orders resolved: {r['all_adjacent_distance_orders_resolved']}." for r in decisions if r['far_to_near_upper_bound'] is not None)
    closure=r6.get('initial_boundary_closure');closure_note=''
    if closure:
        closure_note=f'''A final numerical refinement stayed within the declared 1.2-million-orbit
cap: all {closure['added_initial_rows']:,} still-missing physical rows reached
from the initial numerical source support were added. Its two initial
outgoing residual norms changed from
`{closure['initial_residual_norms_before']}` to `[0,0]`.
This changes the numerical projection, not H, the host or preparation.
The level-8 pilot is retained for interval-overlap and error comparison;
the level-9 refined result is primary. Later-time residuals and unresolved
signals remain explicit; closing the first step is not exact evolution.
'''
    md=f'''# Full-active-host geometry feasibility

**Result: conditional transport on the full host is testable and exhibits
a strong five-site distance hierarchy, but this bounded study does not
establish a dynamically selected spatial phase or a distinctive prediction.**
The six-site result is reported with its actual precision below, including
unresolved receiver signals. Frozen CP-12 is unchanged in every Path A run.

This continues `609db3f`; checkpoint `c2751a0` preserved the protocol and
five-site marginals after a runtime reset. No lost run counts as evidence.
`PROTOCOL.md` fixes the budget, host choice and claim boundaries.

## 1. Fixed edge count is not emergence

For Q>binomial(N-1,2), every host is connected and has diameter at most
two. Thus primary N=5,Q=10 is K5, and **all 455** N=6,Q=12 hosts have
diameter two. Their all-root mean balls `[1,5,6]` and mean pair distance
1.2 are forced for every state and every time. They cannot test extended
distances or establish three-dimensional geometry. This does not exclude
all larger-N phases.

The auxiliary Q=N+1 sequence permits maximum diameters three and four
at N=5,6. It is explicitly different from primary Q=2N scaling.

## 2. Original-neutral dynamics really bias hosts, modestly

Evolve the original equal coherent neutral ray over **all** labelled
fixed-Q hosts to times 0,0.5,1,2,5. The entire active graph is used: no
circulation, matter or electric-state filter. Exact labelled graph-class
counts supply the uniform time-zero baseline.

| N,Q | Connected at time 0 | Connected at time 5 | Maximum-diameter mass at 0 | Maximum-diameter mass at 5 | TV from uniform at 5 |
|---|---:|---:|---:|---:|---:|
{chr(10).join(geometry_rows)}

These shifts are resolved well beyond the full-state numerical errors.
It would be wrong to call them purely combinatorial. However, connectivity
was already high in the supplied ensemble, and ball growth changes very
little. Six-site connected-host mean balls are
`{[round(v,8) for v in dynamics[6]['rows'][0]['geometry']['mean_balls_given_connected'][:5]]}`
at time zero and
`{[round(v,8) for v in six['geometry']['mean_balls_given_connected'][:5]]}`
at time five. Mean connected pair distance changes only from
{dynamics[6]['rows'][0]['geometry']['mean_pair_distance_given_connected']:.8f}
to {six['geometry']['mean_pair_distance_given_connected']:.8f}.

Average degree four is not dimension: supplied N=64 four-regular
next-neighbor-ring and periodic-square controls have balls `[1,5,9,13]`
and `[1,5,13,25]`. No fitted dimension is extracted from these small hosts.

## 3. Host-only source experiment

At waiting time five, the complete instrument records every labelled
host P_g. The feasibility branch is chosen **without inspecting amplitudes**:
maximum connected diameter, greatest labelled multiplicity, minimum
canonical mask; choose the largest peripheral source orbit, then minimum
label. No circulating-core, matter or electric state is selected. The
initial host defines distances for every lag; all later rewiring remains.

| N | Host-class probability | Fixed source in chosen peripheral class | One labelled host | Source variance interval | Exact-host retention at lag 1/4 |
|---|---:|---:|---:|---|---|
{chr(10).join(comparison_rows)}

Use the gauge-invariant phase pulse exp(-i theta n0), with demonstrably
nonzero variance, and measure `chi_0j=2 Im <x,nj y>`. Responses are signed
susceptibilities, not arrival probabilities. Lags 1/8,1/4,1/2 and all fixed
receivers are in the data. Distance-shell means here coincide with fixed
receiver stabilizer orbits; all individual receivers are separately given.

| N | Initial host distance | Raw response interval at lag 1/4 | Nonzero resolved? |
|---|---:|---|---|
{chr(10).join(source_rows)}

Far-below-near decisions and bounds are in `ASSESSMENT.json`. An interval
containing zero remains unresolved, not zero. This is finite-size evidence
conditional on the specified host record. It does not establish an
asymptotic front, a state-independent effective speed or relativistic
propagation. Microscopic hopping on active edges is already part of H.
For six sites, finer distance orders and nonzero distant signals must
not be claimed when their intervals remain unresolved. The neutral-state
preparation bound below is distinct from the larger source-evolution
errors reported for each lag.

At lag 1/4 the direct interval comparisons give:

{source_summary}

These are upper bounds on response magnitude ratios, not estimates of
unresolved distant responses and not an inferred propagation speed.

![Host dynamics and source responses](ACTIVE_HOST.svg)

The figure contains actual model data. Receiver intervals include every
unretained physical row. Downward limits denote unresolved magnitudes;
lines connect finite sampled data without extrapolation.

## 4. Persistence largely reflects a supplied coefficient

An exact full-law identity gives

`B_g^* B_g = (M-Q)/(100 M^2) sum_(e in g) Pi_e`,

where B_g takes a state from one fixed host to all others and Pi_e is its
neutral/removable-edge projector. Therefore

`P(same host at t) >= 1-Q(M-Q)t^2/(100 M^2)`.

At N=6,Q=7, every recorded host, of any shape, has lower retention bound
**0.99937778 at lag 1/2**. For Q=alpha N the cross-host norm is asymptotic
to sqrt(alpha/(50N)); recorded hosts become nearly static at fixed time.
That follows directly from the input rewiring normalization 1/(10M), not
an earned selection of spatial neighborhoods. Sampled exact-host, edge
and source-neighborhood retention are all reported. They are not
continuously monitored first-passage survival probabilities.

A Duhamel comparison with the original frozen-host block gives state
error <=c t and response difference <=4c t+2(c t)^2. All actual source
runs retain full H; the block is a bound/control, not a successful
modified-law experiment. The bound is too loose at these sizes to supply
a sharp propagation front. It cannot be applied to an unrecorded coherent
superposition of many hosts.

## 5. Complete recording is an intervention

Complete-host recording work is `-<R_off>`, evaluated without losing
rewiring that appears on signed-quotient diagonals. All outcomes are
included. At time five:

| N | Unconditional record energy change | Record entropy, bits |
|---|---|---|
{chr(10).join(cost_rows)}

The initial neutral record has work Q(M-Q)/(10M) and entropy
log2 binomial(M,Q). These are not a measuring-device derivation or a heat
dissipation claim. A full labelled N=4 control confirms that recording
also changes later source response: at lag 1/4 its nonselective change
interval is `{interval(checks['complete_record_effects'][0]['response_change_interval'])}`.
No classical host pointer is derived merely by using this instrument.

## 6. Large-N baseline is not large-N quantum dynamics

Seeded fixed-Q uniform graph samples and exact isolated-vertex moments
are reported separately. For Q=alpha N, the uniform ensemble has
E(I)~N exp(-2alpha) isolated vertices and Var(I)=O(N), so its connected
probability tends to zero. The exact finite-N variance bound at N=8192
is about 0.006064 for Q=2N and 0.0005363 for Q=N+1.

Those are **combinatorial** bounds. Full quantum dynamics could bias the
ensemble; this packet does not simulate it beyond six events. Ball means
from the samples are not a derived dimension or continuum limit. The
stated Hoeffding confidence bands cover connected fractions and bounded
ball means; diameter counts are descriptive sample counts.

## 7. A concrete Path B proposal, with an exact obstruction

`PATH_B_DESIGN.md` adds a neutral species with transport coefficient
`J g_e 1_(E_ei E_ej<0)`. This is a **new law**, never used in Path A.
The rule is local in primitive edge fields, gauge preserving and
Hermitian. In probe vacuum the original geometry law is preserved exactly.

But an allowed one-leg charged chord activates in one geometry hop.
The exact specified transfer coefficient is `J^2 T^4/4` at fixed N,
independent of its growing ring distance or ring flux. Its activation
diagonal detuning is -15/8, not a protective gap. Rational controls check
N=4,6,8,12 and fluxes 1,3,8. These supplied-state transfer controls do not
count as original-neutral phase-source experiments or a uniform
fixed-time large-N no-go theorem. The gate alone has not earned locality;
the next design needs a protected low-energy space and uniform
source-to-receiver leakage/virtual-path bounds.

## 8. Quotient-value continuity and decision

The latest review's proposed SCOUT-1 campaign already ran and is frozen
at `a2987fe`. Its actual handoff reports zero distinctive empirical
predictions and retains supplied representation/composition and basin
choices. `REVIEW_CONTINUITY.md` records the checked sources. The cited
`reviews/grut_review_report.md` itself was not found in the checked trees;
no unseen contents are claimed to have been read.

This packet does not reopen that campaign or promote a quotient pattern
to an axiom. H, its coefficients, Q, the preparation, quantum/gauge
structure, external time and the recording instrument remain inputs.
The conditional numerical outputs are not unique fundamental constants
or distinctive empirical predictions.

**Path A remains unestablished as an emergent spatial phase.** It has
earned a better conditional transport candidate than the failed ring
metric and a real, modest finite-size host bias. It has not jointly
demonstrated dynamically selected extended geometry, controlled
large-distance fronts or robust continuum scaling. The budget is closed
without another favourable-host search. Develop Path B's low-energy
protection mathematics before further large new-model runs; preserve all
failed benchmarks and the input ledger. No universal CP-12 no-go or
Milestone A/B/C success is claimed.

## Verification and reproducibility

The six-site neutral calculation has 1,980,546 signed orbits, 133,348
electric-boundary rows and full-state bound {six['state_error']:.9g} at
time five. The host-only source calculation uses {r6['dimension']:,}
retained and {r6['outside']:,} outgoing rows at level {r6['level']}.
Source/error intervals include normalization, numerical initial tails,
all omitted hopping and rewiring, and analytic integration upper sums.

{closure_note}

Independent labelled N=4 source/work/entropy controls, exact raw
cross-host norm identities, all dense primary hosts, equal-degree
counterexamples and the Path B rational controls are executable. Source
checks agree within their intervals. Ordinary floating point and its
reporting allowances are explicit; these are not outward-rounded formal
certificates or external validation. `reproduce.py` verifies saved
evidence; `--full` rebuilds it serially. Canonical GRUT, Stage3, PR12 and
all inherited Hamiltonian sources remain unchanged.
'''
    (P/'REPORT.md').write_text(md)
    assessment=dict(numerical_status='PASS',path_A_spatial_phase='NOT ESTABLISHED',
        conditional_source_tests=decisions,host_selection='Resolved finite-size bias, no continuum/low-dimensional phase established',
        persistence='Universal recorded-host bound is driven by supplied rewiring normalization',
        path_B='Gauge-preserving formation identity, but simple through-field gate fails activation protection',
        distinctive_empirical_prediction=False,earned_absolute_scale=False,SCOUT1_reopened=False,
        no_go_for_every_CP12_phase=False,hamiltonian_unchanged_for_path_A=True,
        source_recording_is_intervention=True,large_N_quantum_dynamics_simulated=False)
    (P/'ASSESSMENT.json').write_text(json.dumps(assessment,indent=2)+'\n')


if __name__=='__main__':main()
