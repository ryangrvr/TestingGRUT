"""Generate the finite-mechanism report from saved numerical records."""
from pathlib import Path
import json,math
from decimal import Decimal,ROUND_CEILING,ROUND_FLOOR

P=Path(__file__).resolve().parent
FILES=['N4_Q4_B256_T20.json','N5_Q5_B128_T10.json','N6_Q6_B128_T10.json']


def load(name):return json.loads((P/name).read_text())


def at(d,t):return next(r for r in d['rows'] if abs(r['T']-t)<1e-8)


def up(x):
    if not x:return '0'
    d=Decimal.from_float(float(x));unit=Decimal(1).scaleb(d.adjusted()-3)
    return f'{d.quantize(unit,rounding=ROUND_CEILING):.3E}'


def endpoint(x,digits=6,upper=False):
    return format(Decimal.from_float(float(x)).quantize(Decimal(1).scaleb(-digits),
                  rounding=ROUND_CEILING if upper else ROUND_FLOOR),'f')


def main():
    models=[load(n) for n in FILES]
    lines=[]
    for d in models:
        r=at(d,10);e=r['bounded_fraction_error_upper']
        keys=['removable_fraction','flux_protected_fraction','low_removability_probability','graph_TV_from_initial']
        digits=13 if e<1e-8 else 8
        vals=[f'{r[k]:.{digits}f}' for k in keys]
        shown_error=up(e+max(abs(float(s)-r[k]) for s,k in zip(vals,keys)))
        lines.append(f"| {d['N']} | {d['budget']} | {vals[0]} | {vals[1]} | {shown_error} | {vals[2]} | {vals[3]} |")
    continuous=[]
    for d in models:
        b=d['continuous_time_bounds'];end=d['rows'][-1]['T']
        continuous.append(f"| {d['N']} | [0, {end:g}] | {endpoint(b['removable_fraction_lower_over_entire_window'])} | {b['persistent_flux_window']} | {endpoint(b['flux_protected_fraction_lower_over_late_window'])} |")
    internal=[]
    for name in ['RECORDS_N4_Q4_B256.json','RECORDS_N5_Q5_B128.json']:
        d=load(name);r=at(d,10)
        S=f"{r['graph_entropy_bits']:.8f}";Se=up(r['graph_entropy_error_upper']+abs(float(S)-r['graph_entropy_bits']))
        internal.append(f"| {d['N']} | {d['graph_count']} | {S} +/- {Se} | {r['graph_classical_entropy_bits']:.8f} | {math.log2(d['graph_count']):.8f} | {r['graph_coherence_l1_normalized']:.6f} |")
    mobility=load('MOBILITY.json');mob=[]
    for r in mobility['rows']:
        if r['lag']==.5:
            value=f"{r['susceptibility']:+.10f}";error=up(r['susceptibility_error_upper']+abs(float(value)-r['susceptibility']))
            mob.append(f"| {r['wait']:g} | {r['lag']:g} | {value} | {error} | {r['before_source_removable_fraction']:.6f} | {r['pulse_work'][0]:+.6f}, {r['pulse_work'][1]:+.6f} |")
    energy=load('ENERGY_AUDIT.json');er=[]
    for m in energy['models'][:3]:
        r=at(m,10)
        F=f"{r['electric_energy']:.12f}";Fe=up(r['electric_energy_error_upper']+abs(float(F)-r['electric_energy']))
        er.append(f"| {m['N']} | {m['exact_total_energy']:.1f} | {F} +/- {Fe} | {r['mass_energy_projected']:.7f} | {r['hopping_energy_projected']:.7f} | {r['rewiring_energy_projected']:.7f} |")
    deletion=[]
    for d in load('DELETION_AUDIT.json')['models'][:3]:
        r=at(d,10);lo,hi=r['coherent_to_dephased_ratio_interval']
        deletion.append(f"| {d['N']} | {r['coherent_to_dephased_ratio']:.6f} | [{endpoint(lo)}, {endpoint(hi,upper=True)}] | {d['maximum_amplification']} |")
    dense=load('N6_Q12_B6_T0.4.json')['rows'][-1]
    more=at(load('N5_Q6_B128_T10.json'),10)
    vacuum=', '.join(f"{at(d,10)['vacuum_flux_probability']:.7f}" for d in models)
    rawdense=load('N6_Q12_B10_T1.json');dense_rows=[];refined_rows=[]
    for t in [.4,.5,1.]:
        r=at(rawdense,t)
        dense_rows.append(f"| {t:g} | {r['removable_fraction']:.8f} | {r['flux_protected_fraction']:.8f} | {up(r['bounded_fraction_error_upper']+5e-9)} |")
    for d in load('REFINED_COUNTS.json')['models']:
        r=next(x for x in d['rows'] if x['detector']=='flux_protected_fraction')
        lo,hi=r['interval']
        refined_rows.append(f"| {d['budget']} | {d['T']:g} | [{endpoint(lo,7)}, {endpoint(hi,7,upper=True)}] | {'Yes' if lo>0 else 'No'} |")
    dense_r_lower=endpoint(rawdense['continuous_time_bounds']['removable_fraction_lower_over_entire_window'])
    text=f'''# CP-12 autonomous flux gate: internal graph records form; spacetime is not derived

**Frozen-law result:** an unforced, topology-neutral pure preparation develops persistent partial gauge-flux protection and internal graph records. Matter remains responsive in the tested four-event state. The unforced state does not reach low mean removability in the tested finite windows, and no extended spatial phase or causal cone is established.

Parent: `4e2b7a876a3a45ed554408ed803b572dd9cab580`. The independent extension is accepted and its previously unexecuted rational checker now passes. No coupling, microscopic law, spatial lattice, phase mask or pinned graph is added.

![Finite unforced protection and graph-record diagnostics](flux_gate.png)

The figure shows the common Q=N comparison through T=10. Shading gives truncation bounds at sampled times; connecting lines guide the eye. Whole-interval statements below use a separate derivative bound.

## 1. The proposed internal process actually occurs

Start with the same equal-amplitude neutral ray over every Q-edge graph. Every edge initially is removable, all vertex modes are empty and all electric fields vanish. Let r count removable edges, k count occupied-edge-fermion edges with nonzero flux, and n_v count vertex fermions. Exactly r+k+n_v=Q.

The main size comparison keeps Q=N at N=4,5,6 and uses full finite-time matrix exponentials, with all outgoing electric-budget transitions retained for a full-law leakage bound. The table is at T=10. The error column bounds each displayed normalized count/event and graph TV, including rounding of its displayed center; it is not a statistical standard error. Displayed interval endpoints are rounded outward relative to the evaluated bounds. This formatting does not certify floating-point roundoff in the underlying calculation.

| N=Q | Electric budget | Mean r/Q | Mean k/Q | Bounded-fraction error | Probability r<=0.1Q | Graph TV from initial |
|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(lines)}

The hole-protected fraction is the remainder 1-r/Q-k/Q. Flux protection and matter occupancy are separated, rather than counting every nonremovable edge as flux protected. These are unrecorded evolutions: no graph-distance apparatus, geometry measurement, low-r postselection or feedback controller enters them.

The probability of n_v=0 with nonzero occupied-edge flux at T=10 is {vacuum} for N=4,5,6 respectively. Each exceeds its bounded-event error in the table. In this sector all charge occupations are neutral and Gauss law makes the nonzero flux a closed circulation. Thus the flux signal is not exclusively open strings accompanying vertex excitations. These small probabilities do not establish an extensive density or persistence of local loops at increasing size.

An additional N=5,Q=6 run at T=10 gives approximately r/Q={more['removable_fraction']:.8f} and k/Q={more['flux_protected_fraction']:.8f}, with normalized-count truncation estimate {up(more['bounded_fraction_error_upper'])}. The denser N=6,Q=12 run is extended to T=1 with budget 10, retaining {rawdense['dimension']:,} orbit states and {rawdense['nnz']:,} Hamiltonian entries:

| Time, N=6,Q=12, budget=10 | Mean r/Q | Mean k/Q | General bounded-fraction error |
|---:|---:|---:|---:|
{chr(10).join(dense_rows)}

The general state bound alone stops resolving the small flux fraction at later times. Detector-specific forward/backward control gives E_O<=2 e d+e^2 instead, without adding a physical recording operation:

| Electric budget | Time | Controlled interval for mean k/Q | Zero excluded? |
|---:|---:|---|---|
{chr(10).join(refined_rows)}

The budget-six general bound at T=0.4 was insufficient, but its refined interval resolves flux formation; the larger budget independently controls the same early signal and its refined endpoint. The higher-density run has the continuous lower bound mean r/Q>={dense_r_lower} throughout [0,1], so it also stays outside the low-mean-removability regime. These are genuine finite higher-density flux tests. They are not a long-time dense-phase test, and six-event twelve-edge graphs have diameter two, so this sector cannot establish an extended causal cone.

## 2. Persistence is bounded over continuous time

Only hopping changes r. A basis row has at most 2Q hopping transitions, each changing r or k by at most one. Therefore both normalized expectation rates have magnitude at most 2|t_hop|=2. At time spacing 0.1, subtracting 0.1 plus the numerical bound from a sampled minimum controls every intervening time.

| N=Q | Whole window | Continuous lower bound on mean r/Q | Flux window | Continuous lower bound on mean k/Q |
|---:|---|---:|---|---:|
{chr(10).join(continuous)}

This establishes sustained partial flux protection and excludes the **specific low-mean-removability regime** throughout these finite intervals. It does not prove a thermodynamic lower bound from three sizes, exclude later behavior, rule out rare protected branches, or establish a no-go for every state. Low mean r is sufficient for an earlier shortcut-suppression bound, not a necessary criterion for all possible emergent locality.

Persistence here concerns the mean amount of protection. It does not measure the lifetime of the same protected edge, the same loop or an intrinsic metric; protection can turn over while its mean stays positive.

The raw low-budget records are retained as regulator controls. Some of their late-time error estimates are vacuous; their apparently converged curves are not promoted to full-law results. Refinement to budgets 128 or 256 is what resolves the claims above.

## 3. Gauge and matter variables write internal graph records

This is a new positive mechanism result. Trace the unforced pure state over its complete fermion/electric configuration in the explicitly stated kinematic embedding. The graph entropy is initially zero. It increases without an external record, while total entropy remains zero.

| N=Q | Labelled graphs | Graph entropy at T=10, bits | Graph-basis Shannon entropy, bits | Maximum graph entropy, bits | Normalized graph l1 coherence |
|---:|---:|---:|---:|---:|---:|
{chr(10).join(internal)}

The five-event calculation unfolds 748,410 physical amplitudes, including the signed fermionic permutation action. Independent full-space controls check the partial trace. Its 6.34-bit graph entropy is graph/environment entanglement, not entropy creation in a closed pure universe.

There is an exact reason protection and records are connected. Given a complete environment configuration mu, every edge with a hole or nonzero flux must be present. These variables reveal the protected core P(mu), leaving at most binom(M-Q+r,r) possible graphs. Thus H(graph | configuration)<=mean log2 binom(M-Q+r,r). If r=0, the environment determines the entire graph, and the graph reduced density matrix is diagonal. The saved r=0 audit verifies this; it does not select that sector in the main run.

Recording a graph and selecting a spatial graph remain separate requirements. The measured graph-basis distribution still has entropy close to its initial maximum, and its TV changes are modest. No three-dimensional volume law, preferred loop family or spatial attractor follows from the internal-record result. A record can preserve an arbitrary geometry.

Moreover, substantial residual coherent-removal capacity remains. Extending the earlier Johnson incidence identity to each physical matter/field fiber gives

    D^dagger D = (M-Q+1)r - (M/Gamma)R.

This evaluates the actual unforced state's coherent deletion norm. Its ratio to the hypothetical fully configuration-dephased value, mean r, is:

| N=Q | Coherent / dephased deletion norm at T=10 | Controlled ratio interval | Maximum ratio at matched mean r |
|---:|---:|---|---:|
{chr(10).join(deletion)}

The dephased comparator is not physically applied and would change energy and entropy. The identity is a bounded-operator resource diagnostic, not a finite-time transport response. The all-graphs preparation has no universally absent edge, so the common-vacant destination formula cannot simply be applied to it. These results nevertheless prevent replacing a measurement of residual coherent-removal capacity with graph entropy alone.

## 4. Matter can still respond after partial protection

At declared waits 5 and 10 in the same four-event unforced preparation, apply the gauge-invariant source exp(-i eta n_0). Detect mean density at the other vertices, independent of graph-distance bins. For the projector n_0, the difference between finite phases +/-pi/2 extracts the linear susceptibility exactly, without finite-difference bias.

| Wait | Lag | Density susceptibility | Absolute truncation bound | Before-source r/Q | Work of + and - pulses |
|---:|---:|---:|---:|---:|---:|
{chr(10).join(mob)}

All six saved source responses exclude zero under their analytic truncation estimates. Removing hopping only after the same preparation makes the response zero. The initial source operation commutes with r and does not instantaneously select protection. The diagnostic has explicit work, and no claim of equal-energy pulse preparations is made. This verifies finite matter mobility alongside partial protection; four events do not test propagation through an extended network or provide a causal speed.

## 5. The stronger resource constraint survives negative hopping

The independent extension proves F>=epsilon_E(Q-r-n_v/2), where F is electric energy. We derive a complementary bound that includes the possibly negative hopping energy. For nu=Q/N and a=2m-epsilon_E/2,

    H >= kappa Q + epsilon_E(Q-r) - N C(nu),
    C(nu) = [sqrt(a^2+8 t_hop^2 nu)-a]/2.

It follows from the singular values of the gauge-covariant vertex/edge hopping matrix, whose square sum is 2Q, and positivity of rewiring. At the exact neutral initial energy E=kappa Q, mean r/N>=nu-C(nu)/epsilon_E at every time. With frozen parameters, this forces a positive removable density for nu>49/2. The inequality is vacuous at nu=1,2,3, so it does not prohibit protection in a candidate degree-six spatial sector or settle the phase diagram.

Electric energy can also be controlled numerically despite the unbounded rotor. Total energy is known exactly, and mass+hopping+rewiring is bounded at fixed N,Q. Subtracting these bounded terms supplies a full-law electric-energy interval; norm error alone would not do so.

| N=Q | Exact total mean energy | Electric energy and inferred error | Mass energy | Hopping energy | Rewiring energy |
|---:|---:|---:|---:|---:|---:|
{chr(10).join(er)}

Complementary component errors are in ENERGY_AUDIT.json. Negative hopping balances positive mass/flux costs. The stronger cost must therefore not be misread as an immediate prohibition based on positive terms alone.

## 6. Protection does not supply a near-destination rule

At fixed matter/field configuration, unprotected neutral connections undergo a Johnson-graph Laplacian: all legal replacement destinations have the same bare coefficient. An exact deletion-map identity proves equal instantaneous edge-creation coupling to any two edges absent throughout a state's support. A supplied flux-triangle comparator has candidate distances 2,3,4 and identical coupling squares 1/7500.

This is an instantaneous statement, not a theorem of equal finite-time responses. The earlier higher-order topology sensitivity remains intact. It shows that the proposed mechanism needs more than protection alone: actual dynamics must organize the protected core spatially and suppress aggregate distant channels while preserving matter motion. Completely freezing an arbitrary graph does not meet that requirement.

## 7. What has and has not passed

| Requirement | Finding |
|---|---|
| Stronger flux/fermion inequality | Confirmed; independent exact checker executed |
| Autonomous flux creation from an unmasked input | Resolved in the tested finite states |
| Persistent partial protection | Controlled over continuous finite windows |
| Internal graph records without external recording | Demonstrated at four and five events |
| Matter-density response after partial protection | Resolved at four events |
| Low mean removability in the tested unforced state | Not attained in the covered windows |
| Spatially selected protected network | Not demonstrated |
| Intrinsic 3D growth and extended causal propagation | Not demonstrated |
| Thermodynamic spatial phase or gravitational continuum | Open |

Q=N gives a connected graph only one independent cycle. These main runs are auxiliary mechanism tests, not the higher-density gauge-capacity/spacetime gate. The supplied pure initial state has energy Q/2, variance 2Q, zero entropy and coherent support over all binom(M,Q) graph configurations. A compact description does not establish its local physical preparation cost or explain why the theory selects it.

Fresh controls include the independent Gauss checker, 16 dense residual and eight dense expectation tests, independent 872- and 6075-state physical matrices, orbit/partial-trace checks, finite source pulses and a fresh rerun of the small unforced trajectory. Larger saved evolutions, dense-sector count refinements and the five-event record unfolding are hash-validated rather than silently rerun by the verifier.

All numeric error evaluations use ordinary floating point, not formal interval arithmetic. The mathematical inequalities concern the unbounded frozen law; the numerical scope remains the stated preparations, sizes and times.

The earned mechanism is now concrete: hopping creates flux/occupancy protection and internal records while matter remains responsive. The missing selection of an extended spatial geometry is an additional scientific problem, not a relabeling of that mechanism as spacetime.
'''
    (P/'REPORT.md').write_text(text)
    print('REPORT.md generated',flush=True)


if __name__=='__main__':main()
