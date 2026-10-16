"""Generate conservative research report and acceptance ledger from saved results."""
from pathlib import Path
import json

P=Path(__file__).resolve().parent


def main():
    r=json.loads((P/'RESULTS.json').read_text());e=json.loads((P/'CONFIGURATION.json').read_text())
    c=r['constants'];form=r['formation'];last=form[-1]
    blocks=sum(len(x['rows']) for x in r['Schur'])
    near=sum(len(x['near_packet_checks']) for x in r['Schur'])
    pairs=sum(x['source_exit_pairs_checked'] for x in e['rows'])
    graphstates=sum(x['complete_reachable_states'] for x in e['rows'])
    prop=sum(len(x['rows']) for x in r['full_propagation'])
    pulse=sum(len(x['rows']) for x in r['sources'])
    ledger=max(abs(x['complete_energy_ledger_error']) for x in form+r['detuned'])
    table='\n'.join(f"| {x['N']} | {x['formation_constants']['observation_window']:.3f} | {x['finite_window_macroscopic_probability']:.9f} | {x['infinite_time_full_probability']:.9f} |" for x in form)
    schur='\n'.join(f"| {x['N']} | {x['dimension']} | {x['P_dimension']} | {x['Q_dimension']} | {len(x['rows'])} |" for x in r['Schur'])
    detuning='\n'.join(f"| {x['N']} | {x['delta']} | {x['all_time_macroscopic_upper']:.9g} |" for x in r['detuning_large_N_bounds'] if x['N'] in [2048,32768] and x['delta'] in [.01,.1,1.])
    (P/'REPORT.md').write_text(f'''# WM2-S2: resonant assembly with protected transport, without a robust phase

S2 resolves S1's closed-seed formation obstruction **at a supplied exact
resonance**. Explicit event storage absorbs released binding energy.
A macroscopic path forms with nonvanishing observation probability,
and a new complete virtual-transport bound survives in the enlarged
protected band. The same calculation rejects a robust untuned-phase
claim: any fixed nonzero storage detuning restores an all-time formation
obstruction at increasing size. A separately tested branching proposal
reopens a constant-length protected shortcut.

This is a constructive conditional result and two architecture failures.
It is not a full WM2 success, a dynamically selected physical dimension,
a theory of gravity or a distinctive prediction about nature.

## What changed from frozen S1

The new Hamiltonian is specified in MODEL.md. It adds one two-level
storage register q_i per event. A head-coupled attachment transfers U
from the free-event counter to that event's storage bit. All adjoints
remain present. There is no bath, irreversible rate or external schedule.

The independent S1 bare append is replaced by storage-assisted append;
unsynchronized append/pop is an explicit probe backreaction, absent in
probe vacuum and retained with all fuel mismatches after probe injection.
This is a separately named law, not an unchanged S1 extension. The low
sector is now **fuel consistency M=0 at every assembly stage**, rather
than S1's full-assembly K=0 sector. No coordinates or fixed event-pair
adjacency enter the primitive law. The path class and normalization rule
remain supplied inputs.

The original CP-12 neutral ray is tensored with the declared bondless
uniform-root storage seed and probe vacuum. The product benchmark keeps
the old Hamiltonian but transfers no energy from it. The new neutral
probe is injected later. Common geometry for original charged CP-12
matter and a low band for the unbounded combined Hamiltonian are not
claimed.

## Formation with a finite, unselected observation window

The complete vacuum component reduces exactly to an N-state radial
Jacobi matrix with links -g sqrt((N-m)/N). Its finite-oscillator
commutator proves equal seed weight 1/N in every eigenstate. Consequently
the infinite-time radial ensemble is I/N; this is an exact reduction,
not a cutoff in physically reachable storage states.

For L=max(2,ceil(N/2)), a single observation at a uniformly sampled time
in [0,T_N], with T_N=32 H_(N-1) sqrt(N(N+1/2))/(pi g), finds at least L
attached events with probability **at least 3/8 for every N>=2**. At
least 1/6 of that window has instantaneous success probability >=1/4.
This follows from the Hermite-root spacing bound and finite dephasing
estimate in MATHEMATICS.md. No favorable instant or label ordering is
selected, but the declared window grows as N log N. This is not an
all-times success probability or continuous-observation hitting law.

| N | Declared window T_N | Exact finite-window p(m>=L) | Infinite-time p(m=N) |
|---:|---:|---:|---:|
{table}

At N=128 the finite-window probability is {last['finite_window_macroscopic_probability']:.9f}.
Full assembly has infinite-time probability 1/N; a nonvanishing
full-N connected phase has not been demonstrated. Unused events remain
isolated in each partial record. Every successful macroscopic relational
record is retained, rather than selecting one exponentially rare word.

## Re-derived complete protection and actual response

Primary supplied parameters: U=16, Lambda=64, g=1, epsilon=J=1/4,
delta=0. The mismatch diagonal begins at 48. The complete off-diagonal
norm is <=3 and PHP norm <=2.5. Thus the auxiliary one-probe low band
is derived with gap **>=42.5** and ||Pi-P0||<={c['spectral_projector_distance_upper']:.9f}.
All 2^(N-1) non-root storage patterns are present in the probe calculation.

Even allowing real assembly and reversal, the minimum protected
configuration operations between attached events equals their path
separation. Every P-to-Q-to-P transfer across distance d requires at
least d+2 primitive operations, with arbitrary final assembly size and
ordering. The complete weighted expansion then gives, for every E<=2.5,

`||R_j^P Sigma(E) P_(i,d;j)|| <= (1/680) (6/91)^d`.

Sources may be arbitrary coherent packets over all qualifying geometry
records. Receivers include every consistent geometry of every size.
For d>1 the direct PHP block vanishes, making this the complete
effective off-metric block bound. Adjacent source packets retain
receiver action >={c['near_packet_lower']:.9f}. The S1 Ud activation
energy barrier does not survive; the new proof pays for resonant
geometry changes through the prefix word bound and full weighted sum.

The full-H Dyson bound is tail_d(3|t|), including inconsistent final
storage states. For an adjacent two-site injected probe packet, source
variance is 1/4 and the exact number-phase response satisfies
|chi_k|<=tail_(d(i,k))(3|t|) tail_(d(j,k))(3|t|). Nearby response at
t=0.005 is bounded below by 0.000795, uniformly in N and in the formed
word. This is a response front in the declared one-dimensional class,
not solely an effective-H element. There is no exponentiation of a
fixed-energy Schur matrix.

Conditional bulk persistence follows from a second full-H estimate.
Destroying a prefix r levels below the current end has probability at
most [tail_r(2(g+epsilon)|t|/sqrt(N))]^2. For a macroscopic interior
prefix this tends to zero even at t proportional to N. The bound uses
the supplied 1/sqrt(N) normalization. Storage is still reversible:
spent charge is not an irreversible lock on one entire word.

## Exact energy and recording account

The initial auxiliary energy is UN+delta, including the supplied root
storage bit. At resonance, free energy U(N-m) plus stored energy Um is
UN for every component. Net attachment stores energy already present
in S1's free-event seed. At detuning, kinetic energy supplies the
remaining delta per attachment. The largest complete-ledger residual
in the finite controls is {ledger:.3g} in the declared energy units.

RESULTS.json reports the complete geometry-record entropy, root included,
and mean projective-record work. Large records have an explicit growing
information cost. Zero mean record work at resonance does not imply a
real apparatus has zero energy, entropy or heat cost. Probe preparation
is an additional controlled operation.

## Robustness failure: local detuning is allowed

The perturbation delta sum q_i is on-event, label symmetric and allowed
by all declared conservation rules. It adds delta(m-1) to the radial
Hamiltonian. For positive denominator the conserved-residual theorem gives

`p(m>=L,t) <= [(g sqrt((N-1)/N)+g sqrt((N-L+1)/N)) / (|delta|(L-1)-2g)]^2`.

For every fixed nonzero delta this is an **all-time O(N^-2)** upper
bound on macroscopic assembly. It is not an extrapolation from a finite
time grid. The exact-resonance equality is supplied, with no earned
symmetry forbidding this detuning. The resonant spectral constants are
not asserted for the extensively detuned band.

| N | Fixed detuning | All-time macroscopic formation upper bound |
|---:|---:|---:|
{detuning}

Therefore S2 fails the frozen robust-phase gate even though its exact
resonant construction has a positive formation theorem.

## Dimension control: naive branching fails before large dynamics

The separately specified two-child/leaf-removal Hamiltonian has no
supplied pair adjacency. However, a leaf can be detached, attached to a
distant leaf, used for one probe hop, detached again and restored. Every
step stays storage-consistent. The whole original tree and all its
storage bits are restored in **five moves**, independent of old distance.

Exact complete fifth-order controls at N=7,13,25,49 give coefficient
-2 with unit primitive amplitudes and zero lower coefficients through
four. With declared normalization this is -2g^4 J/N^2. This coefficient
is not a uniform finite-time large-N theorem. It is a decisive primitive
word obstruction to retaining raw tree distance as the protected metric.

In allowed complete binary-tree records, a radius-five operational ball
at one leaf contains O(N) leaves. Uniform finite-dimensional ball growth
fails on that branching class. No claim that such trees are typically
formed or that all branching grammars fail is made.

## Reproducibility and acceptance

{graphstates:,} finite register/configuration states and {pairs:,} source/exit
tests check independent register rules, label covariance, reverse moves,
energy changes and complete virtual paths. The matrix checks are:

| N | Full one-probe dimension | P0 | Q0 | Schur operator blocks |
|---:|---:|---:|---:|---:|
{schur}

There are {blocks} all-record Schur blocks, {near} adjacent packet checks,
{prop} full-H propagation blocks and {pulse} finite-pulse checks, including
formed partial words. Radial controls reach N=128; analytic results
hold for all N in the declared component. Floating-point controls are
ordinary arithmetic, not outward-rounded numerical certificates.
No large CP-12 trajectory, external validation or empirical fit is reported.

Run `python research/reality_framework/wm2_fuel_lock_20261016/reproduce.py --recompute`
from the repository root. The default command checks saved evidence,
all packet hashes and frozen dependencies without repeating the matrices.

S2 earns resonant macroscopic assembly, complete conditional virtual
protection, actual source-response bounds and bulk persistence. It
does not earn parameter-robust formation, irreversible locking,+dynamically selected dimension, original-matter universality, continuum
gravity or an empirical prediction. ASSESSMENT.json preserves each gate.
''')
    a=dict(candidate='WM2-S2',research_base_commit='7b1d9234deaeaddea3f34c419bdb35b99efd944b',
        explicit_energy_ledger='PASS: released binding energy enters declared reversible storage bits',
        head_only_storage_coupling=True,
        macroscopic_formation='CONDITIONAL PASS: delta=0; uniform-time finite-window probability >=3/8',
        full_assembly_phase='NOT ESTABLISHED: diagonal-ensemble full probability 1/N',
        operational_metric='PROVED in the single-end class including all protected assembly moves',
        complete_virtual_protection=True,source_response_front=True,bulk_prefix_persistence=True,
        rederived_gap_lower=42.5,virtual_word_lower='d+2',complete_Schur_bound='(1/680)*(6/91)^d',
        all_storage_patterns_retained_with_probe=True,
        robust_untuned_formation='FAIL: every fixed nonzero storage detuning has all-time O(N^-2) bound',
        irreversible_lock='FAIL: reverse storage-assisted pop remains active',
        selected_dimension='NOT EARNED: single-end grammar supplies path class',
        naive_branching='FAIL: five protected moves transfer and fully restore across growing raw tree distance',
        common_original_CP12_matter_geometry=False,full_WM2_success=False,
        all_WM2_architectures_ruled_out=False,frozen_CP12_modified=False,S1_modified=False,WM1_modified=False,
        SCOUT1_reopened=False,empirical_prediction=False,earned_quotient_value=False,
        large_CP12_dynamics_run=False,external_independent_validation=False)
    (P/'ASSESSMENT.json').write_text(json.dumps(a,indent=2)+'\n')


if __name__=='__main__':main()
