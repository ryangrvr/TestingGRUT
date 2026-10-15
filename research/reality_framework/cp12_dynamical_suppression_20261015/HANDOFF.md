# GRUT builder/reviewer handoff: dynamical suppression

Parent: `768c46c96eae404a527d9cfee8f1b09abeb61a13`.
Continuation branch: `codex/grut-dynamical-suppression-20261015`.
Read `REPORT.md` for the evidence and `DERIVATION.md` before interpreting
the mechanism ratios. `VERIFICATION.json` and `manifest.json` describe the
checks and exact file provenance.

## What was built

The requested action ratio is now computed under the full frozen
Hamiltonian from the original neutral preparation and its explicit recorded
quarter-threshold charged event. Both integrals, formation weight, outgoing
state-error bounds, event survival and record resources are reported.
No hopping deletion is used in the new evolutions.

The packet adds physical-basis interference references, bounded Pauli and
diagonal-detuning statistics, exact blocking/unlocking controls, an
exploratory neutral-origin filled-ring test and a fixed-ring neighbour
density-response control. The six-site approximation includes two complete
outgoing boundary layers beyond the original numerical basis, constructed
through an intermediate priority extension.

## Main physical findings

At $N=5,Q=10$, wait 2 and delay 0.25, the original primary-event probability
is about 5.963%, but

$$
I_{\rm out}\in[4.65225,4.67287],\quad
I_{\rm in}\in[0.940236,0.948982],\quad
\eta\in[4.90464,4.96759].
$$

Outside hopping action exceeds inside hopping action. Its coherent action
is 2.40842–2.42273 times the physical diagonal reference, indicating
constructive integrated cross terms. About 52.6% of neutral endpoint action
remains Pauli available. The capped diagonal-detuning RMS is about 2.08;
these differences are not dressed spectral gaps.

All nine tested primary five-site ratios have lower bounds above 2.53.
This is negative evidence for the proposed aggregate suppression mechanisms
in these states and time windows, not a no-go theorem for every CP-12 phase.

The new six-site run resolves an increase in the early matched ratio:
at wait 0.5, delay 0.25,

$$
\eta_5\in[4.98907,4.99173],\qquad
\eta_6\in[5.78217,13.04981].
$$

The six-site conditional integrals are
$I_{\rm out}\in[6.67552,8.23764]$ and
$I_{\rm in}\in[0.629002,1.158618]$.
Its original formation probability is only
$[7.17034,7.19592]\times10^{-6}$, so this is not an appreciable phase.

The neutral waiting-state bound improves about 90-fold at wait 0.5 and
17-fold at wait 1. Nevertheless the later charged comparison remains open:
at wait 1, delay 0.25, the six-site $\eta$ interval is
$[0.59687,628.70587]$, overlapping the five-site interval.
At delay 0.5 its core-denominator lower bound is zero, so no finite upper
ratio is resolved. Do not use the later central estimates to claim a
scaling trend. The new reduction improves precision substantially but
does not establish later-time locality or its failure at all sizes.

The rare fully filled five-ring event really does arise from the original
neutral ray, with probability $[5.58411,8.57783]\times10^{-8}$ at wait 2.
It has $\eta(0.25)\in[0.0386366,0.0386617]$ and nonzero core hopping action.
That is a positive finite-time Pauli control, but its formation is far too
rare to count as an appreciable phase. Its immediate density phase source
is a global phase, giving exactly zero response to that particular probe.
Other source operators or later probes remain open.

The auxiliary original-neutral $N=4,Q=5$ fixed-ring record gives a resolved
distance-one neighbour response. At delay 0.25,

$$
\chi_{\rm neighbours}\in[0.0038018298,0.0038019548],\qquad
\Pr(\text{same ring later}\mid\text{record})\simeq0.954664.
$$

The response joint with that same ring still present is also nonzero.
This is a useful probe control, not a $Q=2N$ scaling result and not proof
that ring paths alone carry the signal. The parent already resolves active
diagonal transport across that ring.

## Two exact controls the acceptance rule must respect

One additional exact sector identity narrows the interference search.
If the recorded core has a neutral exterior, first exterior hops have
orthogonal edge-hole/endpoint-field tags. Their squared action is exactly
$\langle\sum_{e=(i,j)\notin C}n_{g_e}(2-n_i-n_j)\rangle$ even in a coherent
superposition. Destructive interference cannot cancel that first action
within the neutral-exterior sector. It can still act later, or in states
with exterior electric dressing. Do not apply this identity to the whole
primary event without checking its support.

For a supplied charged spanning ring with all vertices filled and all
non-core edges neutral,

$$
V_{\rm out}\phi_N=0,\quad
\|V_{\rm out}H\phi_N\|^2=4N,\quad
\|K_{\rm in}\phi_N\|^2=2N.
$$

The full law unlocks blocked neutral edges after a core hop creates a
vacancy. At each fixed $N$,

$$
\eta_N(T)=\frac23T^2+O(T^4).
$$

For the one-vacancy version, the instantaneous ratio is $1/(N-1)$ although
a specified chord spans ring distance $\lfloor(N-1)/2\rfloor$ and has
transfer probability coefficient $T^4/4$. The core denominator is extensive.

These are exact supplied-state mechanism controls, checked for
$N=5,\ldots,12$, with counting proofs in the derivation. They are not
neutral-origin formation results, and their Taylor remainders are not
bounded uniformly in $N$. Do not promote them to an established fixed-time
large-size nonlocality theorem. Their valid conclusion is that the initial
global ratio and initial blocking do not by themselves establish locality.

## Numerical construction worth retaining

The new retained dimensions are 940,985, 1,421,197 and 2,477,611 orbits,
starting from 715,985. Every new-to-old coupling is independently regenerated
and checked against the stored transpose. All outgoing states remain in the
error matrix; the electric fields there are not cut off.

An adaptive outgoing matrix includes rewiring as well as hopping. Use the
change in vertex filling to separate them when constructing hopping action;
count both in the state residual. The small independent test deliberately
contains outgoing rewiring and catches the distinction.

The new filling-sector Jacobi bound proves
$\|K_{\rm in}\|,\|V_{\rm out}\|\le17$ at six sites, versus the older bound 24.
All leading minors of $17I-J$ are strictly positive integers. This bound
holds for all electric fields; it does not assume the retained budget is
the physical Hilbert space. Its five-site counterpart is 14.

Large propagation now reuses the real sparse $H$ through a LinearOperator
with its exact adjoint. It releases Python state labels only after all
physical observables and action matrices are constructed. The independent
small calculation agrees with explicit sparse propagation, action integrals
and residuals to about $10^{-15}$.

## Next concrete physical test

Keep the requested three-part gate. No experiment here passes it:

1. Appreciable original-neutral probability for the candidate geometry.
2. Persistence of the same relational geometry with nonzero intrinsic-distance
   matter response.
3. Suppression of source-resolved propagation that bypasses those distances.

The primary charged event does not store a metric. Its moving-core action
also counts loop creation/destruction, and its two-time event survival can
replace one loop by another. The next transport experiment therefore needs
an explicit fixed relational core record and gauge-invariant source and
receiver operators labelled by that record.

For a density-phase source, first report its variance in the prepared
branch. Zero variance makes the pulse a global phase; a vanishing response
then cannot diagnose an emergent metric. Use an active source and compare
nearby and distant receivers at the same physical delay, retaining the
full-H amplitude and the original formation/record cost. Normalize the
comparison per source rather than against all extensive core activity.

Use the one-vacancy exact chord control as a diagnostic check: the analysis
must identify its direct long-ring-distance matrix element even while the
initial global action ratio decreases.

Five-site $Q=10$ hosts are complete graphs; six-site $Q=12$ hosts have
diameter at most two. Neither is a long-distance or three-dimensional
scaling experiment. Before interpreting larger circulation as space, the
builder must expose how actual matter response depends on the distances of
a geometry that forms under the original law.

GRUT remains the objective. These results justify a sharper locality test,
not a claim that physical spacetime or its continuum dynamics has been
derived.
