# WM2-S7 — positive full-action residuals and measured recovery

**Result: five-event recovery passes; the six-event target remains
computationally unresolved.** The new method avoids cancellation-prone
residual evaluation and uses less peak memory than fresh S5 recovery, but
is slower and returns wider probability bounds. At six events the rank and
error allowances bind before the required observation time; the stop is not
a negative formation result.

S7 implements an exactly closable special case of S6 action enrichment. A
shared receiver frame containing the physical contact orbital captures the
entire emission action, to all Krylov orders, while keeping every ordered
attachment history. A positive structural leakage bound handles imperfect
capture and avoids cancellation in squared-norm subtraction. The unchanged
five-event physical target is recovered with a complete roundoff-aware error
ledger under declared arithmetic assumptions. This is not the general
history-conditioned action-enrichment algorithm, a certified compiler/BLAS
stack or a proof of scalable formation.

## Frozen construction and provenance

The S7 charter was published at `4bf14483a0804980244765906e939e4dbaf29834`,
before controls or primary trajectories, on an isolated branch from S6
`b4c7902`. Its scientific gates are unchanged. Only this folder changes.
S5/S6 saved-record verifiers passed with both manifests preserved.

The original S4 law, neutral seed, all ordered type-word histories, reverse
transitions, original root weights and charge-allowed occupations remain.
Every orbital has all L=64N^2 physical frequency components; there is no
substituted small receiver, irreversible bath or extra photon cutoff. The
method compresses a receiver representation and explicitly pays the omitted
physical-state norm. CP-12 remains a decoupled benchmark factor, as before.

S7 uses the merged-group limit of the S6 atlas. Shared frames retain physical
interference by summing the complete Hamiltonian action in common Fock
coordinates. They deliberately forgo S6's per-history storage advantage.
The method does not solve the general unequal-frame residual calculation.

## Positive action bound

With an exact common isometry W and complete Fock lift V, contact in span W
implies `(1-P)BV=0` to every order. Otherwise,

\[
\|(1-P)BV\|\le {g\sqrt{N-1}\over N}
                 \operatorname{dist}(c,\operatorname{span}W).
\]

An annihilation remains in the frame; an append creates one outside photon.
Every child word has a unique parent, so the outgoing squares sum with the
exact remaining type multiplicities. There is no extra count of histories.
The coefficient-side action still contains every coherent append/pop return.

For a floating frame, directly enclose `||c-W*z||` plus dot-product,
subtraction and constant errors. This upper-bounds distance to its span;
Gram/polar errors separately bound the numerical projected Hamiltonian.
The physical unitary Duhamel identity pays the full residual and future
returns. No norm-subtraction residual is clipped to zero.

In a cancellation control, the true residual is `2^-40` (about 9.095e-13).
Squared-norm subtraction returns zero. S7's positive upper bound is about
9.122e-13. The small physical-matrix controls separately compare the full
unprojected action with the positive structural bound.

The cubic Hermitian exponential action has its unitary Taylor remainder
`(G*h)^4/24` paid. Its three native sparse actions, all basis operations,
normalization, full-H Strang error, physical phases and arithmetic have
separate ledger entries. See MATHEMATICS.md and ARITHMETIC.md.

## Completed five-event recovery

The frozen target has g=.9, common detuning zero, balanced offsets +/-1/4,
L=1600 and T=100. Macroscopic formation means m>=3; full assembly means m=5.
Root components have remaining counts (2,2) and (3,1), with weights 3/5 and
2/5. All final words and receiver fields contribute to each event.

| Component | Complete state-error account | Macroscopic lower | Full-assembly lower |
|---|---:|---:|---:|
| (2,2), weight 3/5 | 0.549138 | 0.698241 | 0.579042 |
| (3,1), weight 2/5 | 0.547993 | 0.699422 | 0.580094 |
| Original root-weighted preparation | component accounts retained | **0.698714** | **0.579463** |

Displayed numbers are rounded summaries; `summary.json` and component JSON
contain the actual outward endpoints. In particular, a displayed lower
rounded upward is not a new certified endpoint. Both upper endpoints are 1.
The gate (completed target, all error budgets, lower >=1/2) passes.

The computed macro populations are nearly one, but the conservative physical
lower bounds above pay the full error account. They are weaker than the
archived S5 bounds, and should not replace them as the sharpest finite
formation estimates. The advance is implementation plus the explicit
arithmetic model, not a claim that formation became more probable.

For the first component, the ledger is approximately .199980 Strang,
.348036 projection, .00110129 polynomial, 1.303e-5 basis, 6.645e-6
normalization, 5.863e-7 phase and 6.403e-8 other action arithmetic. The
positive full-emission residual is at most 1.537e-11 in the accumulated
account. The second component has the same dominant structure. Full raw
ledgers are retained; no original S6 small-tail estimate is retroactively
assigned to these runs.

## Physical and arithmetic verification

The checks compare compact native actions, adjoints, density matrices,
vacancy weights and normalized Fock rotations against immutable S5 Python
operators. Six complete small physical wavefunctions are propagated at
T=.5,2,8. Their full-state discrepancies lie below the complete ledger.
Rational positive-norm and independent rational Taylor tests check the
arithmetic primitives, including the cancellation challenge.

The independent N=3 physical reference uses **all 576 receiver modes** and
a **333,505-dimensional** full sparse physical matrix at T=36. It compares
the full wavefunction, including the complete two-photon tensor. S7's error
is about **0.092928**, below its account **0.540646**. This is a deliberately
loose physical control at the new allocations, not an equal-error comparison
to S5 or S6. It does not prove the larger runs or all-size mathematics.

`controls.json` reports PASS (about 11.4 seconds on this machine). A strict
native rebuild took about .242 seconds; the initial automatic build was not
separately timed. The actual toolchain, flags and arithmetic width are
archived in `environment.json`. No external mathematical review or formal
software/hardware validation is represented.

## Fresh whole-worker comparison

Each component ran sequentially in a fresh process, with imports, setup and
propagation included. Peak RSS is the process high-water mark; sequential
peaks are compared, not added. CPU times closely track wall times. Compiler
rebuild and validation costs above are separate.

| Method / remaining counts | Wall seconds | CPU seconds | Peak MiB | Maximum coefficients |
|---|---:|---:|---:|---:|
| S7 (2,2) | 120.095 | 120.084 | 85.941 | 163,619 |
| S7 (3,1) | 103.883 | 103.874 | 65.875 | 94,626 |
| Immutable S5 P1 (2,2) | 70.579 | 70.571 | 121.387 | 122,149 |
| Immutable S5 P1 (3,1) | 51.040 | 51.028 | 96.441 | 81,549 |

S7 totals about **223.978 seconds**, with peak **85.941 MiB**. Fresh S5
totals about **121.618 seconds**, with peak **121.387 MiB**. Thus this
implementation uses less peak memory but is slower; fewer/more coefficients
alone would not establish that conclusion. S7's maximum coefficient count
is actually larger in both components, while its compact operator maps
need less total resident memory.

The fresh S5 root-weighted macro/full lower bounds reproduce approximately
.771013/.688305, using the P1 complementary-event inequality. Its engine's
legacy `finite_instance_gate` field still tests original S5's .20 state
envelope, not P1's .60 protocol; S7 does not rewrite that engine. The added
P1 probability endpoints and arithmetic status are separate JSON fields.

This is **not an equal-error speed comparison**: rank caps, polynomial
degrees, projection choices, residual allowances and arithmetic accounting
differ. S7's explicit Gram/positive-arithmetic checks and basis operations
cost time. Its wider probability enclosure is reported, not suppressed.

## Joined finite protocol

Using the original preparation's paid formation lower bound, S4's conditional
nearby-response theorem gives a response per original trial at t=.005 of
at least approximately **0.00071824**. `summary.json` contains its downward
endpoint, evaluated using a rational exponential-series enclosure. Source
variance is .25; rejected formation trials contribute zero. All accepted
records retain their actual conditional receiver state. This is an analytic
composition with S4, not a new 1600-mode probe propagation, an improvement on
S5's sharper joined bound or evidence for an all-size spatial phase.

## Six-event result: rank/error stop before the target

After completed controls, five-event recovery and fresh baselines, the
frozen six-event root component ran with counts (2,3), offsets +/-1/4,
g=.9, L=2304 and required T=144. The other equally weighted root is its
proved S5 sign/relabel mirror. It would reuse a completed result, not turn
an incomplete result into a completed one.

The run stopped after **905 of 2340** physical steps at
**t=55.6923076923**, with reason **contact action cannot fit remaining
residual allowance**. It reached rank **28** and **2,358,511** coefficients.
Projection cost was **0.34994075** of .35; complete-emission residual cost
was **0.01889162** of .02. Seven steps paid appreciable uncaptured action
after the rank cap. The next attempted contact action could not fit the
remaining allowance; silently dropping that action is prohibited.

The ledger at the stopped endpoint totals about .446527, including .077324
of incurred splitting error. This does **not** fit the missing evolution to
T=144 into the envelope: no target-time probability interval is reported.
Raw intermediate populations are archived only as diagnostics with the
explicit incomplete-target notice. They are not a formation success or a
low-formation result.

Whole-worker cost was **227.232 wall seconds**, **227.168 CPU seconds**,
and **716.684 MiB** peak RSS. The 600-second and 2048-MiB soft limits were
not reached. The binding constraints were rank and cumulative error, not
the wall budget. No budget extension, easier observation time or additional
favorable family was launched.

At rank 28 the density workspace alone contains 10,029,768 complex entries;
RESOURCE_MODEL.md derives its count and the full coefficient count. These
explain a cost mechanism but do not replace the measured whole-process
peak. Exact action closure is available when the frame contains contact;
at a fixed rank, repeated mode projection can exhaust its allowance. S7
therefore demonstrates neither affordable general action enclosure nor a
new completed six-event formation target.

## Verdict and next boundary

S5 remains the sharpest completed formation frontier. S7 adds a positive
residual oracle and a checked shared-frame action-enclosure implementation,
recovers a solved target under an explicit arithmetic model, and measures
its real tradeoffs. It moves the bottleneck from a cancellation-prone
residual to preserving useful receiver information at affordable rank/error.

The next construction should test a grouped/history-conditioned or other
receiver representation against this positive reference, while keeping all
histories, reverse amplitudes and returns. Its first gate should again be
the solved target, with total error, time and memory measured. Enlarging
rank or silently relaxing projection tolerance is not the earned result
of this campaign. The whole compiler/BLAS arithmetic contract and the
mathematical bounds still invite independent review.

A separate **post-campaign** derivation in NEXT_RESIDUAL.md provides a
concrete refinement before another large run. With R=sqrt(M-k), the exact
shared-frame residual obeys `D*D=(g*sigma/N)^2 R^2`. The resulting step bound
uses the positive remaining-free-event weight and an h^2 commutator term,
rather than charging the worst assembly size throughout the step. It is
checked against two independent small full-Fock matrices, with eight
propagated packets, in WEIGHTED_CONTROLS.json. It was not used in S7's
trajectories or substituted into their accounts. Neither this refinement
nor its small controls proves that rank 28 or the projection budget suffices
at T=144. EXECUTION.md records the separate control-harness correction;
the primary source and outputs are unchanged.

## Scientific scope

This construction keeps the single-end path grammar. It earns neither
selected spatial dimension nor a common metric for original CP-12 charged
matter. A conditional recording/probe consequence may use the earned
formation lower bound and S4's uniform conditional theorem, retaining the
actual receiver and paying rejected trials; recording and probe injection
remain supplied interventions. It is not a direct full receiver-plus-probe
simulation or a permanently locked phase.

Uniform independent-disorder formation, affordable all-size ranks, a
real-energy gap, gravitational continuum dynamics and a distinctive
empirical prediction remain open. No favorable reinterpretation of the
frozen CP-12 locality failure is made.
