# WM2-S5: history-controlled finite formation, with the scaling gate open

S5 supplies the missing computational ingredient: complete attachment
history with an explicit state-error account for the physical receiver.
Under a separately frozen probability-certificate extension, the tested
five-event balanced-disorder families have formation lower bounds above
77% and 80% at the original fixed observation time. The original S5
state-error-envelope stops remain in the record. The six-event target
remains unresolved under the declared compute budget.

This is one unchanged S4 microscopic law and two declared numerical
certificates. It is not a new Hamiltonian or an accumulation of S2's
spectral conclusion with later formation results. Real-energy coupled
poles remain, and all-size independent-disorder formation is still open.

## What is now computed without erasing history

Keep every ordered detuning-type word and its full reverse transition.
Only the receiver's active one-particle frame is compressed. Its columns
retain all L_N=64N^2 physical frequencies; the frame is not a replacement
receiver with fewer physical modes. All charge-allowed active occupations
are kept. No extra photon-number cutoff, subset-only reduction, fitted
decay rate, bath reset or Markov law is introduced.

Free receiver motion is the exact finite-chain unitary. Assisted
emission/absorption takes place in a frame containing contact, up to a
paid 1e-13 threshold. Each subsequent receiver projection contributes its
complete discarded state norm to an accumulated error bound. Unitary
telescoping pays its entire future effect, including every possible
return. The physical calculation therefore does not claim locality or
formation by silently deleting inconvenient paths.

A GRUT-specific commutator estimate removes the extensive receiver norm
from the product-formula constant. For g=.9 and the frozen |delta_i|<=.75,

    C = ||[A,[A,B]]||_upper/24 + ||[B,[B,A]]||_upper/12
      = 0.3666267585889911.

The full-H state error after n Strang steps at T is at most CT^3/n^2,
plus every receiver projection, exponential remainder and recorded
normalization/basis change. Detailed proofs are in MATHEMATICS.md.
Generic active-rank growth and computational complexity remain unproved.

## Original S5 and the separate P1 certificate

The original charter was published before S5 trajectories at
64378573e8251224dd563928667e961a0effb3da. It targets state error <=.20,
with Strang <=.04 and projections <=.15. It completes the four-event
families. Five-event components exhaust the rank/projection envelope;
six-event components reach the per-instance compute limit before T_N.
Their intermediate populations do not count as target-time results.

The P1 charter was separately published before P1 trajectories at
364781b841121e3169041b121e661048fd793524. Its reason is mathematical:
if ||psi-phi||<=e and the computed formation probability is p_hat, then

\[
 p_{\rm physical}\ge
 \max\{0,1-(\sqrt{1-p_{\rm hat}}+e)^2\}.
\]

This follows by bounding the complementary, unformed event. For nearly
certain formation it can be substantially stronger than S5's original
(sqrt(p_hat)-e)_+^2 bound. P1 allows state error <=.60, Strang <=.20 and
projections <=.35, while retaining the same rank/coefficient/compute
ceilings and the same required formation lower bound >=1/2.

The physical Hamiltonian, seed, root weights, disorder families, geometry
event and observation T_N=4N^2 do not change. Completed S5 states are
reused; all unresolved primary components are attempted. A P1 success
is not a success under S5's original .20 full-state envelope. Both
ledgers and every stop are preserved separately.

## Finite-instance results

Each balanced family has offsets +/-1/4 with counts differing by at
most one, common g=.9, and the indicated bar_delta. All label
arrangements with those counts are equivalent by event permutation;
the original root-type weights are retained. There is no chosen root,
preferred full ordering or favorable receiver postselection.

The following are conservative evaluated lower bounds, truncated down
from AUDITED_RESULTS.json. They include the complete numerical reduction
account; they are not outward-rounded arithmetic certificates.

| N | bar_delta | Original S5 envelope | P1 macro formation lower | P1 full-N assembly lower |
|---:|---:|---|---:|---:|
| 4 | 0 | Complete | 0.98479 | 0.96394 |
| 4 | 1/2 | Complete | 0.98033 | 0.95348 |
| 5 | 0 | Rank/projection limit | 0.77101 | 0.68830 |
| 5 | 1/2 | Rank/projection limit | 0.80630 | 0.72173 |
| 6 | 0 | Compute limit | Unresolved | Unresolved |
| 6 | 1/2 | Compute limit | Unresolved | Unresolved |

At N=5 the physical receiver has 1,600 modes and T_N=100 in model units.
The macroscopic threshold is m>=3; the separate full-assembly column
requires all five events. Thus the positive result is not merely the
probability of forming one bond. Every final ordering and field is summed.

The six-event runs do not reach the requested T_N=144 within the
120-second step-boundary budget. Their target-time formation probability
is unknown. A high intermediate population cannot fill that gap, and
these stops are not a physical no-formation theorem.

Secondary original-S5 controls include a four-event fully distinct ramp
and a five-event outlier family. The ramp completes with the original
lower bound above .7536. The outlier has unresolved original root weight
.8; the remaining weight is paid without renormalization, yielding only
a partial lower bound about .1843 under that original protocol. It does
not pass the original complete-family gate.

## Finite continuous neighborhoods, not the all-size disorder gate

A full-H Duhamel comparison covers genuinely distinct detunings near
each solved finite instance. After removing a conserved root phase,
storage perturbations <=h per site and coupling perturbation <=q add
at most T[(N-1)h+q] to the error account. Applying the probability
certificate again gives explicit finite comparison boxes.

For N=5, conservative simultaneous radii are approximately:

| Center | Per-site storage radius | Common-g radius |
|---|---:|---:|
| bar_delta=0, balanced +/-1/4 | 0.0001327 | 0.0005310 |
| bar_delta=1/2, balanced +/-1/4 | 0.0001589 | 0.0006356 |

Intersect these boxes with the frozen physical perturbation class.
Within them the site detunings need not repeat. No broken type-word
isometry is assumed for those new points; the comparison acts in the
complete physical Hilbert space. These small finite-N boxes do not
cover every width-1/4 pattern and do not establish a nonshrinking
all-size stability region.
They are analytic comparisons with solved centers, not claims that
every point was simulated or met the same numerical work budget.

## A single preparation-recording-probe protocol

The formation lower bound can now be joined to S4 for the solved finite
families. Record every macroscopic word at T_N, inject the adjacent
root-pair neutral probe, retain the actual conditional field, and use
the complete one-probe Hamiltonian. S4 is uniform in that field, so
there is no need to know or refresh its detailed state.

For accepted records, the conditional storage-error and response-feedback
bounds remain those of S4. Giving rejected trials zero signal, the
nearby response per original trial at probe time .005 is at least
.00102794836247... times the formation lower bound. The two tested
five-event families therefore give respective signal lower bounds
above 0.0007925 and 0.0008288. Formation acceptance is paid explicitly.

Recording and probe preparation remain supplied interventions, with
S4's work/entropy account. The source variance is 1/4. This is a finite
operational protocol, not a claim of all-N formation or an equilibrium
spatial phase. The path class and its dimension are grammar-supplied.
Original charged-matter universality and a gravitational continuum are
not established.
The joined response bounds are analytical consequences of S4; S5 did
not propagate the enormous complete 1,600-mode one-probe matrix.

## Verification and limitations

Complete small matrices check the bosonic frame transformations,
commutator constants and full wavefunction error. A separate complete
physical-receiver reference at N=3 uses all 576 modes and dimension
333,505. The active calculation's full wavefunction error is about
.01692, below its analytic account .11336. This comparison uses the
entire two-photon tensor, not only a population or norm check.

Uniform full-physical-receiver controls recover the known S3 quadratic
solution. Complete basis populations test the global-sign identity.
There are 10,201 analytic two-dimensional checks of the probability
inequality. Raw trajectories, stop reasons, root weights, error
components and audited interpretations are archived with executable code.
The largest active coefficient dimension observed, 974,051, occurs in
an unfinished six-event run. It is not a completed physical-receiver
matrix diagonalization or a target-time formation result.

The earliest raw records omitted printing the conservative contact
threshold contribution. Those files remain verbatim. AUDITED_RESULTS
adds g T 1e-13 before interpreting any probability, at most 1.296e-11
here. No gate verdict changes. New executions print it directly.

Arithmetic is ordinary Python/NumPy/SciPy double precision. Analytic
reduction bounds exclude an outward-rounded floating-point proof.
Complete-matrix agreement is a computational cross-check, not external
mathematical validation. The two numerical envelopes must not be
collapsed into one stronger full-state result.

The earned advance is controlled finite formation with attachment
history retained, plus its finite joined response protocol. The central
remaining construction problem is a uniform formation lower bound at
increasing size over the frozen independent-disorder class. Coupled
real-energy poles, selected dimension, original-matter universality and
distinctive empirical predictions are not resolved by this computation.
