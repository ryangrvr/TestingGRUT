# WM2-S4: prepared bulk protection, with disordered formation still open

S4's main result is a preparation-conditioned finite-time theorem.
An interior neutral probe in a recorded macroscopic path cannot quickly
create a storage inconsistency or receive appreciable feedback from
the complement, even under the frozen independent storage disorder.
This controls an operational observation window without restoring S2's
spectral gap. It does not establish large-size disordered formation.

The campaign retains S3's path grammar and receiver, changes only
delta sum q_i to sum delta_i q_i, and fixes
delta_i=bar_delta+epsilon_i, |bar_delta|<=1/2, |epsilon_i|<=1/4,
g in [0.9,1.1]. Local charter commit: 4c4d7a8d590b1f8587a571ece32bf0709ff06d9b,
parent S3: 80f5d6ab2e5d74307a9445cbbe9182ba2d4c4131.

## The result that changes the preparation question

For a recorded word of length m, inject the adjacent root-pair probe
(|w_0>+|w_1>)/sqrt(2), retaining the actual receiver state. Its number
variance is 1/4. Let r=max(1,m-2), nu=2.1, epsilon=1/4 and
tau_r(x)=sum_(k>=r)x^k/k!.

Every first storage error costs at least r primitive geometry/probe
operations, including any retreat of the head. Exact full-H Duhamel
then gives

\[
\begin{aligned}
\|Q U(t)\psi\|&\le\min\{1,D_r(t)\},\\
|\chi_j(t)-\chi_j^P(t)|&\le2D_r(t),\\
D_r(t)&=\frac{\epsilon}{\nu}\tau_r(\nu|t|).
\end{aligned}
\]

The comparison chi^P is a mathematical reference; the physical H
retains every Q path. The response difference bounds all repeated
complement visits and returns, including those whose final record
is again consistent. This is stronger than counting final Q records.
For |t|<=r/(4e nu), D_r<=(epsilon/nu)exp(-3r/4).

For every accepted macroscopic record, r>=max(1,ceil(N/2)-2).
The following are conservative analytic bounds, not trajectory data:

| N | Minimum r | Conditional observation time | Q probability upper | Response feedback upper |
|---:|---:|---:|---:|---:|
| 32 | 14 | 0.6131 | 1.08e-11 | 6.56e-6 |
| 64 | 30 | 1.3139 | 4.06e-22 | 4.03e-11 |
| 128 | 62 | 2.7153 | 5.79e-43 | 1.53e-21 |

These bounds hold for all allowed deterministic patterns and all
admissible conditional receiver states. The source must be in the
interior; an end-supported probe has r=1 and lacks this distance gain.
The lifetime is proportional to N while formation takes O(N^2).
The theorem does not supply an equilibrium or indefinite spatial phase,
or a positive nearby response at every time in that entire window.
S3's strictly positive nearby response at t=0.005 remains valid.

## What is learned about coupled poles

For any QHQ spectral window I, the actual P-evolved driving measure

    mu_s(I)=||1_I(QHQ) B U_P(s)psi||^2

is at most epsilon^2 T_(r-1)(nu|s|)^2, where T is the clipped Dyson
tail. Actual Q occupation in that window is also bounded by D_r^2.
This addresses finite-time access, not an absent pole or small total
full-H energy-window weight. Small spectral checks have no near-zero
Q eigenvalues and cannot validate the large-N coupled-pole behavior.

A separate receiver estimate shows that the specific S3 far-field
Weyl witness is not freely available to the original preparation.
Its far-half receiver support begins at 32N^2, beyond the receiver
front during the declared formation window. The probability of any
boson there is bounded by a term of order
(N-1)exp(-26.39356 N^2), divided by the geometry acceptance probability
after conditioning. This does not exclude the entire Q energy-window
subspace: energy projection can alter the witness's mode support.

Indeed the coupled-pole proof extends to every allowed pattern at
sufficiently large N. The real-energy gap remains absent. S4 does not
repair it through terminology or by selecting an inaccessible basis.

## Independent disorder: exact reduction, limited formation evidence

The S3 radial/binomial invariant subspace fails. A uniform one-attachment
state leaks from that subspace with norm exactly the standard deviation
of the non-root detunings. Equal attached subsets are also insufficient:
two different attachment orders already differ at the seed's H^3 term.

An exact reduction survives within classes of equal detuning. It retains
the complete sequence of detuning classes, with append coefficient
-(g/N)sqrt(n_a-k_a) a_0^dagger. This includes every reverse, coherent
history and allowed receiver occupation. With all distinct detunings it
is essentially the original factorial word space, so it is not yet an
efficient all-size solver or a new formation theorem.

The first physical-receiver controls test all binary corner patterns
at N=2,3, bar_delta=-1/2,0,1/2, g=0.9,1.1 and times 4N^2,6N^2,8N^2.
Original uniform root weights are retained. There is no optimized root,
word or receiver selection. Results across the non-root parameter sets:

| N | Minimum macro-assembly probability | Minimum full-assembly probability |
|---:|---:|---:|
| 2 | 0.99829844 | 0.99829844 |
| 3 | 0.99999548 | 0.99575086 |

At these small sizes the macroscopic threshold is only m>=2. These are
finite observations; endpoints do not prove all intermediate patterns,
the whole observation window, or increasing-size formation.

The physical receiver has 256 or 576 modes. The calculations use an
explicitly bounded cut at 128 or 224 modes, without a photon cutoff.
The largest analytic probability error is 2.75e-23, excluding roundoff.
A direct full-576-mode N=3 comparison, dimension 333,505, agrees with
the cut-state evolution to 1.10e-26 in observed state norm at t=72.
The complete energy mean and variance, geometry-record entropy and
mean recording work are saved. Record probabilities include all fields.

Separate two-mode controls at N=3 through 6 include common, balanced,
outlier, ramp and fixed-seed iid patterns. Their pronounced variation
and recurrences illustrate receiver dependence; they do not approximate
the physical L_N model at those late times. The iid sample is not a
high-probability statement. No monotonicity in disorder is claimed.

Only a weak general continuity guarantee is obtained for formation:
w<=1/[160N^2(N-1)] suffices to retain p_macro>=1/2. Its width shrinks
as N^-3. It therefore does not pass the frozen nonshrinking w=1/4 gate.
That gate is OPEN, not disproved by the broken radial reduction.

## Verification and scientific status

The archived S3 Node driver and its full suite ran under Node v24.19.0
in a disposable copy. All 17 archived files stayed unchanged. The new
run compared 5,230 numeric leaves; the largest difference was 4.13e-13.
This is execution consistency, not independent mathematical validation.

Python 3.12 / NumPy 2.3.5 / SciPy 1.17 independently constructed the
new sparse matrices. Checks include 6,874 complete matrix entries
against archived JavaScript, 195 first-error source distances, 12 full
wavefunction type-word isometries, 90 physical-receiver observations,
72 root-averaged corner records, 24 independently propagated phase-pulse
checks, and four large-mode analytic packet checks. Full probe matrices
reach dimension 15,600. Arithmetic is ordinary double precision;
there is no outward-rounded certificate or external referee validation.

| Gate | S4 status |
|---|---|
| Archived S3 Node reproduction | PASS |
| Conditional bulk consistency and full response-feedback bounds | Proved for the frozen independent-detuning class |
| Actual receiver access | Explicit preparation/acceptance-dependent bound |
| Exact reduction with disorder | Type-word history retained; radial/binomial law rejected |
| Nonshrinking-disorder, all-size formation | OPEN |
| Uniform real-energy consistency-sector gap | Absent; pole construction retained |
| Selected dimension / original CP-12 matter universality | Not established |
| Empirical or distinctive fundamental prediction | None |

The next constructive target is large-size formation with independent
detuning in the physical receiver, retaining history and controlling
its numerical reduction. The useful bulk lifetime theorem can be
carried forward as a conditional result; it must not be used to assume
the missing formation probability.
