# GRUT: a working finite-time charged-memory reduction

Both shared chats were recovered and their research handoff checked against
TestingGRUT. The constructive sprint is reproduced at `c2120bc`; its next
research note is `3538c9a`. This packet continues that note with an executed
memory equation, exact short-time coefficients and a controlled reservoir
reduction. The frozen law and original input are unchanged.

The useful new result is a calculable bridge from neutral formation to
charged dynamics that keeps the return amplitudes. It is a finite validation
of the method needed for the larger construction, not a newly demonstrated
spatial phase.

## Formation and return must be calculated together

For a charged-event projector P, the exact amplitude obeys

\[
\dot x=-iPHPx-iPHZ e^{-iZHZ T}\psi_0
-\int_0^T PHZ e^{-iZHZ(T-s)}ZHPx(s)\,ds,
\quad Z=I-P.
\]

This source is supplied by the original neutral ray, not an external drive.
No measurement is performed by this mathematical projection. At N4,Q4,
the event is a charged spanning circulation, the same event used in the
sprint's recorded prototype before its diagnostic measurement. All 15
neutral graph parents remain in the preparation.

| Time | Full charged-event probability | First-transfer approximation | Ratio |
|---:|---:|---:|---:|
| 2 | 0.002036660074 | 0.004563467043 | 2.24 |
| 5 | 0.055022107666 | 0.243061498653 | 4.42 |
| 10 | 0.021947331541 | 0.666110220171 | 30.35 |

The first-transfer approximation includes source and internal propagation
but deletes feedback to the neutral/complement sector. It is a diagnostic
triangular evolution, not a normalized alternative physical model and not
a Markov approximation. These comparisons show why the existing positive
onset coefficients cannot be extrapolated as an irreversible formation law.
No claim of memory decay or thermal rates is made.

## An explicit small interface reproduces the full amplitude

The charged/complement coupling has rank seven in the budget-256 invariant
model: 86 charged coordinates, 491 complement coordinates, 577 total.
The second-kind Volterra history equation therefore uses seven components.
Its source and kernel are computed from the Hamiltonian without rate fitting.
The number of channels depends on the cutoff and is not a new physical
constant. Full spectral preprocessing is still required in this benchmark.

Independent direct exponentials and a four-grid study validate the history
solver. At time 5 its absolute probability errors are:

| Step size | Absolute probability error |
|---:|---:|
| 0.04 | 2.652e-5 |
| 0.02 | 6.620e-6 |
| 0.01 | 1.654e-6 |
| 0.005 | 4.135e-7 |

The amplitude error also falls by approximately four on halving the step.
A degenerate two-level Rabi control and an independent complex rank-two
matrix control pass. This is numerical convergence, not a rigorous bound on
Volterra quadrature error.

## A controlled reservoir realization keeps the memory

Keeping all 86 charged coordinates and 215 complement Krylov coordinates
gives a 301-dimensional Hermitian realization. At time 5 it returns
p=0.055022107665965. Its directly measured finite-model state error is
1.153e-8. The full unbounded-law Duhamel state-error upper estimate is
1.957e-8, before the separate 5e-10 reporting allowance. This includes every
outgoing electric transition, not just the retained matrix residual.

The corresponding conservative probability interval is approximately
[0.05502209825, 0.05502211708]. At this same time direct budget-256 evolution
has a tighter interval, approximately
[0.05502210735, 0.05502210798]. Budget 128 agrees in the central result while
having a wider full-law truncation estimate. These quantities use ordinary
floating-point evaluations of analytical bounds, not formal interval rounding.

This is a modest validated reduction, 577 to 301 in the invariant sector.
An earlier 163-dimensional reservoir already reproduces this one probability
to 3e-11, while its whole-state error is about 0.373. Therefore matching one
observable is insufficient acceptance evidence; the residual criterion must
remain part of the method. No general large-system efficiency claim follows.

## Exact boundary data and the first return sign

For any such sourced projection whose first target amplitude is a*T^d,
the leading difference between full and first-transfer probabilities is

\[
-\frac{2\|C^\dagger a\|^2}{(d+1)(d+2)}T^{2d+2}.
\]

Its sign is nonpositive. This is a short-time consequence of unitary return,
not an all-time suppression theorem. An independent full unregulated integer
calculation at N4,Q4 obtains

\[
p(T)=\frac{1682}{40186125}T^{18}
-\frac{10368151}{128595600000}T^{20}+O(T^{22}),
\]
\[
p(T)-p_{\rm first\ transfer}(T)
=-\frac{3364}{736745625}T^{20}+O(T^{22}).
\]

The first amplitude occurs at order nine; return feedback first changes it
at order eleven. The exact return-norm ratio is six, independently verifying
the sign formula. This provides boundary data for the memory equation while
the earlier N5-N7 coefficients remain the growing-resource short-time targets.

## The theory direction recovered from the chats

GRUT names the research objective. CP-12 is the current candidate microscopic
construction; its survival is not the objective. Its supported mechanism is
matter-to-flux formation, virtual gauge interactions, protection of connections,
and conditional matter propagation. The missing result is still an appreciable
autonomously formed, persistent family with spatial growth and propagation
in the same unfiltered state. Quantum mechanics, the Born rule and Hamiltonian
time remain supplied. No consciousness-based outcome selection is assumed.

One cross-branch result must remain visible: the frozen-coupling ground-state
exclusion at
[`d2fbd05`](https://github.com/ryangrvr/TestingGRUT/commit/d2fbd05c420a209faae618c88fbe4d1902bc999c)
excludes ground concentration on nearly four-regular graphs with subextensive
absolute degree defects at Q=2N. Its report was reread here; its proof and
suite were not rerun in this packet. The controlled heavy-matter organizing
term must not be turned into a frozen-law ground-selection claim. Irregular
spatial candidates and nonequilibrium persistence remain distinct possibilities.

Next, use the sourced memory/reservoir construction for charged loop density
and matter density at Q=2N, without imposing degree-four or a target dimension
as the definition of success. Validate the N5-N7 onset data, demand actual
full-law residual control, and measure formation, return, relational record
persistence, intrinsic growth and density response on those same states.
Compare irregular spatial candidates with random and nonspatial controls.
The present packet makes that next test more concrete; it does not complete it.

## Reproduction

From the repository root, with NumPy and SciPy installed:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/reality_framework/charged_memory_bridge_20261015/exact_onset.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/reality_framework/charged_memory_bridge_20261015/run.py
```

The recorded run used no formation cache. See `RESULTS.json`,
`EXACT_ONSET.json`, `SOURCE_AUDIT.json` and `VERIFICATION.json` for the
executed evidence. Mathematical details are in `DERIVATION.md`.
